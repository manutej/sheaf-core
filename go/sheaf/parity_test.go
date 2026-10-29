package sheaf

// Parity against golden reports produced by the Rust reference (`sheaf analyze` + the go_oracle example), committed
// under testdata/golden. Contract (DECISIONS.md ADR-018):
//   - H0/H1 dims exact, unless the oracle sets near_zero_ambiguous: then the flag must match (and dims are reported).
//   - eigenvalues |dλ| <= 1e-9 * max(1, ||L||_bound); energy, discord, harmonic, potential: per entry
//     |d| <= 1e-9 * max(1, ||w||^2) (scale of the quantity); shares and null statistics: absolute 1e-9.
//   - kernel bases (H0, relative H0) compared through their orthogonal projectors.
//   - error cases: same exit code and error kind; for spec errors the same message.
// Every comparison records its error; TestMain prints the worst observed error per quantity.

import (
	"bytes"
	"encoding/json"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"sync"
	"testing"
)

const tol = 1e-9

var (
	worstMu sync.Mutex
	worst   = map[string]float64{}
	counts  = map[string]int{}
)

func record(q string, e float64) {
	worstMu.Lock()
	defer worstMu.Unlock()
	if w, seen := worst[q]; !seen || e > w || math.IsNaN(e) {
		worst[q] = e
	}
	counts[q]++
}

func TestMain(m *testing.M) {
	code := m.Run()
	if len(worst) > 0 {
		var ks []string
		for k := range worst {
			ks = append(ks, k)
		}
		sort.Strings(ks)
		fmt.Println("parity: worst observed error per quantity (tolerance 1e-9, scaled as in parity_test.go)")
		for _, k := range ks {
			fmt.Printf("  %-28s %-10.3g over %d\n", k, worst[k], counts[k])
		}
	}
	os.Exit(code)
}

func readJSON(t *testing.T, path string) map[string]any {
	t.Helper()
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	var v map[string]any
	if err := json.Unmarshal(b, &v); err != nil {
		t.Fatalf("%s: %v", path, err)
	}
	return v
}

func goldenArgs(t *testing.T) []string {
	b, err := os.ReadFile("../testdata/golden/args.txt")
	if err != nil {
		t.Fatal(err)
	}
	return strings.Fields(string(b))
}

// runCLI runs Main from the go/ directory (paths in goldens are relative to it).
func runCLI(args []string) (int, string, string) {
	wd, _ := os.Getwd()
	os.Chdir("..")
	defer os.Chdir(wd)
	var o, e bytes.Buffer
	code := Main(args, &o, &e)
	return code, o.String(), e.String()
}

func f(v any) float64 {
	if v == nil {
		return math.NaN()
	}
	return v.(float64)
}

func at(v any, path ...string) any {
	for _, p := range path {
		v = v.(map[string]any)[p]
	}
	return v
}

func fs(v any) []float64 {
	a, _ := v.([]any)
	out := make([]float64, len(a))
	for i, x := range a {
		out[i] = f(x)
	}
	return out
}

// relErr = |a-b| / max(1, scale); both null (NaN) counts as equal.
func relErr(a, b, scale float64) float64 {
	if math.IsNaN(a) && math.IsNaN(b) {
		return 0
	}
	return math.Abs(a-b) / math.Max(1, scale)
}

func check(t *testing.T, q string, e float64, what string) {
	t.Helper()
	record(q, e)
	if !(e <= tol) {
		t.Errorf("%s: %s error %.3g > %g", what, q, e, tol)
	}
}

func projector(basis [][]float64, n int) [][]float64 {
	p := make([][]float64, n)
	for i := range p {
		p[i] = make([]float64, n)
	}
	for _, b := range basis {
		for i := 0; i < n; i++ {
			if b[i] == 0 {
				continue
			}
			for j := 0; j < n; j++ {
				p[i][j] += b[i] * b[j]
			}
		}
	}
	return p
}

func projErr(a, b [][]float64, n int) float64 {
	pa, pb := projector(a, n), projector(b, n)
	m := 0.0
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			m = math.Max(m, math.Abs(pa[i][j]-pb[i][j]))
		}
	}
	return m
}

func basisOf(v any) [][]float64 {
	a, _ := v.([]any)
	out := make([][]float64, len(a))
	for i, c := range a {
		out[i] = fs(c)
	}
	return out
}

// sameKeys: the Go report has exactly the reference's key set (recursively, for objects; first element of arrays).
func sameKeys(t *testing.T, path string, a, b any) {
	t.Helper()
	am, aok := a.(map[string]any)
	bm, bok := b.(map[string]any)
	if aok != bok {
		t.Errorf("%s: object vs non-object", path)
		return
	}
	if aok {
		var ka, kb []string
		for k := range am {
			ka = append(ka, k)
		}
		for k := range bm {
			kb = append(kb, k)
		}
		sort.Strings(ka)
		sort.Strings(kb)
		if strings.Join(ka, ",") != strings.Join(kb, ",") {
			t.Errorf("%s: keys %v vs reference %v", path, ka, kb)
			return
		}
		for k := range am {
			if am[k] != nil && bm[k] != nil {
				sameKeys(t, path+"."+k, am[k], bm[k])
			}
		}
		return
	}
	aa, aok := a.([]any)
	ba, bok := b.([]any)
	if aok && bok && len(aa) > 0 && len(ba) > 0 {
		sameKeys(t, path+"[0]", aa[0], ba[0])
	}
}

// fieldsDocumentedToDiffer: solver bookkeeping of the reference's CG / Lanczos routes. The Go port is dense, so
// these are reported but not compared (ADR-018): spectral.solver/rounds/inner_cg/solver_pairs_residuals,
// hodge.cg_iters, extension.route/iters/residual_rel, *.max_mapped_residual, h0 pair residuals.
func TestGoldenReports(t *testing.T) {
	specs, _ := filepath.Glob("../testdata/specs/*.json")
	args := goldenArgs(t)
	n := 0
	for _, spec := range specs {
		name := strings.TrimSuffix(filepath.Base(spec), ".json")
		if strings.HasSuffix(name, ".boundary") || strings.HasPrefix(name, "errf_") {
			continue
		}
		rep := "../testdata/golden/" + name + ".report.json"
		if _, err := os.Stat(rep); err != nil {
			continue
		}
		n++
		t.Run(name, func(t *testing.T) {
			cli := append([]string{"analyze", "testdata/specs/" + name + ".json"}, args...)
			bnd := "testdata/specs/" + name + ".boundary.json"
			if _, err := os.Stat("../" + bnd); err == nil {
				cli = append(cli, "--boundary", bnd)
			}
			code, out, stderr := runCLI(cli)
			if code != 0 {
				t.Fatalf("exit %d: %s", code, stderr)
			}
			var g map[string]any
			if err := json.Unmarshal([]byte(out), &g); err != nil {
				t.Fatal(err)
			}
			r := readJSON(t, rep)
			aux := readJSON(t, "../testdata/golden/"+name+".aux.json")
			compareReport(t, name, g, r, aux)
		})
	}
	if n < 50 {
		t.Fatalf("only %d golden reports found", n)
	}
}

func compareReport(t *testing.T, name string, g, r, aux map[string]any) {
	sameKeys(t, "report", g, r)
	// exact fields
	for _, k := range []string{"format", "method"} {
		if g[k] != r[k] {
			t.Errorf("%s: %v vs %v", k, g[k], r[k])
		}
	}
	for k, v := range r["input"].(map[string]any) {
		if g["input"].(map[string]any)[k] != v {
			t.Errorf("input.%s: %v vs %v", k, g["input"].(map[string]any)[k], v)
		}
	}
	rh0, gh0 := r["h0"].(map[string]any), g["h0"].(map[string]any)
	amb := rh0["near_zero_ambiguous"].(bool)
	if gh0["near_zero_ambiguous"] != rh0["near_zero_ambiguous"] {
		t.Errorf("near_zero_ambiguous: %v vs %v", gh0["near_zero_ambiguous"], rh0["near_zero_ambiguous"])
	}
	dimsMatch := gh0["dim"] == rh0["dim"] && at(g, "h1", "dim") == at(r, "h1", "dim")
	if !dimsMatch && !amb {
		t.Errorf("H0/H1: go %v/%v vs reference %v/%v (not flagged ambiguous)", gh0["dim"], at(g, "h1", "dim"), rh0["dim"], at(r, "h1", "dim"))
	}
	record("h0_h1_dims_mismatch", map[bool]float64{true: 0, false: 1}[dimsMatch])
	for _, k := range []string{"known_kernel_dim", "found_zero_by_solver", "zero_policy", "zero_rule", "cut_consistent_with_zero_tol", "lower_bound_only", "wording"} {
		if dimsMatch && gh0[k] != rh0[k] {
			t.Errorf("h0.%s: %v vs %v", k, gh0[k], rh0[k])
		}
	}
	for _, k := range []string{"zero_tol", "zero_tol_min"} {
		check(t, "h0_cut", relErr(f(gh0[k]), f(rh0[k]), f(rh0[k]))*1, "h0."+k)
	}
	// eigenvalues: |dλ| <= 1e-9 max(1, ||L||_bound)
	lb := f(at(aux, "l0_bound"))
	lscale := math.Max(1, lb)
	full := fs(aux["l0_spectrum"])
	gs := fs(at(g, "spectral", "L0_smallest"))
	rs := fs(at(r, "spectral", "L0_smallest"))
	if len(gs) != len(rs) {
		t.Errorf("L0_smallest: %d vs %d values", len(gs), len(rs))
	} else {
		for i := range gs {
			check(t, "eigenvalue", math.Abs(gs[i]-rs[i])/lscale, fmt.Sprintf("L0_smallest[%d]", i))
		}
	}
	// the full spectrum, Go dense vs Rust dense
	if spec, err := os.ReadFile("../testdata/specs/" + name + ".json"); err == nil {
		s, perr := Parse(spec)
		if perr != nil {
			t.Fatal(perr)
		}
		comps := SpectralComponents(s.Dims, s.Edges, s.Opts.H0, s.Opts.K, s.Opts.DenseMax)
		mine := FullSpectrum(comps)
		if len(mine) != len(full) {
			t.Errorf("full spectrum: %d vs %d", len(mine), len(full))
		} else {
			for i := range mine {
				check(t, "eigenvalue", math.Abs(mine[i]-full[i])/lscale, fmt.Sprintf("spectrum[%d]", i))
			}
		}
		// H0 basis via projectors (dense per-component cohomology(0) vs ours)
		n0 := len(fs(aux["potential"]))
		if dimsMatch {
			check(t, "h0_projector", projErr(H0Basis(s.Dims, comps), basisOf(aux["h0_basis"]), n0), "H0 projector")
		}
		// full Hodge split per entry, scaled by max(1, ||w||^2)
		w := fs(aux["w"])
		d0, _, _ := Coboundary(s.Dims, s.Edges)
		if !amb {
			h := Hodge(s.Dims, s.Edges, comps, d0, w)
			ws := math.Max(1, dot(w, w))
			ex, ha := fs(aux["exact"]), fs(aux["harmonic"])
			for i := range w {
				check(t, "exact_entry", math.Abs(h.Exact[i]-ex[i])/ws, "exact")
				check(t, "harmonic_entry", math.Abs(h.Harmonic[i]-ha[i])/ws, "harmonic")
			}
			pot := fs(aux["potential"])
			ps := 1.0
			for _, x := range pot {
				ps = math.Max(ps, math.Abs(x))
			}
			for i := range pot {
				check(t, "potential_entry", math.Abs(h.Potential[i]-pot[i])/ps, "potential")
			}
		}
	}
	// spectral scalars
	for _, k := range []string{"gap_k", "spectral_gap"} {
		check(t, "eigenvalue", relErr(f(at(g, "spectral", k)), f(at(r, "spectral", k)), lscale*0+1)/lscale*1, "spectral."+k)
	}
	check(t, "lambda_max_bound", relErr(f(at(g, "spectral", "lambda_max_bound")), f(at(r, "spectral", "lambda_max_bound")), lb), "lambda_max_bound")
	if at(g, "spectral", "L1_graph", "kernel_dim") != at(r, "spectral", "L1_graph", "kernel_dim") && !amb {
		t.Errorf("L1 kernel dim differs")
	}
	if amb && !dimsMatch {
		return // ambiguous cut and dims differ: the rest depends on the count (flag matched above)
	}
	// energy / Hodge shares
	wn2 := f(at(r, "hodge", "norm2"))
	ws := math.Max(1, wn2)
	if r["energy_xLx"] != nil || g["energy_xLx"] != nil {
		check(t, "energy", relErr(f(g["energy_xLx"]), f(r["energy_xLx"]), f(r["energy_xLx"])), "energy_xLx")
	}
	check(t, "norm2", relErr(f(at(g, "hodge", "norm2")), wn2, wn2), "hodge.norm2")
	if at(g, "hodge", "cochain") != at(r, "hodge", "cochain") {
		t.Errorf("hodge.cochain differs")
	}
	// go's converged flag is its own certificate (cert_coclosed <= 1e-8); with an unambiguous cut it must agree
	// with the reference's
	if (f(at(g, "hodge", "cert_coclosed")) <= 1e-8) != (at(g, "hodge", "converged") == true) {
		t.Errorf("go hodge.converged inconsistent with its certificate")
	}
	if !amb && at(r, "hodge", "converged") == true && at(g, "hodge", "converged") != true {
		t.Errorf("hodge.converged: go false, reference true")
	}
	// The reference report's split is a CG solve without the zero cut. Its fields are compared only when that
	// solve certified (hodge.converged) and the cut is unambiguous; otherwise the oracle for the split is the dense
	// per-component CochainComplex::hodge in the aux file (already compared entry by entry above), and the
	// shares are checked against it here. Every skip is counted (ADR-018).
	if at(r, "hodge", "converged") != true || amb {
		record("cases_split_vs_dense_oracle", 0)
		ex, ha := fs(aux["exact"]), fs(aux["harmonic"])
		t2 := math.Max(dot(fs(aux["w"]), fs(aux["w"])), 1e-300)
		check(t, "share", math.Abs(f(at(g, "hodge", "exact_share"))-sig(dot(ex, ex)/t2)), "exact_share vs dense oracle")
		check(t, "share", math.Abs(f(at(g, "hodge", "harmonic_share"))-sig(dot(ha, ha)/t2)), "harmonic_share vs dense oracle")
		compareRows(t, g, r, "edges_by_discord", "edge", "discord", ws, "discord")
		compareRows(t, g, r, "nodes_by_incident_discord", "node", "incident_discord", ws, "node_discord")
		return
	}
	for _, k := range []string{"exact_share", "harmonic_share"} {
		check(t, "share", math.Abs(f(at(g, "hodge", k))-f(at(r, "hodge", k))), "hodge."+k)
	}
	for _, nk := range []string{"null_random_gaussian", "null_shuffled_values"} {
		gn, rn := at(g, "hodge", nk).(map[string]any), at(r, "hodge", nk).(map[string]any)
		if gn["reps"] != rn["reps"] {
			t.Errorf("%s.reps: %v vs %v", nk, gn["reps"], rn["reps"])
		}
		for _, k := range []string{"mean", "p05", "p95", "frac_ge_observed", "analytic_expectation_h1_over_dimC1"} {
			if _, ok := rn[k]; ok {
				check(t, "null_stat", relErr(f(gn[k]), f(rn[k]), 0), nk+"."+k)
			}
		}
		if d, ok := rn["degenerate"]; ok && d != gn["degenerate"] {
			t.Errorf("%s.degenerate: %v vs %v", nk, gn["degenerate"], d)
		}
	}
	// rankings: compare per-edge / per-node values by identity (sets above the report's thresholds may tie at ~0)
	compareRows(t, g, r, "edges_by_discord", "edge", "discord", ws, "discord")
	compareRows(t, g, r, "edges_by_harmonic_energy", "edge", "harmonic_energy", ws, "harmonic_energy")
	compareRows(t, g, r, "nodes_by_incident_discord", "node", "incident_discord", ws, "node_discord")
	ps := 1.0
	for _, x := range fs(aux["potential"]) {
		ps = math.Max(ps, math.Abs(x))
	}
	compareRows(t, g, r, "nodes_by_abs_potential", "node", "abs_potential", ps, "abs_potential")
	// extension
	if re, ok := r["extension"].(map[string]any); ok {
		ge, ok := g["extension"].(map[string]any)
		if !ok {
			t.Fatalf("extension missing")
		}
		for _, k := range []string{"boundary_nodes", "boundary_dofs", "h0_relative", "unique", "uniqueness_rule", "note", "method"} {
			if ge[k] != re[k] {
				t.Errorf("extension.%s: %v vs %v", k, ge[k], re[k])
			}
		}
		if ge["converged"] != re["converged"] {
			t.Errorf("extension.converged: %v vs %v", ge["converged"], re["converged"])
		}
		gv, rv := ge["values"].([]any), re["values"].([]any)
		vs := 1.0
		for _, x := range rv {
			for _, y := range fs(x.(map[string]any)["value"]) {
				vs = math.Max(vs, math.Abs(y))
			}
		}
		for i := range rv {
			a, b := fs(gv[i].(map[string]any)["value"]), fs(rv[i].(map[string]any)["value"])
			for j := range b {
				check(t, "extension_value", math.Abs(a[j]-b[j])/vs, "extension.values")
			}
		}
		check(t, "extension_energy", relErr(f(ge["energy"]), f(re["energy"]), vs*vs), "extension.energy")
		if rel, ok := aux["relative"].(map[string]any); ok {
			spec, _ := os.ReadFile("../testdata/specs/" + name + ".json")
			s, _ := Parse(spec)
			bt, _ := os.ReadFile("../testdata/specs/" + name + ".boundary.json")
			bv, err := ParseBoundary(bt, s)
			if err != nil {
				t.Fatal(err)
			}
			d0, on, _ := Coboundary(s.Dims, s.Edges)
			var dofs []int
			for _, b := range bv {
				for a := range b.Value {
					dofs = append(dofs, on[b.Node]+a)
				}
			}
			mine, err := RelativeH0(d0, dofs, s.Opts.H0, s.Opts.DenseMax)
			if err != nil {
				t.Fatal(err)
			}
			check(t, "relative_projector", projErr(mine.Basis, basisOf(rel["basis"]), d0.Cols), "relative H0 projector")
		}
	}
}

func compareRows(t *testing.T, g, r map[string]any, key, id, val string, scale float64, q string) {
	t.Helper()
	index := func(rows any) map[string]float64 {
		m := map[string]float64{}
		for _, x := range rows.([]any) {
			o := x.(map[string]any)
			m[fmt.Sprint(o[id])] = f(o[val])
		}
		return m
	}
	gm, rm := index(g[key]), index(r[key])
	// rows near zero may be on one side of the report's threshold and not the other; values must agree
	threshold := 1e-12 * scale
	for k, rvv := range rm {
		gvv, ok := gm[k]
		if !ok {
			if rvv > threshold {
				t.Errorf("%s: %s=%s missing in go (reference %g)", key, id, k, rvv)
			}
			continue
		}
		check(t, q, math.Abs(gvv-rvv)/scale, key)
	}
	for k, gvv := range gm {
		if _, ok := rm[k]; !ok && gvv > threshold {
			t.Errorf("%s: %s=%s only in go (%g)", key, id, k, gvv)
		}
	}
	// ranking: position by position the values agree (ties within tolerance may swap identities)
	gr, rr := g[key].([]any), r[key].([]any)
	for i := 0; i < len(rr) && i < len(gr); i++ {
		a, b := f(gr[i].(map[string]any)[val]), f(rr[i].(map[string]any)[val])
		if b > threshold || a > threshold {
			check(t, q, math.Abs(a-b)/scale, fmt.Sprintf("%s[%d] (rank order)", key, i))
		}
	}
}

func TestGoldenErrors(t *testing.T) {
	errs, _ := filepath.Glob("../testdata/golden/*.error.json")
	flags := map[string][]string{}
	fc, err := os.ReadFile("../testdata/flag_cases.txt")
	if err != nil {
		t.Fatal(err)
	}
	for _, line := range strings.Split(string(fc), "\n") {
		if line = strings.TrimSpace(line); line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		p := strings.SplitN(line, "|", 2)
		flags[strings.TrimSpace(p[0])] = strings.Fields(p[1])
	}
	if len(errs) < 30 {
		t.Fatalf("only %d error goldens", len(errs))
	}
	for _, path := range errs {
		name := strings.TrimSuffix(filepath.Base(path), ".error.json")
		t.Run(name, func(t *testing.T) {
			want := readJSON(t, path)
			var cli []string
			if tag, ok := strings.CutPrefix(name, "errf_"); ok {
				cli = append([]string{"analyze", "testdata/specs/errf_flags.json"}, flags[tag]...)
			} else {
				cli = []string{"analyze", "testdata/specs/" + name + ".json", "--null-reps", "0"}
				if strings.HasPrefix(name, "errb_") {
					cli = append(cli, "--boundary", "testdata/specs/"+name+".boundary.json")
				}
			}
			code, _, stderr := runCLI(cli)
			wantCode := int(want["exit"].(float64))
			wantErr := strings.TrimSpace(want["stderr"].(string))
			if wantCode == 101 { // the reference panics (k = 0); the Go port refuses it as an input error
				if code != 2 || !strings.HasPrefix(stderr, "P_INPUT_01 flag --k") {
					t.Errorf("k=0: want exit 2 P_INPUT_01, got %d %q", code, stderr)
				}
				return
			}
			if code != wantCode {
				t.Fatalf("exit %d, reference %d (%q vs %q)", code, wantCode, stderr, wantErr)
			}
			got := strings.TrimSpace(stderr)
			kind := func(s string) string { // "P_INPUT_01 json:", "P_INPUT_01 spec:", "P_INPUT_01 flag", "P_INPUT_01 read"
				p := strings.Fields(s)
				if len(p) < 2 {
					return s
				}
				return p[0] + " " + strings.TrimSuffix(p[1], ":")
			}
			if kind(got) != kind(wantErr) {
				t.Errorf("error kind %q vs reference %q (%q vs %q)", kind(got), kind(wantErr), got, wantErr)
			}
			// spec errors carry the same message; JSON syntax / OS messages are library-specific
			if strings.HasPrefix(wantErr, "P_INPUT_01 spec:") && got != wantErr {
				t.Errorf("message %q vs reference %q", got, wantErr)
			}
		})
	}
}
