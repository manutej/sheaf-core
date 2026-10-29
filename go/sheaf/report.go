package sheaf

import (
	"math"
	"sort"
	"strconv"
	"strings"
)

// Report mirrors the "frontier-sheaf/report-v1" document field for field (same key names and nesting).
// Optional values are pointers (nil serialises as null); NaN serialises as null.
type Report struct {
	Input                  Input
	H0                     H0Report
	H1                     H1Report
	Spectral               SpectralReport
	Hodge                  HodgeReport
	EnergyXLX              *float64 // omitted when the spec has no x
	EdgesByHarmonicEnergy  []EdgeRow
	EdgesByDiscord         []EdgeRow
	NodesByAbsPotential    []NodePotential
	NodesByIncidentDiscord []NodeDiscord
	Extension              *ExtensionReport // present only with a boundary
	Assumptions            []string
}

type Input struct {
	Nodes, Edges, DimC0, DimC1, Components, KnownConstantComponents int
	HasX, HasW                                                      bool
}

type Pair struct{ Value, Residual float64 }

type H0Report struct {
	Dim, KnownKernelDim, FoundZeroBySolver int
	ZeroPolicy, ZeroRule                   string
	ZeroTol, ZeroTolMin                    float64
	LargestCountedZero                     *Pair
	SmallestCountedNonzero                 *Pair
	CutConsistentWithZeroTol               bool
	NearZeroAmbiguous                      bool
	CutMarginRatio                         *float64
	LowerBoundOnly                         bool
	Wording                                string
}

type H1Report struct {
	Dim  int
	Rule string
}

type L1Graph struct {
	KernelDim          int
	NonzeroSpectrum    string
	MappedPairsChecked int
	MaxMappedResidual  float64
}

type SpectralReport struct {
	Solver               string
	L0Smallest           []float64
	SolverPairsResiduals []float64
	GapK, SpectralGap    *float64
	LambdaMaxBound       float64
	Rounds               int
	Converged            bool
	L1Graph              L1Graph
}

type NullRandom struct {
	Reps                    int
	Mean                    *float64
	P05, P95                float64
	FracGeObserved          *float64
	AnalyticExpectationH1C1 float64
}

type NullShuffle struct {
	Degenerate     bool
	Note           string
	Reps           int
	Mean           *float64
	P05, P95       float64
	FracGeObserved *float64
}

type HodgeReport struct {
	Cochain            string
	Norm2, ExactShare  float64
	HarmonicShare      float64
	CgIters            int
	CertCoclosed       float64
	Converged          bool
	NullRandomGaussian NullRandom
	NullShuffledValues NullShuffle
}

type EdgeRow struct {
	Edge                                         int
	U, V, Label                                  string
	HarmonicEnergy, Discord, HarmonicFracOfTotal float64
}

type NodePotential struct {
	Node         string
	AbsPotential float64
}

type NodeDiscord struct {
	Node            string
	IncidentDiscord float64
}

type NodeValue struct {
	ID       string
	Boundary bool
	Value    []float64
}

type ExtensionReport struct {
	BoundaryNodes, BoundaryDofs int
	H0Relative                  *int
	Unique                      *bool
	UniquenessRule, Note, Route string
	Converged                   bool
	ResidualRel                 float64
	Iters                       int
	Energy                      float64
	Values                      []NodeValue
}

// ---- JSON encoding (sorted keys, 2-space indent, trailing newline: the reference's byte format) ----

type obj map[string]any

func optF(p *float64) any {
	if p == nil {
		return nil
	}
	return *p
}

func pairV(p *Pair) any {
	if p == nil {
		return nil
	}
	return obj{"value": p.Value, "residual": p.Residual}
}

func edgeRows(rs []EdgeRow) []any {
	out := make([]any, len(rs))
	for i, r := range rs {
		out[i] = obj{"edge": r.Edge, "u": r.U, "v": r.V, "label": r.Label, "harmonic_energy": r.HarmonicEnergy,
			"discord": r.Discord, "harmonic_frac_of_total": r.HarmonicFracOfTotal}
	}
	return out
}

func floats(v []float64) []any {
	out := make([]any, len(v))
	for i, x := range v {
		out[i] = x
	}
	return out
}

func (r *Report) value() obj {
	nodesP := make([]any, len(r.NodesByAbsPotential))
	for i, n := range r.NodesByAbsPotential {
		nodesP[i] = obj{"node": n.Node, "abs_potential": n.AbsPotential}
	}
	nodesI := make([]any, len(r.NodesByIncidentDiscord))
	for i, n := range r.NodesByIncidentDiscord {
		nodesI[i] = obj{"node": n.Node, "incident_discord": n.IncidentDiscord}
	}
	asm := make([]any, len(r.Assumptions))
	for i, a := range r.Assumptions {
		asm[i] = a
	}
	sp := r.Spectral
	hd := r.Hodge
	v := obj{
		"format": ReportFormat,
		"method": "exact",
		"input": obj{"nodes": r.Input.Nodes, "edges": r.Input.Edges, "dim_C0": r.Input.DimC0, "dim_C1": r.Input.DimC1,
			"components": r.Input.Components, "known_constant_components": r.Input.KnownConstantComponents,
			"has_x": r.Input.HasX, "has_w": r.Input.HasW},
		"h0": obj{"method": "exact", "dim": r.H0.Dim, "known_kernel_dim": r.H0.KnownKernelDim, "found_zero_by_solver": r.H0.FoundZeroBySolver,
			"zero_policy": r.H0.ZeroPolicy, "zero_rule": r.H0.ZeroRule, "zero_tol": r.H0.ZeroTol, "zero_tol_min": r.H0.ZeroTolMin,
			"largest_counted_zero": pairV(r.H0.LargestCountedZero), "smallest_counted_nonzero": pairV(r.H0.SmallestCountedNonzero),
			"cut_consistent_with_zero_tol": r.H0.CutConsistentWithZeroTol, "near_zero_ambiguous": r.H0.NearZeroAmbiguous,
			"cut_margin_ratio": optF(r.H0.CutMarginRatio), "lower_bound_only": r.H0.LowerBoundOnly, "wording": r.H0.Wording},
		"h1": obj{"method": "exact", "dim": r.H1.Dim, "rule": r.H1.Rule},
		"spectral": obj{"method": "exact", "solver": sp.Solver, "L0_smallest": floats(sp.L0Smallest),
			"solver_pairs_residuals": floats(sp.SolverPairsResiduals), "gap_k": optF(sp.GapK), "spectral_gap": optF(sp.SpectralGap),
			"lambda_max_bound": sp.LambdaMaxBound, "rounds": sp.Rounds, "converged": sp.Converged, "inner_cg": nil,
			"L1_graph": obj{"kernel_dim": sp.L1Graph.KernelDim, "nonzero_spectrum": sp.L1Graph.NonzeroSpectrum,
				"mapped_pairs_checked": sp.L1Graph.MappedPairsChecked, "max_mapped_residual": sp.L1Graph.MaxMappedResidual}},
		"hodge": obj{"method": "exact", "cochain": hd.Cochain, "norm2": hd.Norm2, "exact_share": hd.ExactShare, "harmonic_share": hd.HarmonicShare,
			"cg_iters": hd.CgIters, "cert_coclosed": hd.CertCoclosed, "converged": hd.Converged,
			"null_random_gaussian": obj{"reps": hd.NullRandomGaussian.Reps, "mean": optF(hd.NullRandomGaussian.Mean),
				"p05": hd.NullRandomGaussian.P05, "p95": hd.NullRandomGaussian.P95, "frac_ge_observed": optF(hd.NullRandomGaussian.FracGeObserved),
				"analytic_expectation_h1_over_dimC1": hd.NullRandomGaussian.AnalyticExpectationH1C1},
			"null_shuffled_values": obj{"degenerate": hd.NullShuffledValues.Degenerate, "note": hd.NullShuffledValues.Note,
				"reps": hd.NullShuffledValues.Reps, "mean": optF(hd.NullShuffledValues.Mean), "p05": hd.NullShuffledValues.P05,
				"p95": hd.NullShuffledValues.P95, "frac_ge_observed": optF(hd.NullShuffledValues.FracGeObserved)}},
		"edges_by_harmonic_energy":  edgeRows(r.EdgesByHarmonicEnergy),
		"edges_by_discord":          edgeRows(r.EdgesByDiscord),
		"nodes_by_abs_potential":    nodesP,
		"nodes_by_incident_discord": nodesI,
		"assumptions":               asm,
	}
	if r.EnergyXLX != nil {
		v["energy_xLx"] = *r.EnergyXLX
	}
	if e := r.Extension; e != nil {
		vals := make([]any, len(e.Values))
		for i, nv := range e.Values {
			vals[i] = obj{"id": nv.ID, "boundary": nv.Boundary, "value": floats(nv.Value)}
		}
		var h0r, uq any
		if e.H0Relative != nil {
			h0r = *e.H0Relative
		}
		if e.Unique != nil {
			uq = *e.Unique
		}
		v["extension"] = obj{"method": "exact", "boundary_nodes": e.BoundaryNodes, "boundary_dofs": e.BoundaryDofs,
			"h0_relative": h0r, "unique": uq, "uniqueness_rule": e.UniquenessRule, "note": e.Note, "route": e.Route,
			"converged": e.Converged, "residual_rel": e.ResidualRel, "iters": e.Iters, "energy": e.Energy, "values": vals}
	}
	return v
}

// JSON renders the report deterministically: sorted keys, 2-space indent, shortest round-trip floats in the
// reference's style ("9.0", "1e-31", "4.5e17"), NaN as null, and a trailing newline.
func (r *Report) JSON() string {
	var b strings.Builder
	writeJSON(&b, r.value(), 0)
	b.WriteByte('\n')
	return b.String()
}

func writeJSON(b *strings.Builder, v any, ind int) {
	pad := func(n int) { b.WriteString(strings.Repeat("  ", n)) }
	switch t := v.(type) {
	case nil:
		b.WriteString("null")
	case bool:
		b.WriteString(strconv.FormatBool(t))
	case int:
		b.WriteString(strconv.Itoa(t))
	case float64:
		b.WriteString(fmtFloat(t))
	case string:
		writeString(b, t)
	case []any:
		if len(t) == 0 {
			b.WriteString("[]")
			return
		}
		b.WriteString("[\n")
		for i, x := range t {
			pad(ind + 1)
			writeJSON(b, x, ind+1)
			if i+1 < len(t) {
				b.WriteByte(',')
			}
			b.WriteByte('\n')
		}
		pad(ind)
		b.WriteByte(']')
	case obj:
		if len(t) == 0 {
			b.WriteString("{}")
			return
		}
		keys := make([]string, 0, len(t))
		for k := range t {
			keys = append(keys, k)
		}
		sort.Strings(keys)
		b.WriteString("{\n")
		for i, k := range keys {
			pad(ind + 1)
			writeString(b, k)
			b.WriteString(": ")
			writeJSON(b, t[k], ind+1)
			if i+1 < len(keys) {
				b.WriteByte(',')
			}
			b.WriteByte('\n')
		}
		pad(ind)
		b.WriteByte('}')
	default:
		panic("sheaf: unsupported JSON value")
	}
}

func writeString(b *strings.Builder, s string) {
	b.WriteByte('"')
	for _, r := range s {
		switch r {
		case '"':
			b.WriteString(`\"`)
		case '\\':
			b.WriteString(`\\`)
		case '\n':
			b.WriteString(`\n`)
		case '\r':
			b.WriteString(`\r`)
		case '\t':
			b.WriteString(`\t`)
		case '\b':
			b.WriteString(`\b`)
		case '\f':
			b.WriteString(`\f`)
		default:
			if r < 0x20 {
				b.WriteString(`\u00`)
				b.WriteByte("0123456789abcdef"[r>>4])
				b.WriteByte("0123456789abcdef"[r&15])
			} else {
				b.WriteRune(r)
			}
		}
	}
	b.WriteByte('"')
}

// fmtFloat formats like serde_json (ryu): shortest round-trip digits; plain notation for decimal exponents in
// (-5, 16], with ".0" for integers; otherwise d[.ddd]e<exp>. NaN/inf become null.
func fmtFloat(x float64) string {
	if math.IsNaN(x) || math.IsInf(x, 0) {
		return "null"
	}
	if x == 0 {
		if math.Signbit(x) {
			return "-0.0"
		}
		return "0.0"
	}
	sign := ""
	if x < 0 {
		sign = "-"
		x = -x
	}
	s := strconv.FormatFloat(x, 'e', -1, 64) // d.ddde±XX
	e := strings.IndexByte(s, 'e')
	digits := strings.Replace(s[:e], ".", "", 1)
	exp, _ := strconv.Atoi(s[e+1:])
	n := len(digits)
	kk := exp + 1 // decimal point position relative to the digit string
	k := kk - n
	switch {
	case k >= 0 && kk <= 16:
		return sign + digits + strings.Repeat("0", k) + ".0"
	case kk > 0 && kk <= 16:
		return sign + digits[:kk] + "." + digits[kk:]
	case kk > -5 && kk <= 0:
		return sign + "0." + strings.Repeat("0", -kk) + digits
	case n == 1:
		return sign + digits + "e" + strconv.Itoa(kk-1)
	default:
		return sign + digits[:1] + "." + digits[1:] + "e" + strconv.Itoa(kk-1)
	}
}
