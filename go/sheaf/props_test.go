package sheaf

// Property tests on seeded random sheaves (mixed stalk dims, multi-edges, zero maps, isolated nodes). They do not
// need the oracle: they check identities that must hold for any correct implementation.

import (
	"fmt"
	"math"
	"math/rand/v2"
	"testing"
)

func randomSpec(r *rand.Rand) *Spec {
	n := 1 + r.IntN(9)
	dims := make([]int, n)
	ids := make([]string, n)
	for i := range dims {
		dims[i] = r.IntN(4)
		ids[i] = fmt.Sprintf("n%d", i)
	}
	var edges []Edge
	if n > 1 {
		for j := 0; j < r.IntN(2*n+2); j++ {
			u := r.IntN(n)
			v := r.IntN(n - 1)
			if v >= u {
				v++
			}
			de := r.IntN(4)
			mk := func(c int) [][]float64 {
				m := make([][]float64, de)
				for i := range m {
					m[i] = make([]float64, c)
					for k := range m[i] {
						if r.IntN(5) > 0 {
							m[i][k] = r.NormFloat64() * math.Pow(10, float64(r.IntN(5)-2))
						}
					}
				}
				return m
			}
			edges = append(edges, Edge{U: u, V: v, De: de, Fu: mk(dims[u]), Fv: mk(dims[v]), Label: fmt.Sprintf("e%d", j)})
		}
	}
	s := &Spec{IDs: ids, Dims: dims, Edges: edges, Opts: DefaultOptions()}
	s.Opts.NullReps = 3
	on := Offsets(dims)
	s.X = make([]float64, on[n])
	for i := range s.X {
		s.X[i] = r.NormFloat64()
	}
	if r.IntN(2) == 0 {
		_, _, oe := Coboundary(dims, edges)
		s.W = make([]float64, oe[len(edges)])
		for i := range s.W {
			s.W[i] = r.NormFloat64()
		}
	}
	return s
}

func forRandom(t *testing.T, count int, fn func(t *testing.T, s *Spec)) {
	for seed := uint64(0); seed < uint64(count); seed++ {
		s := randomSpec(rand.New(rand.NewPCG(99, seed)))
		t.Run(fmt.Sprint(seed), func(t *testing.T) { fn(t, s) })
	}
}

func cochain(s *Spec, d0 *CSR) []float64 {
	if s.W != nil {
		return s.W
	}
	return d0.Mul(s.X)
}

// delta0 applied to every H0 basis vector is ~0 (relative to the component's scale and cut), the basis is
// orthonormal, and its size equals the reported H0 dim.
func TestPropGlobalSectionsAreClosed(t *testing.T) {
	forRandom(t, 300, func(t *testing.T, s *Spec) {
		d0, _, _ := Coboundary(s.Dims, s.Edges)
		comps := SpectralComponents(s.Dims, s.Edges, s.Opts.H0, s.Opts.K, s.Opts.DenseMax)
		sp := Merge(s.Dims, comps, s.Opts.H0, s.Opts.K)
		basis := H0Basis(s.Dims, comps)
		if len(basis) != sp.H0() {
			t.Fatalf("basis %d vs H0 %d", len(basis), sp.H0())
		}
		cutMax := 0.0
		for _, c := range comps {
			cutMax = math.Max(cutMax, c.Cut)
		}
		for i, g := range basis {
			dg := d0.Mul(g)
			if e := dot(dg, dg); e > cutMax*(1+1e-6)+1e-24 {
				t.Errorf("section %d: ||delta g||^2 = %g > cut %g", i, e, cutMax)
			}
			for j, h := range basis[:i+1] {
				want := 0.0
				if i == j {
					want = 1
				}
				if math.Abs(dot(g, h)-want) > 1e-10 {
					t.Errorf("basis not orthonormal at (%d,%d): %g", i, j, dot(g, h))
				}
			}
		}
	})
}

// Hodge: exact ⟂ harmonic, ||w||^2 = ||exact||^2 + ||harmonic||^2, exact = delta0 potential, harmonic is co-closed,
// and potential ⟂ H0 (minimum norm).
func TestPropHodgeOrthogonalAndEnergySplits(t *testing.T) {
	forRandom(t, 300, func(t *testing.T, s *Spec) {
		d0, _, _ := Coboundary(s.Dims, s.Edges)
		comps := SpectralComponents(s.Dims, s.Edges, s.Opts.H0, s.Opts.K, s.Opts.DenseMax)
		w := cochain(s, d0)
		h := Hodge(s.Dims, s.Edges, comps, d0, w)
		ww := math.Max(dot(w, w), 1e-300)
		if e := math.Abs(dot(h.Exact, h.Harmonic)) / ww; e > 1e-12 {
			t.Errorf("exact·harmonic / ||w||^2 = %g", e)
		}
		if e := math.Abs(dot(w, w)-dot(h.Exact, h.Exact)-dot(h.Harmonic, h.Harmonic)) / ww; e > 1e-12 {
			t.Errorf("energy split error %g", e)
		}
		dp := d0.Mul(h.Potential)
		for i := range dp {
			if math.Abs(dp[i]-h.Exact[i]) > 1e-8*math.Sqrt(ww) {
				t.Errorf("delta0 potential != exact at %d: %g vs %g", i, dp[i], h.Exact[i])
				break
			}
		}
		for _, g := range H0Basis(s.Dims, comps) {
			if math.Abs(dot(g, h.Potential)) > 1e-8*math.Max(1, norm(h.Potential)) {
				t.Errorf("potential not orthogonal to H0: %g", dot(g, h.Potential))
			}
		}
		if s.W == nil { // w = delta x is exact by construction (A9), up to the directions the policy counts as zero:
			// each dropped singular direction contributes sigma^2 c^2 <= cut * ||x||^2
			cutMax := 0.0
			for _, c := range comps {
				cutMax = math.Max(cutMax, c.Cut)
			}
			if e := dot(h.Harmonic, h.Harmonic); e > cutMax*dot(s.X, s.X)*(1+1e-9)+1e-24*ww {
				t.Errorf("harmonic energy of an exact cochain %g > cut*||x||^2 = %g", e, cutMax*dot(s.X, s.X))
			}
		}
	})
}

// Energy x^T L0 x = sum of per-edge discord, H1 = dim C1 - rank, and the full spectrum is non-negative with
// exactly H0 entries at or below the cut.
func TestPropEnergyDiscordRankNullity(t *testing.T) {
	forRandom(t, 300, func(t *testing.T, s *Spec) {
		rep, err := Analyze(s)
		if err != nil {
			t.Fatal(err)
		}
		d0, _, oe := Coboundary(s.Dims, s.Edges)
		dx := d0.Mul(s.X)
		sum := 0.0
		for j := range s.Edges {
			sum += blockNorm2(dx, oe, j)
		}
		xlx := dot(s.X, applyL(d0, s.X))
		if math.Abs(sum-xlx) > 1e-10*math.Max(1, xlx) {
			t.Errorf("sum discord %g vs x^T L x %g", sum, xlx)
		}
		if rep.EnergyXLX == nil || math.Abs(*rep.EnergyXLX-sig(sum)) > 1e-10*math.Max(1, sum) {
			t.Errorf("report energy %v vs %g", rep.EnergyXLX, sum)
		}
		if rep.H1.Dim != d0.Rows-(d0.Cols-rep.H0.Dim) {
			t.Errorf("rank-nullity: h1 %d, C1 %d, C0 %d, h0 %d", rep.H1.Dim, d0.Rows, d0.Cols, rep.H0.Dim)
		}
		comps := SpectralComponents(s.Dims, s.Edges, s.Opts.H0, s.Opts.K, s.Opts.DenseMax)
		full := FullSpectrum(comps)
		if len(full) != d0.Cols {
			t.Errorf("spectrum size %d vs dim C0 %d", len(full), d0.Cols)
		}
		lb := d0.NormBound()
		for _, l := range full {
			if l < -1e-12*math.Max(1, lb) || l > lb*(1+1e-12)+1e-12 {
				t.Errorf("eigenvalue %g outside [0, %g]", l, lb)
			}
		}
	})
}

// Relabeling: permuting node order (and reversing edge order) leaves H0, H1, the spectrum, energy, shares and
// every per-edge / per-node value unchanged.
func TestPropInvariantUnderRelabeling(t *testing.T) {
	forRandom(t, 200, func(t *testing.T, s *Spec) {
		r := rand.New(rand.NewPCG(5, uint64(len(s.Dims))))
		n := len(s.Dims)
		perm := r.Perm(n) // new index of old node i
		p := &Spec{Opts: s.Opts, IDs: make([]string, n), Dims: make([]int, n)}
		on := Offsets(s.Dims)
		for i := 0; i < n; i++ {
			p.IDs[perm[i]], p.Dims[perm[i]] = s.IDs[i], s.Dims[i]
		}
		pon := Offsets(p.Dims)
		p.X = make([]float64, len(s.X))
		for i := 0; i < n; i++ {
			copy(p.X[pon[perm[i]]:pon[perm[i]+1]], s.X[on[i]:on[i+1]])
		}
		_, _, oe := Coboundary(s.Dims, s.Edges)
		ne := len(s.Edges)
		for j := ne - 1; j >= 0; j-- {
			e := s.Edges[j]
			p.Edges = append(p.Edges, Edge{U: perm[e.U], V: perm[e.V], De: e.De, Fu: e.Fu, Fv: e.Fv, Label: e.Label})
			if s.W != nil {
				p.W = append(p.W, s.W[oe[j]:oe[j+1]]...)
			}
		}
		a, err1 := Analyze(s)
		b, err2 := Analyze(p)
		if err1 != nil || err2 != nil {
			t.Fatal(err1, err2)
		}
		if a.H0.Dim != b.H0.Dim || a.H1.Dim != b.H1.Dim || a.Input.Components != b.Input.Components {
			t.Fatalf("H0/H1/components changed: %d/%d/%d vs %d/%d/%d", a.H0.Dim, a.H1.Dim, a.Input.Components, b.H0.Dim, b.H1.Dim, b.Input.Components)
		}
		lb := math.Max(1, a.Spectral.LambdaMaxBound)
		for i := range a.Spectral.L0Smallest {
			if math.Abs(a.Spectral.L0Smallest[i]-b.Spectral.L0Smallest[i]) > 1e-9*lb {
				t.Errorf("eigenvalue %d: %g vs %g", i, a.Spectral.L0Smallest[i], b.Spectral.L0Smallest[i])
			}
		}
		if math.Abs(a.Hodge.HarmonicShare-b.Hodge.HarmonicShare) > 1e-9 {
			t.Errorf("harmonic share %g vs %g", a.Hodge.HarmonicShare, b.Hodge.HarmonicShare)
		}
		ws := math.Max(1, a.Hodge.Norm2)
		if a.EnergyXLX != nil && math.Abs(*a.EnergyXLX-*b.EnergyXLX) > 1e-9*ws {
			t.Errorf("energy %g vs %g", *a.EnergyXLX, *b.EnergyXLX)
		}
		byLabel := func(rows []EdgeRow) map[string]float64 {
			m := map[string]float64{}
			for _, r := range rows {
				m[r.Label] = r.Discord
			}
			return m
		}
		am, bm := byLabel(a.EdgesByDiscord), byLabel(b.EdgesByDiscord)
		for k, v := range am {
			if w, ok := bm[k]; ok && math.Abs(v-w) > 1e-9*ws {
				t.Errorf("discord %s: %g vs %g", k, v, w)
			}
		}
		nd := map[string]float64{}
		for _, r := range a.NodesByIncidentDiscord {
			nd[r.Node] = r.IncidentDiscord
		}
		for _, r := range b.NodesByIncidentDiscord {
			if v, ok := nd[r.Node]; ok && math.Abs(v-r.IncidentDiscord) > 1e-9*ws {
				t.Errorf("incident discord %s: %g vs %g", r.Node, v, r.IncidentDiscord)
			}
		}
	})
}

// Harmonic extension: boundary values are kept, the interior satisfies L_II x_I = -L_IB x_B, and with an empty
// interior nothing is solved. Relative H0 with an empty boundary equals H0.
func TestPropExtensionAndRelative(t *testing.T) {
	forRandom(t, 200, func(t *testing.T, s *Spec) {
		d0, on, _ := Coboundary(s.Dims, s.Edges)
		comps := SpectralComponents(s.Dims, s.Edges, s.Opts.H0, s.Opts.K, s.Opts.DenseMax)
		sp := Merge(s.Dims, comps, s.Opts.H0, s.Opts.K)
		rel, err := RelativeH0(d0, nil, s.Opts.H0, s.Opts.DenseMax)
		if err != nil {
			t.Fatal(err)
		}
		// B = {}: relative H0 splits delta0 into coupled blocks (finer than graph components when zero maps or
		// de = 0 edges decouple it), and each block's cut uses its own, smaller, bound, so dim <= H0. The reference
		// behaves identically (pinned by testdata/specs/adv_relative_decoupled.json); its doc's "equals H0" holds
		// only when no singular value lies between the two cuts.
		if rel.Dim > sp.H0() {
			t.Errorf("relative H0 with B = {} is %d > H0 %d", rel.Dim, sp.H0())
		}
		for i, g := range rel.Basis {
			dg := d0.Mul(g)
			if e := dot(dg, dg); e > sp.CutMax*(1+1e-6)+1e-24 {
				t.Errorf("relative section %d not closed: %g", i, e)
			}
		}
		var bd []BoundaryDof
		for i := range s.Dims {
			if i%2 == 0 {
				for a := on[i]; a < on[i+1]; a++ {
					bd = append(bd, BoundaryDof{a, float64(a) - 1.5})
				}
			}
		}
		h := HarmonicExtension(d0, bd)
		isB := map[int]bool{}
		for _, b := range bd {
			isB[b.Dof] = true
			if h.X[b.Dof] != b.Value {
				t.Fatalf("boundary dof %d changed", b.Dof)
			}
		}
		lx := applyL(d0, h.X)
		scale := math.Max(1, d0.NormBound()*math.Max(1, norm(h.X)))
		for i := range lx {
			if !isB[i] && math.Abs(lx[i]) > 1e-8*scale {
				t.Errorf("(L x)_I[%d] = %g", i, lx[i])
			}
		}
		if !h.Converged {
			t.Errorf("extension not converged: residual %g route %s", h.ResidualRel, h.Route)
		}
	})
}

// The report JSON has the reference's formatting rules (spot checks of the float formatter).
func TestFloatFormat(t *testing.T) {
	cases := map[float64]string{0: "0.0", 1: "1.0", 9: "9.0", -2.5: "-2.5", 1e16: "1e16", 1.5e16: "1.5e16", 1e15: "1000000000000000.0",
		123456789012.0: "123456789012.0", 1e-5: "0.00001", 1e-6: "1e-6", 2.24e-5: "0.0000224", 0.0001: "0.0001", 1.97215226305e-31: "1.97215226305e-31",
		0.38196601125: "0.38196601125", 4e-10: "4e-10", 954915028.125: "954915028.125", math.NaN(): "null"}
	for x, want := range cases {
		if got := fmtFloat(x); got != want {
			t.Errorf("fmtFloat(%v) = %s, want %s", x, got, want)
		}
	}
	if fmtExp(1e-10) != "1e-10" || fmtExp(1.5e-20) != "1.5e-20" || fmtExp(0) != "0e0" || fmtExp(1e75) != "1e75" {
		t.Errorf("fmtExp: %s %s %s %s", fmtExp(1e-10), fmtExp(1.5e-20), fmtExp(0), fmtExp(1e75))
	}
}
