package sheaf

import (
	"errors"
	"fmt"
	"math"
	"sort"
	"strconv"
)

// HodgeSplit is the orthogonal split of an edge cochain on a graph (no 2-cells): w = exact + harmonic, with
// exact = delta0 potential and potential the minimum-norm least-squares solution (orthogonal to H0).
type HodgeSplit struct {
	Exact, Harmonic, Potential []float64
	// CertCoclosed = ||delta0^T harmonic|| / ||w|| (the harmonic part must be co-closed).
	CertCoclosed float64
}

// Hodge splits w component by component with a thin SVD of delta0|_c truncated to rank(delta0|_c) as decided
// by THAT component's eigen cut (ADR-015), so H0, H1 and the split all use the one per-component rule.
func Hodge(dims []int, edges []Edge, comps []ComponentSpectrum, d0 *CSR, w []float64) HodgeSplit {
	on := Offsets(dims)
	des := make([]int, len(edges))
	for j, e := range edges {
		des[j] = e.De
	}
	oe := Offsets(des)
	h := HodgeSplit{Exact: make([]float64, len(w)), Potential: make([]float64, d0.Cols)}
	for ci := range comps {
		c := &comps[ci]
		if c.Isolated || c.Dofs == 0 || c.RankAll == 0 {
			continue
		}
		var rows, cols []int
		for _, j := range c.Edges {
			for r := oe[j]; r < oe[j+1]; r++ {
				rows = append(rows, r)
			}
		}
		for _, i := range c.Nodes {
			for q := on[i]; q < on[i+1]; q++ {
				cols = append(cols, q)
			}
		}
		s, u, v := svdThin(c.d.Dense(), len(rows), len(cols))
		r := min(c.RankAll, len(s))
		for k := 0; k < r; k++ {
			coef := 0.0
			for i, row := range rows {
				coef += u[i][k] * w[row]
			}
			for i, row := range rows {
				h.Exact[row] += coef * u[i][k]
			}
			for q, col := range cols {
				h.Potential[col] += v[q][k] * coef / s[k]
			}
		}
	}
	h.Harmonic = make([]float64, len(w))
	for i := range w {
		h.Harmonic[i] = w[i] - h.Exact[i]
	}
	h.CertCoclosed = norm(d0.MulT(h.Harmonic)) / math.Max(norm(w), 1e-300)
	return h
}

type share struct{ exact, harmonic, cert float64 }

func splitShares(dims []int, edges []Edge, comps []ComponentSpectrum, d0 *CSR, w []float64) (share, HodgeSplit) {
	g := Hodge(dims, edges, comps, d0, w)
	t := math.Max(dot(w, w), 1e-300)
	return share{dot(g.Exact, g.Exact) / t, dot(g.Harmonic, g.Harmonic) / t, g.CertCoclosed}, g
}

// sig rounds to 12 significant digits (report precision); NaN/inf become NaN (serialised as null).
func sig(x float64) float64 {
	if x == 0 {
		return 0
	}
	if math.IsNaN(x) || math.IsInf(x, 0) {
		return math.NaN()
	}
	f, _ := strconv.ParseFloat(strconv.FormatFloat(x, 'e', 11, 64), 64)
	return f
}

func sigp(x float64) *float64 { v := sig(x); return &v }

// xorshift is the reference's null-model generator (same seed mixing, same sequence).
type xorshift struct{ s uint64 }

func (r *xorshift) u() uint64 {
	s := r.s
	s ^= s << 13
	s ^= s >> 7
	s ^= s << 17
	r.s = s
	return s
}

func (r *xorshift) gauss() float64 {
	a := (float64(r.u()>>11) + 0.5) / float64(uint64(1)<<53)
	b := float64(r.u()>>11) / float64(uint64(1)<<53)
	return math.Sqrt(-2*math.Log(a)) * math.Cos(2*math.Pi*b)
}

func quant(v []float64, q float64) float64 {
	if len(v) == 0 {
		return math.NaN()
	}
	s := append([]float64(nil), v...)
	sort.Float64s(s)
	return s[int(math.Round(float64(len(s)-1)*q))]
}

func blockNorm2(v []float64, off []int, j int) float64 {
	s := 0.0
	for _, x := range v[off[j]:off[j+1]] {
		s += x * x
	}
	return s
}

const (
	wording     = "certified residual: each eigenvalue lies within +-residual of a true eigenvalue (Weyl). The zero/nonzero labelling is a numerical policy (zero_rule), not a geometric kernel certificate; completeness of the Krylov search (no missed eigenvalue below the cut) is NOT certified; known-kernel dims are exact by construction."
	shuffleNote = "degenerate = every shuffle gives the observed share (e.g. <2 edges per stalk-dim bucket, or a single cycle where shuffling preserves the loop sum): uninformative, not 'not significant'"
	extNote     = "values = the minimum-norm minimiser of ||delta x||^2 with x_B fixed; if not unique it is determined only up to H0(X,B;F)"
)

var assumptions = []string{"A1-A8 of docs/rust-cohomology-design.md", "A9 w := delta x when only x is given (harmonic share 0 by construction)",
	"A10 v1 is graph-only (no 2-cells): no coexact part", "A11 potential = minimum-norm solution of L0 x = d0^T w (orthogonal to H0); it is in node-stalk units, so with raw-unit stalks and scaling maps prefer nodes_by_incident_discord"}

// ErrBadOption is returned by Analyze for k = 0 (only reachable through the --k flag; the spec clamps k to >= 1).
// The reference panics there (exit 101); this package refuses it as an input error instead.
var ErrBadOption = errors.New("P_INPUT_01 flag --k")

// Analyze is the pure core: spec -> report. Dense per-component numerics, no I/O.
func Analyze(spec *Spec) (*Report, error) {
	o := spec.Opts
	if o.K < 1 {
		return nil, fmt.Errorf("%w: must be >= 1", ErrBadOption)
	}
	d0, on, oe := Coboundary(spec.Dims, spec.Edges)
	n0, n1 := d0.Cols, d0.Rows
	_, knownOK, ncomp := Components(spec.Dims, spec.Edges)
	comps := SpectralComponents(spec.Dims, spec.Edges, o.H0, o.K, o.DenseMax)
	sp := Merge(spec.Dims, comps, o.H0, o.K)
	h0 := sp.H0()
	h1 := max(n1+h0-n0, 0)

	// L1 = d0 d0^T check on the mapped nonzero pairs u = d0 v / ||d0 v||
	l1Res, l1Checked := 0.0, 0
	for j := sp.FoundZero; j < len(sp.Pairs); j++ {
		u := d0.Mul(sp.Vectors[j])
		nu := norm(u)
		if nu < 1e-300 {
			continue
		}
		for i := range u {
			u[i] /= nu
		}
		y := d0.Mul(d0.MulT(u))
		s := 0.0
		for i := range y {
			r := y[i] - sp.Pairs[j].Value*u[i]
			s += r * r
		}
		l1Res = math.Max(l1Res, math.Sqrt(s))
		l1Checked++
	}
	k := o.K
	small := make([]float64, min(sp.KnownDim, k+1))
	for _, p := range sp.Pairs {
		small = append(small, p.Value)
	}
	sort.Float64s(small)
	if len(small) > k+1 {
		small = small[:k+1]
	}
	var res []float64
	for _, p := range sp.Pairs[max(sp.FoundZero-(k+1), 0):] {
		res = append(res, sig(p.Residual))
	}
	var gapK *float64
	if len(small) > k {
		gapK = sigp(small[k] - small[k-1])
	}

	// cochain + Hodge split
	var w []float64
	cochain := "w (given edge cochain)"
	if spec.W != nil {
		w = append([]float64(nil), spec.W...)
	} else {
		w = d0.Mul(spec.X)
		cochain = "delta x (derived from node data; exact by construction, A9)"
	}
	hs, g := splitShares(spec.Dims, spec.Edges, comps, d0, w)
	wn2 := dot(w, w)

	// nulls: (a) iid Gaussian edge data of the same shape; (b) edge blocks shuffled within equal stalk dims
	rng := xorshift{o.Seed*0xD1B54A32D192ED03 | 1}
	var nr, ns []float64
	byDe := map[int][]int{}
	for j, e := range spec.Edges {
		byDe[e.De] = append(byDe[e.De], j)
	}
	deKeys := make([]int, 0, len(byDe))
	for d := range byDe {
		deKeys = append(deKeys, d)
	}
	sort.Ints(deKeys)
	for rep := 0; rep < o.NullReps; rep++ {
		r := make([]float64, n1)
		for i := range r {
			r[i] = rng.gauss()
		}
		sh, _ := splitShares(spec.Dims, spec.Edges, comps, d0, r)
		nr = append(nr, sh.harmonic)
		perm := make([]int, len(spec.Edges))
		for i := range perm {
			perm[i] = i
		}
		for _, d := range deKeys {
			grp := byDe[d]
			gg := append([]int(nil), grp...)
			for i := len(gg) - 1; i >= 1; i-- {
				j := int(rng.u() % uint64(i+1))
				gg[i], gg[j] = gg[j], gg[i]
			}
			for t, a := range grp {
				perm[a] = gg[t]
			}
		}
		sw := make([]float64, n1)
		for j := range spec.Edges {
			src := perm[j]
			for q := 0; q < oe[j+1]-oe[j]; q++ {
				sw[oe[j]+q] = w[oe[src]+q]
			}
		}
		if dot(sw, sw) > 0 {
			sh, _ := splitShares(spec.Dims, spec.Edges, comps, d0, sw)
			ns = append(ns, sh.harmonic)
		}
	}
	ge := func(v []float64) *float64 {
		if len(v) == 0 {
			return nil
		}
		c := 0
		for _, x := range v {
			if x >= hs.harmonic-1e-12 {
				c++
			}
		}
		return sigp(float64(c) / float64(len(v)))
	}
	mean := func(v []float64) *float64 {
		if len(v) == 0 {
			return nil
		}
		s := 0.0
		for _, x := range v {
			s += x
		}
		return sigp(s / float64(len(v)))
	}
	analytic := 0.0
	if n1 > 0 {
		analytic = float64(h1) / float64(n1)
	}
	degenerate := len(ns) == 0
	if !degenerate {
		degenerate = true
		for _, x := range ns {
			if !(math.Abs(x-hs.harmonic) <= 1e-12*math.Max(math.Abs(hs.harmonic), 1e-300)+1e-15) {
				degenerate = false
			}
		}
	}

	// rankings
	type er struct {
		j    int
		h, t float64
	}
	ers := make([]er, len(spec.Edges))
	for j := range spec.Edges {
		ers[j] = er{j, blockNorm2(g.Harmonic, oe, j), blockNorm2(w, oe, j)}
	}
	edgeRow := func(x er) EdgeRow {
		e := spec.Edges[x.j]
		frac := 0.0
		if wn2 > 0 {
			frac = x.h / wn2
		}
		return EdgeRow{Edge: x.j, U: spec.IDs[e.U], V: spec.IDs[e.V], Label: e.Label, HarmonicEnergy: sig(x.h), Discord: sig(x.t), HarmonicFracOfTotal: sig(frac)}
	}
	byH := append([]er(nil), ers...)
	sort.SliceStable(byH, func(a, b int) bool {
		if byH[a].h != byH[b].h {
			return byH[a].h > byH[b].h
		}
		return byH[a].j < byH[b].j
	})
	edgesH := []EdgeRow{}
	for _, x := range byH {
		if len(edgesH) >= o.Top {
			break
		}
		if x.h > 1e-24*math.Max(wn2, 1e-300) {
			edgesH = append(edgesH, edgeRow(x))
		}
	}
	byD := append([]er(nil), byH...)
	sort.SliceStable(byD, func(a, b int) bool {
		if byD[a].t != byD[b].t {
			return byD[a].t > byD[b].t
		}
		return byD[a].j < byD[b].j
	})
	edgesD := []EdgeRow{}
	for _, x := range byD {
		if len(edgesD) >= o.Top {
			break
		}
		if x.t > 0 {
			edgesD = append(edgesD, edgeRow(x))
		}
	}
	type nv struct {
		i int
		v float64
	}
	rank := func(v []nv) []nv {
		sort.SliceStable(v, func(a, b int) bool {
			if v[a].v != v[b].v {
				return v[a].v > v[b].v
			}
			return v[a].i < v[b].i
		})
		out := []nv{}
		for _, x := range v {
			if len(out) >= o.Top {
				break
			}
			if x.v > 0 {
				out = append(out, x)
			}
		}
		return out
	}
	pots := make([]nv, len(spec.Dims))
	for i := range spec.Dims {
		pots[i] = nv{i, math.Sqrt(blockNorm2(g.Potential, on, i))}
	}
	nodesP := []NodePotential{}
	for _, x := range rank(pots) {
		nodesP = append(nodesP, NodePotential{Node: spec.IDs[x.i], AbsPotential: sig(x.v)})
	}
	inc := make([]float64, len(spec.Dims))
	for _, x := range byH {
		inc[spec.Edges[x.j].U] += x.t
		inc[spec.Edges[x.j].V] += x.t
	}
	incs := make([]nv, len(spec.Dims))
	for i, v := range inc {
		incs[i] = nv{i, v}
	}
	nodesI := []NodeDiscord{}
	for _, x := range rank(incs) {
		nodesI = append(nodesI, NodeDiscord{Node: spec.IDs[x.i], IncidentDiscord: sig(x.v)})
	}
	var energy *float64
	if spec.X != nil {
		dx := d0.Mul(spec.X)
		energy = sigp(dot(dx, dx))
	}
	var zmax, nzmin *Pair
	if sp.FoundZero > 0 {
		p := sp.Pairs[sp.FoundZero-1]
		zmax = &Pair{sig(p.Value), sig(p.Residual)}
	}
	var gap *float64
	if sp.FoundZero < len(sp.Pairs) {
		p := sp.Pairs[sp.FoundZero]
		nzmin = &Pair{sig(p.Value), sig(p.Residual)}
		gap = sigp(p.Value)
	}
	var margin *float64
	if sp.Margin != nil {
		margin = sigp(*sp.Margin)
	}
	known := 0
	for _, b := range knownOK {
		if b {
			known++
		}
	}
	smallSig := make([]float64, len(small))
	for i, x := range small {
		smallSig[i] = sig(x)
	}
	rep := &Report{
		Input: Input{Nodes: len(spec.Dims), Edges: len(spec.Edges), DimC0: n0, DimC1: n1, Components: ncomp, KnownConstantComponents: known,
			HasX: spec.X != nil, HasW: spec.W != nil},
		H0: H0Report{Dim: h0, KnownKernelDim: sp.KnownDim, FoundZeroBySolver: sp.FoundZero, ZeroPolicy: o.H0.Name(), ZeroRule: o.H0.Describe(),
			ZeroTol: sig(sp.CutMax), ZeroTolMin: sig(sp.CutMin), LargestCountedZero: zmax, SmallestCountedNonzero: nzmin,
			CutConsistentWithZeroTol: sp.Consistent,
			NearZeroAmbiguous:        !sp.Consistent || sp.NearSections || (sp.Margin != nil && *sp.Margin < 10),
			CutMarginRatio:           margin, LowerBoundOnly: sp.LowerBound, Wording: wording},
		H1: H1Report{Dim: h1, Rule: "rank-nullity on a graph: dim C1 - (dim C0 - h0); inherits the h0 cut"},
		Spectral: SpectralReport{Solver: "dense-evd", L0Smallest: smallSig, SolverPairsResiduals: res, GapK: gapK, SpectralGap: gap,
			LambdaMaxBound: sig(d0.NormBound()), Rounds: sp.Rounds, Converged: true,
			L1Graph: L1Graph{KernelDim: h1, NonzeroSpectrum: "equals L0 nonzero spectrum (d0^T d0 vs d0 d0^T)", MappedPairsChecked: l1Checked, MaxMappedResidual: sig(l1Res)}},
		Hodge: HodgeReport{Cochain: cochain, Norm2: sig(wn2), ExactShare: sig(hs.exact), HarmonicShare: sig(hs.harmonic), CgIters: 0,
			CertCoclosed: sig(hs.cert), Converged: !math.IsNaN(hs.cert) && hs.cert <= 1e-8,
			NullRandomGaussian: NullRandom{Reps: len(nr), Mean: mean(nr), P05: sig(quant(nr, 0.05)), P95: sig(quant(nr, 0.95)), FracGeObserved: ge(nr), AnalyticExpectationH1C1: sig(analytic)},
			NullShuffledValues: NullShuffle{Degenerate: degenerate, Note: shuffleNote, Reps: len(ns), Mean: mean(ns), P05: sig(quant(ns, 0.05)), P95: sig(quant(ns, 0.95)), FracGeObserved: ge(ns)}},
		EnergyXLX:              energy,
		EdgesByHarmonicEnergy:  edgesH,
		EdgesByDiscord:         edgesD,
		NodesByAbsPotential:    nodesP,
		NodesByIncidentDiscord: nodesI,
		Assumptions:            assumptions,
	}
	if spec.Boundary != nil {
		rep.Extension = extension(spec, d0, on)
	}
	return rep, nil
}

func extension(spec *Spec, d0 *CSR, on []int) *ExtensionReport {
	o := spec.Opts
	var bd []BoundaryDof
	isb := map[int]bool{}
	for _, b := range spec.Boundary {
		isb[b.Node] = true
		for a, v := range b.Value {
			bd = append(bd, BoundaryDof{on[b.Node] + a, v})
		}
	}
	bdofs := make([]int, len(bd))
	for i, b := range bd {
		bdofs[i] = b.Dof
	}
	rel, relErr := RelativeH0(d0, bdofs, o.H0, o.DenseMax)
	h := HarmonicExtension(d0, bd)
	e := &ExtensionReport{BoundaryNodes: len(spec.Boundary), BoundaryDofs: len(bd),
		UniquenessRule: "unique iff dim H0(X,B;F) = dim ker L_II = 0 (Hansen-Ghrist Prop 4.1); zero cut: " + o.H0.Describe(),
		Note:           extNote, Route: h.Route, Converged: h.Converged, ResidualRel: sig(h.ResidualRel), Iters: 0, Energy: sig(h.Energy)}
	if relErr == nil {
		d := rel.Dim
		u := d == 0
		e.H0Relative, e.Unique = &d, &u
	}
	for i := range spec.Dims {
		v := make([]float64, on[i+1]-on[i])
		for t := range v {
			v[t] = sig(h.X[on[i]+t])
		}
		e.Values = append(e.Values, NodeValue{ID: spec.IDs[i], Boundary: isb[i], Value: v})
	}
	return e
}
