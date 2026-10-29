package sheaf

import (
	"math"
	"sort"
)

// Components labels the connected components of the graph (labels in order of each component's smallest node)
// and flags each component "known constant": every edge has Fu == Fv (same shape, same entries) and injective,
// so the component's constants are exact sections. Isolated nodes form their own (known-constant) components.
func Components(dims []int, edges []Edge) (comp []int, knownConst []bool, count int) {
	n := len(dims)
	p := make([]int, n)
	for i := range p {
		p[i] = i
	}
	find := func(a int) int {
		for p[a] != a {
			p[a] = p[p[a]]
			a = p[a]
		}
		return a
	}
	for _, e := range edges {
		a, b := find(e.U), find(e.V)
		if a != b {
			p[max(a, b)] = min(a, b)
		}
	}
	lab := make([]int, n)
	for i := range lab {
		lab[i] = -1
	}
	comp = make([]int, n)
	for i := 0; i < n; i++ {
		r := find(i)
		if lab[r] < 0 {
			lab[r] = count
			count++
		}
		comp[i] = lab[r]
	}
	knownConst = make([]bool, count)
	for i := range knownConst {
		knownConst[i] = true
	}
	for _, e := range edges {
		c := comp[e.U]
		if !knownConst[c] {
			continue
		}
		same := dims[e.U] == dims[e.V] && len(e.Fu) == len(e.Fv)
		for r := 0; same && r < len(e.Fu); r++ {
			for q := range e.Fu[r] {
				if e.Fu[r][q] != e.Fv[r][q] {
					same = false
					break
				}
			}
		}
		inj := false
		if same {
			s, _ := svdFull(e.Fu, len(e.Fu), dims[e.U])
			inj = numericRank(s, len(e.Fu), dims[e.U]) == dims[e.U]
		}
		if !inj {
			knownConst[c] = false
		}
	}
	return comp, knownConst, count
}

// EigPair is an eigenvalue of L0 with its certified residual ||L v - lambda v|| (v unit): a true eigenvalue lies
// within +-residual of Value (Weyl).
type EigPair struct {
	Value    float64
	Residual float64
}

// ComponentSpectrum is one connected component's dense spectral analysis under the H0 policy.
type ComponentSpectrum struct {
	Nodes []int
	Edges []int
	Dofs  int
	// Isolated: the component has dofs but no edges, so its whole stalk is a section (exact, counted in Known).
	Isolated bool
	// Known is the dim of the deflated known kernel (exact component constants, or the whole stalk if isolated).
	Known int
	// Bound = ||delta_c||_1 ||delta_c||_inf >= ||L0|_c||_2; Cut = the policy's lambda cut for this component.
	Bound, Cut float64
	// Pairs: the reported eigenpairs (ascending, on the complement of the known kernel). For a component with at
	// most DenseMax dofs these are all of them; above it they are the k+1 smallest (widened x4 at most twice while
	// all are zero), exactly the count the reference's iterative route returns.
	Pairs   []EigPair
	Vectors [][]float64 // local (component-dof) eigenvectors of Pairs
	// FoundZero: pairs with Value <= Cut. LowerBound: every reported pair is zero and more exist (H0 is a lower bound).
	FoundZero  int
	LowerBound bool
	// RankAll = rank(delta0|_c) under the cut, counted on the FULL dense spectrum (used by the Hodge split).
	RankAll int
	// Consistent: the cut is consistent with the residuals. NearSect: a counted-zero pair is resolvable from 0.
	Consistent, NearSect bool
	// local data
	d        *CSR
	on       []int
	knownVec [][]float64
	allVals  []float64
	allVecs  [][]float64
}

func subCoboundary(dims []int, edges []Edge, nodes []int, eidx []int) (*CSR, []int, []int) {
	local := make(map[int]int, len(nodes))
	sd := make([]int, len(nodes))
	for li, i := range nodes {
		local[i] = li
		sd[li] = dims[i]
	}
	se := make([]Edge, len(eidx))
	for k, j := range eidx {
		e := edges[j]
		se[k] = Edge{U: local[e.U], V: local[e.V], De: e.De, Fu: e.Fu, Fv: e.Fv}
	}
	return Coboundary(sd, se)
}

func applyL(d *CSR, x []float64) []float64 { return d.MulT(d.Mul(x)) }

func residual(d *CSR, v []float64, lam float64) float64 {
	y := applyL(d, v)
	s := 0.0
	for i := range y {
		r := y[i] - lam*v[i]
		s += r * r
	}
	return math.Sqrt(s)
}

// SpectralComponents runs the per-component analysis. L0 is block diagonal over connected components; each
// component's block gets a dense symmetric EVD (known constants deflated by a shift), and an eigenvalue counts as
// zero iff it is <= the policy's cut for THAT component (ADR-015). No global cut is applied anywhere.
func SpectralComponents(dims []int, edges []Edge, pol H0Policy, k, denseMax int) []ComponentSpectrum {
	comp, known, nc := Components(dims, edges)
	nodesOf := make([][]int, nc)
	for i, c := range comp {
		nodesOf[c] = append(nodesOf[c], i)
	}
	edgesOf := make([][]int, nc)
	for j, e := range edges {
		edgesOf[comp[e.U]] = append(edgesOf[comp[e.U]], j)
	}
	out := make([]ComponentSpectrum, nc)
	for c := 0; c < nc; c++ {
		cs := ComponentSpectrum{Nodes: nodesOf[c], Edges: edgesOf[c], Consistent: true}
		for _, i := range cs.Nodes {
			cs.Dofs += dims[i]
		}
		if cs.Dofs == 0 {
			out[c] = cs
			continue
		}
		if len(cs.Edges) == 0 {
			cs.Isolated = true
			cs.Known = cs.Dofs
			out[c] = cs
			continue
		}
		var son []int
		cs.d, son, _ = subCoboundary(dims, edges, cs.Nodes, cs.Edges)
		cs.on = son
		n := cs.Dofs
		if known[c] {
			nd := len(cs.Nodes)
			f := 1 / math.Sqrt(float64(nd))
			for j := 0; j < dims[cs.Nodes[0]]; j++ {
				q := make([]float64, n)
				for li := 0; li < nd; li++ {
					q[son[li]+j] = f
				}
				cs.knownVec = append(cs.knownVec, q)
			}
		}
		cs.Known = len(cs.knownVec)
		cs.Bound = cs.d.NormBound()
		cs.Cut = pol.ZeroTol(cs.Bound)
		l := make([][]float64, n)
		for i := range l {
			l[i] = make([]float64, n)
		}
		e := make([]float64, n)
		for j := 0; j < n; j++ {
			e[j] = 1
			y := applyL(cs.d, e)
			for i := 0; i < n; i++ {
				l[i][j] = y[i]
			}
			e[j] = 0
		}
		shift := 10*math.Max(cs.Bound, 1) + 10
		for _, q := range cs.knownVec {
			for i := 0; i < n; i++ {
				if q[i] == 0 {
					continue
				}
				for j := 0; j < n; j++ {
					if q[j] != 0 {
						l[i][j] += shift * q[i] * q[j]
					}
				}
			}
		}
		vals, vecs := symEig(l, n)
		free := n - cs.Known
		cs.allVals = vals[:free]
		cs.allVecs = make([][]float64, free)
		for j := 0; j < free; j++ {
			v := make([]float64, n)
			for i := 0; i < n; i++ {
				v[i] = vecs.At(i, j)
			}
			cs.allVecs[j] = v
		}
		for _, x := range cs.allVals {
			if x > cs.Cut {
				cs.RankAll++
			}
		}
		// reported pairs: all of them up to DenseMax dofs; above it the k+1 smallest, widened like the reference
		kl := k
		if n <= denseMax {
			kl = max(n-1-cs.Known, 1)
		}
		for round := 0; ; round++ {
			m := min(kl+1, free)
			cs.Pairs, cs.Vectors = nil, nil
			for j := 0; j < m; j++ {
				cs.Pairs = append(cs.Pairs, EigPair{Value: cs.allVals[j], Residual: residual(cs.d, cs.allVecs[j], cs.allVals[j])})
				cs.Vectors = append(cs.Vectors, cs.allVecs[j])
			}
			cs.FoundZero = 0
			for _, p := range cs.Pairs {
				if p.Value <= cs.Cut {
					cs.FoundZero++
				}
			}
			cs.LowerBound = cs.FoundZero == len(cs.Pairs) && len(cs.Pairs) < free
			if !cs.LowerBound || round == 2 {
				break
			}
			kl *= 4
		}
		if cs.FoundZero > 0 {
			z := cs.Pairs[cs.FoundZero-1]
			cs.Consistent = z.Value+z.Residual <= cs.Cut
		}
		if cs.FoundZero < len(cs.Pairs) {
			nz := cs.Pairs[cs.FoundZero]
			cs.Consistent = cs.Consistent && nz.Value-nz.Residual > cs.Cut
		} else {
			cs.Consistent = cs.Consistent && !cs.LowerBound
		}
		for _, p := range cs.Pairs[:cs.FoundZero] {
			if p.Value > math.Max(64*eps*cs.Bound, 10*p.Residual) {
				cs.NearSect = true
			}
		}
		out[c] = cs
	}
	return out
}

// Spectrum is the merged per-component result (the report's h0/spectral inputs).
type Spectrum struct {
	Comps []ComponentSpectrum
	// KnownDim: exact sections (known constants + isolated stalks). FoundZero: counted-zero pairs over components.
	KnownDim, FoundZero int
	// Pairs: every component's zero pairs (ascending), then the k+1 smallest nonzero pairs overall.
	Pairs []EigPair
	// Vectors: global C0 vectors of the nonzero Pairs (zero pairs have nil vectors).
	Vectors [][]float64
	// CutMin/CutMax: smallest/largest per-component cut; Margin: min over components of smallest nonzero / cut.
	CutMin, CutMax float64
	Margin         *float64
	Consistent     bool
	NearSections   bool
	LowerBound     bool
	Rounds         int
}

// H0 = dim ker L0 under the policy: known sections + counted zeros.
func (s *Spectrum) H0() int { return s.KnownDim + s.FoundZero }

// fminRust mirrors Rust's f64::min (a NaN operand is ignored).
func fminRust(a, b float64) float64 {
	if math.IsNaN(a) {
		return b
	}
	if math.IsNaN(b) {
		return a
	}
	return math.Min(a, b)
}

// Merge combines the components into the report-level spectrum (k = the requested number of nonzero pairs - 1).
func Merge(dims []int, comps []ComponentSpectrum, pol H0Policy, k int) *Spectrum {
	on := Offsets(dims)
	s := &Spectrum{Comps: comps, Consistent: true, CutMin: math.Inf(1)}
	type cand struct {
		p EigPair
		c int
		v []float64
	}
	var zeros []EigPair
	var cands []cand
	for ci := range comps {
		c := &comps[ci]
		if c.Dofs == 0 {
			continue
		}
		if c.Isolated {
			s.KnownDim += c.Dofs
			continue
		}
		s.CutMin = fminRust(s.CutMin, c.Cut)
		s.CutMax = math.Max(s.CutMax, c.Cut)
		s.Consistent = s.Consistent && c.Consistent
		s.NearSections = s.NearSections || c.NearSect
		if c.FoundZero < len(c.Pairs) {
			r := c.Pairs[c.FoundZero].Value / c.Cut
			if s.Margin == nil {
				s.Margin = &r
			} else {
				m := fminRust(*s.Margin, r)
				s.Margin = &m
			}
		}
		s.KnownDim += c.Known
		s.Rounds = 1
		s.LowerBound = s.LowerBound || c.LowerBound
		for j, p := range c.Pairs {
			if p.Value <= c.Cut {
				zeros = append(zeros, p)
			} else {
				cands = append(cands, cand{p, ci, c.Vectors[j]})
			}
		}
	}
	if math.IsInf(s.CutMin, 1) {
		s.CutMin = pol.ZeroTol(0)
		s.CutMax = s.CutMin
	}
	sort.SliceStable(zeros, func(a, b int) bool { return zeros[a].Value < zeros[b].Value })
	sort.SliceStable(cands, func(a, b int) bool {
		if cands[a].p.Value != cands[b].p.Value {
			return cands[a].p.Value < cands[b].p.Value
		}
		return cands[a].c < cands[b].c
	})
	if len(cands) > k+1 {
		cands = cands[:k+1]
	}
	n := on[len(dims)]
	s.FoundZero = len(zeros)
	s.Pairs = append(s.Pairs, zeros...)
	s.Vectors = make([][]float64, len(zeros))
	for _, cd := range cands {
		g := make([]float64, n)
		q := 0
		for _, i := range comps[cd.c].Nodes {
			for t := 0; t < dims[i]; t++ {
				g[on[i]+t] = cd.v[q]
				q++
			}
		}
		s.Pairs = append(s.Pairs, cd.p)
		s.Vectors = append(s.Vectors, g)
	}
	s.Consistent = s.Consistent && !s.LowerBound
	return s
}

// H0Basis returns an orthonormal basis of H0 = ker delta0 under the policy as global C0 vectors: isolated stalks,
// plus, per component, the right singular vectors of delta0|_c with sigma^2 <= that component's cut (SVD of delta0,
// not an EVD of L0, so the basis is not squared-conditioned).
func H0Basis(dims []int, comps []ComponentSpectrum) [][]float64 {
	on := Offsets(dims)
	n := on[len(dims)]
	var out [][]float64
	for ci := range comps {
		c := &comps[ci]
		if c.Dofs == 0 {
			continue
		}
		var cols []int
		for _, i := range c.Nodes {
			for t := on[i]; t < on[i+1]; t++ {
				cols = append(cols, t)
			}
		}
		if c.Isolated {
			for _, col := range cols {
				g := make([]float64, n)
				g[col] = 1
				out = append(out, g)
			}
			continue
		}
		s, v := svdFull(c.d.Dense(), c.d.Rows, c.d.Cols)
		r := 0
		for _, x := range s {
			if x*x > c.Cut {
				r++
			}
		}
		for j := r; j < len(cols); j++ {
			g := make([]float64, n)
			for q, col := range cols {
				g[col] = v.At(q, j)
			}
			out = append(out, g)
		}
	}
	return out
}

// FullSpectrum returns every eigenvalue of L0 (ascending): exact zeros for the known kernel and isolated stalks,
// then each component's full dense spectrum on the complement.
func FullSpectrum(comps []ComponentSpectrum) []float64 {
	var out []float64
	for ci := range comps {
		c := &comps[ci]
		for i := 0; i < c.Known; i++ {
			out = append(out, 0)
		}
		out = append(out, c.allVals...)
	}
	sort.Float64s(out)
	return out
}
