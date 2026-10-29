package sheaf

import (
	"errors"
	"math"
)

// ErrBlockTooLarge: an interior block of relative H0 exceeds DenseMax columns (not computed, never guessed).
var ErrBlockTooLarge = errors.New("relative H0: interior block exceeds dense_max")

// Relative is H0(X,B;F) = ker(delta0 restricted to the interior dofs) = ker L_II: the uniqueness certificate for
// harmonic extension (unique iff Dim == 0).
type Relative struct {
	Dim int
	// Basis: orthonormal basis of H0(X,B;F) as full-length C0 vectors (zero on B).
	Basis [][]float64
	// Cuts: the lambda cut applied to each interior block.
	Cuts []float64
}

// RelativeH0 splits the interior columns of delta0 into blocks linked through shared NONZERO rows and counts, per
// block, singular values with sigma^2 <= policy.ZeroTol(||block||_1 ||block||_inf) as zero (the same single rule as
// absolute H0; with B empty this equals H0).
func RelativeH0(d0 *CSR, boundary []int, pol H0Policy, denseMax int) (*Relative, error) {
	n := d0.Cols
	isB := make([]bool, n)
	for _, b := range boundary {
		if b < 0 || b >= n {
			return nil, errors.New("relative H0: boundary dof out of range")
		}
		isB[b] = true
	}
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
	for r := 0; r < d0.Rows; r++ {
		first := -1
		for k := d0.Ptr[r]; k < d0.Ptr[r+1]; k++ {
			c := d0.Idx[k]
			if isB[c] || d0.Val[k] == 0 {
				continue
			}
			if first < 0 {
				first = c
			} else if a, b := find(first), find(c); a != b {
				p[max(a, b)] = min(a, b)
			}
		}
	}
	blockOf := make([]int, n)
	for i := range blockOf {
		blockOf[i] = -1
	}
	var blocks [][]int
	for i := 0; i < n; i++ {
		if isB[i] {
			continue
		}
		r := find(i)
		if blockOf[r] < 0 {
			blockOf[r] = len(blocks)
			blocks = append(blocks, nil)
		}
		blocks[blockOf[r]] = append(blocks[blockOf[r]], i)
	}
	colPos := make([]int, n)
	for _, cols := range blocks {
		for q, c := range cols {
			colPos[c] = q
		}
	}
	rowsOf := make([][]int, len(blocks))
	for r := 0; r < d0.Rows; r++ {
		for k := d0.Ptr[r]; k < d0.Ptr[r+1]; k++ {
			if c := d0.Idx[k]; d0.Val[k] != 0 && !isB[c] {
				b := blockOf[find(c)]
				rowsOf[b] = append(rowsOf[b], r)
				break
			}
		}
	}
	rel := &Relative{}
	for bi, cols := range blocks {
		if len(cols) > denseMax {
			return nil, ErrBlockTooLarge
		}
		rows := rowsOf[bi]
		b := make([][]float64, len(rows))
		for i, r := range rows {
			b[i] = make([]float64, len(cols))
			for k := d0.Ptr[r]; k < d0.Ptr[r+1]; k++ {
				if c := d0.Idx[k]; !isB[c] && d0.Val[k] != 0 {
					b[i][colPos[c]] += d0.Val[k]
				}
			}
		}
		bound := denseBound(b, len(cols))
		cut := pol.ZeroTol(bound)
		rel.Cuts = append(rel.Cuts, cut)
		s, v := svdFull(b, len(rows), len(cols))
		r := 0
		for _, x := range s {
			if x*x > cut {
				r++
			}
		}
		for j := r; j < len(cols); j++ {
			g := make([]float64, n)
			for q, c := range cols {
				g[c] = v.At(q, j)
			}
			rel.Basis = append(rel.Basis, g)
		}
	}
	rel.Dim = len(rel.Basis)
	return rel, nil
}

func denseBound(b [][]float64, cols int) float64 {
	col := make([]float64, cols)
	rmax := 0.0
	for _, row := range b {
		r := 0.0
		for j, v := range row {
			r += math.Abs(v)
			col[j] += math.Abs(v)
		}
		rmax = math.Max(rmax, r)
	}
	cmax := 0.0
	for _, c := range col {
		cmax = math.Max(cmax, c)
	}
	return cmax * rmax
}

// BoundaryDof fixes C0 dof Dof to Value.
type BoundaryDof struct {
	Dof   int
	Value float64
}

// Extension is the harmonic extension: x_I minimising ||delta x||^2 with x_B fixed (L_II x_I = -L_IB x_B).
type Extension struct {
	X []float64
	// Route: "trivial" (no interior), "dense" (Cholesky of L_II), or "dense-pinv" (L_II singular: minimum-norm
	// solution through its eigendecomposition; the reference takes its CG route here and reports "cg").
	Route       string
	ResidualRel float64
	Converged   bool
	Energy      float64
}

// HarmonicExtension solves the extension densely. Residual scale and pivot tolerance use ||L_II||'s bound
// (delta restricted to the interior columns), as in the reference.
func HarmonicExtension(d0 *CSR, bd []BoundaryDof) Extension {
	n := d0.Cols
	isB := make([]bool, n)
	x := make([]float64, n)
	xb := make([]float64, n)
	for _, b := range bd {
		isB[b.Dof] = true
		x[b.Dof] = b.Value
		xb[b.Dof] = b.Value
	}
	var interior []int
	for i := 0; i < n; i++ {
		if !isB[i] {
			interior = append(interior, i)
		}
	}
	ni := len(interior)
	lxb := applyL(d0, xb)
	rhs := make([]float64, ni)
	for k, i := range interior {
		rhs[k] = -lxb[i]
	}
	// ||delta_I||_1 ||delta_I||_inf
	c1 := make([]float64, n)
	rinf := 0.0
	for r := 0; r < d0.Rows; r++ {
		s := 0.0
		for k := d0.Ptr[r]; k < d0.Ptr[r+1]; k++ {
			if c := d0.Idx[k]; !isB[c] {
				s += math.Abs(d0.Val[k])
				c1[c] += math.Abs(d0.Val[k])
			}
		}
		rinf = math.Max(rinf, s)
	}
	cmax := 0.0
	for _, v := range c1 {
		cmax = math.Max(cmax, v)
	}
	lb := rinf * cmax
	e := Extension{Route: "trivial"}
	if ni > 0 {
		a := make([][]float64, ni)
		f := make([]float64, n)
		for j, col := range interior {
			f[col] = 1
			y := applyL(d0, f)
			f[col] = 0
			for i, row := range interior {
				if a[i] == nil {
					a[i] = make([]float64, ni)
				}
				a[i][j] = y[row]
			}
		}
		if xi, ok := cholSolve(a, rhs, 1e-12*math.Max(lb, 1e-300)); ok {
			e.Route = "dense"
			for k, i := range interior {
				x[i] = xi[k]
			}
		} else {
			e.Route = "dense-pinv"
			// minimum-norm solution of the consistent system: every eigen-direction above working precision
			// (numpy pinv rule n*eps*lambda_max) is solved, as the reference's CG does; only the numerical kernel
			// of L_II is dropped.
			vals, vecs := symEig(a, ni)
			cut := float64(ni) * eps * math.Max(vals[ni-1], 0)
			for j := 0; j < ni; j++ {
				if vals[j] <= cut {
					continue
				}
				c := 0.0
				for i := 0; i < ni; i++ {
					c += vecs.At(i, j) * rhs[i]
				}
				c /= vals[j]
				for k, i := range interior {
					x[i] += c * vecs.At(k, j)
				}
			}
		}
	}
	lx := applyL(d0, x)
	resAbs, xin := 0.0, 0.0
	for _, i := range interior {
		resAbs += lx[i] * lx[i]
		xin += x[i] * x[i]
	}
	scale := math.Max(math.Max(norm(rhs), lb*math.Sqrt(xin)), 1e-300)
	if ni > 0 {
		e.ResidualRel = math.Sqrt(resAbs) / scale
	}
	e.Converged = !math.IsNaN(e.ResidualRel) && !math.IsInf(e.ResidualRel, 0) && e.ResidualRel <= 1e-10
	dx := d0.Mul(x)
	e.Energy = dot(dx, dx)
	e.X = x
	return e
}

// cholSolve: Cholesky of the symmetric a and solve a x = b; false if a pivot <= tol (not numerically SPD).
func cholSolve(a [][]float64, b []float64, tol float64) ([]float64, bool) {
	n := len(b)
	l := make([][]float64, n)
	for i := range l {
		l[i] = append([]float64(nil), a[i]...)
	}
	for j := 0; j < n; j++ {
		d := l[j][j]
		for k := 0; k < j; k++ {
			d -= l[j][k] * l[j][k]
		}
		if !(d > tol) {
			return nil, false
		}
		d = math.Sqrt(d)
		l[j][j] = d
		for i := j + 1; i < n; i++ {
			s := l[i][j]
			for k := 0; k < j; k++ {
				s -= l[i][k] * l[j][k]
			}
			l[i][j] = s / d
		}
	}
	y := append([]float64(nil), b...)
	for i := 0; i < n; i++ {
		for k := 0; k < i; k++ {
			y[i] -= l[i][k] * y[k]
		}
		y[i] /= l[i][i]
	}
	for i := n - 1; i >= 0; i-- {
		for k := i + 1; k < n; k++ {
			y[i] -= l[k][i] * y[k]
		}
		y[i] /= l[i][i]
	}
	return y, true
}
