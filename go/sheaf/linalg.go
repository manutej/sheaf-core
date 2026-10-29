package sheaf

import (
	"math"
	"sort"

	"gonum.org/v1/gonum/mat"
)

// CSR is a row-major sparse matrix. Rows keep their entries sorted by column and store nonzeros only, so
// products sum in the same order as the reference (bitwise-equal delta0 x).
type CSR struct {
	Rows, Cols int
	Ptr        []int
	Idx        []int
	Val        []float64
}

// Mul returns A x.
func (a *CSR) Mul(x []float64) []float64 {
	y := make([]float64, a.Rows)
	for i := 0; i < a.Rows; i++ {
		s := 0.0
		for p := a.Ptr[i]; p < a.Ptr[i+1]; p++ {
			s += a.Val[p] * x[a.Idx[p]]
		}
		y[i] = s
	}
	return y
}

// MulT returns A^T y.
func (a *CSR) MulT(y []float64) []float64 {
	out := make([]float64, a.Cols)
	for i := 0; i < a.Rows; i++ {
		for p := a.Ptr[i]; p < a.Ptr[i+1]; p++ {
			out[a.Idx[p]] += a.Val[p] * y[i]
		}
	}
	return out
}

// Dense returns A as a row-major dense matrix (rows x cols).
func (a *CSR) Dense() [][]float64 {
	m := make([][]float64, a.Rows)
	for i := range m {
		m[i] = make([]float64, a.Cols)
		for p := a.Ptr[i]; p < a.Ptr[i+1]; p++ {
			m[i][a.Idx[p]] += a.Val[p]
		}
	}
	return m
}

// NormBound returns ||A||_1 * ||A||_inf, an upper bound on ||A^T A||_2.
func (a *CSR) NormBound() float64 {
	col := make([]float64, a.Cols)
	rmax := 0.0
	for i := 0; i < a.Rows; i++ {
		r := 0.0
		for p := a.Ptr[i]; p < a.Ptr[i+1]; p++ {
			r += math.Abs(a.Val[p])
			col[a.Idx[p]] += math.Abs(a.Val[p])
		}
		rmax = math.Max(rmax, r)
	}
	cmax := 0.0
	for _, c := range col {
		cmax = math.Max(cmax, c)
	}
	return cmax * rmax
}

// Offsets returns prefix sums [0, d0, d0+d1, ...].
func Offsets(d []int) []int {
	o := make([]int, len(d)+1)
	for i, x := range d {
		o[i+1] = o[i] + x
	}
	return o
}

// Coboundary assembles delta0 (dim C1 x dim C0) with (delta x)_j = Fv x_v - Fu x_u, plus node and edge offsets.
// Edge order, multi-edges and orientation are kept; exact zeros are not stored.
func Coboundary(dims []int, edges []Edge) (*CSR, []int, []int) {
	on := Offsets(dims)
	des := make([]int, len(edges))
	for j, e := range edges {
		des[j] = e.De
	}
	oe := Offsets(des)
	a := &CSR{Rows: oe[len(edges)], Cols: on[len(dims)], Ptr: make([]int, oe[len(edges)]+1)}
	type ent struct {
		c int
		v float64
	}
	for j, e := range edges {
		for r := 0; r < e.De; r++ {
			row := make([]ent, 0, dims[e.U]+dims[e.V])
			for c, x := range e.Fu[r] {
				if v := -x; v != 0 {
					row = append(row, ent{on[e.U] + c, v})
				}
			}
			for c, x := range e.Fv[r] {
				if x != 0 {
					row = append(row, ent{on[e.V] + c, x})
				}
			}
			sort.SliceStable(row, func(p, q int) bool { return row[p].c < row[q].c })
			for _, t := range row {
				a.Idx = append(a.Idx, t.c)
				a.Val = append(a.Val, t.v)
			}
			a.Ptr[oe[j]+r+1] = len(a.Idx)
		}
	}
	return a, on, oe
}

func dot(a, b []float64) float64 {
	s := 0.0
	for i := range a {
		s += a[i] * b[i]
	}
	return s
}

func norm(a []float64) float64 { return math.Sqrt(dot(a, a)) }

// svdValues returns the singular values (descending) and V (cols x cols, columns = right singular vectors) of
// a rows x cols dense matrix. A matrix with no rows or no cols has no singular values and V = I.
func svdFull(b [][]float64, rows, cols int) ([]float64, *mat.Dense) {
	if rows == 0 || cols == 0 {
		v := mat.NewDense(max(cols, 1), max(cols, 1), nil)
		for i := 0; i < cols; i++ {
			v.Set(i, i, 1)
		}
		return nil, v
	}
	m := mat.NewDense(rows, cols, nil)
	for i := range b {
		for j, x := range b[i] {
			m.Set(i, j, x)
		}
	}
	var s mat.SVD
	if !s.Factorize(m, mat.SVDFullV) {
		panic("sheaf: SVD did not converge")
	}
	var v mat.Dense
	s.VTo(&v)
	return s.Values(nil), &v
}

// numericRank is the numpy/LAPACK default rank: #{s > max(m,n) * eps * s_max}.
func numericRank(s []float64, rows, cols int) int {
	smax := 0.0
	for _, x := range s {
		smax = math.Max(smax, x)
	}
	tol := float64(max(rows, cols)) * eps * smax
	r := 0
	for _, x := range s {
		if x > tol {
			r++
		}
	}
	return r
}

const eps = 2.220446049250313e-16

// symEig returns ascending eigenvalues and eigenvectors (columns) of the symmetric n x n matrix a (row-major).
func symEig(a [][]float64, n int) ([]float64, *mat.Dense) {
	s := mat.NewSymDense(n, nil)
	for i := 0; i < n; i++ {
		for j := 0; j <= i; j++ {
			s.SetSym(i, j, a[i][j])
		}
	}
	var e mat.EigenSym
	if !e.Factorize(s, true) {
		panic("sheaf: symmetric eigendecomposition did not converge")
	}
	var v mat.Dense
	e.VectorsTo(&v)
	return e.Values(nil), &v
}

// svdThin returns s (descending, k = min(rows, cols)), U (rows x k) and V (cols x k) as row-major slices.
func svdThin(b [][]float64, rows, cols int) ([]float64, [][]float64, [][]float64) {
	if rows == 0 || cols == 0 {
		return nil, make([][]float64, rows), make([][]float64, cols)
	}
	m := mat.NewDense(rows, cols, nil)
	for i := range b {
		for j, x := range b[i] {
			m.Set(i, j, x)
		}
	}
	var s mat.SVD
	if !s.Factorize(m, mat.SVDThin) {
		panic("sheaf: SVD did not converge")
	}
	var u, v mat.Dense
	s.UTo(&u)
	s.VTo(&v)
	k := min(rows, cols)
	us := make([][]float64, rows)
	for i := range us {
		us[i] = make([]float64, k)
		for j := 0; j < k; j++ {
			us[i][j] = u.At(i, j)
		}
	}
	vs := make([][]float64, cols)
	for i := range vs {
		vs[i] = make([]float64, k)
		for j := 0; j < k; j++ {
			vs[i][j] = v.At(i, j)
		}
	}
	return s.Values(nil), us, vs
}
