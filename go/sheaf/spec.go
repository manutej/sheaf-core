package sheaf

import (
	"bytes"
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"math"
	"math/bits"
	"strconv"
	"strings"
)

// Format identifiers (wire format ids shared with the reference implementation).
const (
	SpecFormat     = "frontier-sheaf/v1"
	ReportFormat   = "frontier-sheaf/report-v1"
	BoundaryFormat = "frontier-sheaf/boundary-v1"
)

// Admission limits, applied before any allocation (same values and order as the reference parser).
const (
	// MaxTotalDim caps the total node stalk dimension and every edge stalk dimension.
	MaxTotalDim = (1 << 31) - 1
	// MaxMapBytes caps the cumulative dense restriction-map storage (8 bytes per entry).
	MaxMapBytes = 4 << 30
	// MaxAbsEntry is the largest admitted |entry| in maps, x, w, boundary values and h0_absolute_tol.
	MaxAbsEntry = 1e75
)

// SpecError is an input error. Kind is "json" (malformed JSON) or "spec" (valid JSON, invalid spec).
// Error() renders "P_INPUT_01 <kind>: <msg>", the reference CLI's stderr line.
type SpecError struct {
	Kind string
	Msg  string
}

func (e *SpecError) Error() string { return "P_INPUT_01 " + e.Kind + ": " + e.Msg }

func fieldErr(format string, a ...any) error {
	return &SpecError{Kind: "spec", Msg: fmt.Sprintf(format, a...)}
}
func jsonErr(err error) error { return &SpecError{Kind: "json", Msg: err.Error()} }

// ErrorKind returns "json" or "spec" for a *SpecError and "" otherwise.
func ErrorKind(err error) string {
	var se *SpecError
	if errors.As(err, &se) {
		return se.Kind
	}
	return ""
}

// ZeroCutKind selects the H0 zero-cut rule.
type ZeroCutKind int

const (
	// PerComponent: lambda counts as zero iff lambda <= Rel * max(1, ||L0|_c||_bound) on its own component c.
	PerComponent ZeroCutKind = iota
	// Absolute: lambda counts as zero iff lambda <= Abs on every component (caller-chosen permitted delta).
	Absolute
)

// DefaultRel is the default relative zero cut.
const DefaultRel = 1e-10

// H0Policy is the single H0 zero-cut rule used by H0, rank(delta0) in H1 and the Hodge split, and relative H0.
type H0Policy struct {
	Kind ZeroCutKind
	Rel  float64
	Abs  float64
}

// DefaultPolicy is the per-component cut with rel = 1e-10.
func DefaultPolicy() H0Policy { return H0Policy{Kind: PerComponent, Rel: DefaultRel} }

// ZeroTol is the lambda cut for one component whose L0 block has norm bound componentBound.
func (p H0Policy) ZeroTol(componentBound float64) float64 {
	if p.Kind == Absolute {
		return p.Abs
	}
	return p.Rel * math.Max(componentBound, 1)
}

// Name is the machine-readable policy name echoed in reports.
func (p H0Policy) Name() string {
	if p.Kind == Absolute {
		return "absolute"
	}
	return "per_component_eig"
}

// Describe is the human-readable rule echoed in reports.
func (p H0Policy) Describe() string {
	if p.Kind == Absolute {
		return fmt.Sprintf("absolute: lambda <= %s on every component (caller-chosen permitted delta, ADR-015)", fmtExp(p.Abs))
	}
	return fmt.Sprintf("per component c: lambda <= %s * max(1, ||L0|_c||_bound) (ADR-012/ADR-015)", fmtExp(p.Rel))
}

// Options mirror the spec "options" object and the CLI flags.
type Options struct {
	K        int
	NullReps int
	Seed     uint64
	Top      int
	// Method is accepted for compatibility ("shift-invert" or "polynomial"); this package always uses dense EVD.
	Method   string
	DenseMax int
	H0       H0Policy
}

// DefaultOptions are the reference defaults.
func DefaultOptions() Options {
	return Options{K: 6, NullReps: 20, Seed: 1, Top: 10, Method: "shift-invert", DenseMax: 1500, H0: DefaultPolicy()}
}

// Edge j of a graph sheaf: (delta x)_j = Fv x_v - Fu x_u, Fu is De x dim(u), Fv is De x dim(v), row-major.
type Edge struct {
	U, V   int
	De     int
	Fu, Fv [][]float64
	Label  string
}

// BoundaryValue fixes node Node to Value (length = its stalk dim) for harmonic extension.
type BoundaryValue struct {
	Node  int
	Value []float64
}

// Spec is a parsed "frontier-sheaf/v1" spec.
type Spec struct {
	IDs   []string
	Dims  []int
	Edges []Edge
	X     []float64 // nil if absent
	W     []float64 // nil if absent
	Opts  Options
	// Boundary is set by ParseBoundary / --boundary; nil means no extension section.
	Boundary []BoundaryValue
}

// ---- JSON value helpers (serde_json::Value semantics) ----

func decodeJSON(text []byte) (any, error) {
	dec := json.NewDecoder(bytes.NewReader(text))
	dec.UseNumber()
	var v any
	if err := dec.Decode(&v); err != nil {
		return nil, jsonErr(err)
	}
	if _, err := dec.Token(); err != io.EOF {
		return nil, jsonErr(errors.New("trailing characters"))
	}
	if err := checkNumbers(v); err != nil {
		return nil, err
	}
	return v, nil
}

// checkNumbers rejects number literals that overflow f64 (serde_json: "number out of range").
func checkNumbers(v any) error {
	switch t := v.(type) {
	case json.Number:
		if f, err := strconv.ParseFloat(string(t), 64); err != nil && math.IsInf(f, 0) {
			return jsonErr(fmt.Errorf("number out of range: %s", string(t)))
		}
	case []any:
		for _, x := range t {
			if err := checkNumbers(x); err != nil {
				return err
			}
		}
	case map[string]any:
		for _, x := range t {
			if err := checkNumbers(x); err != nil {
				return err
			}
		}
	}
	return nil
}

func get(v any, key string) (any, bool) {
	m, ok := v.(map[string]any)
	if !ok {
		return nil, false
	}
	x, ok := m[key]
	return x, ok
}

// asU64: a non-negative integer literal that fits u64 (serde_json as_u64).
func asU64(v any) (uint64, bool) {
	n, ok := v.(json.Number)
	if !ok {
		return 0, false
	}
	s := string(n)
	if s == "" || strings.ContainsAny(s, ".eE-+") {
		return 0, false
	}
	u, err := strconv.ParseUint(s, 10, 64)
	return u, err == nil
}

func asF64(v any) (float64, bool) {
	n, ok := v.(json.Number)
	if !ok {
		return 0, false
	}
	f, err := strconv.ParseFloat(string(n), 64)
	if err != nil && math.IsInf(f, 0) {
		return 0, false
	}
	return f, true
}

func asStr(v any) (string, bool) { s, ok := v.(string); return s, ok }
func asArr(v any) ([]any, bool)  { a, ok := v.([]any); return a, ok }

func clampInt(u uint64) int {
	if u > math.MaxInt {
		return math.MaxInt
	}
	return int(u)
}

func num(v any, what string) (float64, error) {
	x, ok := asF64(v)
	if !ok {
		return 0, fieldErr("%s: not a number", what)
	}
	if math.IsNaN(x) || math.IsInf(x, 0) {
		return 0, fieldErr("%s: not finite", what)
	}
	if math.Abs(x) > MaxAbsEntry {
		return 0, fieldErr("%s: |value| %s exceeds %s (L0 / energy would overflow f64)", what, fmtExp(x), fmtExp(MaxAbsEntry))
	}
	return x, nil
}

func matrix(v any, rows, cols int, what string) ([][]float64, error) {
	r, ok := asArr(v)
	if !ok {
		return nil, fieldErr("%s: must be a list of rows", what)
	}
	if len(r) != rows {
		return nil, fieldErr("%s: %d rows, expected %d", what, len(r), rows)
	}
	m := make([][]float64, len(r))
	for i, rowv := range r {
		row, ok := asArr(rowv)
		if !ok {
			return nil, fieldErr("%s: row %d not a list", what, i)
		}
		if len(row) != cols {
			return nil, fieldErr("%s: row %d has %d cols, expected %d", what, i, len(row), cols)
		}
		m[i] = make([]float64, cols)
		for j, x := range row {
			f, err := num(x, what)
			if err != nil {
				return nil, err
			}
			m[i][j] = f
		}
	}
	return m, nil
}

func vectors(v any, lens []int, what string) ([]float64, error) {
	r, ok := asArr(v)
	if !ok {
		return nil, fieldErr("%s: must be a list", what)
	}
	if len(r) != len(lens) {
		return nil, fieldErr("%s: %d entries, expected %d", what, len(r), len(lens))
	}
	out := []float64{}
	for i, rowv := range r {
		row, ok := asArr(rowv)
		if !ok {
			return nil, fieldErr("%s[%d]: not a list", what, i)
		}
		if len(row) != lens[i] {
			return nil, fieldErr("%s[%d]: length %d, expected %d", what, i, len(row), lens[i])
		}
		for _, x := range row {
			f, err := num(x, what)
			if err != nil {
				return nil, err
			}
			out = append(out, f)
		}
	}
	return out, nil
}

func identity(rows, cols int) [][]float64 {
	m := make([][]float64, rows)
	for i := range m {
		m[i] = make([]float64, cols)
		if i < cols {
			m[i][i] = 1
		}
	}
	return m
}

// Parse parses a "frontier-sheaf/v1" spec. Checks run in the reference order with the reference messages.
func Parse(text []byte) (*Spec, error) {
	v, err := decodeJSON(text)
	if err != nil {
		return nil, err
	}
	if f, _ := get(v, "format"); f != SpecFormat {
		return nil, fieldErr("format must be \"%s\"", SpecFormat)
	}
	nv, _ := get(v, "nodes")
	nodes, ok := asArr(nv)
	if !ok {
		return nil, fieldErr("nodes: missing list")
	}
	ids := make([]string, 0, len(nodes))
	dims := make([]int, 0, len(nodes))
	for i, n := range nodes {
		idv, _ := get(n, "id")
		id, ok := asStr(idv)
		if !ok {
			id = strconv.Itoa(i)
		}
		dv, _ := get(n, "dim")
		d, ok := asU64(dv)
		if !ok {
			return nil, fieldErr("nodes[%d].dim: missing non-negative integer", i)
		}
		ids = append(ids, id)
		dims = append(dims, clampInt(d))
	}
	total := uint64(0)
	for _, d := range dims {
		var carry uint64
		total, carry = bits.Add64(total, uint64(d), 0)
		if carry != 0 || total > MaxTotalDim {
			return nil, fieldErr("nodes: total stalk dimension exceeds %d", MaxTotalDim)
		}
	}
	index := make(map[string]int, len(ids))
	for i, s := range ids {
		index[s] = i
	}
	if len(index) != len(ids) {
		return nil, fieldErr("nodes: duplicate id")
	}
	endpoint := func(x any, present bool, j int, w string) (int, error) {
		if present {
			switch t := x.(type) {
			case json.Number:
				if u, ok := asU64(t); ok && u < uint64(len(ids)) {
					return int(u), nil
				}
				return 0, fieldErr("edges[%d].%s: index out of range", j, w)
			case string:
				if i, ok := index[t]; ok {
					return i, nil
				}
				return 0, fieldErr("edges[%d].%s: unknown node id %s", j, w, t)
			}
		}
		return 0, fieldErr("edges[%d].%s: missing", j, w)
	}
	ev, _ := get(v, "edges")
	es, ok := asArr(ev)
	if !ok {
		return nil, fieldErr("edges: missing list")
	}
	edges := make([]Edge, 0, len(es))
	mapTotal := uint64(0)
	for j, e := range es {
		uv, up := get(e, "u")
		u, err := endpoint(uv, up, j, "u")
		if err != nil {
			return nil, err
		}
		vv, vp := get(e, "v")
		w, err := endpoint(vv, vp, j, "v")
		if err != nil {
			return nil, err
		}
		if u == w {
			return nil, fieldErr("edges[%d]: self-loop", j)
		}
		fuv, fuPresent := get(e, "Fu")
		fvv, fvPresent := get(e, "Fv")
		// de = the "de" field, else rows of Fu, else rows of Fv, else dim(u)
		de := dims[u]
		dev, _ := get(e, "de")
		if d, ok := asU64(dev); ok {
			de = clampInt(d)
		} else if a, ok := asArr(fuv); ok {
			de = len(a)
		} else if a, ok := asArr(fvv); ok {
			de = len(a)
		}
		// two dense maps (de x dim u, de x dim v), measured before they are allocated
		hi, lo := bits.Mul64(uint64(de), uint64(dims[u])+uint64(dims[w]))
		ok8 := hi == 0
		var b uint64
		if ok8 {
			var hi2 uint64
			hi2, b = bits.Mul64(lo, 8)
			ok8 = hi2 == 0
		}
		if ok8 {
			var carry uint64
			mapTotal, carry = bits.Add64(mapTotal, b, 0)
			if carry != 0 {
				mapTotal = math.MaxUint64
			}
		} else {
			mapTotal = math.MaxUint64
		}
		if de > MaxTotalDim || mapTotal > MaxMapBytes {
			return nil, fieldErr("edges[%d]: restriction maps exceed %d bytes (dense storage)", j, MaxMapBytes)
		}
		mk := func(f any, present bool, d int, name string) ([][]float64, error) {
			if present {
				return matrix(f, de, d, fmt.Sprintf("edges[%d].%s", j, name))
			}
			if d == de {
				return identity(de, d), nil
			}
			return nil, fieldErr("edges[%d].%s: omitted but dim %d != de %d", j, name, d, de)
		}
		fu, err := mk(fuv, fuPresent, dims[u], "Fu")
		if err != nil {
			return nil, err
		}
		fv, err := mk(fvv, fvPresent, dims[w], "Fv")
		if err != nil {
			return nil, err
		}
		lv, _ := get(e, "label")
		label, ok := asStr(lv)
		if !ok {
			label = ids[u] + "->" + ids[w]
		}
		edges = append(edges, Edge{U: u, V: w, De: de, Fu: fu, Fv: fv, Label: label})
	}
	des := make([]int, len(edges))
	for j, e := range edges {
		des[j] = e.De
	}
	s := &Spec{IDs: ids, Dims: dims, Edges: edges, Opts: DefaultOptions()}
	if xv, ok := get(v, "x"); ok && xv != nil {
		if s.X, err = vectors(xv, dims, "x"); err != nil {
			return nil, err
		}
	}
	if wv, ok := get(v, "w"); ok && wv != nil {
		if s.W, err = vectors(wv, des, "w"); err != nil {
			return nil, err
		}
	}
	if s.X == nil && s.W == nil {
		return nil, fieldErr("need x (node data) or w (edge cochain)")
	}
	if o, ok := get(v, "options"); ok {
		if k, ok := get(o, "k"); ok {
			if u, ok := asU64(k); ok {
				s.Opts.K = max(clampInt(u), 1)
			}
		}
		if k, ok := get(o, "null_reps"); ok {
			if u, ok := asU64(k); ok {
				s.Opts.NullReps = clampInt(u)
			}
		}
		if k, ok := get(o, "seed"); ok {
			if u, ok := asU64(k); ok {
				s.Opts.Seed = u
			}
		}
		if k, ok := get(o, "top"); ok {
			if u, ok := asU64(k); ok {
				s.Opts.Top = clampInt(u)
			}
		}
		if k, ok := get(o, "dense_max"); ok {
			if u, ok := asU64(k); ok {
				s.Opts.DenseMax = clampInt(u)
			}
		}
		if t, ok := get(o, "h0_absolute_tol"); ok {
			f, err := num(t, "options.h0_absolute_tol")
			if err != nil {
				return nil, err
			}
			if f < 0 {
				return nil, fieldErr("options.h0_absolute_tol: must be >= 0")
			}
			s.Opts.H0 = H0Policy{Kind: Absolute, Abs: f}
		}
		if m, ok := get(o, "method"); ok {
			if ms, ok := asStr(m); ok {
				switch ms {
				case "shift-invert", "polynomial":
					s.Opts.Method = ms
				default:
					return nil, fieldErr("options.method: unknown %s", ms)
				}
			}
		}
	}
	return s, nil
}

// ParseBoundary parses a boundary file {"format": "frontier-sheaf/boundary-v1" (optional),
// "boundary": [{"node": "a" | 0, "value": [...]}, ...]} against spec s.
func ParseBoundary(text []byte, s *Spec) ([]BoundaryValue, error) {
	v, err := decodeJSON(text)
	if err != nil {
		return nil, err
	}
	if f, ok := get(v, "format"); ok && f != BoundaryFormat {
		return nil, fieldErr("boundary format must be \"%s\"", BoundaryFormat)
	}
	lv, _ := get(v, "boundary")
	list, ok := asArr(lv)
	if !ok {
		return nil, fieldErr("boundary: missing list")
	}
	seen := make([]bool, len(s.Dims))
	out := []BoundaryValue{}
	for j, b := range list {
		node := -1
		if nv, ok := get(b, "node"); ok {
			switch t := nv.(type) {
			case json.Number:
				if u, ok := asU64(t); ok && u < uint64(len(s.Dims)) {
					node = int(u)
				}
			case string:
				for i, id := range s.IDs {
					if id == t {
						node = i
						break
					}
				}
			}
		}
		if node < 0 {
			return nil, fieldErr("boundary[%d].node: missing or unknown", j)
		}
		if seen[node] {
			return nil, fieldErr("boundary[%d].node: duplicate", j)
		}
		seen[node] = true
		vv, _ := get(b, "value")
		val, ok := asArr(vv)
		if !ok {
			return nil, fieldErr("boundary[%d].value: not a list", j)
		}
		if len(val) != s.Dims[node] {
			return nil, fieldErr("boundary[%d].value: length %d, expected %d", j, len(val), s.Dims[node])
		}
		vals := make([]float64, len(val))
		for i, x := range val {
			if vals[i], err = num(x, fmt.Sprintf("boundary[%d].value", j)); err != nil {
				return nil, err
			}
		}
		out = append(out, BoundaryValue{Node: node, Value: vals})
	}
	return out, nil
}

// fmtExp renders x like Rust's `{:e}` (shortest round-trip mantissa, exponent without '+' or zero padding).
func fmtExp(x float64) string {
	s := strconv.FormatFloat(x, 'e', -1, 64)
	i := strings.IndexByte(s, 'e')
	if i < 0 {
		return s
	}
	mant, exp := s[:i], s[i+1:]
	neg := strings.HasPrefix(exp, "-")
	exp = strings.TrimLeft(exp, "+-")
	exp = strings.TrimLeft(exp, "0")
	if exp == "" {
		exp = "0"
	}
	if neg {
		exp = "-" + exp
	}
	return mant + "e" + exp
}
