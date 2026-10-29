// Command gen writes the generated parity fixtures into testdata/specs (deterministic: fixed seeds, PCG).
// Run from go/: `go run ./internal/gen`. Every fixture is synthetic: seeded random sheaves or hand-built
// adversarial cases. Shared fixtures (python/parity/fixtures) and the labelled toys are copied, not generated.
package main

import (
	"encoding/json"
	"fmt"
	"math"
	"math/rand/v2"
	"os"
	"path/filepath"
	"sort"
)

type obj = map[string]any

var outDir = "testdata/specs"

func write(name string, v any) {
	b, err := json.MarshalIndent(v, "", " ")
	if err != nil {
		panic(err)
	}
	if err := os.WriteFile(filepath.Join(outDir, name), append(b, '\n'), 0o644); err != nil {
		panic(err)
	}
}

// writeRaw writes bytes verbatim (malformed-JSON cases).
func writeRaw(name, s string) {
	if err := os.WriteFile(filepath.Join(outDir, name), []byte(s), 0o644); err != nil {
		panic(err)
	}
}

func node(id string, dim int) obj { return obj{"id": id, "dim": dim} }

func mat(r *rand.Rand, rows, cols int) [][]float64 {
	m := make([][]float64, rows)
	for i := range m {
		m[i] = make([]float64, cols)
		for j := range m[i] {
			m[i][j] = round(r.NormFloat64())
		}
	}
	return m
}

func vec(r *rand.Rand, n int) []float64 {
	v := make([]float64, n)
	for i := range v {
		v[i] = round(r.NormFloat64())
	}
	return v
}

// round keeps fixtures short and exactly representable in JSON (12 significant digits).
func round(x float64) float64 {
	var y float64
	fmt.Sscanf(fmt.Sprintf("%.11e", x), "%g", &y)
	return y
}

func spec(nodes []obj, edges []obj, extra obj) obj {
	s := obj{"format": "frontier-sheaf/v1", "nodes": nodes, "edges": edges}
	for k, v := range extra {
		s[k] = v
	}
	return s
}

func eye(n int) [][]float64 {
	m := make([][]float64, n)
	for i := range m {
		m[i] = make([]float64, n)
		m[i][i] = 1
	}
	return m
}

func rot(t float64) [][]float64 {
	return [][]float64{{math.Cos(t), -math.Sin(t)}, {math.Sin(t), math.Cos(t)}}
}

func randomSheaves() {
	for seed := uint64(1); seed <= 24; seed++ {
		r := rand.New(rand.NewPCG(20261029, seed))
		n := 2 + r.IntN(9)
		dims := make([]int, n)
		nodes := make([]obj, n)
		for i := range dims {
			dims[i] = r.IntN(4) // 0..3, mixed stalk dims including 0
			if seed%3 != 0 && dims[i] == 0 {
				dims[i] = 1 + r.IntN(3)
			}
			nodes[i] = node(fmt.Sprintf("n%d", i), dims[i])
		}
		ne := r.IntN(2*n + 1)
		var edges []obj
		for j := 0; j < ne; j++ {
			u := r.IntN(n)
			v := r.IntN(n - 1)
			if v >= u {
				v++
			}
			de := r.IntN(4)
			if seed%4 != 0 && de == 0 {
				de = 1 + r.IntN(2)
			}
			e := obj{"u": u, "v": v, "Fu": mat(r, de, dims[u]), "Fv": mat(r, de, dims[v])}
			if de == 0 {
				e["de"] = 0
			}
			if r.IntN(3) == 0 {
				e["u"], e["v"] = fmt.Sprintf("n%d", u), fmt.Sprintf("n%d", v)
			}
			if r.IntN(2) == 0 {
				e["label"] = fmt.Sprintf("e%d", j)
			}
			if dims[u] == dims[v] && de == dims[u] && r.IntN(4) == 0 { // identity default (omitted maps)
				delete(e, "Fu")
				delete(e, "Fv")
			}
			edges = append(edges, e)
		}
		extra := obj{}
		mode := seed % 3
		if mode != 1 {
			xs := make([][]float64, n)
			for i := range xs {
				xs[i] = vec(r, dims[i])
			}
			extra["x"] = xs
		}
		if mode != 0 {
			ws := make([][]float64, len(edges))
			for j, e := range edges {
				de := 0
				if fu, ok := e["Fu"].([][]float64); ok {
					de = len(fu)
				} else {
					de = dims[e2i(e["u"])]
				}
				ws[j] = vec(r, de)
			}
			extra["w"] = ws
		}
		extra["_comment"] = fmt.Sprintf("RANDOM synthetic sheaf (go/internal/gen, PCG(20261029, %d)); mixed stalk dims", seed)
		if edges == nil {
			edges = []obj{}
		}
		write(fmt.Sprintf("rand_%02d.json", seed), spec(nodes, edges, extra))
	}
}

func e2i(v any) int {
	switch t := v.(type) {
	case int:
		return t
	case string:
		var i int
		fmt.Sscanf(t, "n%d", &i)
		return i
	}
	panic("endpoint")
}

// cycleNearCut: a 1-dim cycle whose holonomy is 1+eps, so the smallest L0 eigenvalue is ~eps^2/len.
func cycleNearCut(name string, eps float64, comment string) {
	var nodes []obj
	var edges []obj
	m := 4
	for i := 0; i < m; i++ {
		nodes = append(nodes, node(fmt.Sprintf("c%d", i), 1))
	}
	for i := 0; i < m; i++ {
		fv := 1.0
		if i == m-1 {
			fv = 1 + eps
		}
		edges = append(edges, obj{"u": i, "v": (i + 1) % m, "Fu": [][]float64{{1}}, "Fv": [][]float64{{fv}}})
	}
	write(name, spec(nodes, edges, obj{"x": [][]float64{{1}, {1}, {1}, {1}}, "_comment": comment}))
}

func adversarial() {
	c := "ADVERSARIAL synthetic case (go/internal/gen): "
	// two-scale components: global cut would give H0 = 3, per-component cut gives 2
	write("adv_two_scale.json", spec([]obj{node("a", 1), node("b", 1), node("c", 1), node("d", 1)},
		[]obj{{"u": "a", "v": "b", "Fu": [][]float64{{1e5}}, "Fv": [][]float64{{1e5}}}, {"u": "c", "v": "d", "Fu": [][]float64{{0.1}}, "Fv": [][]float64{{0.1}}}},
		obj{"x": [][]float64{{1}, {2}, {3}, {4}}, "_comment": c + "two components at scales 1e5 and 1e-1 (per-component cut)"}))
	write("adv_two_scale_generic.json", spec([]obj{node("a", 2), node("b", 2), node("c", 2), node("d", 2)},
		[]obj{{"u": "a", "v": "b", "Fu": [][]float64{{1e5, 0}}, "Fv": [][]float64{{2e5, 0}}}, {"u": "c", "v": "d", "Fu": [][]float64{{0.1, 0.05}}, "Fv": [][]float64{{0, 0.1}}}},
		obj{"x": [][]float64{{1, 0}, {2, 1}, {3, 0}, {4, 1}}, "_comment": c + "two scales, non-constant maps (no known kernel)"}))
	// near-cut eigenvalues: smallest eigenvalue placed around the per-component cut (oracle flags ambiguity)
	for i, eps := range []float64{3e-6, 1e-5, 2.5e-5, 6e-5, 2e-4} {
		cycleNearCut(fmt.Sprintf("adv_near_cut_%d.json", i), eps, fmt.Sprintf("%s1-dim 4-cycle with holonomy 1+%g: smallest eigenvalue near the zero cut", c, eps))
	}
	// multi-edges in both orientations with different edge stalk dims
	r := rand.New(rand.NewPCG(7, 7))
	write("adv_multi_edges.json", spec([]obj{node("p", 2), node("q", 2), node("s", 1)},
		[]obj{{"u": "p", "v": "q"}, {"u": "q", "v": "p"}, {"u": "p", "v": "q", "Fu": mat(r, 1, 2), "Fv": mat(r, 1, 2)}, {"u": "q", "v": "p", "Fu": mat(r, 3, 2), "Fv": mat(r, 3, 2)},
			{"u": "s", "v": "p", "Fu": [][]float64{{1}}, "Fv": [][]float64{{1, 0}}}, {"u": "p", "v": "s", "Fu": [][]float64{{0, 1}}, "Fv": [][]float64{{1}}}},
		obj{"x": [][]float64{vec(r, 2), vec(r, 2), vec(r, 1)}, "_comment": c + "parallel edges, both orientations, de 1..3"}))
	// zero maps: an edge whose maps are all zero still connects the component but constrains nothing
	write("adv_zero_maps.json", spec([]obj{node("a", 2), node("b", 2), node("c", 2)},
		[]obj{{"u": "a", "v": "b", "Fu": [][]float64{{0, 0}, {0, 0}}, "Fv": [][]float64{{0, 0}, {0, 0}}}, {"u": "b", "v": "c", "Fu": [][]float64{{0, 0}}, "Fv": [][]float64{{1, 2}}}, {"u": "a", "v": "c"}},
		obj{"x": [][]float64{{1, 2}, {3, 4}, {5, 6}}, "w": [][]float64{{1, -1}, {2}, {0.5, 0.25}}, "_comment": c + "all-zero and one-sided zero restriction maps"}))
	// isolated vertices, a dim-0 vertex, and a de = 0 edge
	write("adv_isolated.json", spec([]obj{node("a", 2), node("b", 2), node("lone", 3), node("empty", 0), node("c", 1)},
		[]obj{{"u": "a", "v": "b", "Fu": rot(0.3), "Fv": eye(2)}, {"u": "b", "v": "c", "Fu": [][]float64{{1, 1}}, "Fv": [][]float64{{2}}}, {"u": "empty", "v": "c", "de": 0, "Fu": [][]float64{}, "Fv": [][]float64{}}},
		obj{"x": [][]float64{{1, 0}, {0, 1}, {7, 8, 9}, {}, {3}}, "_comment": c + "isolated vertex, dim-0 vertex, de = 0 edge"}))
	write("adv_no_edges.json", spec([]obj{node("a", 2), node("b", 1)}, []obj{}, obj{"x": [][]float64{{1, 2}, {3}}, "_comment": c + "nodes but no edges"}))
	write("adv_empty_graph.json", spec([]obj{}, []obj{}, obj{"x": [][]float64{}, "_comment": c + "empty graph"}))
	// dense_max below the component size: only the k+1 smallest pairs are reported (widened while all are zero)
	var nodes []obj
	var edges []obj
	for i := 0; i < 8; i++ {
		nodes = append(nodes, node(fmt.Sprintf("v%d", i), 2))
	}
	for i := 0; i < 8; i++ {
		edges = append(edges, obj{"u": i, "v": (i + 1) % 8})
		if i%2 == 0 {
			edges = append(edges, obj{"u": i, "v": (i + 3) % 8})
		}
	}
	xs := make([][]float64, 8)
	for i := range xs {
		xs[i] = vec(r, 2)
	}
	write("adv_dense_max_small.json", spec(nodes, edges, obj{"x": xs, "options": obj{"dense_max": 4, "k": 2}, "_comment": c + "known-constant component above dense_max (k+1 truncation)"}))
	var zedges []obj
	for i := 0; i < 7; i++ {
		zedges = append(zedges, obj{"u": i, "v": i + 1, "Fu": [][]float64{{0, 0}}, "Fv": [][]float64{{0, 0}}})
	}
	write("adv_widening.json", spec(nodes, zedges, obj{"x": xs, "options": obj{"dense_max": 4, "k": 1}, "_comment": c + "all-zero maps above dense_max: k+1 pairs all zero, so the reference widens k x4 (twice); H0 complete"}))
	// 60 dofs of all-zero maps with dense_max below that and k = 1: widening stops at k = 16 with all 17 pairs zero,
	// so H0 is only a lower bound (lower_bound_only, ambiguous)
	var wn []obj
	var we []obj
	for i := 0; i < 30; i++ {
		wn = append(wn, node(fmt.Sprintf("z%d", i), 2))
		if i > 0 {
			we = append(we, obj{"u": i - 1, "v": i, "Fu": [][]float64{{0, 0}}, "Fv": [][]float64{{0, 0}}})
		}
	}
	wx := make([][]float64, 30)
	for i := range wx {
		wx[i] = vec(r, 2)
	}
	write("adv_lower_bound.json", spec(wn, we, obj{"x": wx, "options": obj{"dense_max": 8, "k": 1}, "_comment": c + "60 zero-map dofs above dense_max: widening exhausted, H0 lower bound only"}))
	// absolute tolerance and option handling
	write("adv_absolute_tol.json", spec([]obj{node("a", 1), node("b", 1), node("c", 1)},
		[]obj{{"u": 0, "v": 1, "Fu": [][]float64{{1}}, "Fv": [][]float64{{1e-4}}}, {"u": 1, "v": 2, "Fu": [][]float64{{1e-4}}, "Fv": [][]float64{{1}}}, {"u": 0, "v": 2, "Fu": [][]float64{{1}}, "Fv": [][]float64{{1.0000001}}}},
		obj{"x": [][]float64{{1}, {1e4}, {1}}, "options": obj{"h0_absolute_tol": 1e-6, "top": 2, "k": 1, "method": "polynomial"}, "_comment": c + "absolute zero cut, top/k/method options"}))
	// harmonic content on a w-only cochain with two components (one unpinned by the boundary)
	write("adv_harmonic_two_comp.json", spec([]obj{node("a", 1), node("b", 1), node("c", 1), node("d", 2), node("e", 2)},
		[]obj{{"u": "a", "v": "b"}, {"u": "b", "v": "c"}, {"u": "c", "v": "a"}, {"u": "d", "v": "e"}, {"u": "e", "v": "d", "Fu": rot(0.7), "Fv": eye(2)}},
		obj{"w": [][]float64{{1}, {2}, {-0.5}, {0.3, -0.2}, {1, 1}}, "_comment": c + "w-only, harmonic content, two components"}))
	write("adv_harmonic_two_comp.boundary.json", obj{"boundary": []obj{{"node": "a", "value": []float64{1}}, {"node": 2, "value": []float64{4}}}})
	write("adv_boundary_all.json", spec([]obj{node("a", 1), node("b", 2)}, []obj{{"u": "a", "v": "b", "Fu": [][]float64{{1}, {2}}, "Fv": eye(2)}},
		obj{"x": [][]float64{{1}, {1, 2}}, "_comment": c + "boundary covers every node (trivial extension)"}))
	write("adv_boundary_all.boundary.json", obj{"format": "frontier-sheaf/boundary-v1", "boundary": []obj{{"node": "a", "value": []float64{1}}, {"node": "b", "value": []float64{2, 2}}}})
	write("adv_boundary_path.json", spec([]obj{node("a", 1), node("b", 1), node("c", 1), node("d", 1), node("e", 1)},
		[]obj{{"u": "a", "v": "b"}, {"u": "b", "v": "c"}, {"u": "c", "v": "d"}, {"u": "d", "v": "e"}}, obj{"x": [][]float64{{0}, {0}, {0}, {0}, {0}}, "_comment": c + "path, unique linear extension"}))
	write("adv_boundary_path.boundary.json", obj{"boundary": []obj{{"node": "a", "value": []float64{1}}, {"node": 4, "value": []float64{5}}}})
	// a component large enough that the reference uses its iterative (Lanczos) route
	var bn []obj
	var be []obj
	m := 70
	for i := 0; i < m; i++ {
		bn = append(bn, node(fmt.Sprintf("L%d", i), 2))
	}
	for i := 0; i < m; i++ {
		be = append(be, obj{"u": i, "v": (i + 1) % m, "Fu": mat(r, 2, 2), "Fv": mat(r, 2, 2)})
		if i%5 == 0 {
			be = append(be, obj{"u": i, "v": (i + 17) % m, "Fu": mat(r, 1, 2), "Fv": mat(r, 1, 2)})
		}
	}
	bx := make([][]float64, m)
	for i := range bx {
		bx[i] = vec(r, 2)
	}
	write("adv_lanczos_route.json", spec(bn, be, obj{"x": bx, "options": obj{"dense_max": 16, "k": 4}, "_comment": c + "140-dof component above dense_max: the reference runs shift-invert Lanczos"}))

	// ---- error cases: must fail like the reference (same exit code and error kind; same message for spec errors)
	ok := func() obj {
		return spec([]obj{node("a", 2), node("b", 2)}, []obj{{"u": "a", "v": "b"}}, obj{"x": [][]float64{{1, 2}, {3, 4}}})
	}
	errCase := func(name string, mutate func(s obj)) {
		s := ok()
		mutate(s)
		s["_comment"] = c + "invalid input, must be rejected"
		write("err_"+name+".json", s)
	}
	errCase("mismatched_map_rows", func(s obj) {
		s["edges"] = []obj{{"u": "a", "v": "b", "Fu": [][]float64{{1, 0}}, "Fv": [][]float64{{1, 0}, {0, 1}}}}
	})
	errCase("mismatched_map_cols", func(s obj) {
		s["edges"] = []obj{{"u": "a", "v": "b", "Fu": [][]float64{{1, 0, 0}}, "Fv": [][]float64{{1, 0}}}}
	})
	errCase("omitted_map_dim", func(s obj) {
		s["nodes"] = []obj{node("a", 2), node("b", 3)}
		s["x"] = [][]float64{{1, 2}, {3, 4, 5}}
	})
	errCase("self_loop", func(s obj) { s["edges"] = []obj{{"u": "a", "v": "a"}} })
	errCase("unknown_node", func(s obj) { s["edges"] = []obj{{"u": "a", "v": "zz"}} })
	errCase("index_out_of_range", func(s obj) { s["edges"] = []obj{{"u": 0, "v": 2}} })
	errCase("missing_endpoint", func(s obj) { s["edges"] = []obj{{"u": 0}} })
	errCase("duplicate_id", func(s obj) { s["nodes"] = []obj{node("a", 2), node("a", 2)} })
	errCase("bad_format", func(s obj) { s["format"] = "frontier-sheaf/v2" })
	errCase("no_cochain", func(s obj) { delete(s, "x") })
	errCase("x_wrong_length", func(s obj) { s["x"] = [][]float64{{1, 2}, {3}} })
	errCase("w_wrong_count", func(s obj) { s["w"] = [][]float64{{1, 2}, {3, 4}} })
	errCase("dim_not_integer", func(s obj) { s["nodes"] = []obj{node("a", 2), {"id": "b", "dim": 2.5}} })
	errCase("dim_negative", func(s obj) { s["nodes"] = []obj{node("a", 2), {"id": "b", "dim": -1}} })
	errCase("huge_entry", func(s obj) { s["x"] = [][]float64{{1e76, 2}, {3, 4}} })
	errCase("huge_dims", func(s obj) { s["nodes"] = []obj{node("a", 200000), node("b", 200000)} })
	errCase("negative_abs_tol", func(s obj) { s["options"] = obj{"h0_absolute_tol": -1} })
	errCase("unknown_method", func(s obj) { s["options"] = obj{"method": "qr"} })
	errCase("map_not_list", func(s obj) { s["edges"] = []obj{{"u": "a", "v": "b", "Fu": "I", "Fv": [][]float64{{1, 0}, {0, 1}}}} })
	errCase("entry_not_number", func(s obj) { s["x"] = []any{[]any{1, "2"}, []any{3, 4}} })
	writeRaw("err_bad_json.json", `{"format": "frontier-sheaf/v1", "nodes": [`)
	writeRaw("err_trailing.json", `{"format":"frontier-sheaf/v1","nodes":[{"id":"a","dim":1}],"edges":[],"x":[[1]]} x`)
	writeRaw("err_number_overflow.json", `{"format":"frontier-sheaf/v1","nodes":[{"id":"a","dim":1},{"id":"b","dim":1}],"edges":[{"u":"a","v":"b"}],"x":[[1e400],[1]]}`)
	// boundary errors (valid spec, invalid boundary file)
	for _, b := range []struct {
		name string
		v    obj
	}{
		{"wrong_length", obj{"boundary": []obj{{"node": "a", "value": []float64{1}}}}},
		{"unknown_node", obj{"boundary": []obj{{"node": "zz", "value": []float64{1, 2}}}}},
		{"duplicate", obj{"boundary": []obj{{"node": "a", "value": []float64{1, 2}}, {"node": 0, "value": []float64{1, 2}}}}},
		{"bad_format", obj{"format": "x", "boundary": []obj{}}},
	} {
		s := ok()
		s["_comment"] = c + "valid spec with an invalid boundary file, must be rejected"
		write("errb_"+b.name+".json", s)
		write("errb_"+b.name+".boundary.json", b.v)
	}
	// flag errors (valid spec, bad flag): recorded in the manifest
	s := ok()
	s["_comment"] = c + "valid spec; the manifest runs it with bad flags"
	write("errf_flags.json", s)
}

func main() {
	if len(os.Args) > 1 {
		outDir = os.Args[1]
	}
	randomSheaves()
	adversarial()
	ents, _ := os.ReadDir(outDir)
	var names []string
	for _, e := range ents {
		names = append(names, e.Name())
	}
	sort.Strings(names)
	fmt.Printf("%d files in %s\n", len(names), outDir)
}
