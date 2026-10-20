//! Review F6: `CellComplex::graph` / `Sheaf::from_graph` keep edge order, multi-edges and u > v orientation, so the
//! sheaf-cohom coboundary equals v1's own assembled coboundary entry for entry (and H0 agrees under the one policy).
use faer::Mat;
use sheaf_cli::v1::{self, Edge};
use sheaf_cell::{complex::CellComplex, sheaf::Sheaf};
use sheaf_cohom::cohomology::CochainComplex;

fn lcg(s: &mut u64) -> f64 { *s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); ((*s >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0 }
fn rand_mat(r: usize, c: usize, s: &mut u64) -> Mat<f64> { Mat::from_fn(r, c, |_, _| lcg(s)) }

fn check(dims: &[usize], edges: &[Edge]) {
    let (d0, _, _) = v1::coboundary(dims, edges);
    let tuples: Vec<_> = edges.iter().map(|e| (e.u, e.v, e.fu.clone(), e.fv.clone())).collect();
    let sh = Sheaf::from_graph(dims, &tuples).unwrap();
    let ours = sh.coboundary_dense(0);
    assert_eq!((ours.nrows(), ours.ncols()), (d0.nrows, d0.ncols));
    let mut theirs = Mat::<f64>::zeros(d0.nrows, d0.ncols);
    for r in 0..d0.nrows { for k in d0.ptr[r]..d0.ptr[r + 1] { theirs[(r, d0.idx[k])] += d0.val[k]; } }
    for i in 0..d0.nrows { for j in 0..d0.ncols { assert_eq!(ours[(i, j)], theirs[(i, j)], "entry ({i},{j})"); } }
    // same H0 under the one policy: dense cohomology(0) vs the v1 per-component spectral driver
    let h_dense = CochainComplex::from_sheaf(&sh).unwrap().cohomology(0).unwrap().dim;
    let h_v1 = v1::spectral_components(dims, edges, &Default::default()).spectral.h0_dim();
    assert_eq!(h_dense, h_v1);
}

#[test]
fn graph_builder_keeps_order_multi_edges_and_orientation() {
    let mut s = 7u64;
    let dims = [2usize, 1, 3, 2, 2];
    let mk = |u: usize, v: usize, de: usize, s: &mut u64| Edge { u, v, fu: rand_mat(de, dims[u], s), fv: rand_mat(de, dims[v], s), label: String::new() };
    // u > v edges, a triple multi-edge (0,1) in both orientations, and an unsorted edge order
    let edges = vec![mk(3, 0, 2, &mut s), mk(0, 1, 1, &mut s), mk(1, 0, 2, &mut s), mk(0, 1, 3, &mut s), mk(4, 2, 2, &mut s), mk(2, 1, 1, &mut s)];
    check(&dims, &edges);
    // simplicial() would merge the multi-edge and re-orient: it must NOT be used for v1/m1 data
    let pairs: Vec<(usize, usize)> = edges.iter().map(|e| (e.u, e.v)).collect();
    assert!(CellComplex::simplicial(dims.len(), &pairs, &[]).n_cells[1] < edges.len());
    assert_eq!(CellComplex::graph(dims.len(), &pairs).unwrap().n_cells[1], edges.len());
}

#[test]
fn graph_builder_parity_random() {
    for seed in 1..30u64 {
        let mut s = seed;
        let n = 3 + (seed as usize % 6);
        let dims: Vec<usize> = (0..n).map(|i| 1 + (i + seed as usize) % 3).collect();
        let ne = n + (seed as usize % 5);
        let mut edges = vec![];
        for j in 0..ne {
            let u = (j * 7 + seed as usize) % n; let mut v = (j * 3 + 1 + seed as usize / 2) % n; if v == u { v = (u + 1) % n; }
            let de = 1 + (j % 3);
            edges.push(Edge { u, v, fu: rand_mat(de, dims[u], &mut s), fv: rand_mat(de, dims[v], &mut s), label: String::new() });
        }
        check(&dims, &edges);
    }
}

#[test]
fn graph_builder_rejects_self_loops_and_bad_indices() {
    assert!(CellComplex::graph(2, &[(0, 0)]).is_err());
    assert!(CellComplex::graph(2, &[(0, 2)]).is_err());
}
