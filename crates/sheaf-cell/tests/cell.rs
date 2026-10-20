//! sheaf-cell: the single graph δ⁰ assembly, energy / discord, CG minimum-norm, family labels.
use faer::Mat;
use sheaf_cell::{graph::{self, GraphEdge}, num::{family::*, krylov::{cg, SolveTol}, norm, sparse::Csr}, ops, sheaf::Sheaf};

fn lcg(s: &mut u64) -> f64 { *s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); ((*s >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0 }
fn edges(seed: u64) -> (Vec<usize>, Vec<GraphEdge>) {
    let dims = vec![2, 3, 1, 2, 2]; let mut s = seed;
    // includes u > v, parallel edges and a non-square map
    let pairs = [(0, 1), (2, 1), (1, 2), (3, 0), (4, 3), (3, 4), (4, 2)];
    let de = [2, 1, 2, 3, 2, 1, 2];
    let es = pairs.iter().zip(de).map(|(&(u, v), d)| GraphEdge { u, v, fu: Mat::from_fn(d, dims[u], |_, _| lcg(&mut s)), fv: Mat::from_fn(d, dims[v], |_, _| lcg(&mut s)), label: String::new() }).collect();
    (dims, es)
}

#[test]
fn graph_coboundary_is_the_sheaf_coboundary_bitwise() {
    let (dims, es) = edges(3);
    let (d0, on, oe) = graph::coboundary(&dims, &es);
    let t: Vec<_> = es.iter().map(|e| (e.u, e.v, e.fu.clone(), e.fv.clone())).collect();
    let s = Sheaf::from_graph(&dims, &t).unwrap();
    let dense = s.coboundary_dense(0);
    assert_eq!((d0.nrows, d0.ncols), (dense.nrows(), dense.ncols()));
    assert_eq!(on, s.offsets(0)); assert_eq!(oe, s.offsets(1));
    // independent per-edge formula (δx)_e = F_v x_v − F_u x_u
    let mut sd = 11; let x: Vec<f64> = (0..d0.ncols).map(|_| lcg(&mut sd)).collect();
    let dx = d0.mul(&x);
    for (j, e) in es.iter().enumerate() { for r in 0..e.fu.nrows() {
        let want: f64 = (0..e.fv.ncols()).map(|c| e.fv[(r, c)] * x[on[e.v] + c]).sum::<f64>() - (0..e.fu.ncols()).map(|c| e.fu[(r, c)] * x[on[e.u] + c]).sum::<f64>();
        assert!((dx[oe[j] + r] - want).abs() < 1e-14, "edge {j} row {r}");
    } }
    let disc = ops::edge_discord(&d0, &x, &oe);
    assert!((disc.iter().sum::<f64>() - ops::energy(&d0, &x)).abs() < 1e-12);
    let (comp, _, nc) = graph::components(dims.len(), &es);
    assert_eq!(nc, 1); assert!(comp.iter().all(|&c| c == 0));
}

#[test]
fn cg_from_zero_is_minimum_norm_on_a_singular_system() {
    // constant sheaf on a path: L0 singular (constants); b in range; CG from 0 must have zero mean
    let n = 9; let t: Vec<_> = (0..n - 1).map(|i| (i, i + 1, Mat::<f64>::identity(1, 1), Mat::<f64>::identity(1, 1))).collect();
    let s = Sheaf::from_graph(&vec![1; n], &t).unwrap();
    let (r, c, tr) = s.coboundary_triplets(0); let d = Csr::from_triplets(r, c, tr);
    let mut sd = 5; let xt: Vec<f64> = (0..n).map(|_| lcg(&mut sd)).collect();
    let b = d.mul_t(&d.mul(&xt));
    let ap = |v: &[f64]| d.mul_t(&d.mul(v));
    let out = cg(&ap, &b, vec![0.0; n], norm(&b), 4.0, SolveTol::default());
    assert!(!out.breakdown && out.residual_norm <= 1e-12 * norm(&b));
    assert!(out.x.iter().sum::<f64>().abs() < 1e-10, "min-norm solution is orthogonal to the constants");
    let mean = xt.iter().sum::<f64>() / n as f64;
    for i in 0..n { assert!((out.x[i] - (xt[i] - mean)).abs() < 1e-9); }
}

#[test]
fn family_and_arithmetic_labels_are_typed() {
    assert_eq!(serde_json::to_string(&Numerical::f64(3usize)).unwrap(), r#"{"value":3,"method":"exact","arithmetic":"f64"}"#);
    let e = ExactFp::new(vec![1usize, 1], 2);
    assert_eq!(e.p(), 2);
    assert_eq!(serde_json::to_string(&e).unwrap(), r#"{"value":[1,1],"method":"exact","arithmetic":"F_2"}"#);
    let a = Approximate::new(0.5f64, vec!["interval box".into()]);
    assert!(serde_json::to_string(&a).unwrap().contains(r#""method":"approximate""#));
    assert_eq!(serde_json::to_string(&Family::Approximate).unwrap(), r#""approximate""#);
}
