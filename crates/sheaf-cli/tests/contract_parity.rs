//! Bridge sheaf-contract (m1, gated, byte-identical) <-> sheaf-cell through the adapter `sheaf_cli::m1`: on the m1
//! fixture spec, m1's per-edge `discord()` must equal ||(delta_0 x)_e||^2 from sheaf-cell's coboundary (and
//! `ops::edge_discord`), and the total must equal x^T L0 x (L0 = delta_0^T delta_0).
use sheaf_cell::{graph, num::sparse::Csr, ops, sheaf::Sheaf};
use sheaf_contract::SheafSpec;

fn load() -> SheafSpec {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/m1/spec.json");
    serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap()
}
/// m1 SheafSpec -> sheaf-cell Sheaf via the cli adapter (keeps m1's edge order and multi-edges; edge j = 1-cell j).
fn to_cohom(s: &SheafSpec) -> Sheaf { let (dims, edges) = sheaf_cli::m1::to_graph(s); graph::to_sheaf(&dims, &edges).unwrap() }

#[test]
fn m1_discord_equals_cohom_coboundary_per_edge() {
    let s = load();
    assert_eq!(s.validate(), Ok(()));
    let sh = to_cohom(&s);
    let (r, c, t) = sh.coboundary_triplets(0);
    let d0 = Csr::from_triplets(r, c, t);
    let x: Vec<f64> = s.x.iter().flatten().copied().collect();
    let dx = d0.mul(&x);
    let m1 = s.discord();
    let mut total = 0.0;
    for j in 0..s.edges.len() {
        let ours: f64 = dx[j * s.d..(j + 1) * s.d].iter().map(|v| v * v).sum();
        assert!((ours - m1[j]).abs() <= 1e-9 * m1[j].abs().max(1.0), "edge {j}: cohom {ours} vs m1 {}", m1[j]);
        total += ours;
    }
    let lx = d0.mul_t(&d0.mul(&x));
    let xlx: f64 = x.iter().zip(&lx).map(|(a, b)| a * b).sum();
    assert!((xlx - total).abs() <= 1e-9 * total.max(1.0));
    let (dims, edges) = sheaf_cli::m1::to_graph(&s);
    let (d0g, _, oe) = graph::coboundary(&dims, &edges);
    for (j, e) in ops::edge_discord(&d0g, &x, &oe).iter().enumerate() { assert!((e - m1[j]).abs() <= 1e-9 * m1[j].abs().max(1.0), "edge_discord {j}"); }
    assert!((ops::energy(&d0g, &x) - total).abs() <= 1e-9 * total.max(1.0));
    let m1_total: f64 = m1.iter().sum();
    assert!((m1_total - total).abs() <= 1e-9 * total.max(1.0));
}

/// m1 allows parallel edges (different types on the same pair); the bridge keeps them as separate 1-cells.
#[test]
fn m1_multi_edge_discord_parity() {
    use sheaf_contract::{Edge, EdgeType};
    let mut s = load();
    let e0 = s.edges.iter().find(|e| e.ty == EdgeType::About).map(|e| (e.u, e.v)).unwrap();
    s.edges.push(Edge { u: e0.0, v: e0.1, ty: EdgeType::Related });
    s.edges.push(Edge { u: e0.0, v: e0.1, ty: EdgeType::SameSlug });
    assert_eq!(s.validate(), Ok(()));
    let sh = to_cohom(&s);
    assert_eq!(sh.complex.n_cells[1], s.edges.len());
    let (r, c, t) = sh.coboundary_triplets(0);
    let dx = Csr::from_triplets(r, c, t).mul(&s.x.iter().flatten().copied().collect::<Vec<_>>());
    for (j, want) in s.discord().iter().enumerate() {
        let ours: f64 = dx[j * s.d..(j + 1) * s.d].iter().map(|v| v * v).sum();
        assert!((ours - want).abs() <= 1e-9 * want.abs().max(1.0), "edge {j}");
    }
}
