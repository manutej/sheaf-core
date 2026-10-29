//! Oracle extras for the Go port's golden tests (go/testdata/golden/*.aux.json); the report itself comes from the
//! `sheaf` binary. For one v1 spec (+ optional boundary file) prints one JSON object with the quantities a report
//! only shows rounded or truncated: the full dense L0 spectrum, the full per-edge Hodge split of the report's cochain,
//! the minimum-norm potential, and orthonormal bases of H0 (dense per-component `cohomology(0)`) and of the relative
//! H0(X,B;F), all under the spec's own H0 policy. Bases are compared through their projectors, never raw vectors.
//! usage: cargo run --release -p sheaf-cli --example go_oracle -- spec.json [boundary.json]
use faer::Side;
use sheaf_cell::graph::to_sheaf;
use sheaf_cli::v1;
use sheaf_cohom::cohomology::CochainComplex;
use sheaf_cohom::hodge_graph::graph_hodge_cg_cert;
use sheaf_cohom::spectral::HodgeOp;
use serde_json::json;

fn cols(m: &faer::Mat<f64>) -> Vec<Vec<f64>> { (0..m.ncols()).map(|j| (0..m.nrows()).map(|i| m[(i, j)]).collect()).collect() }

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let spec = v1::parse(&std::fs::read_to_string(&a[1]).expect("read spec")).expect("valid spec");
    let (d0, on, oe) = v1::coboundary(&spec.dims, &spec.edges);
    let op = HodgeOp::l0(&d0);
    let mut spectrum: Vec<f64> = if op.n == 0 { vec![] } else {
        let e = op.dense().self_adjoint_eigen(Side::Lower).expect("evd");
        (0..op.n).map(|i| e.S().column_vector()[i]).collect()
    };
    spectrum.sort_by(|x, y| x.partial_cmp(y).unwrap());
    let w = match (&spec.w, &spec.x) { (Some(w), _) => w.clone(), (None, Some(x)) => d0.mul(x), _ => unreachable!() };
    let g = graph_hodge_cg_cert(&d0, &w, 1e-12, 50_000, 1e-8);
    let cc = CochainComplex::from_sheaf(&to_sheaf(&spec.dims, &spec.edges).expect("sheaf")).expect("complex").with_h0_policy(spec.opts.h0);
    let h0_basis = if op.n == 0 { vec![] } else { cols(&cc.cohomology(0).expect("h0").harmonic_basis) };
    let (exact, harmonic, potential) = if d0.nrows == 0 { (vec![], vec![], vec![0.0; op.n]) } else { let h = cc.hodge(1, &w); (h.exact, h.harmonic, h.potential) };
    let relative = a.get(2).map(|p| {
        let b = v1::parse_boundary(&std::fs::read_to_string(p).expect("read boundary"), &spec).expect("valid boundary");
        let bd = sheaf_solve::node_boundary(&on, &b);
        let dofs: Vec<usize> = bd.iter().map(|x| x.0).collect();
        match sheaf_cohom::relative::h0_relative(&d0, &dofs, spec.opts.h0, spec.opts.dense_max) {
            Ok(r) => json!({"dim": r.value.dim, "basis": cols(&r.value.basis)}),
            Err(_) => serde_json::Value::Null,
        }
    });
    println!("{}", json!({
        "node_offsets": on, "edge_offsets": oe, "l0_bound": op.norm_bound(), "l0_spectrum": spectrum,
        "w": w, "exact": exact, "harmonic": harmonic, "potential": potential, "cg_harmonic": g.harmonic, "cg_potential": g.potential,
        "h0_basis": h0_basis, "relative": relative,
    }));
}
