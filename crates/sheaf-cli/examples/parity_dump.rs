//! For each v1 spec path: dense coboundary, dense L0, energy, per-edge discord, H0 (v1 spectral path). One JSON line each.
use sheaf_cohom::spectral::*; use sheaf_cli::v1;
fn main() {
    for p in std::env::args().skip(1) {
        let s = v1::parse(&std::fs::read_to_string(&p).unwrap()).unwrap();
        let (d0, on, oe) = v1::coboundary(&s.dims, &s.edges);
        let mut delta = vec![0.0; d0.nrows * d0.ncols];
        for i in 0..d0.nrows { for q in d0.ptr[i]..d0.ptr[i + 1] { delta[i * d0.ncols + d0.idx[q]] += d0.val[q]; } }
        let op = HodgeOp::l0(&d0); let l = op.dense();
        let lap: Vec<f64> = (0..op.n).flat_map(|i| (0..op.n).map(move |j| (i, j))).map(|(i, j)| l[(i, j)]).collect();
        let (rep, _) = v1::analyze(&s, None);
        let (energy, disc) = match &s.x { Some(x) => { let r = d0.mul(x); let mut y = vec![0.0; op.n]; op.apply(x, &mut y);
            (Some(dot(x, &y)), Some((0..s.edges.len()).map(|j| r[oe[j]..oe[j + 1]].iter().map(|z| z * z).sum::<f64>()).collect::<Vec<_>>())) }, None => (None, None) };
        let _ = on;
        let name = std::path::Path::new(&p).file_stem().unwrap().to_string_lossy().to_string();
        println!("{}", serde_json::json!({"name": name, "delta": delta, "laplacian": lap, "energy": energy, "edge_discord": disc, "h0": rep.h0.dim}));
    }
}
