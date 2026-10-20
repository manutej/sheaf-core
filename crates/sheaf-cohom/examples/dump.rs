//! Writes fixtures (dims, coboundary triplets, Rust betti + rank gaps) as JSON for the Python/Julia cross-check.
use sheaf_cell::fixtures::catalogue;
use sheaf_cohom::cohomology::CochainComplex;
fn main() {
    let mut out = String::from("[\n");
    let cat = catalogue();
    for (n_i, n) in cat.iter().enumerate() {
        let cc = CochainComplex::from_sheaf(&n.sheaf).unwrap();
        let betti: Vec<String> = (0..cc.dims.len()).map(|k| cc.cohomology(k).unwrap().dim.to_string()).collect();
        let gaps: Vec<String> = (0..cc.dims.len()).map(|k| format!("{:.3e}", cc.cohomology(k).unwrap().min_gap_ratio().min(1e300))).collect();
        let mut deltas = Vec::new();
        for k in 0..cc.delta.len() {
            let (r, c, t) = n.sheaf.coboundary_triplets(k);
            let ts: Vec<String> = t.iter().map(|(i, j, v)| format!("[{},{},{:.17e}]", i, j, v)).collect();
            deltas.push(format!("{{\"rows\":{},\"cols\":{},\"t\":[{}]}}", r, c, ts.join(",")));
        }
        let edges: Vec<String> = n.sheaf.complex.verts[1].iter().map(|e| format!("[{},{}]", e[0], e[1])).collect();
        out += &format!("{{\"name\":\"{}\",\"dims\":{:?},\"betti_rust\":[{}],\"min_gap\":[{}],\"edges\":[{}],\"delta\":[{}]}}{}\n",
            n.name, cc.dims, betti.join(","), gaps.join(","), edges.join(","), deltas.join(","), if n_i + 1 < cat.len() { "," } else { "" });
    }
    out += "]\n";
    std::fs::write("out/fixtures.json", out).unwrap();
    println!("wrote out/fixtures.json ({} fixtures)", cat.len());
}
