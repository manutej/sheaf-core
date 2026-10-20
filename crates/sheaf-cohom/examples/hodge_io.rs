//! Plain-text driver: graph / 2-complex + edge cochains in, Betti numbers + Hodge split out (JSON).
//! Input (stdin): `V n`, `E a b` lines, `F a b c` lines, then blocks `C name` followed by one line of
//! edge values in sorted-edge order (a<b, lexicographic), orientation a->b (value = x_b - x_a convention).
//! Identity (constant R^1) sheaf, per attribute. Used by product/metaprompt-audit.
use sheaf_cell::{complex::CellComplex, sheaf::Sheaf};
use sheaf_cohom::cohomology::CochainComplex;
use std::io::Read;
fn main() {
    let mut s = String::new(); std::io::stdin().read_to_string(&mut s).unwrap();
    let (mut n, mut edges, mut faces, mut cochains) = (0usize, vec![], vec![], vec![]);
    let mut lines = s.lines().peekable();
    while let Some(l) = lines.next() {
        let t: Vec<&str> = l.split_whitespace().collect();
        if t.is_empty() { continue; }
        match t[0] {
            "V" => n = t[1].parse().unwrap(),
            "E" => edges.push((t[1].parse::<usize>().unwrap(), t[2].parse::<usize>().unwrap())),
            "F" => faces.push([t[1].parse::<usize>().unwrap(), t[2].parse().unwrap(), t[3].parse().unwrap()]),
            "C" => { let v: Vec<f64> = lines.next().unwrap().split_whitespace().map(|x| x.parse().unwrap()).collect(); cochains.push((t[1].to_string(), v)); }
            _ => panic!("bad line {l}"),
        }
    }
    let cx = CellComplex::simplicial(n, &edges, &faces);
    let ne = cx.n_cells[1];
    let sh = Sheaf::constant(cx, 1); sh.validate().unwrap();
    let cc = CochainComplex::from_sheaf(&sh).unwrap();
    let betti: Vec<String> = (0..cc.dims.len()).map(|k| cc.cohomology(k).unwrap().dim.to_string()).collect();
    let f = |v: &Vec<f64>| v.iter().map(|x| format!("{:.10e}", x)).collect::<Vec<_>>().join(",");
    let mut out = format!("{{\"dims\":{:?},\"betti\":[{}],\"cochains\":[", cc.dims, betti.join(","));
    for (i, (name, w)) in cochains.iter().enumerate() {
        assert_eq!(w.len(), ne, "cochain {name} has {} values, complex has {ne} edges", w.len());
        let h = cc.hodge(1, w); let (ee, eh, ec) = h.energies();
        out += &format!("{}{{\"name\":\"{}\",\"E_exact\":{:.10e},\"E_harm\":{:.10e},\"E_coex\":{:.10e},\"cert_max\":{:.3e},\"potential_cond\":{:.3e},\"harmonic\":[{}],\"coexact\":[{}],\"potential\":[{}]}}",
            if i > 0 { "," } else { "" }, name, ee, eh, ec, h.cert.max(), h.potential_cond, f(&h.harmonic), f(&h.coexact), f(&h.potential));
    }
    out += "]}";
    println!("{out}");
}
