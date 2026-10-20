//! Scale benchmark for the v1 spectral path (single run, box-local numbers).
//! usage: spectral_scale synthetic N EXTRA | arxiv EDGE_CSV | fb15k SPEC_TXT   [--k K] [--analyze NULL_REPS] [--method polynomial] [--out FILE]
use sheaf_cohom::spectral::*; use sheaf_cli::v1::{self, Edge};
use faer::Mat;
use std::time::Instant;
fn rss_peak_mb() -> f64 {
    let s = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
    s.lines().find(|l| l.starts_with("VmHWM")).and_then(|l| l.split_whitespace().nth(1)).and_then(|x| x.parse::<f64>().ok()).unwrap_or(0.0) / 1024.0
}
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let flag = |f: &str| a.iter().position(|x| x == f).map(|i| a[i + 1].clone());
    let k: usize = flag("--k").map(|x| x.parse().unwrap()).unwrap_or(6);
    let t0 = Instant::now();
    let (dims, edges, x): (Vec<usize>, Vec<Edge>, Option<Vec<f64>>) = match a[1].as_str() {
        "synthetic" => {
            let n: usize = a[2].parse().unwrap(); let extra: usize = a[3].parse().unwrap();
            let mut e: Vec<(usize, usize)> = (0..n).map(|i| (i.min((i + 1) % n), i.max((i + 1) % n))).collect();
            let mut st: u64 = 0x2545F4914F6CDD1D;
            for _ in 0..extra { st ^= st << 13; st ^= st >> 7; st ^= st << 17; let p = (st % n as u64) as usize; st ^= st << 13; st ^= st >> 7; st ^= st << 17; let q = (st % n as u64) as usize; if p != q { e.push((p.min(q), p.max(q))); } }
            e.sort(); e.dedup();
            let one = Mat::<f64>::identity(1, 1);
            (vec![1; n], e.into_iter().map(|(u, v)| Edge { u, v, fu: one.clone(), fv: one.clone(), label: String::new() }).collect(), None)
        }
        "arxiv" => {
            let s = std::fs::read_to_string(&a[2]).unwrap();
            let mut e: Vec<(usize, usize)> = s.lines().filter_map(|l| { let mut t = l.split(','); let p: usize = t.next()?.trim().parse().ok()?; let q: usize = t.next()?.trim().parse().ok()?; if p == q { None } else { Some((p.min(q), p.max(q))) } }).collect();
            e.sort(); e.dedup();
            let n = e.iter().map(|&(_, v)| v).max().unwrap() + 1;
            let one = Mat::<f64>::identity(1, 1);
            (vec![1; n], e.into_iter().map(|(u, v)| Edge { u, v, fu: one.clone(), fv: one.clone(), label: String::new() }).collect(), None)
        }
        "fb15k" => {
            let s = std::fs::read_to_string(&a[2]).unwrap(); let mut ls = s.lines();
            let h: Vec<usize> = ls.next().unwrap().split_whitespace().map(|x| x.parse().unwrap()).collect(); let (n, m, d) = (h[0], h[1], h[2]);
            let mut x = Vec::with_capacity(n * d); for _ in 0..n { x.extend(ls.next().unwrap().split_whitespace().map(|v| v.parse::<f64>().unwrap())); }
            let mut es = Vec::with_capacity(m);
            for _ in 0..m { let t: Vec<&str> = ls.next().unwrap().split_whitespace().collect();
                let f: Vec<f64> = t[3..].iter().map(|z| z.parse().unwrap()).collect();
                es.push(Edge { u: t[0].parse().unwrap(), v: t[1].parse().unwrap(), fu: Mat::from_fn(d, d, |i, j| f[i * d + j]), fv: Mat::from_fn(d, d, |i, j| f[d * d + i * d + j]), label: t[2].to_string() }); }
            (vec![d; n], es, Some(x))
        }
        _ => panic!("mode"),
    };
    let t_load = t0.elapsed().as_secs_f64();
    let t1 = Instant::now();
    let (d0, on, _oe) = v1::coboundary(&dims, &edges);
    let (_comp, ok, ncomp) = v1::components(dims.len(), &edges);
    let _ = &on;
    let t_build = t1.elapsed().as_secs_f64();
    let method = if flag("--method").as_deref() == Some("polynomial") { Method::Polynomial } else { Method::ShiftInvert };
    let o = SpectralOpts { k, method, ..Default::default() };
    let t2 = Instant::now();
    let hc = v1::spectral_components(&dims, &edges, &o); let sp = &hc.spectral;
    let t_spec = t2.elapsed().as_secs_f64();
    let rss_spec = rss_peak_mb();
    let mut out = serde_json::json!({"mode": a[1], "nodes": dims.len(), "edges": edges.len(), "dim_C0": d0.ncols, "dim_C1": d0.nrows, "components": ncomp,
        "known_constant_components": ok.iter().filter(|&&b| b).count(), "k": k, "method": sp.method, "load_s": t_load, "build_s": t_build, "spectral_s": t_spec,
        "rounds": sp.rounds, "converged": sp.converged, "inner_cg": sp.inner.map(|(s, i, r)| serde_json::json!([s, i, r])), "sigma": sp.sigma,
        "known_kernel_dim": sp.known_kernel_dim, "found_zero": sp.found_zero, "h0_dim": sp.h0_dim(), "zero_tol_max": hc.cut_max, "zero_tol_min": hc.cut_min, "cut_consistent": hc.cut_consistent, "near_sections": hc.near_sections, "lower_bound_only": sp.h0_lower_bound_only,
        "eigenvalues": sp.pairs.iter().map(|p| p.value).collect::<Vec<_>>(), "residuals": sp.pairs.iter().map(|p| p.residual).collect::<Vec<_>>(),
        "spectral_gap": sp.spectral_gap(), "lambda_max_bound": sp.lambda_max_bound, "peak_rss_mb_after_spectral": rss_spec});
    if let Some(r) = flag("--analyze") {
        let w = if x.is_none() { Some((0..d0.nrows).map(|i| ((i as f64) * 0.618).sin()).collect()) } else { None };
        let spec = v1::Spec { ids: (0..dims.len()).map(|i| i.to_string()).collect(), dims: dims.clone(), edges, x, w, opts: v1::Options { k, null_reps: r.parse().unwrap(), method, ..Default::default() }, boundary: None };
        let t3 = Instant::now(); let (rep, _) = v1::analyze(&spec, None); let rep = rep.to_value(); let ta = t3.elapsed().as_secs_f64();
        out["analyze_s"] = serde_json::json!(ta); out["analyze_hodge"] = rep["hodge"].clone(); out["analyze_h0"] = rep["h0"]["dim"].clone();
        out["peak_rss_mb_total"] = serde_json::json!(rss_peak_mb());
    }
    let s = serde_json::to_string_pretty(&out).unwrap();
    if let Some(p) = flag("--out") { std::fs::write(p, &s).unwrap(); }
    println!("{s}");
}
