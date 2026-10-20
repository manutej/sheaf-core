//! v1 pipeline core: JSON spec with PER-EDGE restriction maps -> deterministic Report + Effects (ADR-005).
//!
//! Spec (format "frontier-sheaf/v1"):
//! { "format": "frontier-sheaf/v1",
//!   "nodes": [ {"id": "a", "dim": 2}, ... ],
//!   "edges": [ {"u": 0 | "a", "v": 1 | "b", "Fu": [[..]], "Fv": [[..]], "label": "..."}, ... ],
//!   "x": [[..per node..], ...]            (optional node data; then w := delta x)
//!   "w": [[..per edge..], ...]            (optional edge cochain, e.g. observed per-hop drift)
//!   "options": {"k": 6, "null_reps": 20, "seed": 1, "top": 10} }
//! Convention (Hansen-Ghrist, = sheafcore.py): (delta x)_e = Fv x_v - Fu x_u, Fu: R^{dim u} -> R^{de}.
//! Fu/Fv omitted => identity (requires dim u == dim v). Edge stalk dim de = rows of Fu (or "de" field).
//! Graph only (no 2-cells) in v1: the Hodge split of w is exact + harmonic (no coexact part).
//! A9: if only x is given, w = delta x is exact by construction, so its harmonic share is 0 (reported, not a finding).
use sheaf_cell::num::sparse::Csr;
use sheaf_cohom::hodge_graph::graph_hodge_cg_cert;
use sheaf_cohom::spectral::{self, dot, H0Policy, HodgeOp, KnownKernel, Method, SpectralOpts};
use faer::Mat;
use crate::report::*;
use serde_json::Value;

#[derive(Debug, Clone, PartialEq)]
pub enum SpecError { Json(String), Field(String) }
impl std::fmt::Display for SpecError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { match self { SpecError::Json(s) => write!(f, "P_INPUT_01 json: {s}"), SpecError::Field(s) => write!(f, "P_INPUT_01 spec: {s}") } }
}

// The graph sheaf, δ⁰ assembly, components and the per-component H0 driver live in sheaf-cell / sheaf-cohom
// (one implementation each); re-exported here for callers of the v1 API.
pub use sheaf_cell::graph::{components, coboundary, offsets, GraphEdge as Edge};
pub use sheaf_cohom::spectral::{spectral_components, H0Components};

#[derive(Debug, Clone)]
pub struct Options { pub k: usize, pub null_reps: usize, pub seed: u64, pub top: usize, pub method: Method, pub dense_max: usize,
                     /// the one H0 zero-cut policy (ADR-015): default per-component ADR-012; `options.h0_absolute_tol` = permitted delta
                     pub h0: H0Policy }
impl Default for Options { fn default() -> Self { Options { k: 6, null_reps: 20, seed: 1, top: 10, method: Method::ShiftInvert, dense_max: 1500, h0: H0Policy::adr012() } } }

#[derive(Debug, Clone)]
pub struct Spec { pub ids: Vec<String>, pub dims: Vec<usize>, pub edges: Vec<Edge>, pub x: Option<Vec<f64>>, pub w: Option<Vec<f64>>, pub opts: Options,
                  /// boundary data for harmonic extension: (node, value in its stalk); set by `parse_boundary` / `--boundary`
                  pub boundary: Option<Vec<(usize, Vec<f64>)>> }

fn num(v: &Value, what: &str) -> Result<f64, SpecError> {
    let x = v.as_f64().ok_or_else(|| SpecError::Field(format!("{what}: not a number")))?;
    if !x.is_finite() { return Err(SpecError::Field(format!("{what}: not finite"))); }
    // ADR-017 P-NUM-01: finite but huge entries overflowed L0 = δᵀδ / xᵀLx to inf/NaN and the run PANICKED (faer evd
    // NoConvergence at 1e300, NaN sort at 1e200). |entry| <= 1e75 keeps every product chain (x·F·F·x) inside f64.
    if x.abs() > MAX_ABS_ENTRY { return Err(SpecError::Field(format!("{what}: |value| {x:e} exceeds {MAX_ABS_ENTRY:e} (L0 / energy would overflow f64)"))); }
    Ok(x)
}
fn matrix(v: &Value, rows: Option<usize>, cols: usize, what: &str) -> Result<Mat<f64>, SpecError> {
    let r = v.as_array().ok_or_else(|| SpecError::Field(format!("{what}: must be a list of rows")))?;
    if let Some(n) = rows { if r.len() != n { return Err(SpecError::Field(format!("{what}: {} rows, expected {n}", r.len()))); } }
    let mut m = Mat::<f64>::zeros(r.len(), cols);
    for (i, row) in r.iter().enumerate() {
        let row = row.as_array().ok_or_else(|| SpecError::Field(format!("{what}: row {i} not a list")))?;
        if row.len() != cols { return Err(SpecError::Field(format!("{what}: row {i} has {} cols, expected {cols}", row.len()))); }
        for (j, x) in row.iter().enumerate() { m[(i, j)] = num(x, what)?; }
    }
    Ok(m)
}
fn vectors(v: &Value, lens: &[usize], what: &str) -> Result<Vec<f64>, SpecError> {
    let r = v.as_array().ok_or_else(|| SpecError::Field(format!("{what}: must be a list")))?;
    if r.len() != lens.len() { return Err(SpecError::Field(format!("{what}: {} entries, expected {}", r.len(), lens.len()))); }
    let mut out = Vec::new();
    for (i, (row, &d)) in r.iter().zip(lens).enumerate() {
        let row = row.as_array().ok_or_else(|| SpecError::Field(format!("{what}[{i}]: not a list")))?;
        if row.len() != d { return Err(SpecError::Field(format!("{what}[{i}]: length {}, expected {d}", row.len()))); }
        for x in row { out.push(num(x, what)?); }
    }
    Ok(out)
}

/// Total node (and per-edge) stalk dimension cap, as m1's R-DIM-02 MAX_TOTAL_DIM.
pub const MAX_TOTAL_DIM: usize = (1 << 31) - 1;
/// Total dense restriction-map storage admitted by `parse` (m1's MAX_SPARSE_BYTES, 4 GiB).
pub const MAX_MAP_BYTES: usize = 4 << 30;
/// Largest admitted |entry| in maps, x, w and boundary values (ADR-017 P-NUM-01).
pub const MAX_ABS_ENTRY: f64 = 1e75;

pub fn parse(text: &str) -> Result<Spec, SpecError> {
    let v: Value = serde_json::from_str(text).map_err(|e| SpecError::Json(e.to_string()))?;
    if v.get("format").and_then(|f| f.as_str()) != Some("frontier-sheaf/v1") { return Err(SpecError::Field("format must be \"frontier-sheaf/v1\"".into())); }
    let nodes = v.get("nodes").and_then(|n| n.as_array()).ok_or_else(|| SpecError::Field("nodes: missing list".into()))?;
    let mut ids = vec![]; let mut dims = vec![];
    for (i, n) in nodes.iter().enumerate() {
        let id = n.get("id").and_then(|x| x.as_str()).map(String::from).unwrap_or_else(|| i.to_string());
        let d = n.get("dim").and_then(|x| x.as_u64()).ok_or_else(|| SpecError::Field(format!("nodes[{i}].dim: missing non-negative integer")))? as usize;
        ids.push(id); dims.push(d);
    }
    // ADR-017 P-RES-01: admission BEFORE any allocation. A 120-byte spec with two nodes of dim 200000 and an edge with
    // omitted maps used to allocate a 200000x200000 identity (320 GB) and abort; dim 2^64-1 wrapped the offsets.
    if dims.iter().try_fold(0usize, |a, &d| a.checked_add(d)).map_or(true, |t| t > MAX_TOTAL_DIM) {
        return Err(SpecError::Field(format!("nodes: total stalk dimension exceeds {MAX_TOTAL_DIM}")));
    }
    let index: std::collections::HashMap<&str, usize> = ids.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
    if index.len() != ids.len() { return Err(SpecError::Field("nodes: duplicate id".into())); }
    let endpoint = |x: Option<&Value>, j: usize, w: &str| -> Result<usize, SpecError> {
        match x {
            Some(Value::Number(n)) => n.as_u64().map(|u| u as usize).filter(|&u| u < ids.len()).ok_or_else(|| SpecError::Field(format!("edges[{j}].{w}: index out of range"))),
            Some(Value::String(s)) => index.get(s.as_str()).cloned().ok_or_else(|| SpecError::Field(format!("edges[{j}].{w}: unknown node id {s}"))),
            _ => Err(SpecError::Field(format!("edges[{j}].{w}: missing"))),
        }
    };
    let es = v.get("edges").and_then(|n| n.as_array()).ok_or_else(|| SpecError::Field("edges: missing list".into()))?;
    let mut edges = vec![]; let mut map_total = 0usize;
    for (j, e) in es.iter().enumerate() {
        let u = endpoint(e.get("u"), j, "u")?; let w_ = endpoint(e.get("v"), j, "v")?;
        if u == w_ { return Err(SpecError::Field(format!("edges[{j}]: self-loop"))); }
        let de_field = e.get("de").and_then(|x| x.as_u64()).map(|x| x as usize);
        let infer = |f: Option<&Value>| f.and_then(|m| m.as_array()).map(|r| r.len());
        let de = de_field.or(infer(e.get("Fu"))).or(infer(e.get("Fv"))).unwrap_or(dims[u]);
        let mk = |f: Option<&Value>, d: usize, w: &str| -> Result<Mat<f64>, SpecError> {
            match f { Some(m) => matrix(m, Some(de), d, &format!("edges[{j}].{w}")),
                      None => if d == de { Ok(Mat::identity(de, d)) } else { Err(SpecError::Field(format!("edges[{j}].{w}: omitted but dim {d} != de {de}"))) } }
        };
        // ADR-017 P-RES-01: the two dense maps (de x dim u, de x dim v) are measured before they are allocated
        let map_bytes = de.checked_mul(dims[u].saturating_add(dims[w_])).and_then(|k| k.checked_mul(8));
        map_total = map_bytes.and_then(|b| map_total.checked_add(b)).unwrap_or(usize::MAX);
        if de > MAX_TOTAL_DIM || map_total > MAX_MAP_BYTES { return Err(SpecError::Field(format!("edges[{j}]: restriction maps exceed {MAX_MAP_BYTES} bytes (dense storage)"))); }
        let fu = mk(e.get("Fu"), dims[u], "Fu")?; let fv = mk(e.get("Fv"), dims[w_], "Fv")?;
        let label = e.get("label").and_then(|x| x.as_str()).map(String::from).unwrap_or_else(|| format!("{}->{}", ids[u], ids[w_]));
        edges.push(Edge { u, v: w_, fu, fv, label });
    }
    let des: Vec<usize> = edges.iter().map(|e| e.fu.nrows()).collect();
    let x = match v.get("x") { Some(Value::Null) | None => None, Some(xv) => Some(vectors(xv, &dims, "x")?) };
    let w = match v.get("w") { Some(Value::Null) | None => None, Some(wv) => Some(vectors(wv, &des, "w")?) };
    if x.is_none() && w.is_none() { return Err(SpecError::Field("need x (node data) or w (edge cochain)".into())); }
    let mut opts = Options::default();
    if let Some(o) = v.get("options") {
        if let Some(k) = o.get("k").and_then(|x| x.as_u64()) { opts.k = (k as usize).max(1); }
        if let Some(k) = o.get("null_reps").and_then(|x| x.as_u64()) { opts.null_reps = k as usize; }
        if let Some(k) = o.get("seed").and_then(|x| x.as_u64()) { opts.seed = k; }
        if let Some(k) = o.get("top").and_then(|x| x.as_u64()) { opts.top = k as usize; }
        if let Some(k) = o.get("dense_max").and_then(|x| x.as_u64()) { opts.dense_max = k as usize; }
        if let Some(t) = o.get("h0_absolute_tol") {
            let t = num(t, "options.h0_absolute_tol")?;
            if t < 0.0 { return Err(SpecError::Field("options.h0_absolute_tol: must be >= 0".into())); }
            opts.h0 = H0Policy::absolute(t);
        }
        match o.get("method").and_then(|x| x.as_str()) { None | Some("shift-invert") => {}, Some("polynomial") => opts.method = Method::Polynomial,
            Some(m) => return Err(SpecError::Field(format!("options.method: unknown {m}"))) }
    }
    Ok(Spec { ids, dims, edges, x, w, opts, boundary: None })
}

/// Boundary file (`--boundary`): {"format": "frontier-sheaf/boundary-v1" (optional),
///   "boundary": [ {"node": "a" | 0, "value": [..stalk dim..]}, ... ]}
pub fn parse_boundary(text: &str, spec: &Spec) -> Result<Vec<(usize, Vec<f64>)>, SpecError> {
    let v: Value = serde_json::from_str(text).map_err(|e| SpecError::Json(e.to_string()))?;
    if let Some(f) = v.get("format") { if f.as_str() != Some("frontier-sheaf/boundary-v1") { return Err(SpecError::Field("boundary format must be \"frontier-sheaf/boundary-v1\"".into())); } }
    let list = v.get("boundary").and_then(|b| b.as_array()).ok_or_else(|| SpecError::Field("boundary: missing list".into()))?;
    let mut seen = vec![false; spec.dims.len()]; let mut out = vec![];
    for (j, b) in list.iter().enumerate() {
        let node = match b.get("node") {
            Some(Value::Number(n)) => n.as_u64().map(|u| u as usize).filter(|&u| u < spec.dims.len()),
            Some(Value::String(s)) => spec.ids.iter().position(|i| i == s),
            _ => None,
        }.ok_or_else(|| SpecError::Field(format!("boundary[{j}].node: missing or unknown")))?;
        if seen[node] { return Err(SpecError::Field(format!("boundary[{j}].node: duplicate"))); }
        seen[node] = true;
        let val = b.get("value").and_then(|x| x.as_array()).ok_or_else(|| SpecError::Field(format!("boundary[{j}].value: not a list")))?;
        if val.len() != spec.dims[node] { return Err(SpecError::Field(format!("boundary[{j}].value: length {}, expected {}", val.len(), spec.dims[node]))); }
        out.push((node, val.iter().map(|x| num(x, &format!("boundary[{j}].value"))).collect::<Result<Vec<_>, _>>()?));
    }
    Ok(out)
}

/// Harmonic extension (sheaf-solve) + uniqueness H⁰(X,B;F) (sheaf-cohom::relative), composed here.
fn extension(spec: &Spec, d0: &Csr, on: &[usize], bnodes: &[(usize, Vec<f64>)]) -> ExtensionReport {
    let o = &spec.opts;
    let bd = sheaf_solve::node_boundary(on, bnodes);
    let bdofs: Vec<usize> = bd.iter().map(|p| p.0).collect();
    let rel = sheaf_cohom::relative::h0_relative(d0, &bdofs, o.h0, o.dense_max).ok();
    let h = sheaf_solve::harmonic_extension(d0, &bd, None, &sheaf_solve::Options { dense_max: o.dense_max.min(400), ..Default::default() })
        .expect("boundary validated by parse_boundary");
    let isb: std::collections::HashSet<usize> = bnodes.iter().map(|b| b.0).collect();
    ExtensionReport {
        method: Exact, boundary_nodes: bnodes.len(), boundary_dofs: bd.len(),
        h0_relative: rel.as_ref().map(|r| r.value.dim), unique: rel.as_ref().map(|r| r.value.unique()),
        uniqueness_rule: format!("unique iff dim H0(X,B;F) = dim ker L_II = 0 (Hansen-Ghrist Prop 4.1); zero cut: {}", o.h0.describe()),
        note: "values = the minimum-norm minimiser of ||delta x||^2 with x_B fixed; if not unique it is determined only up to H0(X,B;F)",
        route: h.cert.route, converged: h.cert.converged, residual_rel: sig(h.cert.residual_rel), iters: h.cert.iters, energy: sig(h.cert.energy),
        values: (0..spec.dims.len()).map(|i| NodeValue { id: spec.ids[i].clone(), boundary: isb.contains(&i), value: h.x[on[i]..on[i + 1]].iter().map(|&v| sig(v)).collect() }).collect(),
    }
}

fn sig(x: f64) -> f64 { // 12 significant digits: reproducible across thread counts / platforms at report level (NaN -> null)
    if x == 0.0 || !x.is_finite() { return if x.is_finite() { 0.0 } else { f64::NAN }; }
    format!("{:.11e}", x).parse::<f64>().unwrap()
}
struct Rng(u64);
impl Rng {
    fn u(&mut self) -> u64 { let mut s = self.0; s ^= s << 13; s ^= s >> 7; s ^= s << 17; self.0 = s; s }
    fn gauss(&mut self) -> f64 { let a = ((self.u() >> 11) as f64 + 0.5) / (1u64 << 53) as f64; let b = ((self.u() >> 11) as f64) / (1u64 << 53) as f64; (-2.0 * a.ln()).sqrt() * (2.0 * std::f64::consts::PI * b).cos() }
}
fn quant(v: &mut Vec<f64>, q: f64) -> f64 { v.sort_by(|a, b| a.partial_cmp(b).unwrap()); if v.is_empty() { return f64::NAN; } v[((v.len() - 1) as f64 * q).round() as usize] }

#[derive(Debug, Clone, PartialEq)]
pub enum Effect { EmitReport { path: Option<String>, contents: String } }

pub struct HodgeShare { pub exact: f64, pub harmonic: f64, pub cert: f64, pub converged: bool, pub iters: usize }
fn split(d0: &Csr, w: &[f64]) -> (HodgeShare, sheaf_cohom::hodge_graph::GraphHodge) {
    let g = graph_hodge_cg_cert(d0, w, 1e-12, 50_000, 1e-8);
    let t = dot(w, w).max(1e-300);
    (HodgeShare { exact: dot(&g.exact, &g.exact) / t, harmonic: dot(&g.harmonic, &g.harmonic) / t, cert: g.cert_coclosed, converged: g.converged, iters: g.iters }, g)
}

/// Pure core: spec -> (typed Report, effects). No I/O.
pub fn analyze(spec: &Spec, out_path: Option<String>) -> (Report, Vec<Effect>) {
    let o = &spec.opts;
    let (d0, on, oe) = coboundary(&spec.dims, &spec.edges);
    let (n0, n1) = (d0.ncols, d0.nrows);
    let (comp, known_ok, ncomp) = components(spec.dims.len(), &spec.edges);
    let node_comp: Vec<Option<u32>> = comp.iter().map(|&c| if known_ok[c] { Some(c as u32) } else { None }).collect();
    let known = KnownKernel::from_components(&node_comp, &on);
    let so = SpectralOpts { k: o.k, method: o.method, dense_max: o.dense_max, h0: o.h0, seed: o.seed.wrapping_mul(0x9E3779B97F4A7C15) | 1, ..Default::default() };
    let _ = &known;
    let hc = spectral_components(&spec.dims, &spec.edges, &so);
    let sp = &hc.spectral;
    let h0 = sp.h0_dim();
    let h1 = (n1 + h0).saturating_sub(n0); // rank-nullity on a graph: dim C1 - rank d0, rank d0 = dim C0 - h0
    // L1 = d0 d0^T: its nonzero spectrum equals L0's; verify by mapping u = d0 v / sqrt(lambda) and checking the residual on L1
    let l1 = HodgeOp::l1_graph(&d0);
    let mut l1_res = 0.0f64; let mut l1_checked = 0;
    for (p, v) in sp.pairs.iter().zip(&sp.vectors).skip(sp.found_zero) {
        let mut u = d0.mul(v); let nu = spectral::norm(&u); if nu < 1e-300 { continue; }
        u.iter_mut().for_each(|x| *x /= nu); l1_res = l1_res.max(l1.residual(&u, p.value)); l1_checked += 1;
    }
    let k = o.k;
    let mut small: Vec<f64> = vec![0.0; sp.known_kernel_dim.min(k + 1)];
    small.extend(sp.pairs.iter().map(|p| p.value)); small.sort_by(|a, b| a.partial_cmp(b).unwrap()); small.truncate(k + 1);
    let n_zero_pairs = sp.found_zero;
    let res: Vec<f64> = sp.pairs.iter().skip(n_zero_pairs.saturating_sub(k + 1)).map(|p| sig(p.residual)).collect();
    let (zmax, nzmin) = sp.cut();
    let gap_k = if small.len() > k { Some(small[k] - small[k - 1]) } else { None };
    // ---- cochain + Hodge split ----
    let (w, source) = match (&spec.w, &spec.x) { (Some(w), _) => (w.clone(), "w (given edge cochain)"), (None, Some(x)) => (d0.mul(x), "delta x (derived from node data; exact by construction, A9)") , _ => unreachable!() };
    let (hs, g) = split(&d0, &w);
    let wn2 = dot(&w, &w);
    // ---- nulls: (a) iid Gaussian edge data, same shape; (b) shuffled edge blocks (within equal stalk dim) ----
    let mut rng = Rng(o.seed.wrapping_mul(0xD1B54A32D192ED03) | 1);
    let mut nr = vec![]; let mut ns = vec![];
    for _ in 0..o.null_reps {
        let r: Vec<f64> = (0..n1).map(|_| rng.gauss()).collect(); nr.push(split(&d0, &r).0.harmonic);
        let mut perm: Vec<usize> = (0..spec.edges.len()).collect();
        let mut by_de: std::collections::BTreeMap<usize, Vec<usize>> = Default::default();
        for (j, e) in spec.edges.iter().enumerate() { by_de.entry(e.fu.nrows()).or_default().push(j); }
        for g_ in by_de.values() { let mut gg = g_.clone(); for i in (1..gg.len()).rev() { let j = (rng.u() % (i as u64 + 1)) as usize; gg.swap(i, j); } for (a, b) in g_.iter().zip(gg) { perm[*a] = b; } }
        let mut sw = vec![0.0; n1];
        for j in 0..spec.edges.len() { let src = perm[j]; for r in 0..(oe[j + 1] - oe[j]) { sw[oe[j] + r] = w[oe[src] + r]; } }
        if dot(&sw, &sw) > 0.0 { ns.push(split(&d0, &sw).0.harmonic); }
    }
    let ge = |v: &Vec<f64>| if v.is_empty() { None } else { Some(sig(v.iter().filter(|&&x| x >= hs.harmonic - 1e-12).count() as f64 / v.len() as f64)) };
    let mean = |v: &Vec<f64>| if v.is_empty() { None } else { Some(sig(v.iter().sum::<f64>() / v.len() as f64)) };
    let null_random = NullRandom { reps: nr.len(), mean: mean(&nr), p05: sig(quant(&mut nr.clone(), 0.05)), p95: sig(quant(&mut nr.clone(), 0.95)), frac_ge_observed: ge(&nr),
                             analytic_expectation_h1_over_dim_c1: sig(if n1 > 0 { h1 as f64 / n1 as f64 } else { 0.0 }) };
    let degenerate = ns.is_empty() || ns.iter().all(|&v| (v - hs.harmonic).abs() <= 1e-12 * hs.harmonic.abs().max(1e-300) + 1e-15);
    let null_shuffle = NullShuffle { degenerate, note: "degenerate = every shuffle gives the observed share (e.g. <2 edges per stalk-dim bucket, or a single cycle where shuffling preserves the loop sum): uninformative, not 'not significant'", reps: ns.len(), mean: mean(&ns), p05: sig(quant(&mut ns.clone(), 0.05)), p95: sig(quant(&mut ns.clone(), 0.95)), frac_ge_observed: ge(&ns) };
    // ---- rankings ----
    let blk = |v: &[f64], o_: &[usize], j: usize| -> f64 { v[o_[j]..o_[j + 1]].iter().map(|x| x * x).sum() };
    let mut er: Vec<(usize, f64, f64)> = (0..spec.edges.len()).map(|j| (j, blk(&g.harmonic, &oe, j), blk(&w, &oe, j))).collect();
    er.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap().then(a.0.cmp(&b.0)));
    let edge_row = |&(j, h, t): &(usize, f64, f64)| { let e = &spec.edges[j]; EdgeRow { edge: j, u: spec.ids[e.u].clone(), v: spec.ids[e.v].clone(), label: e.label.clone(), harmonic_energy: sig(h), discord: sig(t), harmonic_frac_of_total: sig(if wn2 > 0.0 { h / wn2 } else { 0.0 }) } };
    let edges_h: Vec<EdgeRow> = er.iter().filter(|r| r.1 > 1e-24 * wn2.max(1e-300)).take(o.top).map(edge_row).collect();
    let mut ed = er.clone(); ed.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap().then(a.0.cmp(&b.0)));
    let edges_d: Vec<EdgeRow> = ed.iter().filter(|r| r.2 > 0.0).take(o.top).map(edge_row).collect();
    let mut nrk: Vec<(usize, f64)> = (0..spec.dims.len()).map(|i| (i, blk(&g.potential, &on, i).sqrt())).collect();
    nrk.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap().then(a.0.cmp(&b.0)));
    let nodes_p: Vec<NodePotential> = nrk.iter().filter(|r| r.1 > 0.0).take(o.top).map(|&(i, p)| NodePotential { node: spec.ids[i].clone(), abs_potential: sig(p) }).collect();
    // incident edge discord per node: edge-stalk units (the maps' normalisation), so unit-independent unlike |potential| (A8)
    let mut inc = vec![0.0f64; spec.dims.len()];
    for &(j, _, t) in &er { inc[spec.edges[j].u] += t; inc[spec.edges[j].v] += t; }
    let mut irk: Vec<(usize, f64)> = inc.into_iter().enumerate().collect();
    irk.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap().then(a.0.cmp(&b.0)));
    let nodes_i: Vec<NodeDiscord> = irk.iter().filter(|r| r.1 > 0.0).take(o.top).map(|&(i, d)| NodeDiscord { node: spec.ids[i].clone(), incident_discord: sig(d) }).collect();
    let energy: Option<f64> = spec.x.as_ref().map(|x| sig(dot(&d0.mul(x), &d0.mul(x))));
    let pair = |p: &spectral::EigPair| Pair { value: sig(p.value), residual: sig(p.residual) };
    let report = Report {
        format: "frontier-sheaf/report-v1",
        // ADR-013: core cohomology outputs are the exact family (numerical linear algebra on the exact complex, with
        // certified residuals); an ordered/interval layer must carry "approximate" and never write into these fields.
        method: Exact,
        input: Input { nodes: spec.dims.len(), edges: spec.edges.len(), dim_c0: n0, dim_c1: n1, components: ncomp,
            known_constant_components: known_ok.iter().filter(|&&b| b).count(), has_x: spec.x.is_some(), has_w: spec.w.is_some() },
        h0: H0Report { method: Exact, dim: h0, known_kernel_dim: sp.known_kernel_dim, found_zero_by_solver: sp.found_zero,
            zero_policy: so.h0.name(), zero_rule: so.h0.describe(), zero_tol: sig(hc.cut_max), zero_tol_min: sig(hc.cut_min),
            largest_counted_zero: zmax.map(pair), smallest_counted_nonzero: nzmin.map(pair),
            cut_consistent_with_zero_tol: hc.cut_consistent,
            // one rule, so "ambiguous" now means: the cut is not consistent with the residuals, or a counted-zero pair is
            // numerically resolvable from 0 (near section), or the nearest nonzero eigenvalue is within 10x of its cut
            near_zero_ambiguous: !hc.cut_consistent || hc.near_sections || hc.margin_ratio.is_some_and(|m| m < 10.0),
            cut_margin_ratio: hc.margin_ratio.map(sig), lower_bound_only: sp.h0_lower_bound_only,
            wording: "certified residual: each eigenvalue lies within +-residual of a true eigenvalue (Weyl). The zero/nonzero labelling is a numerical policy (zero_rule), not a geometric kernel certificate; completeness of the Krylov search (no missed eigenvalue below the cut) is NOT certified; known-kernel dims are exact by construction." },
        h1: H1Report { method: Exact, dim: h1, rule: "rank-nullity on a graph: dim C1 - (dim C0 - h0); inherits the h0 cut" },
        spectral: SpectralReport { method: Exact, solver: sp.method, l0_smallest: small.iter().map(|&x| sig(x)).collect(),
            solver_pairs_residuals: res, gap_k: gap_k.map(sig), spectral_gap: sp.spectral_gap().map(sig), lambda_max_bound: sig(sp.lambda_max_bound),
            rounds: sp.rounds, converged: sp.converged,
            inner_cg: sp.inner.map(|(s_, i, r)| InnerCg { solves: s_, iters: i, worst_rel_residual: sig(r) }),
            l1_graph: L1Graph { kernel_dim: h1, nonzero_spectrum: "equals L0 nonzero spectrum (d0^T d0 vs d0 d0^T)", mapped_pairs_checked: l1_checked, max_mapped_residual: sig(l1_res) } },
        hodge: HodgeReport { method: Exact, cochain: source, norm2: sig(wn2), exact_share: sig(hs.exact), harmonic_share: sig(hs.harmonic),
            cg_iters: hs.iters, cert_coclosed: sig(hs.cert), converged: hs.converged, null_random_gaussian: null_random, null_shuffled_values: null_shuffle },
        energy_xlx: energy,
        edges_by_harmonic_energy: edges_h,
        edges_by_discord: edges_d,
        nodes_by_abs_potential: nodes_p,
        nodes_by_incident_discord: nodes_i,
        extension: spec.boundary.as_ref().map(|b| extension(spec, &d0, &on, b)),
        assumptions: vec!["A1-A8 of docs/rust-cohomology-design.md", "A9 w := delta x when only x is given (harmonic share 0 by construction)",
            "A10 v1 is graph-only (no 2-cells): no coexact part", "A11 potential = minimum-norm solution of L0 x = d0^T w (orthogonal to H0); it is in node-stalk units, so with raw-unit stalks and scaling maps prefer nodes_by_incident_discord"],
    };
    let contents = report.to_json();
    (report, vec![Effect::EmitReport { path: out_path, contents }])
}
