//! v1 pipeline tests (moved from the cohomology spike tests/spectral.rs; bodies unchanged).
use faer::Side;
use sheaf_cohom::spectral::*;
use sheaf_cli::v1;

const TRI: &str = r#"{"format":"frontier-sheaf/v1","nodes":[{"id":"a","dim":1},{"id":"b","dim":1},{"id":"c","dim":1}],
 "edges":[{"u":"a","v":"b"},{"u":"b","v":"c"},{"u":"a","v":"c"}],"w":[[0],[0],[1]],"options":{"null_reps":8}}"#;

#[test]
fn v1_triangle_harmonic_share_is_one_third_and_report_is_deterministic() {
    let s = v1::parse(TRI).unwrap();
    let (r, e1) = v1::analyze(&s, None);
    let r = r.to_value(); // typed Report -> JSON view for field checks
    let (_, e2) = v1::analyze(&s, None);
    assert_eq!(e1, e2, "report must be byte-identical across runs");
    let h = r["hodge"]["harmonic_share"].as_f64().unwrap();
    assert!((h - 1.0 / 3.0).abs() < 1e-10, "{h}");
    assert_eq!(r["h0"]["dim"], 1); assert_eq!(r["h1"]["dim"], 1);
    assert_eq!(r["edges_by_harmonic_energy"].as_array().unwrap().len(), 3);
}

#[test]
fn v1_x_only_is_exact_and_rotation_loop_kills_h0() {
    // EN->ES->FR with +90 deg hops and a +90 closing map: holonomy 270 deg => no global section (H0 = 0)
    let spec = r#"{"format":"frontier-sheaf/v1","nodes":[{"id":"en","dim":2},{"id":"es","dim":2},{"id":"fr","dim":2}],
     "edges":[{"u":"en","v":"es","Fu":[[0,-1],[1,0]],"Fv":[[1,0],[0,1]]},{"u":"es","v":"fr","Fu":[[0,-1],[1,0]],"Fv":[[1,0],[0,1]]},
              {"u":"en","v":"fr","Fu":[[1,0],[0,1]],"Fv":[[0,-1],[1,0]]}],
     "x":[[1,0],[0,1],[-1,0]],"options":{"null_reps":4}}"#;
    let s = v1::parse(spec).unwrap();
    let (r, _) = v1::analyze(&s, None);
    let r = r.to_value(); // typed Report -> JSON view for field checks
    assert_eq!(r["h0"]["dim"], 0);
    assert!(r["hodge"]["harmonic_share"].as_f64().unwrap() < 1e-12);
    assert!(r["energy_xLx"].as_f64().unwrap() > 1.0);
}

#[test]
fn v1_rejects_bad_specs() {
    assert!(v1::parse("{").is_err());
    assert!(v1::parse(r#"{"format":"frontier-sheaf/v1","nodes":[{"id":"a","dim":1}],"edges":[],"w":null}"#).is_err());
    let bad_shape = r#"{"format":"frontier-sheaf/v1","nodes":[{"id":"a","dim":2},{"id":"b","dim":1}],"edges":[{"u":0,"v":1,"Fu":[[1,0]],"Fv":[[1,2]]}],"x":[[0,0],[1]]}"#;
    assert!(v1::parse(bad_shape).is_err());
    let omitted_mismatch = r#"{"format":"frontier-sheaf/v1","nodes":[{"id":"a","dim":2},{"id":"b","dim":1}],"edges":[{"u":0,"v":1}],"x":[[0,0],[1]]}"#;
    assert!(v1::parse(omitted_mismatch).is_err());
}

#[test]
fn v1_multi_component_h0_matches_dense_count() {
    // components: (a) 80-cycle with a 60deg twist (R^2), (b) random-map tree on 70 nodes (R^2; tree => H0 = 2),
    // (c) constant 90-cycle (R^2, known kernel), (d) one isolated node of dim 3. Expected H0 = 0 + 2 + 2 + 3 = 7.
    let r = |t: f64| vec![vec![t.cos(), -t.sin()], vec![t.sin(), t.cos()]];
    let mut nodes = vec![]; let mut edges = vec![];
    let mut st: u64 = 99; let mut rnd = || { st ^= st << 13; st ^= st >> 7; st ^= st << 17; ((st % 2001) as f64) / 1000.0 - 1.0 };
    for i in 0..80 { nodes.push(serde_json::json!({"id": format!("a{i}"), "dim": 2})); }
    for i in 0..80usize { let fv = if i == 0 { r(std::f64::consts::PI / 3.0) } else { vec![vec![1.0, 0.0], vec![0.0, 1.0]] }; edges.push(serde_json::json!({"u": i, "v": (i + 1) % 80, "Fv": fv})); }
    for i in 0..70 { nodes.push(serde_json::json!({"id": format!("b{i}"), "dim": 2})); }
    for i in 1..70usize { let m: Vec<Vec<f64>> = (0..2).map(|_| (0..2).map(|_| rnd()).collect()).collect(); let m2: Vec<Vec<f64>> = (0..2).map(|_| (0..2).map(|_| rnd()).collect()).collect(); edges.push(serde_json::json!({"u": 80 + i, "v": 80 + (i - 1) / 2, "Fu": m, "Fv": m2})); }
    for i in 0..90 { nodes.push(serde_json::json!({"id": format!("c{i}"), "dim": 2})); }
    for i in 0..90usize { edges.push(serde_json::json!({"u": 150 + i, "v": 150 + (i + 1) % 90})); }
    nodes.push(serde_json::json!({"id": "iso", "dim": 3}));
    let ne = edges.len();
    let spec = serde_json::json!({"format": "frontier-sheaf/v1", "nodes": nodes, "edges": edges, "w": vec![vec![0.5, -0.25]; ne], "options": {"null_reps": 2, "k": 4}});
    let mut s = v1::parse(&spec.to_string()).unwrap();
    s.opts.dense_max = 100; // all three non-trivial components (160, 140, 180 dofs) go through Lanczos
    let (rep, _) = v1::analyze(&s, None);
    let rep = rep.to_value();
    let (d0, _, _) = v1::coboundary(&s.dims, &s.edges);
    let op = HodgeOp::l0(&d0);
    let ev = op.dense().self_adjoint_eigenvalues(Side::Lower).unwrap();
    let tol = EIG_TOL * op.norm_bound().max(1.0);
    let dense_h0 = ev.iter().filter(|&&x| x <= tol).count();
    assert_eq!(dense_h0, 7);
    assert_eq!(rep["h0"]["dim"], 7, "{}", serde_json::to_string_pretty(&rep["h0"]).unwrap());
    let small: Vec<f64> = rep["spectral"]["L0_smallest"].as_array().unwrap().iter().map(|x| x.as_f64().unwrap()).collect();
    let nz: Vec<f64> = ev.iter().cloned().filter(|&x| x > tol).collect();
    // L0_smallest = k+1 smallest overall (here all zeros since H0 = 7 > k+1); the first nonzero must match the spectral gap
    assert!(small.iter().all(|&x| x <= tol));
    let gap = rep["spectral"]["spectral_gap"].as_f64().unwrap();
    assert!((gap - nz[0]).abs() <= 1e-8 * ev.last().unwrap(), "gap {gap} vs dense {}", nz[0]);
}

#[test]
fn v1_two_scale_components_use_local_zero_cut() {
    // Grok v1 #1 repro: a-b at scale 1e5, c-d at scale 1e-1. True H0 = 2; the global sheafcore rule says 3.
    let spec = r#"{"format":"frontier-sheaf/v1","nodes":[{"id":"a","dim":1},{"id":"b","dim":1},{"id":"c","dim":1},{"id":"d","dim":1}],
      "edges":[{"u":"a","v":"b","Fu":[[1e5]],"Fv":[[1e5]]},{"u":"c","v":"d","Fu":[[0.1]],"Fv":[[0.1]]}],"x":[[1],[2],[3],[4]],"options":{"null_reps":2}}"#;
    let (r, _) = v1::analyze(&v1::parse(spec).unwrap(), None);
    // ADR-015: one rule (per-component); the old global R-H0-08 recount (which said 3) is gone, so no ambiguity
    assert_eq!(r.h0.dim, 2);
    assert_eq!(r.h0.zero_policy, "per_component_eig");
    assert!(!r.h0.near_zero_ambiguous && r.h0.cut_consistent_with_zero_tol);
    assert!(r.h0.zero_tol > r.h0.zero_tol_min, "two components, two cuts");
    // a context-chosen absolute permitted delta is the same rule family, echoed in the report: 0.03 > lambda(c-d) = 0.02
    let abs = spec.replace(r#""options":{"null_reps":2}"#, r#""options":{"null_reps":2,"h0_absolute_tol":0.03}"#);
    let (r2, _) = v1::analyze(&v1::parse(&abs).unwrap(), None);
    assert_eq!(r2.h0.zero_policy, "absolute"); assert_eq!(r2.h0.dim, 3);
    assert!(v1::parse(&spec.replace(r#""options":{"null_reps":2}"#, r#""options":{"h0_absolute_tol":-1}"#)).is_err());
}

#[test]
fn v1_options_method_parsed_and_unknown_rejected() {
    let base = |m: &str| format!(r#"{{"format":"frontier-sheaf/v1","nodes":[{{"id":"a","dim":1}},{{"id":"b","dim":1}}],"edges":[{{"u":0,"v":1}}],"w":[[1]],"options":{{"method":"{m}"}}}}"#);
    assert_eq!(v1::parse(&base("polynomial")).unwrap().opts.method, Method::Polynomial);
    assert!(v1::parse(&base("magic")).is_err());
}

/// ADR-017 P-RES-01: a ~120-byte spec must be refused before any allocation (it used to try a 320 GB identity map and
/// abort), and a stalk dim near 2^64 must not wrap the offsets (it panicked with index out of bounds).
#[test]
fn adversarial_tiny_spec_huge_dims_refused_before_alloc() {
    let huge = r#"{"format":"frontier-sheaf/v1","nodes":[{"id":"a","dim":200000},{"id":"b","dim":200000}],"edges":[{"u":0,"v":1}],"w":[]}"#;
    let e = sheaf_cli::v1::parse(huge).err().expect("must be refused");
    assert!(e.to_string().contains("restriction maps exceed"), "{e}");
    let wrap = r#"{"format":"frontier-sheaf/v1","nodes":[{"id":"a","dim":18446744073709551615},{"id":"b","dim":1}],"edges":[],"w":[]}"#;
    let e = sheaf_cli::v1::parse(wrap).err().expect("must be refused");
    assert!(e.to_string().contains("total stalk dimension"), "{e}");
    let de = r#"{"format":"frontier-sheaf/v1","nodes":[{"id":"a","dim":1},{"id":"b","dim":1}],"edges":[{"u":0,"v":1,"de":4294967296}],"w":[[0]]}"#;
    assert!(sheaf_cli::v1::parse(de).is_err());
}

/// ADR-017 P-NUM-01: finite-but-huge entries used to crash `analyze` (evd NoConvergence panic / NaN sort panic);
/// they are now refused at parse with a P_INPUT_01 error.
#[test]
fn adversarial_overflowing_entries_refused_not_panicking() {
    for t in [r#"{"format":"frontier-sheaf/v1","nodes":[{"id":"a","dim":1},{"id":"b","dim":1},{"id":"c","dim":1}],"edges":[{"u":0,"v":1,"Fu":[[1e300]],"Fv":[[1e300]]},{"u":1,"v":2}],"x":[[1e300],[-1e300],[0]]}"#,
              r#"{"format":"frontier-sheaf/v1","nodes":[{"id":"a","dim":1},{"id":"b","dim":1}],"edges":[{"u":0,"v":1,"Fu":[[1e200]],"Fv":[[1e200]]}],"w":[[1]]}"#] {
        let e = sheaf_cli::v1::parse(t).err().expect("must be refused");
        assert!(e.to_string().contains("exceeds 1e75"), "{e}");
    }
    // the largest admitted magnitude runs to completion without panicking
    let ok = r#"{"format":"frontier-sheaf/v1","nodes":[{"id":"a","dim":1},{"id":"b","dim":1},{"id":"c","dim":1}],"edges":[{"u":0,"v":1,"Fu":[[1e75]],"Fv":[[-1e75]]},{"u":1,"v":2,"Fu":[[1e75]],"Fv":[[1e75]]}],"x":[[1e75],[-1e75],[1e75]]}"#;
    let s = sheaf_cli::v1::parse(ok).unwrap();
    let _ = sheaf_cli::v1::analyze(&s, None);
}

#[test]
fn v1_truncated_lists_report_totals_sheaf_glue_graph_001() {
    // Bend sheaf-glue-graph-001 :267 counterexample: 3 conflicting edges, top 2 -> 2 listed, 1 hidden with no count.
    let spec = r#"{"format":"frontier-sheaf/v1","nodes":[{"id":"a","dim":1},{"id":"b","dim":1},{"id":"c","dim":1}],
     "edges":[{"u":"a","v":"b"},{"u":"b","v":"c"},{"u":"a","v":"c"}],"w":[[1],[2],[3]],"options":{"null_reps":2,"top":2}}"#;
    let (r, _) = v1::analyze(&v1::parse(spec).unwrap(), None);
    assert_eq!(r.edges_by_discord.len(), 2);
    assert_eq!(r.edges_by_discord_total, 3, "hidden conflicting edge must be counted");
    assert_eq!(r.edges_by_harmonic_energy.len(), r.edges_by_harmonic_energy_total.min(2));
    assert_eq!(r.nodes_by_incident_discord.len(), 2);
    assert_eq!(r.nodes_by_incident_discord_total, 3);
    assert_eq!(r.nodes_by_abs_potential.len(), r.nodes_by_abs_potential_total.min(2));
    let v = r.to_value();
    assert_eq!(v["format"], "frontier-sheaf/report-v1");
    assert_eq!(v["edges_by_discord_total"], 3);
}
