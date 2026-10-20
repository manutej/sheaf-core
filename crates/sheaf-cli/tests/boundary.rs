//! `--boundary`: harmonic extension (sheaf-solve) + uniqueness (sheaf-cohom::relative) in the v1 report.
use sheaf_cli::v1;

const PATH: &str = r#"{"format":"frontier-sheaf/v1","nodes":[{"id":"a","dim":1},{"id":"b","dim":1},{"id":"c","dim":1},{"id":"d","dim":1},{"id":"e","dim":1}],
 "edges":[{"u":"a","v":"b"},{"u":"b","v":"c"},{"u":"c","v":"d"},{"u":"d","v":"e"}],"x":[[0],[0],[0],[0],[0]],"options":{"null_reps":2}}"#;

#[test]
fn path_extension_is_linear_and_unique() {
    let mut s = v1::parse(PATH).unwrap();
    let (r0, _) = v1::analyze(&s, None);
    assert!(r0.to_value().get("extension").is_none(), "no --boundary => no extension section (report unchanged)");
    s.boundary = Some(v1::parse_boundary(r#"{"format":"frontier-sheaf/boundary-v1","boundary":[{"node":"a","value":[1]},{"node":4,"value":[5]}]}"#, &s).unwrap());
    let (r, _) = v1::analyze(&s, None);
    let e = r.extension.clone().unwrap();
    assert_eq!((e.h0_relative, e.unique, e.converged), (Some(0), Some(true), true));
    for (i, nv) in e.values.iter().enumerate() { assert!((nv.value[0] - (1.0 + i as f64)).abs() < 1e-10, "{i}: {:?}", nv); }
    assert!((e.energy - 4.0).abs() < 1e-10);
    let v = r.to_value();
    assert_eq!(v["extension"]["method"], "exact"); assert_eq!(v["extension"]["route"], "dense");
}

#[test]
fn unpinned_component_is_reported_not_unique_and_bad_boundary_rejected() {
    let spec = r#"{"format":"frontier-sheaf/v1","nodes":[{"id":"a","dim":2},{"id":"b","dim":2},{"id":"c","dim":2},{"id":"d","dim":2}],
     "edges":[{"u":"a","v":"b"},{"u":"c","v":"d"}],"x":[[0,0],[0,0],[0,0],[0,0]],"options":{"null_reps":2}}"#;
    let mut s = v1::parse(spec).unwrap();
    s.boundary = Some(v1::parse_boundary(r#"{"boundary":[{"node":"a","value":[1,2]}]}"#, &s).unwrap());
    let e = v1::analyze(&s, None).0.extension.unwrap();
    assert_eq!((e.h0_relative, e.unique), (Some(2), Some(false)));
    assert_eq!(e.values[1].value, vec![1.0, 2.0]);
    assert_eq!(e.values[2].value, vec![0.0, 0.0], "minimum-norm representative on the unpinned component");
    assert!(v1::parse_boundary(r#"{"boundary":[{"node":"a","value":[1]}]}"#, &s).is_err());
    assert!(v1::parse_boundary(r#"{"boundary":[{"node":"zz","value":[1,2]}]}"#, &s).is_err());
    assert!(v1::parse_boundary(r#"{"boundary":[{"node":"a","value":[1,2]},{"node":0,"value":[1,2]}]}"#, &s).is_err());
}
