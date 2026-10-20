//! ADR-017 adversarial panel: tests that kill surviving cargo-mutants mutants of the m1 admission code.
use sheaf_cell::m1::{check_edges, M1Edge, M1Error};

/// Mutant `check_edges -> Ok(())` survived: with it, the C ABI admitted out-of-range endpoints, wrong map shapes and
/// NaN maps (the last one silently, into δ). Each must be refused with the right typed error.
#[test]
fn check_edges_refuses_each_bad_edge() {
    let dv = [2usize, 1];
    let ok = [1.0, 2.0]; let one = [3.0]; let nan = [f64::NAN, 1.0];
    assert_eq!(check_edges(&dv, &[M1Edge { u: 0, v: 1, de: 1, fu: &ok, fv: &one }]), Ok(()));
    assert!(matches!(check_edges(&dv, &[M1Edge { u: 0, v: 2, de: 1, fu: &ok, fv: &one }]), Err(M1Error::Dim(_))));
    assert!(matches!(check_edges(&dv, &[M1Edge { u: 2, v: 0, de: 1, fu: &one, fv: &ok }]), Err(M1Error::Dim(_))));
    assert!(matches!(check_edges(&dv, &[M1Edge { u: 0, v: 1, de: 1, fu: &one, fv: &one }]), Err(M1Error::Dim(_))));
    assert!(matches!(check_edges(&dv, &[M1Edge { u: 0, v: 1, de: 1, fu: &ok, fv: &ok }]), Err(M1Error::Dim(_))));
    assert!(matches!(check_edges(&dv, &[M1Edge { u: 0, v: 1, de: 1, fu: &nan, fv: &one }]), Err(M1Error::Finite(_))));
    assert!(matches!(check_edges(&dv, &[M1Edge { u: 0, v: 1, de: 1, fu: &ok, fv: &[f64::INFINITY] }]), Err(M1Error::Finite(_))));
    // de = 2 onto a 2-dim stalk: the map must hold de*dim = 4 entries, not de/dim = 1 (kills `* -> /` on both sides)
    let four = [1.0, 0.0, 0.0, 1.0];
    assert_eq!(check_edges(&[2, 2], &[M1Edge { u: 0, v: 1, de: 2, fu: &four, fv: &four }]), Ok(()));
    assert!(matches!(check_edges(&[2, 2], &[M1Edge { u: 0, v: 1, de: 2, fu: &four, fv: &one }]), Err(M1Error::Dim(_))));
    assert!(matches!(check_edges(&[2, 2], &[M1Edge { u: 0, v: 1, de: 2, fu: &one, fv: &four }]), Err(M1Error::Dim(_))));
}

/// Mutant `Display for M1Error -> Ok(())` survived: callers match on the SPEC prefixes, so the text is contract.
#[test]
fn m1_error_messages_carry_spec_prefixes() {
    assert_eq!(M1Error::Dim("x".into()).to_string(), "R-DIM-02: x");
    assert_eq!(M1Error::Finite("y".into()).to_string(), "R-FINITE-01: y");
}
