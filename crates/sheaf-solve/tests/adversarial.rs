//! ADR-017 round 2 regression (P-NUM-05): an ill-scaled harmonic extension was reported converged when it was not a
//! solution (δ = [2^20, 2^-20], x_0 = 1 on the boundary: returned x_1 = 0, energy 2^40, converged = true).
use sheaf_solve::{harmonic_extension, Csr, Options};

#[test]
fn ill_scaled_extension_is_solved_or_not_certified() {
    let d = Csr { nrows: 1, ncols: 2, ptr: vec![0, 2], idx: vec![0, 1], val: vec![2f64.powi(20), 2f64.powi(-20)] };
    let h = harmonic_extension(&d, &[(0, 1.0)], None, &Options::default()).unwrap();
    // the true minimiser: 2^20·1 ... (δx) = 2^-20 x1 - 2^20·1 ... sign convention aside, energy 0 is attainable
    if h.cert.converged { assert!(h.cert.energy <= 1e-6, "converged=true but energy {}", h.cert.energy); assert!((h.x[1].abs() - 2f64.powi(40)).abs() <= 1e-6 * 2f64.powi(40)); }
    assert!(h.cert.converged, "dense route should solve this exactly");
    // CG route (dense_max = 0) must not claim convergence on a non-solution either
    let c = harmonic_extension(&d, &[(0, 1.0)], None, &Options { dense_max: 0, ..Default::default() }).unwrap();
    assert!(!c.cert.converged || c.cert.energy <= 1e-6, "CG: converged=true with energy {}", c.cert.energy);
}

/// An overflowed residual scale is not a certificate. δ = [2^-270, ·] with the second column empty, B = ∅, start
/// (2^540, 0): ‖x_start‖² = 2^1080 = +inf in f64, so residual/scale was 1/inf = 0 and the untouched start
/// (L x = (1, 0), exact relative residual 1, energy 2^540; the minimiser is x = 0) came back converged = true.
#[test]
fn overflowed_scale_is_not_certified() {
    let d = Csr { nrows: 1, ncols: 2, ptr: vec![0, 1], idx: vec![0], val: vec![2f64.powi(-270)] };
    for dense_max in [400, 0] {
        let h = harmonic_extension(&d, &[], Some(&[2f64.powi(540), 0.0]), &Options { dense_max, ..Default::default() }).unwrap();
        assert!(!h.cert.converged || h.cert.energy <= 1e-6, "converged=true with ‖(Lx)_I‖ = {} and energy {}", h.cert.residual_abs, h.cert.energy);
    }
}

/// An underflowed residual norm is not a certificate either. δ = [2^-600, 2^-600], x_0 = 2^600 on the boundary: the
/// returned x_1 = 0 has (Lx)_1 = 2^-600 and energy 1 (the minimiser x_1 = -2^600 has energy 0), but (2^-600)² is 0 in
/// f64, so ‖(Lx)_I‖ = sqrt(0) = 0 against the 1e-300 floor came back residual_rel = 0, converged = true (exact value: 1).
#[test]
fn underflowed_residual_norm_is_not_certified() {
    let d = Csr { nrows: 1, ncols: 2, ptr: vec![0, 2], idx: vec![0, 1], val: vec![2f64.powi(-600); 2] };
    for dense_max in [400, 0] {
        let h = harmonic_extension(&d, &[(0, 2f64.powi(600))], None, &Options { dense_max, ..Default::default() }).unwrap();
        assert!(h.cert.residual_abs >= 2f64.powi(-600) || h.x[1] != 0.0, "‖(Lx)_I‖ reported as {}", h.cert.residual_abs);
        assert!(!h.cert.converged || h.cert.energy <= 1e-6, "converged=true with x_1 = {} and energy {}", h.x[1], h.cert.energy);
    }
}
