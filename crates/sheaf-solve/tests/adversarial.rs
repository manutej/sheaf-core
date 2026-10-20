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
