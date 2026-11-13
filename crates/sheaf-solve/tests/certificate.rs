//! Certificate soundness on arbitrary δ (not only sheaf coboundaries): the certificate fields are the documented
//! functions of the RETURNED x, recomputed here from a dense copy of δ, and `converged` is never claimed against a
//! scale that is not a finite number. Entries span 2^±20 (the ADR-017 P-NUM-05 regime), interiors may be singular.
use proptest::prelude::*;
use sheaf_solve::{harmonic_extension, Csr, Options};

fn nrm(v: &[f64]) -> f64 { v.iter().map(|a| a * a).sum::<f64>().sqrt() }
/// (δᵀδx restricted to `keep`), with δ dense row-major; same summation order as `Csr::mul` / `Csr::mul_t`.
fn l_keep(d: &[Vec<f64>], x: &[f64], keep: &[bool]) -> Vec<f64> {
    let y: Vec<f64> = d.iter().map(|r| r.iter().zip(x).map(|(a, b)| a * b).sum()).collect();
    (0..x.len()).filter(|&c| keep[c]).map(|c| { let mut s = 0.0; for (r, yr) in d.iter().zip(&y) { s += r[c] * yr; } s }).collect()
}
fn close(a: f64, b: f64) -> bool { (a - b).abs() <= 1e-12 * a.abs().max(b.abs()) }

/// a value m·2^e with m ∈ −4..=4 (0 allowed), e ∈ −20..=20
fn scaled() -> impl Strategy<Value = f64> { (-4i32..=4, -20i32..=20).prop_map(|(m, e)| m as f64 * 2f64.powi(e)) }

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]
    #[test]
    fn certificate_is_the_documented_function_of_x(
        (m, n, vals, mask, bvals, x0, use_x0) in (1usize..5, 1usize..6).prop_flat_map(|(m, n)| (
            Just(m), Just(n), prop::collection::vec(prop::option::weighted(0.6, scaled()), m * n), prop::collection::vec(any::<bool>(), n),
            prop::collection::vec(scaled(), n), prop::collection::vec(scaled(), n), any::<bool>())),
        cg in any::<bool>(),
    ) {
        let d: Vec<Vec<f64>> = (0..m).map(|r| (0..n).map(|c| vals[r * n + c].unwrap_or(0.0)).collect()).collect();
        let trip: Vec<(usize, usize, f64)> = (0..m * n).filter_map(|k| vals[k].map(|v| (k / n, k % n, v))).collect();
        let delta = Csr::from_triplets(m, n, trip);
        let bd: Vec<(usize, f64)> = (0..n).filter(|&c| mask[c]).map(|c| (c, bvals[c])).collect();
        let o = Options { dense_max: if cg { 0 } else { 400 }, max_iter: 200, ..Default::default() };
        let h = harmonic_extension(&delta, &bd, use_x0.then_some(&x0[..]), &o).unwrap();
        let (x, c) = (&h.x, &h.cert);

        // boundary values are returned bit for bit
        for &(i, v) in &bd { prop_assert_eq!(x[i].to_bits(), v.to_bits()); }
        let int: Vec<bool> = mask.iter().map(|b| !b).collect();
        let pick = |v: &[f64]| -> Vec<f64> { (0..n).filter(|&i| int[i]).map(|i| v[i]).collect() };
        // residual_abs = ‖(Lx)_I‖ and energy = ‖δx‖², both on the returned x
        let res = nrm(&l_keep(&d, x, &int));
        prop_assert!(close(c.residual_abs, res), "residual_abs {} vs recomputed {}", c.residual_abs, res);
        let e: f64 = d.iter().map(|r| { let y: f64 = r.iter().zip(x).map(|(a, b)| a * b).sum(); y * y }).sum();
        prop_assert!(close(c.energy, e), "energy {} vs ‖δx‖² {}", c.energy, e);
        // scale = max(‖L_IB x_B‖, ‖δ_I‖₁‖δ_I‖∞·max(‖x_I‖, ‖x_start,I‖), 1e-300)
        let xb: Vec<f64> = (0..n).map(|i| if mask[i] { bvals[i] } else { 0.0 }).collect();
        let rinf = d.iter().map(|r| pick(&r.iter().map(|a| a.abs()).collect::<Vec<_>>()).iter().sum::<f64>()).fold(0.0, f64::max);
        let c1 = pick(&(0..n).map(|j| d.iter().map(|r| r[j].abs()).sum::<f64>()).collect::<Vec<_>>()).into_iter().fold(0.0, f64::max);
        let start = if use_x0 { nrm(&pick(&x0)) } else { 0.0 };
        let scale = nrm(&l_keep(&d, &xb, &int)).max(rinf * c1 * nrm(&pick(x)).max(start)).max(1e-300);
        prop_assert!(scale.is_finite());
        if int.contains(&true) { prop_assert!(close(c.residual_rel, res / scale), "residual_rel {} vs {}/{}", c.residual_rel, res, scale); }
        else { prop_assert_eq!(c.residual_rel, 0.0); }
        // converged ⇔ residual within max(rtol, 1e-10) of that scale: a certified x satisfies the documented bound
        prop_assert_eq!(c.converged, c.residual_rel <= 1e-10, "{:?}", c);
        if c.converged { prop_assert!(res <= 1e-10 * scale * (1.0 + 1e-9), "converged but ‖(Lx)_I‖ = {} > 1e-10·{}", res, scale); }
    }
}
