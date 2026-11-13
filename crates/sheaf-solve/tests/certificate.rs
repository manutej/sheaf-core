//! Certificate soundness on arbitrary δ (not only sheaf coboundaries): the certificate fields are the documented
//! functions of the RETURNED x, recomputed here from a dense copy of δ, and `converged` is never claimed against a
//! scale that is not a finite number. Entries span 2^±20 (the ADR-017 P-NUM-05 regime), interiors may be singular.
//! Third test: which x comes back when the interior is singular by construction (an unpinned dof, a duplicated column).
//! Last test: no absolute tolerance, a change of units (δ·2^k, data·2^j) changes the answer by exactly those units.
use proptest::prelude::*;
use sheaf_solve::{harmonic_extension, Csr, Options, Route};

fn nrm(v: &[f64]) -> f64 { v.iter().map(|a| a * a).sum::<f64>().sqrt() }
/// (δᵀδx restricted to `keep`), with δ dense row-major; same summation order as `Csr::mul` / `Csr::mul_t`.
fn l_keep(d: &[Vec<f64>], x: &[f64], keep: &[bool]) -> Vec<f64> {
    let y: Vec<f64> = d.iter().map(|r| r.iter().zip(x).map(|(a, b)| a * b).sum()).collect();
    (0..x.len()).filter(|&c| keep[c]).map(|c| { let mut s = 0.0; for (r, yr) in d.iter().zip(&y) { s += r[c] * yr; } s }).collect()
}
fn close(a: f64, b: f64) -> bool { (a - b).abs() <= 1e-12 * a.abs().max(b.abs()) }

/// a value m·2^e with m ∈ −4..=4 (0 allowed), e ∈ −20..=20
fn scaled() -> impl Strategy<Value = f64> { (-4i32..=4, -20i32..=20).prop_map(|(m, e)| m as f64 * 2f64.powi(e)) }
/// a nonzero value ±m·2^e with m ∈ 1..=4, e ∈ −600..=600: squares and products leave the f64 range
fn extreme() -> impl Strategy<Value = f64> { (1i32..=4, any::<bool>(), -600i32..=600).prop_map(|(m, neg, e)| if neg { -1.0 } else { 1.0 } * m as f64 * 2f64.powi(e)) }
/// (m, n, optional entries of an m × n δ, boundary mask, boundary values, start, use the start?)
type Problem = (usize, usize, Vec<Option<f64>>, Vec<bool>, Vec<f64>, Vec<f64>, bool);
fn problem<S: Strategy<Value = f64>>(v: fn() -> S) -> impl Strategy<Value = Problem> {
    (1usize..5, 1usize..6).prop_flat_map(move |(m, n)| (
        Just(m), Just(n), prop::collection::vec(prop::option::weighted(0.6, v()), m * n), prop::collection::vec(any::<bool>(), n),
        prop::collection::vec(v(), n), prop::collection::vec(v(), n), any::<bool>()))
}
/// the documented scale max(‖L_IB x_B‖, ‖δ_I‖₁‖δ_I‖∞·max(‖x_I‖, ‖x_start,I‖), 1e-300); `int` marks the interior columns
fn scale(d: &[Vec<f64>], int: &[bool], l_ib_xb: f64, xi: f64, start: f64) -> f64 {
    let rinf = d.iter().map(|r| r.iter().zip(int).filter(|p| *p.1).map(|p| p.0.abs()).sum::<f64>()).fold(0.0, f64::max);
    let c1 = (0..int.len()).filter(|&j| int[j]).map(|j| d.iter().map(|r| r[j].abs()).sum::<f64>()).fold(0.0, f64::max);
    l_ib_xb.max(rinf * c1 * xi.max(start)).max(1e-300)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]
    #[test]
    fn certificate_is_the_documented_function_of_x(
        (m, n, vals, mask, bvals, x0, use_x0) in problem(scaled),
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
        let start = if use_x0 { nrm(&pick(&x0)) } else { 0.0 };
        let scale = scale(&d, &int, nrm(&l_keep(&d, &xb, &int)), nrm(&pick(x)), start);
        prop_assert!(scale.is_finite());
        if int.contains(&true) { prop_assert!(close(c.residual_rel, res / scale), "residual_rel {} vs {}/{}", c.residual_rel, res, scale); }
        else { prop_assert_eq!(c.residual_rel, 0.0); }
        // converged ⇔ residual within max(rtol, 1e-10) of that scale: a certified x satisfies the documented bound
        prop_assert_eq!(c.converged, c.residual_rel <= 1e-10, "{:?}", c);
        if c.converged { prop_assert!(res <= 1e-10 * scale * (1.0 + 1e-9), "converged but ‖(Lx)_I‖ = {} > 1e-10·{}", res, scale); }
    }

    /// The same bound where f64 overflows: entries, boundary values and starts up to 4·2^±600 make ‖x‖, ‖L_IB x_B‖ or
    /// the L_II bound +inf. `converged` is then refused: it always means ‖(Lx)_I‖ ≤ 1e-10·scale with a FINITE scale.
    #[test]
    fn converged_means_a_finite_scale_bounds_the_residual((m, n, vals, mask, bvals, x0, use_x0) in problem(extreme), cg in any::<bool>()) {
        let d: Vec<Vec<f64>> = (0..m).map(|r| (0..n).map(|c| vals[r * n + c].unwrap_or(0.0)).collect()).collect();
        let delta = Csr::from_triplets(m, n, (0..m * n).filter_map(|k| vals[k].map(|v| (k / n, k % n, v))).collect());
        let bd: Vec<(usize, f64)> = (0..n).filter(|&c| mask[c]).map(|c| (c, bvals[c])).collect();
        let o = Options { dense_max: if cg { 0 } else { 400 }, max_iter: 200, ..Default::default() };
        let c = harmonic_extension(&delta, &bd, use_x0.then_some(&x0[..]), &o).unwrap();
        prop_assert_eq!(c.cert.converged, c.cert.residual_rel <= 1e-10, "{:?}", c.cert);
        if !c.cert.converged { return Ok(()); }
        let int: Vec<bool> = mask.iter().map(|b| !b).collect();
        // (δᵀδv)_I through the sparse products: a structural zero times an overflowed δv is skipped, not NaN
        let l_i = |v: &[f64]| -> Vec<f64> { delta.mul_t(&delta.mul(v)).into_iter().zip(&int).filter(|p| *p.1).map(|p| p.0).collect() };
        let pick = |v: &[f64]| -> Vec<f64> { (0..n).filter(|&i| int[i]).map(|i| v[i]).collect() };
        let xb: Vec<f64> = (0..n).map(|i| if mask[i] { bvals[i] } else { 0.0 }).collect();
        let (res, s) = (nrm(&l_i(&c.x)), scale(&d, &int, nrm(&l_i(&xb)), nrm(&pick(&c.x)), if use_x0 { nrm(&pick(&x0)) } else { 0.0 }));
        prop_assert!(s.is_finite() && res <= 1e-10 * s * (1.0 + 1e-9), "converged but ‖(Lx)_I‖ = {} against scale {} ({:?})", res, s, c.cert);
        // no entry of (Lx)_I exceeds the reported norm or the bound: |entry| is exact, a square below 2^-1074 is 0
        let top = l_i(&c.x).into_iter().fold(0.0f64, |a, v| a.max(v.abs()));
        prop_assert!(top <= c.cert.residual_abs && top <= 1e-10 * s * (1.0 + 1e-9), "converged but an entry of (Lx)_I is {:e} against scale {:e} ({:?})", top, s, c.cert);
    }

    /// A non-unique extension is the start's component in ker L_II plus the minimum-norm solution (module docs), here
    /// on the two kernel vectors that can be built by hand. δ gets two more interior columns: `t`, a copy of interior
    /// column `s` with the same start (e_s − e_t ∈ ker L_II), and `z`, with no entries (e_z ∈ ker L_II). Cholesky must
    /// refuse this L_II on the default route, x_z must come back untouched and x_s, x_t equal bit for bit.
    #[test]
    fn unpinned_dof_keeps_its_start_and_twin_columns_stay_equal(
        (m, n, vals, mut mask, bvals, mut x0, use_x0) in problem(scaled), cg in any::<bool>(), pick in 0usize..5,
    ) {
        let (s, t, z) = (pick % n, n, n + 1);
        mask[s] = false;
        let mut trip: Vec<(usize, usize, f64)> = (0..m * n).filter_map(|k| vals[k].map(|v| (k / n, k % n, v))).collect();
        let twin: Vec<_> = trip.iter().filter(|e| e.1 == s).map(|e| (e.0, t, e.2)).collect();
        trip.extend(twin);
        let delta = Csr::from_triplets(m, n + 2, trip);
        let bd: Vec<(usize, f64)> = (0..n).filter(|&c| mask[c]).map(|c| (c, bvals[c])).collect();
        x0.extend([x0[s], bvals[s]]);   // starts of t (= start of s) and of z (any value: s is interior, bvals[s] is free)
        let o = Options { dense_max: if cg { 0 } else { 400 }, max_iter: 200, ..Default::default() };
        let h = harmonic_extension(&delta, &bd, use_x0.then_some(&x0[..]), &o).unwrap();
        prop_assert_eq!(h.cert.route, Route::Cg, "a singular L_II was solved by Cholesky: {:?}", h.cert);
        prop_assert_eq!(h.x[z], if use_x0 { x0[z] } else { 0.0 }, "the unpinned dof moved ({:?})", h.cert);
        prop_assert_eq!(h.x[s].to_bits(), h.x[t].to_bits(), "twin columns: x_s = {}, x_t = {} ({:?})", h.x[s], h.x[t], h.cert);
    }

    /// No absolute tolerance anywhere: a change of units changes only the units. δ·2^k (edge units) and boundary values
    /// and start ·2^j (stalk units) are exact in f64 and commute with every rounded +, −, ×, ÷ and (even powers) sqrt,
    /// so both routes must take the same decisions and return x·2^j, the same residual_rel and verdict, bit for bit.
    #[test]
    fn a_change_of_units_changes_only_the_units(
        (m, n, vals, mask, bvals, x0, use_x0) in problem(scaled), cg in any::<bool>(), k in -30i32..=30, j in -30i32..=30,
    ) {
        let (c, a) = (2f64.powi(k), 2f64.powi(j));
        let delta = |s: f64| Csr::from_triplets(m, n, (0..m * n).filter_map(|i| vals[i].map(|v| (i / n, i % n, s * v))).collect());
        let bd = |s: f64| -> Vec<(usize, f64)> { (0..n).filter(|&i| mask[i]).map(|i| (i, s * bvals[i])).collect() };
        let x0a: Vec<f64> = x0.iter().map(|v| a * v).collect();
        let o = Options { dense_max: if cg { 0 } else { 400 }, max_iter: 200, ..Default::default() };
        let h = harmonic_extension(&delta(1.0), &bd(1.0), use_x0.then_some(&x0[..]), &o).unwrap();
        let g = harmonic_extension(&delta(c), &bd(a), use_x0.then_some(&x0a[..]), &o).unwrap();
        let (p, q) = (&h.cert, &g.cert);
        prop_assert_eq!((q.route, q.iters, q.breakdown, q.converged), (p.route, p.iters, p.breakdown, p.converged), "{:?} vs {:?}", q, p);
        for i in 0..n { prop_assert_eq!(g.x[i].to_bits(), (a * h.x[i]).to_bits(), "x[{}]: {:e} vs 2^{}·{:e}", i, g.x[i], j, h.x[i]); }
        prop_assert_eq!(q.residual_rel.to_bits(), p.residual_rel.to_bits(), "residual_rel {:e} vs {:e}", q.residual_rel, p.residual_rel);
        prop_assert_eq!(q.residual_abs.to_bits(), (c * c * a * p.residual_abs).to_bits(), "residual_abs {:e} vs {:e}", q.residual_abs, p.residual_abs);
        prop_assert_eq!(q.energy.to_bits(), (c * c * a * a * p.energy).to_bits(), "energy {:e} vs {:e}", q.energy, p.energy);
    }
}
