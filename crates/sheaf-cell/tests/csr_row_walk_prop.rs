//! `Csr::mul_t` and `HodgeOp::diagonal` against the plain index loops they replace, on well-formed CSR
//! (`ptr.len() == nrows + 1`, as `Csr::from_triplets` builds it). The outcome is compared whole: a panic on one side
//! must be a panic on the other, and a returned vector must match bit for bit (so NaN, ±0 and ±inf count too).
use proptest::prelude::*;
use sheaf_cell::{num::sparse::Csr, ops::HodgeOp};
use std::panic::{catch_unwind, AssertUnwindSafe};

// The oracle is the index loop itself, kept verbatim.
#[allow(clippy::needless_range_loop)]
fn mul_t_ref(a: &Csr, y: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0; a.ncols];
    for i in 0..a.nrows { for p in a.ptr[i]..a.ptr[i + 1] { out[a.idx[p]] += a.val[p] * y[i]; } }
    out
}

#[allow(clippy::needless_range_loop)]
fn diagonal_ref(n: usize, down: Option<&Csr>, up: Option<&Csr>) -> Vec<f64> {
    let mut d = vec![0.0; n];
    if let Some(a) = down { for i in 0..a.nrows { for p in a.ptr[i]..a.ptr[i + 1] { d[i] += a.val[p] * a.val[p]; } } }
    if let Some(b) = up { for i in 0..b.nrows { for p in b.ptr[i]..b.ptr[i + 1] { d[b.idx[p]] += b.val[p] * b.val[p]; } } }
    d
}

fn any_f64() -> impl Strategy<Value = f64> {
    prop_oneof![any::<u64>().prop_map(f64::from_bits), (-3i8..=3).prop_map(f64::from)]
}

/// Up to 4 x 4, up to 8 entries (repeats allowed, they stay separate entries), any `f64` values.
fn arb_csr() -> impl Strategy<Value = Csr> {
    (0usize..5, 1usize..5).prop_flat_map(|(nrows, ncols)| {
        let len = if nrows == 0 { 0..=0 } else { 0..=8 };
        prop::collection::vec((0..nrows.max(1), 0..ncols, any_f64()), len)
            .prop_map(move |t| Csr::from_triplets(nrows, ncols, t))
    })
}

/// `Some(bits)` for a returned vector, `None` for a panic.
fn outcome(f: impl FnOnce() -> Vec<f64>) -> Option<Vec<u64>> {
    catch_unwind(AssertUnwindSafe(f)).ok().map(|v| v.iter().map(|x| x.to_bits()).collect())
}

proptest! {
    /// ∀ well-formed `a`, ∀ `y` of any length (shorter than, equal to or longer than `a.nrows`):
    /// `a.mul_t(y)` and the index loop have the same outcome.
    #[test]
    fn mul_t_matches_index_loop(a in arb_csr(), y in prop::collection::vec(any_f64(), 0..7)) {
        prop_assert_eq!(outcome(|| a.mul_t(&y)), outcome(|| mul_t_ref(&a, &y)));
    }

    /// ∀ `n < 7`, ∀ well-formed optional `down` and `up` (their sizes need not agree with `n`):
    /// `HodgeOp::general(n, down, up).diagonal()` and the index loops have the same outcome.
    #[test]
    fn diagonal_matches_index_loops(
        n in 0usize..7,
        down in prop::option::of(arb_csr()),
        up in prop::option::of(arb_csr()),
    ) {
        let op = HodgeOp::general(n, down.as_ref(), up.as_ref());
        prop_assert_eq!(outcome(|| op.diagonal()), outcome(|| diagonal_ref(n, down.as_ref(), up.as_ref())));
    }
}
