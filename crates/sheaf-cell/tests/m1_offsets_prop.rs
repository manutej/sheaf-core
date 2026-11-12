//! Prefix-sum law for `m1::offsets`: block `i` owns the index range `o[i]..o[i] + d[i]`.
use proptest::prelude::*;
use sheaf_cell::m1::offsets;

proptest! {
    /// ∀ dims `d` (length below 24, each entry below 1000, so no sum overflows):
    /// `o = offsets(d)` has `d.len() + 1` entries, starts at 0, and each step
    /// `o[i + 1] - o[i]` equals `d[i]`. Hence the last entry is the total dimension
    /// and the blocks tile `0..total` in order with no gap or overlap.
    #[test]
    fn m1_offsets_are_prefix_sums(d in prop::collection::vec(0usize..1000, 0..24)) {
        let o = offsets(&d);
        prop_assert_eq!(o.len(), d.len() + 1);
        prop_assert_eq!(o[0], 0);
        for (i, &w) in d.iter().enumerate() {
            prop_assert_eq!(o[i + 1] - o[i], w);
        }
        prop_assert_eq!(o[d.len()], d.iter().sum::<usize>());
    }
}
