//! Entry-magnitude gate of `v1::parse` (ADR-017 P-NUM-01), the refusal that
//! `adversarial_overflowing_entries_refused_not_panicking` in tests/v1.rs reads through `expect_err`.
use proptest::prelude::*;
use sheaf_cli::v1;

/// Small integers, values within 32 ulp of the cap on either side, or any normal f64; either sign.
fn arb_entry() -> impl Strategy<Value = f64> {
    let cap = v1::MAX_ABS_ENTRY.to_bits();
    prop_oneof![
        3 => (-4i32..=4).prop_map(f64::from),
        2 => (0u64..=64, any::<bool>()).prop_map(move |(k, neg)| { let m = f64::from_bits(cap - 32 + k); if neg { -m } else { m } }),
        1 => prop::num::f64::NORMAL,
    ]
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(300))]
    /// ∀ five finite entries placed in Fu, Fv, x[0], x[1], w of a two-node, one-edge spec (all stalks of dim 1),
    /// each written as a JSON number token and read back by serde_json as the reference value:
    /// if some reference value has |value| > MAX_ABS_ENTRY, `parse` is an error naming the cap;
    /// otherwise `parse` succeeds and stores each of the five reference values bit for bit.
    #[test]
    fn parse_admits_entries_iff_within_cap(entries in prop::array::uniform5(arb_entry())) {
        let tok = entries.map(|v| format!("{v:e}"));
        let want = tok.clone().map(|t| serde_json::from_str::<f64>(&t).unwrap());
        let text = format!(
            r#"{{"format":"frontier-sheaf/v1","nodes":[{{"id":"a","dim":1}},{{"id":"b","dim":1}}],"edges":[{{"u":0,"v":1,"Fu":[[{}]],"Fv":[[{}]]}}],"x":[[{}],[{}]],"w":[[{}]]}}"#,
            tok[0], tok[1], tok[2], tok[3], tok[4]);
        let got = v1::parse(&text);
        if want.iter().any(|v| v.abs() > v1::MAX_ABS_ENTRY) {
            let e = got.expect_err("an entry above the cap must be refused");
            prop_assert!(e.to_string().contains("exceeds 1e75"), "{}", e);
        } else {
            let s = got.unwrap();
            let read = [s.edges[0].fu[(0, 0)], s.edges[0].fv[(0, 0)], s.x.as_ref().unwrap()[0], s.x.as_ref().unwrap()[1], s.w.as_ref().unwrap()[0]];
            prop_assert_eq!(read.map(f64::to_bits), want.map(f64::to_bits));
        }
    }
}
