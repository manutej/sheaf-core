//! Bridge B: the Rust implementation satisfies every law id of LAWS.bend / Laws.hs, exhaustively over the
//! finite types and over ticks 0..=6 (the Bend proofs cover all Nat; this checks the transcription).
use sheaf_kit::*;
use Receipt::*;

const R: [Receipt; 4] = [Holds, Broken, Blocked, Hole];
const V: [Verdict; 5] = [Verdict::Consistent, Verdict::FixableNodes, Verdict::LoopConflict, Verdict::Underdetermined, Verdict::Inconclusive];
const B: [bool; 2] = [false, true];
const N: std::ops::RangeInclusive<u64> = 0..=6;

fn ev(bounded: bool, xonly: bool, c: u64, tc: u64, d: u64, td: u64) -> Evidence {
    Evidence { bounded, xonly, circ: Ticks::new(c), tol_c: Ticks::new(tc), disc: Ticks::new(d), tol_d: Ticks::new(td) }
}
fn each6(mut f: impl FnMut(u64, u64, u64, u64)) { for a in N { for b in N { for c in N { for d in N { f(a, b, c, d) } } } } }
const T: Cert = Cert { certified: true };
const F: Cert = Cert { certified: false };

#[test] fn sound_uncertified() { for b in B { for x in B { each6(|c, tc, d, td| assert_eq!(classify(F, ev(b, x, c, tc, d, td)), Verdict::Inconclusive)) } } }
#[test] fn undetermined() { for x in B { each6(|c, tc, d, td| assert_eq!(classify(T, ev(false, x, c, tc, d, td)), Verdict::Underdetermined)) } }
#[test] fn xonly_no_loop() { each6(|c, tc, d, td| assert_ne!(classify(T, ev(true, true, c, tc, d, td)), Verdict::LoopConflict)) }
#[test] fn circ_only() {
    for x in B { for c2 in N { each6(|c1, tc, d, td| if (c1 >= tc) == (c2 >= tc) {
        assert_eq!(classify(T, ev(true, x, c1, tc, d, td)), classify(T, ev(true, x, c2, tc, d, td)))
    }) } }
}
#[test] fn mono_tolerance() {
    for k in N { each6(|c, tc, d, td| if classify(T, ev(true, false, c, tc, d, td)) == Verdict::Consistent {
        assert_eq!(classify(T, ev(true, false, c, tc + k, d, td)), Verdict::Consistent)
    }) }
}
#[test] fn receipt_map_hole() { for v in V { assert_eq!(to_receipt(false, v), Hole) } }
#[test] fn receipt_map_uncertified() { for b in B { for x in B { each6(|c, tc, d, td| assert_eq!(to_receipt(true, classify(F, ev(b, x, c, tc, d, td))), Blocked)) } } }
#[test] fn join_assoc() { for a in R { for b in R { for c in R { assert_eq!(join(a, join(b, c)), join(join(a, b), c)) } } } }
#[test] fn join_comm() { for a in R { for b in R { assert_eq!(join(a, b), join(b, a)) } } }
#[test] fn join_idem() { for a in R { assert_eq!(join(a, a), a) } }
#[test] fn join_unit_l() { for a in R { assert_eq!(join(Holds, a), a) } }
#[test] fn join_unit_r() { for a in R { assert_eq!(join(a, Holds), a) } }
#[test] fn join_hole_absorbs() { for a in R { assert_eq!(join(a, Hole), Hole) } }
#[test] fn join_is_max() { for a in R { for b in R { assert_eq!(join(a, b), a.max(b)) } } }
#[test] fn join_cover_missing_is_hole() {
    // every cover of length <= 4 over {Some(r), None} that contains a None joins to HOLE
    let alpha: Vec<Option<Receipt>> = R.iter().copied().map(Some).chain(std::iter::once(None)).collect();
    let mut covers: Vec<Vec<Option<Receipt>>> = vec![vec![]];
    for _ in 0..4 { let next: Vec<_> = covers.iter().flat_map(|c| alpha.iter().map(move |a| { let mut v = c.clone(); v.push(*a); v })).collect(); covers.extend(next); covers.dedup(); }
    for c in covers.iter().filter(|c| c.contains(&None)) { assert_eq!(join_all(c), Hole, "{c:?}") }
}
#[test] fn bridge_table_matches_golden() {
    let golden = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/bridge.golden")).unwrap();
    assert_eq!(bridge_table().join("\n") + "\n", golden);
}
#[test] fn loop_iff_circ() {
    for x in B { each6(|c, tc, d, td| assert_eq!(classify(T, ev(true, x, c, tc, d, td)) == Verdict::LoopConflict, !x && c >= tc)) }
}
#[test] fn fixable_iff_disc() {
    for x in B { each6(|c, tc, d, td| assert_eq!(classify(T, ev(true, x, c, tc, d, td)) == Verdict::FixableNodes, d >= td && !(!x && c >= tc))) }
}
#[test] fn consistent_iff() {
    for x in B { each6(|c, tc, d, td| assert_eq!(classify(T, ev(true, x, c, tc, d, td)) == Verdict::Consistent, !(!x && c >= tc) && !(d >= td))) }
}
#[test] fn holds_iff_consistent() { for v in V { assert_eq!(to_receipt(true, v) == Holds, v == Verdict::Consistent) } }
#[test] fn full_table_matches_golden() {
    let golden = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/bridge_full.golden")).unwrap();
    assert_eq!(full_table().join("\n") + "\n", golden);
}
