//! sheaf-kit: the discrete verdict layer. Transcribed from products/bend/main.bend (UNCHECKED until the atom gate
//! records ALL PROOFS CHECK for it) and
//! products/haskell/Laws.hs (enforced there); tests here carry the same law ids (bridge B).
//!
//! Types-tier laws live in the signatures:
//! - I.quantize.classify: loop and discord evidence are different unit types (`Ticks<Loop>`, `Ticks<Disc>`),
//!   so a circulation can never be compared against a discord tolerance.
//! - S.quantize.only-crossing: nothing here accepts an f64.
//! - S.classify.total: `classify` is a total function into a closed enum.
use std::marker::PhantomData;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Loop;
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Disc;

/// A quantity measured in ticks of one domain tolerance. Only Quantize may construct it from a float.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Ticks<U>(u64, PhantomData<U>);
impl<U> Ticks<U> {
    pub const fn new(n: u64) -> Self { Ticks(n, PhantomData) }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Verdict { Consistent, FixableNodes, LoopConflict, Underdetermined, Inconclusive }

/// Touchstone receipt verdicts, in their order HOLDS < BROKEN < BLOCKED < HOLE.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Receipt { Holds, Broken, Blocked, Hole }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Cert { pub certified: bool }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Evidence {
    /// false iff a boundary was given and the extension is not unique
    pub bounded: bool,
    /// no independent edge observations (w := delta x)
    pub xonly: bool,
    pub circ: Ticks<Loop>,
    pub tol_c: Ticks<Loop>,
    pub disc: Ticks<Disc>,
    pub tol_d: Ticks<Disc>,
}

fn nodes(disc_over: bool) -> Verdict { if disc_over { Verdict::FixableNodes } else { Verdict::Consistent } }

fn loops(xonly: bool, circ_over: bool, disc_over: bool) -> Verdict {
    match (xonly, circ_over) {
        (true, _) => nodes(disc_over),
        (false, true) => Verdict::LoopConflict,
        (false, false) => nodes(disc_over),
    }
}

pub fn classify(cert: Cert, ev: Evidence) -> Verdict {
    if !cert.certified { return Verdict::Inconclusive; }
    if !ev.bounded { return Verdict::Underdetermined; }
    loops(ev.xonly, ev.circ >= ev.tol_c, ev.disc >= ev.tol_d)
}

/// `agree`: the core's Hodge split and the product's Circulate agree. Disagreement is a HOLE, never averaged.
pub fn to_receipt(agree: bool, v: Verdict) -> Receipt {
    if !agree { return Receipt::Hole; }
    match v {
        Verdict::Consistent => Receipt::Holds,
        Verdict::FixableNodes | Verdict::LoopConflict => Receipt::Broken,
        Verdict::Underdetermined | Verdict::Inconclusive => Receipt::Blocked,
    }
}

/// Same case shape as the Bend `join` (not `max`), so mutants line up across the three layers.
pub fn join(a: Receipt, b: Receipt) -> Receipt {
    use Receipt::*;
    match (a, b) {
        (Holds, b) => b,
        (Hole, _) => Hole,
        (Broken, Holds) | (Broken, Broken) => Broken,
        (Broken, Blocked) => Blocked,
        (Broken, Hole) => Hole,
        (Blocked, Hole) => Hole,
        (Blocked, _) => Blocked,
    }
}

/// One entry per declared component of the cover; `None` is a component that produced no receipt.
pub fn join_all(cover: &[Option<Receipt>]) -> Receipt {
    let mut acc = Receipt::Holds;
    for c in cover.iter().rev() {
        acc = match c { None => Receipt::Hole, Some(r) => join(*r, acc) };
    }
    acc
}

impl Verdict {
    pub fn show(self) -> &'static str {
        match self {
            Verdict::Consistent => "consistent", Verdict::FixableNodes => "fixable", Verdict::LoopConflict => "loop",
            Verdict::Underdetermined => "underdetermined", Verdict::Inconclusive => "inconclusive",
        }
    }
}
impl Receipt {
    pub fn show(self) -> &'static str {
        match self { Receipt::Holds => "HOLDS", Receipt::Broken => "BROKEN", Receipt::Blocked => "BLOCKED", Receipt::Hole => "HOLE" }
    }
}

/// The bridge table: the same ten rows as main.bend and Laws.hs, in the same order.
pub fn bridge_table() -> Vec<&'static str> {
    let ev = |bounded, xonly, c, tc, d, td| Evidence { bounded, xonly, circ: Ticks::new(c), tol_c: Ticks::new(tc), disc: Ticks::new(d), tol_d: Ticks::new(td) };
    let t = Cert { certified: true }; let f = Cert { certified: false };
    vec![
        classify(f, ev(true, false, 9, 1, 9, 1)).show(),
        classify(t, ev(false, false, 9, 1, 9, 1)).show(),
        classify(t, ev(true, true, 9, 1, 0, 1)).show(),
        classify(t, ev(true, true, 9, 1, 2, 1)).show(),
        classify(t, ev(true, false, 3, 3, 0, 1)).show(),
        classify(t, ev(true, false, 2, 3, 5, 5)).show(),
        classify(t, ev(true, false, 2, 3, 4, 5)).show(),
        to_receipt(false, Verdict::LoopConflict).show(),
        join(Receipt::Broken, Receipt::Blocked).show(),
        join_all(&[Some(Receipt::Holds), None, Some(Receipt::Holds)]).show(),
    ]
}

/// Bridge B, full domain: the same rows, in the same order, as `runghc Laws.hs fulltable`.
pub fn full_table() -> Vec<String> {
    let bools = [false, true];
    let b = |t: bool| if t { "1" } else { "0" };
    let vs = [Verdict::Consistent, Verdict::FixableNodes, Verdict::LoopConflict, Verdict::Underdetermined, Verdict::Inconclusive];
    let rs = [Receipt::Holds, Receipt::Broken, Receipt::Blocked, Receipt::Hole];
    let mut out = vec![];
    for c in bools { for bo in bools { for x in bools { for ci in 0..=6u64 { for tc in 0..=6u64 { for d in 0..=6u64 { for td in 0..=6u64 {
        let ev = Evidence { bounded: bo, xonly: x, circ: Ticks::new(ci), tol_c: Ticks::new(tc), disc: Ticks::new(d), tol_d: Ticks::new(td) };
        out.push(format!("{} {} {} {ci} {tc} {d} {td} {}", b(c), b(bo), b(x), classify(Cert { certified: c }, ev).show()));
    } } } } } } }
    for a in bools { for v in vs { out.push(format!("receipt {} {} {}", b(a), v.show(), to_receipt(a, v).show())); } }
    for p in rs { for q in rs { out.push(format!("join {} {} {}", p.show(), q.show(), join(p, q).show())); } }
    out
}
