//! The root composite: which stage output feeds which stage input slot. External inputs: the product's
//! Corpus, CorePolicy, DomainPolicy, an optional Boundary, the G0 Manifest, and (for Diff) a prior run's Finding.
use crate::{Color, Op};

pub enum Source { From(Op), Ext(Color) }

/// Stages whose output leaves the pipeline (product outputs). Every other stage's output must be consumed.
pub const TERMINAL: [Op; 4] = [Op::Bind, Op::Explain, Op::Diff, Op::Fill];

/// (consumer, slot, source)
pub fn root() -> Vec<(Op, usize, Source)> {
    use crate::Op::*; use Source::{From, Ext};
    vec![
        (Ingest, 0, Ext(Color::Corpus)),
        (Encode, 0, From(Ingest)),
        (Validate, 0, From(Encode)),
        (Analyze, 0, From(Validate)),
        (Analyze, 1, Ext(Color::CorePolicy)),
        (Analyze, 2, Ext(Color::Boundary)),
        (Certify, 0, From(Analyze)),
        (Circulate, 0, From(Validate)),
        (Quantize, 0, From(Analyze)),
        (Quantize, 1, From(Circulate)),
        (Quantize, 2, Ext(Color::DomainPolicy)),
        (Classify, 0, From(Certify)),
        (Classify, 1, From(Quantize)),
        (Attribute, 0, From(Analyze)),
        (Attribute, 1, From(Ingest)),
        (Attribute, 2, From(Classify)),
        (Attribute, 3, From(Circulate)),
        (Judge, 0, From(Attribute)),
        (Judge, 1, From(Certify)),
        (Join, 0, From(Judge)),
        (Join, 1, From(Judge)),
        (Bind, 0, From(Join)),
        (Bind, 1, Ext(Color::Manifest)),
        (Bind, 2, From(Attribute)),
        (Explain, 0, From(Attribute)),
        (Diff, 0, From(Attribute)),
        (Diff, 1, Ext(Color::Finding)),
        (Fill, 0, From(Attribute)),
        (Fill, 1, From(Ingest)),
    ]
}

pub fn adjacent(a: Op, b: Op) -> bool {
    root().iter().any(|(c, _, s)| *c == b && matches!(s, Source::From(p) if *p == a))
}

/// Distinct producer -> consumer pairs between stages.
pub fn edges() -> Vec<(Op, Op)> {
    let mut v: Vec<(Op, Op)> = vec![];
    for (c, _, s) in root() { if let Source::From(p) = s { if !v.contains(&(p, c)) { v.push((p, c)); } } }
    v
}

/// Two stages overlap when they read or write a common color (the sheaf "overlap" of a glue law).
pub fn overlap(a: Op, b: Op) -> bool {
    let (ia, oa) = a.signature(); let (ib, ob) = b.signature();
    let ca: Vec<Color> = ia.iter().copied().chain(std::iter::once(oa)).collect();
    ib.iter().copied().chain(std::iter::once(ob)).any(|c| ca.contains(&c))
}
