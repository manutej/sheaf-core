//! The root composite: which stage output feeds which stage input slot. External inputs are the
//! product's DomainGraph, an optional Boundary and the Policy.
use crate::{Color, Op};

pub enum Source { From(Op), Ext(Color) }

/// (consumer, slot, source)
pub fn root() -> Vec<(Op, usize, Source)> {
    use crate::Op::*; use Source::{From, Ext};
    vec![
        (Encode, 0, Ext(Color::DomainGraph)),
        (Validate, 0, From(Encode)),
        (Analyze, 0, From(Validate)),
        (Extend, 0, From(Validate)),
        (Extend, 1, Ext(Color::Boundary)),
        (Certify, 0, From(Analyze)),
        (Circulate, 0, From(Validate)),
        (Quantize, 0, From(Analyze)),
        (Quantize, 1, From(Circulate)),
        (Quantize, 2, Ext(Color::Policy)),
        (Classify, 0, From(Certify)),
        (Classify, 1, From(Quantize)),
        (Attribute, 0, From(Analyze)),
        (Attribute, 1, Ext(Color::DomainGraph)),
        (Attribute, 2, From(Classify)),
        (Receipt, 0, From(Attribute)),
        (Receipt, 1, From(Certify)),
        (Join, 0, From(Receipt)),
        (Join, 1, From(Receipt)),
    ]
}

pub fn adjacent(a: Op, b: Op) -> bool {
    root().iter().any(|(c, _, s)| *c == b && matches!(s, Source::From(p) if *p == a))
}

/// Two stages overlap when they read or write a common color (the sheaf "overlap" of a glue law).
pub fn overlap(a: Op, b: Op) -> bool {
    let (ia, oa) = a.signature(); let (ib, ob) = b.signature();
    let ca: Vec<Color> = ia.iter().copied().chain(std::iter::once(oa)).collect();
    ib.iter().copied().chain(std::iter::once(ob)).any(|c| ca.contains(&c))
}
