//! Well-formedness of the law ontology. These are the spec's own laws: if one fails, the spec is ill-formed
//! before any code exists. v0.2 adds the checks the three adversarial seats showed were missing.
use sheaf_laws::*;
use sheaf_laws::capabilities::{capabilities, Cover};
use sheaf_laws::wiring::{self, Source, TERMINAL};
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

fn by_id() -> HashMap<&'static str, Law> { laws().into_iter().map(|l| (l.id, l)).collect() }
fn repo(p: &str) -> PathBuf { Path::new(env!("CARGO_MANIFEST_DIR")).join("..").join(p) }
fn read(p: &str) -> String { std::fs::read_to_string(repo(p)).unwrap_or_else(|e| panic!("{p}: {e}")) }

#[test]
fn ids_are_unique() {
    let mut seen = HashSet::new();
    for l in laws() { assert!(seen.insert(l.id), "duplicate law id {}", l.id); }
}

#[test]
fn depth_is_capped_and_parents_are_strictly_shallower() {
    let all = laws(); let m = by_id();
    for l in &all {
        let d = depth(&all, l.id);
        assert!(d <= 3, "{} has depth {d} > 3", l.id);
        assert!(d <= l.layer as usize, "{} depth {d} exceeds its layer {:?}", l.id, l.layer);
        match (l.layer, l.parent) {
            (Layer::Global, None) => {}
            (Layer::Global, Some(p)) => panic!("global law {} has parent {p}", l.id),
            (_, None) => panic!("non-global law {} has no parent (orphan)", l.id),
            (_, Some(p)) => {
                let pl = m.get(p).unwrap_or_else(|| panic!("{} parent {p} does not exist", l.id));
                assert!(pl.layer < l.layer, "{} must refine a shallower law, not {p}", l.id);
            }
        }
    }
}

#[test]
fn every_global_law_has_a_derivation_and_derivations_resolve() {
    let m = by_id();
    for l in laws() {
        for d in l.derives_from { assert!(m.contains_key(d), "{} derives from unknown law {d}", l.id); }
        if l.layer == Layer::Global { assert!(!l.derives_from.is_empty(), "global law {} is stated but not derived", l.id); }
    }
}

#[test]
fn derivation_graph_is_acyclic() {
    let m = by_id();
    fn visit(id: &'static str, m: &HashMap<&'static str, Law>, stack: &mut Vec<&'static str>, done: &mut HashSet<&'static str>) {
        if done.contains(id) { return; }
        assert!(!stack.contains(&id), "derivation cycle through {id}: {stack:?}");
        stack.push(id);
        for d in m[id].derives_from { visit(d, m, stack, done); }
        stack.pop(); done.insert(id);
    }
    let mut done = HashSet::new();
    for l in laws() { visit(l.id, &m, &mut vec![], &mut done); }
}

/// Holes propagate: a law that derives from a Hole must itself be declared a Hole, so no open question can sit
/// under a green parent. (v0.1 had a test that only forbade deriving through a hole; Status::Hole was never
/// used, so that test was vacuous: seat 3 §4.)
#[test]
fn holes_propagate_upward() {
    let all = laws();
    for l in &all {
        if effective_status(&all, l.id) == Status::Hole {
            assert_eq!(l.status, Status::Hole, "{} derives from a Hole but is declared Proposed", l.id);
        }
    }
    assert!(all.iter().any(|l| l.status == Status::Hole), "no Hole at all: either the spec is complete or the check is vacuous");
}

#[test]
fn root_wiring_is_well_typed_and_every_slot_is_fed_once() {
    let w = wiring::root();
    for op in Op::ALL {
        let (ins, _) = op.signature();
        for (slot, color) in ins.iter().enumerate() {
            let feeds: Vec<_> = w.iter().filter(|(c, s, _)| *c == op && *s == slot).collect();
            assert_eq!(feeds.len(), 1, "{op:?} slot {slot} ({color:?}) fed {} times", feeds.len());
            let got = match &feeds[0].2 { Source::From(p) => p.signature().1, Source::Ext(c) => *c };
            assert_eq!(got, *color, "{op:?} slot {slot} expects {color:?}, wired {got:?}");
        }
    }
}

/// COMPOSE witness: every stage output reaches a consumer, except the product outputs. (v0.1: Extend dangled.)
#[test]
fn no_dangling_outputs() {
    for op in Op::ALL {
        if TERMINAL.contains(&op) { continue; }
        let consumed = wiring::root().iter().any(|(_, _, s)| matches!(s, Source::From(p) if *p == op));
        assert!(consumed, "{op:?}'s output is consumed by nothing (dangling)");
    }
}

#[test]
fn every_wired_edge_has_a_glue_law() {
    let all = laws();
    let missing: Vec<_> = wiring::edges().into_iter()
        .filter(|(a, b)| !all.iter().any(|l| l.scope == Scope::Edge(*a, *b)))
        .map(|(a, b)| format!("{a:?} -> {b:?}")).collect();
    assert!(missing.is_empty(), "wired edges with no glue law: {missing:?}");
}

#[test]
fn edge_laws_sit_on_real_interfaces() {
    for l in laws() {
        if let Scope::Edge(a, b) = l.scope {
            assert!(wiring::adjacent(a, b) || wiring::overlap(a, b), "{} glues {a:?} and {b:?}, which neither connect nor overlap", l.id);
        }
    }
}

#[test]
fn every_stage_has_a_law() {
    let all = laws();
    for op in Op::ALL {
        let has = all.iter().any(|l| matches!(l.scope, Scope::Stage(o) if o == op));
        assert!(has, "stage {op:?} has no stage law");
    }
}

#[test]
fn tiers_are_honest_about_ownership() {
    for l in laws() {
        if let Scope::Stage(op) | Scope::Product(_, op) = l.scope {
            if op.owner() == Owner::Core { assert!(l.tier != Tier::Formal, "{} claims Formal on core-owned {op:?}", l.id); }
        }
        if l.derives_from.is_empty() {
            if l.tier == Tier::Formal { assert!(!l.evidence.is_empty(), "{} is Formal but names no Bend law", l.id); }
            if l.tier == Tier::CoreCited { assert!(!l.evidence.is_empty(), "{} is CoreCited but names no core test", l.id); }
        }
    }
}

/// Weakest-link valuation: a composed law is no stronger than its weakest part.
#[test]
fn composed_laws_do_not_overclaim_their_tier() {
    let m = by_id();
    fn eff(m: &HashMap<&'static str, Law>, id: &str) -> Tier {
        let l = &m[id];
        if l.derives_from.is_empty() { return l.tier; }
        l.derives_from.iter().map(|d| eff(m, d)).max().unwrap().max(l.tier)
    }
    let bad: Vec<_> = laws().into_iter().filter(|l| !l.derives_from.is_empty())
        .filter_map(|l| { let w = l.derives_from.iter().map(|d| eff(&m, d)).max().unwrap(); (l.tier < w).then(|| format!("{} declares {:?}, weakest part {:?}", l.id, l.tier, w)) })
        .collect();
    assert!(bad.is_empty(), "overclaimed tiers:\n{}", bad.join("\n"));
}

#[test]
fn core_cited_tests_exist_in_sheaf_core() {
    let mut src = String::new();
    fn walk(p: &Path, out: &mut String) {
        for e in std::fs::read_dir(p).unwrap().flatten() {
            let p = e.path();
            if p.is_dir() { walk(&p, out) } else if p.extension().is_some_and(|x| x == "rs") { out.push_str(&std::fs::read_to_string(&p).unwrap()); }
        }
    }
    walk(&repo("../crates"), &mut src);
    for l in laws().iter().filter(|l| l.tier == Tier::CoreCited) {
        for t in l.evidence { assert!(src.contains(&format!("fn {t}(")), "{} cites core test {t}, which does not exist in crates/", l.id); }
    }
}

fn stated(file: &str, prefix: &str) -> Vec<String> {
    read(file).lines().filter_map(|l| l.strip_prefix(prefix)).map(|l| l.split([':', '(', '"']).next().unwrap().trim().to_string()).filter(|s| !s.is_empty()).collect()
}

/// Same ids at every layer: each law stated in LAWS.bend has a def in PROOF.bend, a Haskell property and a Rust
/// test of the same name; no Bend law is an orphan of the ontology. (seat 2 S2-16, seat 1 S1-10)
#[test]
fn formal_ids_match_across_bend_haskell_rust() {
    let bend = stated("bend/LAWS.bend", "law ");
    let proof = read("bend/PROOF.bend"); let hs = read("haskell/Laws.hs"); let rs = read("kit/tests/laws.rs");
    let registered: Vec<&str> = laws().iter().flat_map(|l| l.evidence.iter().copied()).collect();
    for b in &bend {
        assert!(registered.contains(&b.as_str()), "LAWS.bend states `{b}`, which no ontology law cites");
        assert!(proof.contains(&format!("def Laws.{b}(")), "LAWS.bend `{b}` has no def in PROOF.bend");
        assert!(hs.contains(&format!("(\"{b}\"")), "LAWS.bend `{b}` has no Haskell property of the same id");
        assert!(rs.contains(&format!("fn {b}(")), "LAWS.bend `{b}` has no Rust test of the same id");
    }
    let all = laws();
    let leaves: Vec<_> = all.iter().filter(|l| l.tier == Tier::Formal && l.derives_from.is_empty()).collect();
    let backed = leaves.iter().filter(|l| { let (b, n) = formal_backing(l, &read("bend/LAWS.bend"), &proof, &hs, &rs); n > 0 && b == n }).count();
    println!("Formal leaves fully backed (Bend law + proof + Haskell + Rust): {backed} of {}", leaves.len());
}

#[test]
fn every_spec_capability_has_stage_and_law() {
    let m = by_id();
    let mut ids = HashSet::new();
    for c in capabilities() {
        assert!(ids.insert(c.id), "duplicate capability {}", c.id);
        match c.cover {
            Cover::Deferred(why) => assert!(why.len() > 20, "{} deferred without a real reason", c.id),
            Cover::Laws(ops, ls) => {
                assert!(!ops.is_empty() && !ls.is_empty(), "{} ({}) has no stage or no law: GAP", c.id, c.what);
                for l in ls { assert!(m.contains_key(l), "{} cites unknown law {l}", c.id); }
            }
        }
    }
}

#[test]
fn no_python_on_the_gate_path() {
    for dir in ["laws", "kit", "haskell", "bend", "gates"] {
        fn walk(p: &Path, out: &mut Vec<PathBuf>) {
            if let Ok(rd) = std::fs::read_dir(p) { for e in rd.flatten() { let p = e.path(); if p.ends_with("target") { continue; } if p.is_dir() { walk(&p, out) } else { out.push(p) } } }
        }
        let mut files = vec![]; walk(&repo(dir), &mut files);
        for f in &files {
            assert!(f.extension().map_or(true, |x| x != "py"), "Python file on the gate path: {}", f.display());
            if dir == "gates" && f.extension().is_some_and(|x| x == "sh") {
                let s = std::fs::read_to_string(f).unwrap();
                assert!(!s.lines().any(|l| !l.trim_start().starts_with('#') && l.contains("python")), "{} invokes python", f.display());
            }
        }
    }
}

#[test]
fn report_holes_and_counts() {
    let all = laws();
    let holes: Vec<_> = all.iter().filter(|l| l.status == Status::Hole).map(|l| l.id).collect();
    println!("{} laws; Holes: {:?}", all.len(), holes);
}
