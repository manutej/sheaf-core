//! Well-formedness of the law ontology. These are the spec's own laws: if one fails, the spec is
//! ill-formed before any code exists.
use sheaf_laws::*;
use sheaf_laws::wiring::{self, Source};
use std::collections::{HashMap, HashSet};

fn by_id() -> HashMap<&'static str, Law> { laws().into_iter().map(|l| (l.id, l)).collect() }

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
                assert!(pl.layer < l.layer, "{} (layer {:?}) must refine a shallower law, not {p} ({:?})", l.id, l.layer, pl.layer);
            }
        }
    }
}

#[test]
fn every_global_law_has_a_derivation_and_derivations_resolve() {
    let m = by_id();
    for l in laws() {
        for d in l.derives_from { assert!(m.contains_key(d), "{} derives from unknown law {d}", l.id); }
        if l.layer == Layer::Global {
            assert!(!l.derives_from.is_empty(), "global law {} is stated but not derived: spec HOLE", l.id);
        }
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

#[test]
fn no_law_is_derived_through_a_hole() {
    let m = by_id();
    for l in laws() {
        for d in l.derives_from { assert!(m[d].status != Status::Hole, "{} derives through hole {d}", l.id); }
    }
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
        let has = all.iter().any(|l| matches!(l.scope, Scope::Stage(o) if o == op) || matches!(l.scope, Scope::Edge(a, b) if a == op || b == op));
        assert!(has, "stage {op:?} has no law");
    }
}

#[test]
fn tiers_are_honest_about_ownership() {
    for l in laws() {
        let op = match l.scope { Scope::Stage(o) | Scope::Product(_, o) => Some(o), _ => None };
        if let Some(op) = op {
            // the core's numerics are never claimed Formal on the products side
            if op.owner() == Owner::Core { assert!(l.tier != Tier::Formal, "{} claims Formal on core-owned stage {op:?}", l.id); }
        }
        // a leaf is discharged by its own evidence; a composed law by its derivation
        if l.derives_from.is_empty() {
            if l.tier == Tier::Formal { assert!(!l.evidence.is_empty(), "{} is Formal but names no Bend law", l.id); }
            if l.tier == Tier::CoreCited { assert!(!l.evidence.is_empty(), "{} is CoreCited but names no core test", l.id); }
        }
    }
}

/// Weakest-link valuation of the derivation tree: a composed law is no stronger than its weakest part.
#[test]
fn composed_laws_do_not_overclaim_their_tier() {
    let m = by_id();
    let mut bad = vec![];
    for l in laws() {
        if l.derives_from.is_empty() { continue; }
        let weakest = l.derives_from.iter().map(|d| effective_tier(&m, d)).max().unwrap();
        if l.tier < weakest { bad.push(format!("{} declares {:?} but its weakest part is {:?}", l.id, l.tier, weakest)); }
    }
    assert!(bad.is_empty(), "overclaimed tiers:\n{}", bad.join("\n"));
}

fn effective_tier(m: &HashMap<&'static str, Law>, id: &str) -> Tier {
    let l = &m[id];
    if l.derives_from.is_empty() { return l.tier; }
    l.derives_from.iter().map(|d| effective_tier(m, d)).max().unwrap().max(l.tier)
}

#[test]
fn core_cited_tests_exist_in_sheaf_core() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../crates");
    let mut src = String::new();
    fn walk(p: &std::path::Path, out: &mut String) {
        for e in std::fs::read_dir(p).unwrap().flatten() {
            let p = e.path();
            if p.is_dir() { walk(&p, out) } else if p.extension().is_some_and(|x| x == "rs") { out.push_str(&std::fs::read_to_string(&p).unwrap()); }
        }
    }
    walk(&root, &mut src);
    for l in laws() {
        if l.tier == Tier::CoreCited || l.tier == Tier::Differential || l.id.starts_with("G") {
            for t in l.evidence {
                if t.contains('/') { continue; } // a path, not a test name
                if l.tier == Tier::Formal { continue; }
                assert!(src.contains(&format!("fn {t}(")), "{} cites core test {t}, which does not exist in crates/", l.id);
            }
        }
    }
}

/// LAWS.bend may state only laws the ontology knows (no orphan laws); report how many Formal leaves are
/// stated in Bend yet versus planned.
#[test]
fn bend_laws_are_registered() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../bend/LAWS.bend");
    let text = std::fs::read_to_string(path).unwrap();
    let stated: Vec<&str> = text.lines().filter_map(|l| l.strip_prefix("law ")).map(|l| l.trim_end_matches(':')).collect();
    let registered: Vec<&str> = laws().iter().flat_map(|l| l.evidence.iter().copied()).collect();
    for s in &stated { assert!(registered.contains(s), "LAWS.bend states `{s}`, which no ontology law cites"); }
    let formal_leaves: Vec<_> = laws().into_iter().filter(|l| l.tier == Tier::Formal && l.derives_from.is_empty()).collect();
    let covered = formal_leaves.iter().filter(|l| l.evidence.iter().any(|e| stated.contains(e))).count();
    println!("Formal leaves stated in LAWS.bend: {covered} of {}", formal_leaves.len());
}
