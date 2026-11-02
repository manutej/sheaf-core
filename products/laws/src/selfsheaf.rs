//! The ontology as a sheaf, so sheaf-core can check the spec that specifies it.
//!
//! coverage(): nodes = ROOT + laws + stages + colors, constant R^1 sheaf, edges = refines / derives / scoped-to /
//!   reads-writes / global-to-root. Boundary = leaf laws pinned to their evidence level
//!   (1 executed here, 0.5 implemented but proof pending, 0 planned). H0 counts islands (expect 1);
//!   the harmonic extension gives every interior node an evidence-support score; unique = false means some part
//!   of the ontology rests on no pinned leaf.
//! contradiction(): our laws plus Spec v1's, joined by agrees (+1 = identity) / contradicts (Fu = -1) maps.
//!   Boundary = our global laws pinned adopted (+1). The extension forces a value on each v1 law: near +1 it
//!   survives, near -1 it must be retracted; the energy is the irreducible disagreement.
use crate::{laws, Law, Op, Color, Scope, Tier};

pub struct Graph { pub nodes: Vec<String>, pub edges: Vec<(String, String, f64, String)> }

impl Graph {
    fn node(&mut self, id: &str) { if !self.nodes.iter().any(|n| n == id) { self.nodes.push(id.to_string()); } }
    fn edge(&mut self, u: &str, v: &str, sign: f64, label: &str) { self.node(u); self.node(v); self.edges.push((u.into(), v.into(), sign, label.into())); }
    pub fn spec_json(&self, x: Option<&dyn Fn(&str) -> f64>) -> String {
        let mut s = String::from("{\"format\": \"frontier-sheaf/v1\",\n \"nodes\": [");
        s += &self.nodes.iter().map(|n| format!("{{\"id\": \"{n}\", \"dim\": 1}}")).collect::<Vec<_>>().join(", ");
        s += "],\n \"edges\": [\n";
        s += &self.edges.iter().map(|(u, v, sg, l)| format!("  {{\"u\": \"{u}\", \"v\": \"{v}\", \"Fu\": [[{sg}]], \"Fv\": [[1]], \"label\": \"{l}\"}}")).collect::<Vec<_>>().join(",\n");
        s += "\n ],\n";
        match x {
            Some(f) => s += &format!(" \"x\": [{}],\n", self.nodes.iter().map(|n| format!("[{}]", f(n))).collect::<Vec<_>>().join(", ")),
            None => s += &format!(" \"x\": [{}],\n", self.nodes.iter().map(|_| "[0]").collect::<Vec<_>>().join(", ")),
        }
        s += &format!(" \"options\": {{\"null_reps\": 0, \"top\": {}}}\n}}\n", self.edges.len().max(self.nodes.len()));
        s
    }
}

pub fn boundary_json(pins: &[(String, f64)]) -> String {
    format!("{{\"format\": \"frontier-sheaf/boundary-v1\", \"boundary\": [{}]}}\n",
        pins.iter().map(|(n, v)| format!("{{\"node\": \"{n}\", \"value\": [{v}]}}")).collect::<Vec<_>>().join(", "))
}

fn op_id(o: Op) -> String { format!("op:{o:?}") }
fn color_id(c: Color) -> String { format!("color:{c:?}") }

pub fn is_leaf(all: &[Law], l: &Law) -> bool { l.derives_from.is_empty() && !all.iter().any(|c| c.parent == Some(l.id)) }

/// Evidence level of a leaf, from what exists in the repo now (not from claims).
pub fn evidence_level(l: &Law, rust_tests: &str, laws_bend: &str) -> f64 {
    let in_rust = l.evidence.iter().any(|e| rust_tests.contains(&format!("fn {e}(")));
    let in_bend = l.evidence.iter().any(|e| laws_bend.contains(&format!("law {e}:")));
    match l.tier {
        Tier::CoreCited => 1.0,                                   // sheaf-core tests ran 70/70 on this box
        Tier::Formal if in_rust && in_bend => 0.5,                // Rust-checked; Bend proof not yet run
        Tier::Formal if in_bend => 0.25,
        _ if in_rust => 1.0,
        _ => 0.0,                                                 // planned
    }
}

pub fn coverage() -> Graph {
    let all = laws();
    let mut g = Graph { nodes: vec!["ROOT".into()], edges: vec![] };
    for o in Op::ALL { let (ins, out) = o.signature(); for c in ins { g.edge(&op_id(o), &color_id(*c), 1.0, "reads"); } g.edge(&op_id(o), &color_id(out), 1.0, "writes"); }
    for l in &all {
        g.node(l.id);
        if let Some(p) = l.parent { g.edge(l.id, p, 1.0, "refines"); }
        for d in l.derives_from { g.edge(l.id, d, 1.0, "derives"); }
        match l.scope {
            Scope::Root => g.edge(l.id, "ROOT", 1.0, "global"),
            Scope::Stage(o) | Scope::Product(_, o) => g.edge(l.id, &op_id(o), 1.0, "scoped"),
            Scope::Edge(a, b) => { g.edge(l.id, &op_id(a), 1.0, "glues"); g.edge(l.id, &op_id(b), 1.0, "glues"); }
        }
    }
    g
}

/// Spec v1's laws (from the "Sheaf Products — Complete Spec v1" doc, Bend layer section) and how they relate to ours.
pub fn contradiction() -> (Graph, Vec<(String, f64)>) {
    let all = laws();
    let mut g = Graph { nodes: vec![], edges: vec![] };
    // our global laws and the ones v1 touches
    for l in all.iter().filter(|l| l.layer == crate::Layer::Global) { g.node(l.id); }
    let rel: &[(&str, &str, f64)] = &[
        ("v1.L1-memory_encode_sound", "S.encode.valid", 1.0),
        ("v1.L2-identity_encode_sound", "S.encode.valid", 1.0),
        ("v1.L3-encode_preserves_edge_count", "S.encode.index", 1.0),
        ("v1.L4-consistent_when_energy_zero", "G5", 1.0),
        ("v1.L5-cyclic_if_harmonic_share_material", "G4", -1.0),
        ("v1.L6-local_if_exact_share_material", "G4", -1.0),
        ("v1.L7-mixed_if_both_shares_material", "G4", -1.0),
        ("v1.L8-flag_iff_status_differs", "P.memory.flag-iff-differs", 1.0),
        ("v1.L9-claim_flag_rule", "P.memory.claim-rule", 1.0),
        ("v1.L10-equal_status_never_flags", "P.memory.claim-rule", 1.0),
        ("v1.verdict-share-threshold-5pct", "A.classify.circ-only", -1.0),
        ("v1.verdict-share-threshold-5pct", "v1.L5-cyclic_if_harmonic_share_material", 1.0),
        ("v1.null-p95-reported-not-gated", "S.quantize.floor", 1.0),
        ("v1.python-tool-layer", "G2", -1.0),
    ];
    for (a, b, s) in rel { g.edge(a, b, *s, if *s > 0.0 { "agrees" } else { "contradicts" }); }
    // connect our laws to their globals so the pins propagate
    for (_a, b, _) in rel { if !b.starts_with("v1.") { let mut cur = *b; while let Some(p) = all.iter().find(|l| l.id == cur).and_then(|l| l.parent) { g.edge(cur, p, 1.0, "refines"); cur = p; } } }
    let pins = all.iter().filter(|l| l.layer == crate::Layer::Global).map(|l| (l.id.to_string(), 1.0)).collect();
    (g, pins)
}
