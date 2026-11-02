//! Every capability the product specs promise (SPEC.md v0.1 + Spec v1), traced to the stages that carry it and
//! the laws that govern it. Seeded from seat 3's trace (review/seat3-coverage.md §1). The test
//! `every_spec_capability_has_stage_and_law` fails on any row with no stage or no law unless it is Deferred
//! with a reason. Adding a capability to a spec means adding a row here.
use crate::Op::{self, *};

pub enum Cover { Laws(&'static [Op], &'static [&'static str]), Deferred(&'static str) }

pub struct Capability { pub id: &'static str, pub what: &'static str, pub spec: &'static str, pub cover: Cover }

macro_rules! cap {
    ($id:literal, $what:literal, $spec:literal, [$($op:ident),*], [$($l:literal),*]) => {
        Capability { id: $id, what: $what, spec: $spec, cover: Cover::Laws(&[$($op),*], &[$($l),*]) }
    };
    ($id:literal, $what:literal, $spec:literal, deferred $why:literal) => {
        Capability { id: $id, what: $what, spec: $spec, cover: Cover::Deferred($why) }
    };
}

pub fn capabilities() -> Vec<Capability> {
    vec![
    cap!("C01", "sheaf_analyze: report + verdict", "specv1 tool contract", [Validate, Analyze, Certify, Quantize, Classify], ["S.validate.mirror", "S.analyze.outcome", "S.certify.wf", "S.classify.total"]),
    cap!("C02", "sheaf_extend: fill with uniqueness", "specv1; SPEC §0", [Analyze, Certify, Attribute], ["S.analyze.extension-unique", "A.certify.unique-null", "A.attribute.fill"]),
    cap!("C03", "sheaf_explain / kb_explain / recon_explain", "specv1; SPEC P1, P2", [Explain], ["S.explain.pure", "S.explain.wording", "A.explain.loop-blame", "A.explain.resolve"]),
    cap!("C04", "kb_diff and --baseline (exit 2 on a new flag)", "SPEC P1; specv1 P1", [Diff], ["S.diff.keyed", "S.diff.partition", "S.diff.reflexive", "S.diff.exit"]),
    cap!("C05", "--selftest: rerun is byte-identical", "specv1 P1", [Analyze, Bind], ["G3", "S.analyze.thread-pin", "S.bind.binding"]),
    cap!("C06", "exit codes 3 / 2 / 0", "specv1 P1", [Ingest, Diff], ["S.ingest.total-or-reject", "S.diff.exit"]),
    cap!("C07", "kb_fill write-back (only when unique; ask otherwise)", "SPEC P1; specv1 P1 UI", [Fill], ["S.fill.determined-only", "S.fill.ask", "P.memory.writeback"]),
    cap!("C08", "CoreCrash as a typed finding", "SPEC §2, §3.3", [Analyze, Certify], ["S.analyze.outcome", "I.analyze.certify-error"]),
    cap!("C09", "P_INPUT_01 text returned unchanged", "specv1 tool contract", [Analyze], ["S.analyze.outcome"]),
    cap!("C10", "domain-v1 -> frontier-sheaf/v1 compiler", "SPEC §3.1", [Encode], ["S.encode.valid", "S.encode.index", "S.encode.mode"]),
    cap!("C11", "adapters: markdown wiki, JSON, Notion, Obsidian, CSV + recon.yaml", "SPEC P1, P2", [Ingest], ["S.ingest.total-or-reject", "S.ingest.ids-stable", "I.ingest.encode", "P.memory.status-table"]),
    cap!("C12", "schema / ontology-alignment adapter", "specv1 uc5; SPEC P5", [Ingest, Encode, Analyze], ["P.schema.holonomy-drop", "S.encode.valid"]),
    cap!("C13", "translation-loop adapter", "specv1 uc6; SPEC P5", [Encode, Circulate], ["P.translation.loop-is-harmonic"]),
    cap!("C14", "recon.yaml static or dated rates", "SPEC P2", [Encode], ["P.recon.rate-loop", "P.recon.dated-rate"]),
    cap!("C15", "batch K keys as components; shard if slow", "SPEC P2", [Join], ["P.recon.per-key-join", "S.join.cover", "I.judge.join"]),
    cap!("C16", "exceptions CSV / dated exception report", "SPEC P2", [Explain], ["S.explain.pure", "P.studio.traceable"]),
    cap!("C17", "scheduled nightly runs", "SPEC P2; specv1 v3", deferred "iteration of the whole composite plus Diff, run by a scheduled task outside the operad; no new law, C04 covers the comparison"),
    cap!("C18", "consensus_check(beliefs, frames?)", "SPEC P3", [Encode, Classify], ["P.consensus.pairwise-for-loops", "P.consensus.frame-holonomy"]),
    cap!("C19", "multi-agent-eval hook and per-round timeline", "SPEC P3; specv1 P3", [Diff], ["P.consensus.round-diff"]),
    cap!("C20", "outlier ranked #1 by incident discord", "SPEC E1a, E2b, E3a", [Attribute], ["A.attribute.ranking", "A.attribute.fixable", "P.consensus.outlier-rank"]),
    cap!("C21", "loop attribution with blame = loop", "SPEC §3.2", [Circulate, Attribute], ["I.circulate.attribute", "A.attribute.loop-blame"]),
    cap!("C22", "Studio verdict card, certificate chip, heat, math drawer", "SPEC P4", [Attribute], ["P.studio.traceable"]),
    cap!("C23", "Studio boundary mode", "SPEC P4 view 3", [Analyze], ["P.studio.pin-is-boundary", "S.analyze.extension-pins"]),
    cap!("C24", "Studio numbered vertical steps with hover highlight", "SPEC P4; specv1 P1 UI", [Explain], ["P.studio.steps"]),
    cap!("C25", "Studio what-if slider", "specv1 P4 v1", deferred "needs incremental recompute in the core (CR-4); out of the ontology until CR-4 lands"),
    cap!("C26", "same Finding over CLI, MCP and HTTP", "SPEC §5; specv1 versions", [Bind], ["G8", "S.surface.envelope"]),
    cap!("C27", "report bindings generated from report.rs", "SPEC §3.2; T3", [Certify], ["I.analyze.pin"]),
    cap!("C28", "Finding certificate hashes", "SPEC §3.2", [Bind], ["S.bind.binding"]),
    cap!("C29", "Finding context shown only as context", "SPEC §1 finding 1", [Attribute], ["A.attribute.context-only"]),
    cap!("C30", "nulls reported, never a gate", "specv1; SPEC §1", [Quantize, Attribute], ["I.analyze.quantize", "A.attribute.context-only"]),
    cap!("C31", "verdict rules in order, first match wins", "SPEC §3.3", [Classify], ["S.classify.total", "A.classify.consistent-iff"]),
    cap!("C32", "per-loop harmonic energy vs Circulate", "SPEC CR-5", [Circulate, Quantize], ["I.circulate.quantize"]),
    cap!("C33", "x-only flag; no 'no loops' badge for x-only", "SPEC §1 finding 4", [Classify, Attribute], ["A.classify.xonly", "A.attribute.x-only"]),
    cap!("C34", "tolerances in domain units; h0_absolute_tol to the core", "SPEC §3.1", [Analyze, Quantize], ["S.certify.policy", "S.quantize.floor"]),
    cap!("C35", "core options incl. top (complete rankings)", "core main.rs", [Analyze], ["S.analyze.top-complete", "S.bind.binding"]),
    cap!("C36", "h0_relative null above dense_max", "report.rs", [Certify], ["A.certify.unique-null"]),
    cap!("C37", "Spec v1 L1-L3 encoders sound", "specv1 Bend layer", [Encode], ["S.encode.valid", "S.encode.index"]),
    cap!("C38", "Spec v1 L4-L7 share verdicts (retracted)", "specv1 Bend layer", [Classify], ["G4", "A.classify.circ-only"]),
    cap!("C39", "Spec v1 L8-L10 memory flag rules", "specv1 Bend layer", [Encode], ["P.memory.flag-iff-differs", "P.memory.claim-rule"]),
    cap!("C40", "encoder differential across implementations", "specv1 differential", [Ingest, Encode], ["I.ingest.encode", "I.encode.validate"]),
    cap!("C41", "map conflicts (rate triangle, twisted frames, odd contradiction)", "seat1 S1-09", [Classify], ["A.classify.h0-policy", "P.recon.rate-triangle", "P.consensus.frame-holonomy", "P.memory.contradiction-parity"]),
    ]
}
