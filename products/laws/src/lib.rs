//! The law ontology of sheaf-products, as typed data (v0.2, after three independent adversarial seats).
//!
//! Shape (operadic): the product is ONE composite of many-in/one-out stages (`Op`), wired by typed edges
//! (`Color`, see wiring.rs). Laws live at three depths below the root:
//!   L0 Global   — about the root composite (the whole pipeline, any product)
//!   L1 Stage    — about one stage (local) or one wired edge between stages (interface / glue)
//!   L2 Atom     — about one function inside a stage
//!   L3 Product  — a product's refinement of a stage or atom law (leaf; depth cap 3)
//! A global law must be DERIVABLE from lower laws (`derives_from`). A law derived from a `Hole` is itself a
//! `Hole` (holes propagate up), so an open question can never hide under a green parent.
//! Tiers follow the Touchstone oracle ladder; "proved" is only ever said of `Formal` with a checked proof.
//! Review record: products/laws/review/ (seat1-math, seat2-touchstone, seat3-coverage, meta-operad, CONSOLIDATED).

pub mod wiring;
pub mod selfsheaf;
pub mod capabilities;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Layer { Global = 0, Stage = 1, Atom = 2, Product = 3 }

/// Typed edges of the pipeline (what flows between stages).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Color {
    Corpus,       // markdown notes / JSON / CSV + recon.yaml / Notion export
    DomainGraph,  // sheaf-kit/domain-v1 (entities, relations, values, anchors) with stable ids
    Spec,         // frontier-sheaf/v1 bytes
    ValidSpec,    // spec the core parser accepted
    Boundary,     // frontier-sheaf/boundary-v1, possibly empty (empty = plain analysis)
    CorePolicy,   // k, null_reps, seed, method, dense_max, h0_absolute_tol, top, thread count
    DomainPolicy, // tolerances in domain units (loop, discord), expected H0, kappa
    Outcome,      // Report (frontier-sheaf/report-v1, f64, core-owned) | CoreError{exit, stderr}
    Cert,         // certificate projection of an Outcome (bools, residuals, policy echo, reason)
    Circ,         // per-fundamental-cycle circulations of w (products side; scalar/invertible maps only)
    Evidence,     // quantised evidence: comparison bits + saturated ticks; the only f64 -> Nat crossing
    Verdict,      // domain verdict: consistent | fixable | loop | map-conflict | underdetermined | inconclusive
    Finding,      // finding-v1: verdict + attribution (fixable, loops, fills, context)
    Explanation,  // plain-language text for agents and people
    Delta,        // run-over-run change set (new / resolved / unchanged / drifted)
    TsVerdict,    // Touchstone verdict HOLDS | BROKEN | BLOCKED | HOLE
    Manifest,     // G0 core manifest + laws.lock + thread count
    Receipt,      // TsVerdict bound to sha256 of manifest, spec, outcome, policies, finding
}

/// The generators of the operad: each stage is a many-in, one-out typed operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Op { Ingest, Encode, Validate, Analyze, Certify, Circulate, Quantize, Classify, Attribute, Judge, Join, Bind, Explain, Diff, Fill }

impl Op {
    pub fn signature(self) -> (&'static [Color], Color) {
        use Color::*;
        match self {
            Op::Ingest    => (&[Corpus], DomainGraph),
            Op::Encode    => (&[DomainGraph], Spec),
            Op::Validate  => (&[Spec], ValidSpec),
            Op::Analyze   => (&[ValidSpec, CorePolicy, Boundary], Outcome),
            Op::Certify   => (&[Outcome], Cert),
            Op::Circulate => (&[ValidSpec], Circ),
            Op::Quantize  => (&[Outcome, Circ, DomainPolicy], Evidence),
            Op::Classify  => (&[Cert, Evidence], Verdict),
            Op::Attribute => (&[Outcome, DomainGraph, Verdict, Circ], Finding),
            Op::Judge     => (&[Finding, Cert], TsVerdict),
            Op::Join      => (&[TsVerdict, TsVerdict], TsVerdict),
            Op::Bind      => (&[TsVerdict, Manifest, Finding], Receipt),
            Op::Explain   => (&[Finding], Explanation),
            Op::Diff      => (&[Finding, Finding], Delta),
            Op::Fill      => (&[Finding, DomainGraph], DomainGraph),
        }
    }
    pub fn owner(self) -> Owner {
        match self { Op::Validate | Op::Analyze => Owner::Core, _ => Owner::Products }
    }
    /// Effect signature (interface.lock).
    pub fn effect(self) -> &'static str {
        match self {
            Op::Analyze => "subprocess(sheaf analyze), total via Outcome, deterministic under pinned threads",
            Op::Validate => "subprocess parse, total via Outcome",
            _ => "pure, total, deterministic",
        }
    }
    pub const ALL: [Op; 15] = [Op::Ingest, Op::Encode, Op::Validate, Op::Analyze, Op::Certify, Op::Circulate, Op::Quantize,
        Op::Classify, Op::Attribute, Op::Judge, Op::Join, Op::Bind, Op::Explain, Op::Diff, Op::Fill];
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Owner { Core, Products }

/// Where a law lives in the tree.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Scope {
    Root,                         // the whole composite, or the verification process around it
    Stage(Op),                    // inside one stage
    Edge(Op, Op),                 // the interface from producer to consumer (glue law)
    Product(Prod, Op),            // a product's refinement of a stage
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Prod { Memory, Recon, Consensus, Studio, Schema, Translation }

/// Touchstone oracle tier: how the law is enforced. Strongest first.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Tier {
    Formal,       // Bend proof (+ --verdict) AND Haskell properties under the same id; Rust conforms (bridge B)
    Types,        // enforced by Rust types (newtype / enum / no constructor) or a mechanical repo check
    CoreCited,    // already tested in sheaf-core; pinned by test name + G0 manifest; not re-proved here
    Differential, // two independent implementations must agree (core binary vs typed Rust); disagreement = HOLE
    Property,     // Haskell + Rust property checks, same id
}

/// Spec status. `Proposed`: stated, believed true, not necessarily built. `Hole`: an open question or a known
/// violation; the law cannot be discharged as things stand. Holes propagate to every law derived from them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status { Proposed, Hole }

#[derive(Debug, Clone)]
pub struct Law {
    pub id: &'static str,
    pub layer: Layer,
    pub scope: Scope,
    pub tier: Tier,
    pub status: Status,
    pub statement: &'static str,
    /// For CoreCited: sheaf-core test names. For Formal: Bend law names. Otherwise artifacts or review ids.
    pub evidence: &'static [&'static str],
    /// Laws (lower layer, or glue at the same layer) whose conjunction yields this one.
    pub derives_from: &'static [&'static str],
    /// The law this one refines (its parent in the ontology tree). Root laws have none.
    pub parent: Option<&'static str>,
}

macro_rules! law {
    ($id:literal, $layer:ident, $scope:expr, $tier:ident, $status:ident, $stmt:literal,
     ev[$($ev:literal),*], from[$($d:literal),*], parent $p:expr) => {
        Law { id: $id, layer: Layer::$layer, scope: $scope, tier: Tier::$tier, status: Status::$status,
              statement: $stmt, evidence: &[$($ev),*], derives_from: &[$($d),*], parent: $p }
    };
}

pub fn laws() -> Vec<Law> {
    use Op::*; use Prod::*; use Scope::*;
    vec![
    // ============================================================ L0 GLOBAL
    law!("G0", Global, Root, Types, Proposed,
        "Core invariance: crates/**, Cargo.lock, Cargo.toml and rust-toolchain.toml are byte-identical to the pinned manifest before every run; no product law is discharged by editing the core.",
        ev["products/gates/core.sha256", "products/gates/g0_invariance.sh"], from["S.bind.manifest", "I.analyze.pin"], parent None),
    law!("G1", Global, Root, Property, Proposed,
        "Claim soundness: a receipt other than BLOCKED is issued only when the core produced a Report (no CoreError), the run is certified (spectral, Hodge and extension converged, cut consistent, not ambiguous, relative H0 computed when a boundary was given, report well-formed), and the policy floor is cleared.",
        ev["seat1 S1-08", "seat3 S3-04"], from["S.certify.wf", "A.certify.converged", "A.certify.cut", "A.certify.unique-null", "A.classify.uncertified", "S.judge.map", "I.analyze.certify", "I.analyze.certify-error"], parent None),
    law!("G2", Global, Root, Property, Proposed,
        "Traceability: every number in a Finding is a Report field verbatim, a Circulate value, or a tick produced by Quantize from one; rankings are complete (no truncation); no f64 arithmetic happens outside the core, Circulate and Quantize.",
        ev["seat3 S3-11"], from["S.attribute.projection", "S.quantize.only-crossing", "I.analyze.attribute"], parent None),
    law!("G3", Global, Root, Property, Proposed,
        "Determinism: equal (spec, CorePolicy, Boundary) bytes under a pinned thread count give equal verdicts, findings and receipts (excluding the timestamp). Report bytes may differ in the last digits of spectral fields across thread counts, which is why the count is pinned and bound.",
        ev["seat1 S1-04 (RAYON_NUM_THREADS 1/2/8 gave three report md5s)"], from["S.analyze.thread-pin", "S.analyze.determinism", "S.join.monoid", "A.circulate.forest", "S.bind.binding"], parent None),
    law!("G4", Global, Root, Property, Proposed,
        "Gauge invariance: on the Circulate domain (scalar or square-invertible maps), the loop verdict is unchanged by adding any exact cochain delta-p to w, provided the run stays certified and the policy floor is still cleared after the shift. Consequence: no verdict reads harmonic SHARE (Spec v1 L5-L7 retracted; E2a arbitrage at share 2e-5).",
        ev["E2a", "seat1 S1-02, S1-15"], from["A.circulate.exact-zero", "A.circulate.linear", "A.classify.circ-only", "S.circulate.domain", "S.quantize.floor"], parent None),
    law!("G5", Global, Root, Property, Proposed,
        "Clean-twin completeness: x a global section (x-mode) or w = delta p (w-mode) is reported consistent / HOLDS. In w-mode node discord is NotApplicable (an exact w is consistent pairwise evidence, not a node error), so only circulations can fire.",
        ev["seat1 S1-01 (w = delta p ranked nodes by |w|^2)"], from["G4", "A.quantize.disc-source", "S.quantize.floor", "S.analyze.x-only-exact"], parent None),
    law!("G6", Global, Root, Differential, Hole,
        "Compositionality (parallel associativity): the receipt of a multi-component spec equals the join of the receipts of its components analysed separately; findings partition by component. HOLE until S.analyze.hodge-local is shown: batching can mask a small component's loop under the core's |w|^2-relative row filter.",
        ev["seat1 S1-07 (d_batch_mask vs d_tri_alone)"], from["S.analyze.component-local", "S.analyze.hodge-local", "S.join.monoid", "S.join.cover", "I.judge.join"], parent None),
    law!("G7", Global, Root, Types, Proposed,
        "Tier honesty: 'proved' is attached only to Formal laws whose Bend proof and --verdict checked; numeric fields are worded 'certified residual'; approximate layers never write exact fields (ADR-008, ADR-013).",
        ev["seat2 S2-12"], from["S.bind.tier", "S.explain.wording"], parent None),
    law!("G8", Global, Root, Property, Proposed,
        "Surface equivalence: for the same (Corpus, policies, Boundary) the CLI, the MCP tool and the HTTP handler return byte-identical Findings; a transport only wraps.",
        ev["seat3 S3-09"], from["S.surface.envelope", "S.explain.pure", "G3"], parent None),
    law!("G9", Global, Root, Differential, Hole,
        "Able to fail: every blocking oracle kills its planted mutants at Touchstone layers 0 (Bend), 1 (Haskell), 2 (Rust), 3 (interface) and 4 (verifier, planted-bug corpus with a hidden holdout). HOLE until layers 3 and 4 exist.",
        ev["seat2 S2-05"], from["V.mut.L0", "V.mut.L1", "V.mut.L2", "V.mut.L3", "V.mut.L4"], parent None),
    law!("G10", Global, Root, Property, Hole,
        "Separation of duties: laws are signed by a human before they gate anything; verifier-owned paths (laws, tests, golden tables, gates) are not writable by builder agents; verdicts come from a different model family than the builder. HOLE: all three are unbuilt today.",
        ev["seat2 S2-01, S2-03, S2-17"], from["V.sep.signed", "V.sep.fence", "V.sep.family"], parent None),
    law!("G11", Global, Root, Types, Proposed,
        "Language and ship path: no Python on any gate, oracle or receipt path; shipped verdict code is the Rust transcription, never Bend-compiled output (a Bend proof never approves code alone).",
        ev["seat2 S2-13"], from["V.lang.no-python", "V.lang.ship-rust"], parent None),

    // ============================================================ L1 STAGES (local)
    law!("S.ingest.total-or-reject", Stage, Stage(Ingest), Types, Proposed,
        "An unknown status, kind or endpoint, or a self-loop, yields Rejected(reason) (exit 3), never a default value.", ev["specv1 memory_audit exit 3"], from[], parent Some("G1")),
    law!("S.ingest.ids-stable", Stage, Stage(Ingest), Formal, Proposed,
        "Entity ids are the note slug or (system, key): equal corpora give equal ids, so Diff can match across runs.", ev["ingest_ids_stable"], from[], parent Some("G3")),
    law!("S.encode.valid", Stage, Stage(Encode), Formal, Proposed,
        "Encode is valid by construction: for every well-formed domain graph the spec passes Validate (endpoints exist, u != v, map shapes de x dim, identity only when dims agree, finite entries <= 1e75).",
        ev["encode_valid"], from["A.encode.endpoints", "A.encode.shapes"], parent Some("G1")),
    law!("S.encode.index", Stage, Stage(Encode), Formal, Proposed,
        "Encode fixes an order-preserving bijection between domain relations and spec edges; Attribute and Explain use its inverse.", ev["encode_index_bij"], from[], parent Some("G2")),
    law!("S.encode.mode", Stage, Stage(Encode), Formal, Proposed,
        "Mode honesty: w is sent iff every relation has an observation; a partial observation set is rejected.", ev["encode_mode"], from[], parent Some("G5")),
    law!("S.encode.one-cochain", Stage, Stage(Encode), Types, Proposed,
        "Exactly one of x, w is sent; both is rejected (the core silently prefers w, v1.rs:238).", ev["seat1 S1-18"], from[], parent Some("G5")),
    law!("S.validate.mirror", Stage, Stage(Validate), Differential, Proposed,
        "The products-side validator accepts exactly what the core parser accepts, on generated and mutated specs.",
        ev["v1_rejects_bad_specs", "check_edges_refuses_each_bad_edge"], from[], parent Some("G1")),
    law!("S.analyze.outcome", Stage, Stage(Analyze), Types, Proposed,
        "A core run yields exactly one of Report or CoreError{exit 1|2|signal|timeout, stderr}; a CoreError is never coerced into a Report; the P_INPUT_01 text is kept byte for byte.", ev["seat3 S3-04"], from[], parent Some("G1")),
    law!("S.analyze.complex", Stage, Stage(Analyze), CoreCited, Proposed,
        "delta delta = 0; a twisted map on a filled face is rejected as not a complex.",
        ev["d_squared_is_zero", "twisted_map_on_filled_face_is_rejected_as_not_a_complex"], from[], parent Some("G1")),
    law!("S.analyze.hodge", Stage, Stage(Analyze), CoreCited, Proposed,
        "Hodge split: exact and harmonic parts are orthogonal and sum to w; well defined for any map shape.", ev["hodge_components_orthogonal_and_sum"], from[], parent Some("G4")),
    law!("S.analyze.x-only-exact", Stage, Stage(Analyze), CoreCited, Proposed,
        "With node data only, w := delta x is exact: harmonic share is zero.", ev["v1_x_only_is_exact_and_rotation_loop_kills_h0"], from[], parent Some("G5")),
    law!("S.analyze.holonomy", Stage, Stage(Analyze), CoreCited, Proposed,
        "Holonomy around an unfilled loop removes global sections; K4 keeps its loop space.", ev["holonomy_kills_h0_and_k4_keeps_loop_space"], from[], parent Some("G4")),
    law!("S.analyze.component-local", Stage, Stage(Analyze), CoreCited, Proposed,
        "H0 cut and spectra are per connected component (L0 block diagonal).", ev["v1_multi_component_h0_matches_dense_count", "v1_two_scale_components_use_local_zero_cut"], from[], parent Some("G6")),
    law!("S.analyze.hodge-local", Stage, Stage(Analyze), Differential, Hole,
        "Per-edge harmonic energies of a multi-component spec equal those of each component analysed alone up to the certified residual, and no row a policy-visible tolerance would count is dropped by the |w|^2-relative filter or by top. OPEN: the Hodge CG is a global solve with a global stop rule.",
        ev["seat1 S1-07"], from[], parent Some("G6")),
    law!("S.analyze.determinism", Stage, Stage(Analyze), CoreCited, Proposed,
        "Two runs of the v1 path in one process give an identical report (same threads).", ev["v1_triangle_harmonic_share_is_one_third_and_report_is_deterministic"], from[], parent Some("G3")),
    law!("S.analyze.thread-pin", Stage, Stage(Analyze), Types, Proposed,
        "The products always run the core with RAYON_NUM_THREADS pinned (default 1) and record the count in CorePolicy.", ev["seat1 S1-04"], from[], parent Some("G3")),
    law!("S.analyze.top-complete", Stage, Stage(Analyze), Types, Proposed,
        "The core is always invoked with top >= |V| + |E| of the spec, so every ranking is complete; a truncated ranking is BLOCKED(incomplete).", ev["seat3 S3-11"], from[], parent Some("G2")),
    law!("S.analyze.extension-pins", Stage, Stage(Analyze), Differential, Proposed,
        "Pinned boundary values come back equal after the report's 12-significant-digit rounding; products compare sig(pin), never the raw pin.", ev["seat1 S1-05"], from[], parent Some("G2")),
    law!("S.analyze.extension-unique", Stage, Stage(Analyze), CoreCited, Proposed,
        "unique iff dim H0(X,B;F) = 0 under the run's H0 policy; an unpinned component is reported not unique; an empty boundary equals plain H0.",
        ev["unpinned_component_is_reported_not_unique_and_bad_boundary_rejected", "empty_boundary_equals_absolute_h0_on_the_catalogue"], from[], parent Some("G1")),
    law!("S.analyze.extension-minimal", Stage, Stage(Analyze), CoreCited, Proposed,
        "The extension minimises energy with the boundary fixed; when not unique the values are the minimum-norm minimiser.",
        ev["extension_minimises_energy", "path_extension_is_linear_and_unique", "no_boundary_reduces_to_h0_projection"], from[], parent Some("G5")),
    law!("S.certify.wf", Stage, Stage(Certify), Property, Proposed,
        "Report well-formedness: h1 = dim_C1 + h0 - dim_C0 (graphs); if norm2 > 0 then |exact_share + harmonic_share - 1| <= 2 cert_coclosed + 1e-11, else both shares are 0 and w-mode evidence is NotApplicable; unique = (h0_relative == 0) when present. Any violation is BLOCKED(report).",
        ev["rank_nullity_and_euler", "seat1 S1-13"], from[], parent Some("G1")),
    law!("S.certify.policy", Stage, Stage(Certify), Types, Proposed,
        "Cert carries zero_policy, zero_rule and h0_absolute_tol; verdicts are never compared across runs with different policies.", ev["seat1 S1-08"], from[], parent Some("G3")),
    law!("S.circulate.domain", Stage, Stage(Circulate), Types, Proposed,
        "Circulate is defined only when every edge map is a scalar or a square invertible matrix on uniform stalks (de = dim u = dim v); otherwise it returns NotApplicable and the loop decision is single-oracle (core harmonic energy), recorded as such in the receipt.",
        ev["seat1 S1-02 (a_tree_nonsurj, b_rate_triangle_exact, c_zero_map_edge)"], from[], parent Some("G4")),
    law!("S.circulate.exact", Stage, Stage(Circulate), Differential, Proposed,
        "On its domain, Circulate returns for each fundamental cycle of a fixed spanning forest the transported sum of w around it; exact for rational scalar rates.",
        ev[], from["A.circulate.exact-zero", "A.circulate.linear", "A.circulate.basis-size", "A.circulate.forest"], parent Some("G4")),
    law!("S.quantize.only-crossing", Stage, Stage(Quantize), Types, Proposed,
        "Quantize is the only place an f64 becomes discrete; Classify accepts comparison bits and Ticks<Unit> only, never f64.", ev[], from[], parent Some("G2")),
    law!("S.quantize.threshold", Stage, Stage(Quantize), Types, Proposed,
        "Threshold faithfulness by construction: each comparison v >= tol is made once in f64 inside Quantize and carried as a Bool; ticks only carry magnitude.",
        ev["seat1 S1-12"], from["A.quantize.compare-once", "A.quantize.saturate"], parent Some("G4")),
    law!("S.quantize.floor", Stage, Stage(Quantize), Property, Proposed,
        "Fail-before-freeze for policies: tol >= kappa * max(certified residuals, cert_coclosed * |w|, 1e-12 |v| (report rounding), 1e-24 |w|^2 for energies (row filter), ulp); otherwise BLOCKED(policy).",
        ev["seat1 S1-03, S1-12"], from[], parent Some("G5")),
    law!("S.quantize.nonfinite", Stage, Stage(Quantize), Property, Proposed,
        "NaN, inf and null report values quantise to BLOCKED, never to 0.", ev[], from[], parent Some("G1")),
    law!("S.classify.total", Stage, Stage(Classify), Types, Proposed,
        "Classify is a total, terminating function to a closed verdict enum (Bend totality checker + Rust closed enum).",
        ev["seat1 S1-16"], from["A.classify.uncertified", "A.classify.undetermined", "A.classify.xonly", "A.classify.mono", "A.classify.circ-only", "A.classify.loop-iff", "A.classify.disc-iff", "A.classify.consistent-iff", "A.classify.h0-policy"], parent Some("G1")),
    law!("S.attribute.projection", Stage, Stage(Attribute), Property, Proposed,
        "Attribute only projects: Finding numbers are built through FromReport<T> / FromCirc<T> constructors; there is no public constructor from f64.",
        ev[], from["A.attribute.loop-blame", "A.attribute.fixable", "A.attribute.fill", "A.attribute.x-only", "A.attribute.context-only", "A.attribute.ranking"], parent Some("G2")),
    law!("S.judge.map", Stage, Stage(Judge), Formal, Proposed,
        "Verdict map: consistent -> HOLDS (and only consistent); fixable, loop or map-conflict -> BROKEN; inconclusive or underdetermined -> BLOCKED; two oracles disagreeing -> HOLE. Never averaged.",
        ev["receipt_map_hole", "receipt_map_uncertified", "holds_iff_consistent"], from[], parent Some("G1")),
    law!("S.join.monoid", Stage, Stage(Join), Formal, Proposed,
        "Join is a bounded join-semilattice on HOLDS < BROKEN < BLOCKED < HOLE: associative, commutative, idempotent, HOLDS the unit, HOLE absorbing.",
        ev["join_assoc", "join_comm", "join_idem", "join_unit_l", "join_unit_r", "join_hole_absorbs"], from[], parent Some("G6")),
    law!("S.join.cover", Stage, Stage(Join), Formal, Proposed,
        "A system verdict is a join over a declared cover; a missing component contributes HOLE (a child can never claim system-level HOLDS).", ev["join_cover_missing_is_hole"], from[], parent Some("G6")),
    law!("S.join.is-max", Stage, Stage(Join), Property, Proposed,
        "join(a, b) = max(a, b) in the order HOLDS < BROKEN < BLOCKED < HOLE.", ev["join_is_max"], from[], parent Some("G6")),
    law!("S.join.no-drop", Stage, Stage(Join), Types, Proposed,
        "The system receipt lists every non-HOLDS component verdict with its finding or question; joining never discards a BROKEN under a HOLE.", ev["seat2 S2-14"], from[], parent Some("G6")),
    law!("S.bind.binding", Stage, Stage(Bind), Property, Proposed,
        "A receipt binds sha256 of the core manifest, laws.lock, spec bytes, outcome bytes, CorePolicy (incl. thread count), DomainPolicy and finding; changing any one changes the receipt hash.", ev["seat1 S1-04", "seat3 S3-12"], from[], parent Some("G3")),
    law!("S.bind.tier", Stage, Stage(Bind), Types, Proposed,
        "Receipt.tier is an enum; proved() is Some only for Tier::Formal with a non-null, checked proof hash.", ev["seat2 S2-12"], from[], parent Some("G7")),
    law!("S.bind.manifest", Stage, Stage(Bind), Types, Proposed,
        "Bind refuses unless the G0 manifest check passed in this run; the manifest covers the core only (laws are frozen separately in laws.lock).", ev["seat2 S2-01", "seat3 S3-18"], from[], parent Some("G0")),
    law!("S.explain.pure", Stage, Stage(Explain), Types, Proposed,
        "An Explanation is a function of the Finding alone: no Report access, no model call.", ev["seat3 S3-02"], from[], parent Some("G8")),
    law!("S.explain.wording", Stage, Stage(Explain), Types, Proposed,
        "Every Explanation carries 'certified residual, not a proof that any fact is true' verbatim and never the word 'proved' for a numeric result.", ev["seat3 S3-02"], from[], parent Some("G7")),
    law!("S.diff.keyed", Stage, Stage(Diff), Formal, Proposed,
        "Findings are matched by entity / edge id (stable Ingest ids), never by rank position.", ev["diff_keyed"], from[], parent Some("G3")),
    law!("S.diff.partition", Stage, Stage(Diff), Formal, Proposed,
        "Every entry of either Finding lands in exactly one of {new, resolved, unchanged, drifted}.", ev["diff_partition"], from[], parent Some("G3")),
    law!("S.diff.reflexive", Stage, Stage(Diff), Formal, Proposed,
        "Diff(f, f) is empty.", ev["diff_reflexive"], from[], parent Some("G3")),
    law!("S.diff.exit", Stage, Stage(Diff), Types, Proposed,
        "Exit code 2 iff new is non-empty; 3 on Rejected; 0 otherwise (the --baseline CI gate).", ev["specv1 memory_audit --baseline"], from[], parent Some("G3")),
    law!("S.fill.determined-only", Stage, Stage(Fill), Formal, Proposed,
        "A domain value changes only if its fill is marked determined (extension unique).", ev["fill_determined"], from[], parent Some("G1")),
    law!("S.fill.ask", Stage, Stage(Fill), Types, Proposed,
        "When any fill is not determined the output is Ask(free_entities), not a graph.", ev["seat3 S3-05"], from[], parent Some("G1")),
    law!("S.fill.idempotent", Stage, Stage(Fill), Property, Proposed,
        "Re-running the pipeline on written-back data makes no further change: the remaining flags are irreducible.", ev["specv1 memory_audit step 3"], from[], parent Some("G5")),
    law!("S.surface.envelope", Stage, Root, Types, Proposed,
        "Transports (CLI, MCP stdio, HTTP) wrap a Finding in an envelope and never alter it.", ev["seat3 S3-09"], from[], parent Some("G8")),
    law!("V.mut.L0", Stage, Root, Differential, Proposed,
        "Bend layer: every planted mutant of the atoms makes at least one closed law fail to check (judged per law id, never by text diff).", ev["products/gates/atom_gate.sh step 7"], from[], parent Some("G9")),
    law!("V.mut.L1", Stage, Root, Differential, Proposed,
        "Haskell layer: the same mutants applied to Laws.hs each make a law BROKEN.", ev["products/gates/atom_gate.sh step 7b"], from[], parent Some("G9")),
    law!("V.mut.L2", Stage, Root, Differential, Proposed,
        "Rust layer: cargo-mutants plus constructor-swap mutants are each killed by a law test (not only by the golden table).", ev["products/kit/receipts"], from[], parent Some("G9")),
    law!("V.mut.L3", Stage, Root, Differential, Hole,
        "Interface layer: mutating interface.lock (swap unit types, reorder fields, flip defaults) breaks a build or a consumer test.", ev["seat2 S2-05"], from[], parent Some("G9")),
    law!("V.mut.L4", Stage, Root, Differential, Hole,
        "Verifier layer: a planted-bug corpus kept outside the repo by hash, re-planted weekly, with a hidden holdout scored separately.", ev["seat2 S2-05"], from[], parent Some("G9")),
    law!("V.sep.signed", Stage, Root, Types, Hole,
        "A law id enters laws.lock only with a human signed_by and a checker receipt (HOLDS or BROKEN, never UNCHECKED).", ev["seat2 S2-01, S2-03"], from[], parent Some("G10")),
    law!("V.sep.fence", Stage, Root, Types, Hole,
        "A hook denies builder writes to products/bend/LAWS.bend, laws.lock, kit/tests, kit/bridge*, haskell/, gates/.", ev["seat2 S2-03"], from[], parent Some("G10")),
    law!("V.sep.family", Stage, Root, Property, Hole,
        "At least one verifying seat per review round is from a different model family than the builder (today: all Claude).", ev["seat2 S2-17"], from[], parent Some("G10")),
    law!("V.lang.no-python", Stage, Root, Types, Proposed,
        "No .py file under products/{laws,kit,haskell,bend,gates} and no python invocation in any gate script (checked by a test).", ev["tests/ontology.rs no_python_on_the_gate_path"], from[], parent Some("G11")),
    law!("V.lang.ship-rust", Stage, Root, Types, Proposed,
        "The shipped classifier, judge and join are sheaf-kit (Rust); Bend output is never linked into a product.", ev["seat2 S2-13"], from[], parent Some("G11")),

    // ============================================================ L1 EDGES (glue: one per wired producer -> consumer pair)
    law!("I.ingest.encode", Stage, Edge(Ingest, Encode), Differential, Proposed,
        "The Rust ingester and the reference adapter emit byte-identical DomainGraphs on the fixture set.", ev["seat3 S3-06, C40"], from[], parent Some("G3")),
    law!("I.encode.validate", Stage, Edge(Encode, Validate), Differential, Proposed,
        "Every encoder output is accepted by the core binary (exit 0) on fixtures and generated domain graphs.", ev[], from["S.encode.valid", "S.validate.mirror"], parent Some("G1")),
    law!("I.validate.analyze", Stage, Edge(Validate, Analyze), Types, Proposed,
        "Analyze receives exactly the bytes Validate accepted, plus CorePolicy and Boundary as separate typed inputs.", ev["seat3 S3-10"], from[], parent Some("G3")),
    law!("I.validate.circulate", Stage, Edge(Validate, Circulate), Types, Proposed,
        "Circulate's edge list is ValidSpec's, in spec index order.", ev["seat3 S3-10"], from[], parent Some("G2")),
    law!("I.analyze.pin", Stage, Edge(Analyze, Certify), Types, Proposed,
        "The report binding is generated from report.rs at the pinned commit; a missing or renamed field fails to parse and yields BLOCKED(report-schema).", ev[], from[], parent Some("G0")),
    law!("I.analyze.certify", Stage, Edge(Analyze, Certify), Types, Proposed,
        "Certify reads only certificate fields and copies them unchanged.", ev[], from["I.analyze.pin"], parent Some("G1")),
    law!("I.analyze.certify-error", Stage, Edge(Analyze, Certify), Types, Proposed,
        "A CoreError maps to Cert{certified: false, reason: CoreCrash | InputRejected | Timeout}.", ev["seat3 S3-04"], from[], parent Some("G1")),
    law!("I.analyze.quantize", Stage, Edge(Analyze, Quantize), Types, Proposed,
        "Quantize reads exactly edges_by_harmonic_energy[].harmonic_energy, nodes_by_incident_discord[].incident_discord, hodge.cert_coclosed, h0.zero_tol and the residuals; never harmonic_share and never the nulls.", ev["seat3 S3-10, S3-13"], from[], parent Some("G4")),
    law!("I.circulate.quantize", Stage, Edge(Circulate, Quantize), Differential, Proposed,
        "Two-oracle agreement (constant R^1 sheaf, C the cycle x edge incidence of the basis, c = C w): the core's unfiltered harmonic energy E_h = c^T (C C^T)^-1 c lies in [|c|^2 / lambda_max(C C^T), |c|^2 / lambda_min(C C^T)]; outside the interval -> HOLE.",
        ev["seat1 S1-03 (j_cycle100: E_h = c^2 / L)"], from["S.circulate.exact", "S.analyze.hodge"], parent Some("G4")),
    law!("I.certify.classify", Stage, Edge(Certify, Classify), Types, Proposed,
        "certified = converged AND cut AND not unique-null AND wf; nothing else reaches Classify from Certify.", ev["seat3 S3-10"], from[], parent Some("G1")),
    law!("I.quantize.classify", Stage, Edge(Quantize, Classify), Types, Proposed,
        "Evidence and tolerance share one unit: Ticks<Loop> and Ticks<Disc> do not mix.", ev[], from[], parent Some("G2")),
    law!("I.analyze.attribute", Stage, Edge(Analyze, Attribute), Types, Proposed,
        "Attribute reads only the complete rankings and extension values (S.analyze.top-complete holds).", ev["seat3 S3-11"], from[], parent Some("G2")),
    law!("I.ingest.attribute", Stage, Edge(Ingest, Attribute), Types, Proposed,
        "Finding entity ids are Ingest ids; Attribute never invents an id.", ev["seat3 S3-06"], from[], parent Some("G2")),
    law!("I.classify.attribute", Stage, Edge(Classify, Attribute), Types, Proposed,
        "Finding.verdict is the Classify verdict by construction (newtype, no setter).", ev["seat3 S3-10"], from[], parent Some("G1")),
    law!("I.circulate.attribute", Stage, Edge(Circulate, Attribute), Formal, Proposed,
        "Every loops[i].edges is a fundamental cycle of Circulate's basis whose circulation cleared tolerance; edge ids are Encode indices.", ev["loops_are_basis_cycles"], from[], parent Some("G2")),
    law!("I.attribute.judge", Stage, Edge(Attribute, Judge), Formal, Proposed,
        "The Touchstone verdict is a function of the Finding's verdict and the oracle-agreement bit only.", ev["receipt_map_hole"], from["S.judge.map"], parent Some("G1")),
    law!("I.certify.judge", Stage, Edge(Certify, Judge), Types, Proposed,
        "Judge sees certification only through Cert (the single-oracle flag from S.circulate.domain travels here).", ev["seat3 S3-10"], from[], parent Some("G1")),
    law!("I.judge.join", Stage, Edge(Judge, Join), Differential, Proposed,
        "The cover is built from the spec's connected components (union-find on the spec sent, count checked against report.input.components); each Join input carries a distinct component id; mismatch -> HOLE.",
        ev["seat1 S1-14", "seat2 S2-15", "seat3 S3-10"], from[], parent Some("G6")),
    law!("I.join.bind", Stage, Edge(Join, Bind), Types, Proposed,
        "Bind receives the system verdict together with every non-HOLDS component (S.join.no-drop).", ev["seat2 S2-14"], from[], parent Some("G6")),
    law!("I.attribute.bind", Stage, Edge(Attribute, Bind), Property, Proposed,
        "The receipt hash covers the Finding bytes.", ev[], from[], parent Some("G3")),
    law!("I.attribute.explain", Stage, Edge(Attribute, Explain), Types, Proposed,
        "Explain receives the Finding only.", ev["seat3 S3-02"], from[], parent Some("G8")),
    law!("I.attribute.diff", Stage, Edge(Attribute, Diff), Types, Proposed,
        "Both Diff inputs are full Findings of runs with equal policies (S.certify.policy); otherwise Diff refuses.", ev["seat3 S3-03"], from[], parent Some("G3")),
    law!("I.attribute.fill", Stage, Edge(Attribute, Fill), Types, Proposed,
        "Fill reads fills[] only; determined comes from Cert-backed uniqueness.", ev["seat3 S3-05"], from[], parent Some("G1")),
    law!("I.ingest.fill", Stage, Edge(Ingest, Fill), Types, Proposed,
        "Fill writes back to exactly the ids Ingest produced, in the corpus's own format.", ev["seat3 S3-05"], from[], parent Some("G2")),
    law!("I.encode.attribute", Stage, Edge(Encode, Attribute), Formal, Proposed,
        "Round trip: Attribute's edge map composed with Encode's index is the identity on relations.", ev["index_roundtrip"], from["S.encode.index"], parent Some("G2")),

    // ============================================================ L2 ATOMS
    law!("A.encode.endpoints", Atom, Stage(Encode), Formal, Proposed,
        "Every encoded edge has u < n, v < n and u != v.", ev["encode_ends_ok"], from[], parent Some("S.encode.valid")),
    law!("A.encode.shapes", Atom, Stage(Encode), Formal, Proposed,
        "Every encoded map has de rows and dim(endpoint) columns; omitted maps only when dims agree.", ev["encode_shapes_ok"], from[], parent Some("S.encode.valid")),
    law!("A.certify.converged", Atom, Stage(Certify), Types, Proposed,
        "spectral.converged AND hodge.converged AND (no boundary OR extension.converged).", ev["seat3 S3-16"], from[], parent Some("S.certify.wf")),
    law!("A.certify.cut", Atom, Stage(Certify), Types, Proposed,
        "cut_consistent_with_zero_tol AND NOT near_zero_ambiguous.", ev["seat3 S3-16"], from[], parent Some("S.certify.wf")),
    law!("A.certify.unique-null", Atom, Stage(Certify), Types, Proposed,
        "unique = null or lower_bound_only = true -> not certified, reason NotComputed; never read as unique = true (three-valued Bounded | Free | Unknown).", ev["seat1 S1-08", "seat3 S3-07"], from[], parent Some("S.certify.wf")),
    law!("A.circulate.exact-zero", Atom, Stage(Circulate), Formal, Proposed,
        "Telescoping with maps: for scalar rational rates (num, den), the transported circulation of delta p around any cycle is zero (over Nat: product-scaled heads equal product-scaled tails).", ev["circ_exact_zero"], from[], parent Some("S.circulate.exact")),
    law!("A.circulate.linear", Atom, Stage(Circulate), Formal, Proposed,
        "circ(w1 + w2) = circ(w1) + circ(w2) on every cycle (signed values as Bool + Nat pairs).", ev["circ_linear"], from[], parent Some("S.circulate.exact")),
    law!("A.circulate.basis-size", Atom, Stage(Circulate), Differential, Proposed,
        "The fundamental cycle count |E| - |V| + #components equals the core's h1 only for the constant R^1 sheaf with identity maps; otherwise no relation is claimed.", ev["seat1 S1-06"], from[], parent Some("S.circulate.exact")),
    law!("A.circulate.forest", Atom, Stage(Circulate), Formal, Proposed,
        "The spanning forest is BFS in spec edge order from the lowest node index, so the cycle basis is a function of the spec bytes.", ev["forest_deterministic"], from[], parent Some("S.circulate.exact")),
    law!("A.quantize.compare-once", Atom, Stage(Quantize), Types, Proposed,
        "The only f64 comparisons are v >= tol inside Quantize; Classify reads the resulting Bool.", ev["seat1 S1-12"], from[], parent Some("S.quantize.threshold")),
    law!("A.quantize.saturate", Atom, Stage(Quantize), Formal, Proposed,
        "Ticks saturate at 2^53; a saturated tick never changes a comparison bit.", ev["ticks_saturate"], from[], parent Some("S.quantize.threshold")),
    law!("A.quantize.disc-source", Atom, Stage(Quantize), Types, Proposed,
        "In x-mode disc is the incident discord of delta x; in w-mode disc is NotApplicable (an exact w is consistent evidence) and only circulations can fire.", ev["seat1 S1-01"], from[], parent Some("S.quantize.floor")),
    law!("A.classify.uncertified", Atom, Stage(Classify), Formal, Proposed,
        "Not certified -> inconclusive, whatever the evidence.", ev["sound_uncertified"], from[], parent Some("S.classify.total")),
    law!("A.classify.undetermined", Atom, Stage(Classify), Formal, Proposed,
        "Certified, boundary given, not unique -> underdetermined.", ev["undetermined"], from[], parent Some("S.classify.total")),
    law!("A.classify.xonly", Atom, Stage(Classify), Formal, Proposed,
        "x-only evidence never yields a loop conflict.", ev["xonly_no_loop"], from[], parent Some("S.classify.total")),
    law!("A.classify.mono", Atom, Stage(Classify), Formal, Proposed,
        "Raising the loop tolerance never turns consistent into loop conflict (a theorem needing a Nat.is_ge monotonicity lemma; in OPEN.bend until proved).", ev["mono_tolerance"], from[], parent Some("S.classify.total")),
    law!("A.classify.circ-only", Atom, Stage(Classify), Formal, Proposed,
        "The loop decision reads only whether the max circulation clears its tolerance.", ev["circ_only"], from[], parent Some("S.classify.total")),
    law!("A.classify.loop-iff", Atom, Stage(Classify), Formal, Proposed,
        "Certified and bounded: LoopConflict exactly when the input is not x-only and the max circulation clears its tolerance.", ev["loop_iff_circ"], from[], parent Some("S.classify.total")),
    law!("A.classify.disc-iff", Atom, Stage(Classify), Formal, Proposed,
        "Certified and bounded: FixableNodes exactly when discord clears its tolerance and no loop fires (added after mutation testing).", ev["fixable_iff_disc"], from[], parent Some("S.classify.total")),
    law!("A.classify.consistent-iff", Atom, Stage(Classify), Formal, Proposed,
        "Certified and bounded: Consistent exactly when no loop fires and discord is below tolerance (seat 2: the Consistent branch was unpinned).", ev["consistent_iff"], from[], parent Some("S.classify.total")),
    law!("A.classify.h0-policy", Atom, Stage(Classify), Formal, Proposed,
        "Certified and h0 < h0_expected -> MapConflict (the maps themselves admit no section), which outranks fixable and is never reported as a node error (seat 1: E3c twisted frames read HOLDS).", ev["h0_map_conflict"], from[], parent Some("S.classify.total")),
    law!("A.attribute.loop-blame", Atom, Stage(Attribute), Types, Proposed,
        "A loop finding names a cycle (two or more edges, closed walk), never a single edge as the cause.", ev[], from[], parent Some("S.attribute.projection")),
    law!("A.attribute.fixable", Atom, Stage(Attribute), Formal, Proposed,
        "The fixable set equals {v : incident-discord ticks >= tolerance} over all nodes (complete rankings).", ev["fixable_above_tol"], from[], parent Some("S.attribute.projection")),
    law!("A.attribute.ranking", Atom, Stage(Attribute), Property, Proposed,
        "Fixable nodes are listed by incident discord descending, ties by entity id.", ev["seat3 C20"], from[], parent Some("S.attribute.projection")),
    law!("A.attribute.fill", Atom, Stage(Attribute), Formal, Proposed,
        "A fill is marked determined only when the extension is certified unique.", ev["fill_determined"], from[], parent Some("S.attribute.projection")),
    law!("A.attribute.x-only", Atom, Stage(Attribute), Types, Proposed,
        "For x-only input the loops field is NotApplicable, never an empty list.", ev[], from[], parent Some("S.attribute.projection")),
    law!("A.attribute.context-only", Atom, Stage(Attribute), Types, Proposed,
        "Finding.context (h0, h1, harmonic_share, nulls) has no reader inside the pipeline; Explain may render it.", ev["seat3 S3-13"], from[], parent Some("S.attribute.projection")),
    law!("A.explain.loop-blame", Atom, Stage(Explain), Formal, Proposed,
        "For a loop entry the text names at least two edges and no single edge as the cause (L-WORD-1).", ev["explain_loop_blame"], from[], parent Some("S.explain.wording")),
    law!("A.explain.resolve", Atom, Stage(Explain), Property, Proposed,
        "what_would_resolve names the node for a fixable entry and the cycle edges for a loop entry.", ev["seat3 S3-02"], from[], parent Some("S.explain.pure")),

    // ============================================================ L3 PRODUCT REFINEMENTS (leaves)
    law!("P.memory.status-table", Product, Product(Memory, Ingest), Formal, Proposed,
        "The status string -> quarter-unit level map is the fixed 7-row table of Spec v1 (verified 1 ... refuted -0.5).", ev["status_table"], from[], parent Some("S.ingest.total-or-reject")),
    law!("P.memory.flag-iff-differs", Product, Product(Memory, Encode), Formal, Proposed,
        "Quarter units with gain 4: discord is (level_u - level_v)^2 exactly; a supports-edge flags iff the statuses differ (Spec v1 L8).", ev["flag_iff_status_differs"], from[], parent Some("A.encode.shapes")),
    law!("P.memory.claim-rule", Product, Product(Memory, Encode), Formal, Proposed,
        "claim_of edge flags iff assertiveness exceeds the source's status (L9); equal status never flags (L10).", ev["claim_flag_rule", "equal_status_never_flags"], from[], parent Some("A.encode.shapes")),
    law!("P.memory.contradiction-parity", Product, Product(Memory, Classify), Formal, Proposed,
        "A cycle of contradicts (-1) maps has a nonzero section iff it has an even number of contradicts edges; odd -> MapConflict (E1b: H0 = 0).", ev["contradiction_parity"], from[], parent Some("A.classify.h0-policy")),
    law!("P.memory.writeback", Product, Product(Memory, Fill), Formal, Proposed,
        "Status write-back happens only for fills marked determined.", ev["fill_determined"], from[], parent Some("S.fill.determined-only")),
    law!("P.recon.rate-loop", Product, Product(Recon, Circulate), Formal, Proposed,
        "Rates as exact rationals (num, den): a conversion loop closes iff the product of numerators equals the product of denominators.", ev["rate_loop_closes"], from[], parent Some("A.circulate.exact-zero")),
    law!("P.recon.rate-triangle", Product, Product(Recon, Classify), Formal, Proposed,
        "A conversion loop that does not close is MapConflict (blames the rates), never a fixable record (seat 1: b_rate_triangle_exact blamed a record).", ev["h0_map_conflict"], from[], parent Some("A.classify.h0-policy")),
    law!("P.recon.unit-slip-is-exact", Product, Product(Recon, Classify), Formal, Proposed,
        "A single system's unit slip is a node error (exact), never a loop conflict: circulations are unchanged.", ev["slip_is_exact"], from[], parent Some("A.classify.circ-only")),
    law!("P.recon.minor-units", Product, Product(Recon, Quantize), Property, Proposed,
        "Tolerances are in the currency's minor unit (cents).", ev[], from[], parent Some("S.quantize.floor")),
    law!("P.recon.dated-rate", Product, Product(Recon, Encode), Types, Proposed,
        "Every relation carries its rate's as-of date and one spec uses one as-of date (otherwise a 'loop' is a rate change).", ev["seat3 S3-14"], from[], parent Some("S.encode.valid")),
    law!("P.recon.per-key-join", Product, Product(Recon, Join), Differential, Proposed,
        "K keys batched as one spec give the join of the K per-key verdicts, provided keys of very different scale are not batched together (S1-07).", ev["S1 scale run: 5000 keys 0.49 s"], from[], parent Some("S.join.monoid")),
    law!("P.consensus.pairwise-for-loops", Product, Product(Consensus, Encode), Formal, Proposed,
        "Loop findings require the pairwise-claim encoding (w); absolute beliefs (x) can only yield outlier findings.", ev["encode_mode"], from[], parent Some("S.encode.mode")),
    law!("P.consensus.frame-holonomy", Product, Product(Consensus, Classify), Formal, Proposed,
        "With rotation frames, no consensus in any frame (h0 below expected) is MapConflict, never HOLDS (E3c).", ev["h0_map_conflict"], from[], parent Some("A.classify.h0-policy")),
    law!("P.consensus.outlier-rank", Product, Product(Consensus, Attribute), Property, Proposed,
        "Under identity maps the outlier agent is ranked first by incident discord (E3a).", ev["E3a"], from[], parent Some("A.attribute.ranking")),
    law!("P.consensus.round-diff", Product, Product(Consensus, Diff), Formal, Proposed,
        "The per-round timeline is Diff over consecutive rounds keyed by agent id.", ev["diff_keyed"], from[], parent Some("S.diff.keyed")),
    law!("P.studio.traceable", Product, Product(Studio, Attribute), Property, Proposed,
        "Every rendered number maps to a Finding or Report field path; the page performs no arithmetic.", ev[], from[], parent Some("S.attribute.projection")),
    law!("P.studio.steps", Product, Product(Studio, Explain), Property, Proposed,
        "Classify's rules appear as a numbered vertical list in evaluation order; hovering step i highlights exactly the entries rule i produced.", ev["seat3 S3-19"], from[], parent Some("S.explain.pure")),
    law!("P.studio.pin-is-boundary", Product, Product(Studio, Analyze), Types, Proposed,
        "A pinned node renders extension.values[].boundary = true; its badge reads determined iff unique = Some(true).", ev["seat3 S3-19"], from[], parent Some("S.analyze.extension-unique")),
    law!("P.schema.holonomy-drop", Product, Product(Schema, Analyze), CoreCited, Proposed,
        "A mapping cycle that does not compose to the identity drops H0 (uc5).", ev["holonomy_kills_h0_and_k4_keeps_loop_space"], from[], parent Some("S.analyze.holonomy")),
    law!("P.translation.loop-is-harmonic", Product, Product(Translation, Circulate), Formal, Proposed,
        "Per-hop drift that closes around a translation loop has zero circulation; a loop drift shows only as circulation (uc6).", ev["circ_exact_zero"], from[], parent Some("A.circulate.exact-zero")),
    ]
}

/// Effective status: a law derived (transitively) from a Hole is a Hole.
pub fn effective_status(all: &[Law], id: &str) -> Status {
    let l = all.iter().find(|l| l.id == id).unwrap();
    if l.status == Status::Hole { return Status::Hole; }
    if l.derives_from.iter().any(|d| effective_status(all, d) == Status::Hole) { Status::Hole } else { Status::Proposed }
}

pub fn depth(laws: &[Law], id: &str) -> usize {
    let mut d = 0; let mut cur = id;
    while let Some(p) = laws.iter().find(|l| l.id == cur).and_then(|l| l.parent) { d += 1; cur = p; if d > 8 { break; } }
    d
}

/// For a Formal law: how many of its evidence names are backed at all four layers (LAWS.bend law, PROOF.bend def,
/// Haskell property, Rust test with the same name), out of how many it cites.
pub fn formal_backing(l: &Law, laws_bend: &str, proof_bend: &str, hs: &str, rs: &str) -> (usize, usize) {
    let backed = l.evidence.iter().filter(|e| laws_bend.contains(&format!("law {e}:")) && proof_bend.contains(&format!("def Laws.{e}("))
        && hs.contains(&format!("(\"{e}\"")) && rs.contains(&format!("fn {e}("))).count();
    (backed, l.evidence.len())
}
