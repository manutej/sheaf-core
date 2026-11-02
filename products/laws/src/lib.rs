//! The law ontology of sheaf-products, as typed data.
//!
//! Shape (operadic): the product is ONE composite of many-in/one-out stages (`Op`), wired by typed
//! edges (`Color`). Laws live at three depths:
//!   L0 Global   — about the root composite (the whole pipeline, any product)
//!   L1 Stage    — about one stage (local) or one edge between stages (interface / glue)
//!   L2 Atom     — about one function inside a stage
//!   L3 Product  — a product's refinement of a stage law (leaf; depth cap)
//! A global law must be DERIVABLE from lower laws (`derives_from`); a global law with no derivation, or a
//! derivation through a law marked `Hole`, is a spec gap. That is the operadic consistency check on the
//! spec itself: the law stated directly at the root must agree with the law composed from its parts.
//!
//! Tiers follow the Touchstone oracle ladder; "proved" is only ever said of `Formal`.

pub mod wiring;
pub mod selfsheaf;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Layer { Global = 0, Stage = 1, Atom = 2, Product = 3 }

/// Typed edges of the pipeline (what flows between stages).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Color {
    DomainGraph, // product input (sheaf-kit/domain-v1)
    Spec,        // frontier-sheaf/v1 bytes
    ValidSpec,   // spec the core accepted
    Boundary,    // frontier-sheaf/boundary-v1
    Report,      // frontier-sheaf/report-v1 (f64; core-owned)
    Cert,        // certificate projection of a Report (bools, residuals, integer dims)
    Policy,      // tolerances in domain units
    Circ,        // per-fundamental-cycle circulations of w (product layer, not core)
    Evidence,    // quantised Nat evidence: the only f64 -> Nat crossing
    Verdict,     // domain verdict (consistent | fixable | loop | underdetermined | inconclusive)
    Finding,     // finding-v1 (verdict + attribution)
    Receipt,     // Touchstone receipt (HOLDS | BROKEN | BLOCKED | HOLE, bound to hashes)
}

/// The generators of the operad: each stage is a many-in, one-out typed operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Op { Encode, Validate, Analyze, Extend, Certify, Circulate, Quantize, Classify, Attribute, Receipt, Join }

impl Op {
    pub fn signature(self) -> (&'static [Color], Color) {
        use Color::*;
        match self {
            Op::Encode    => (&[DomainGraph], Spec),
            Op::Validate  => (&[Spec], ValidSpec),
            Op::Analyze   => (&[ValidSpec], Report),
            Op::Extend    => (&[ValidSpec, Boundary], Report),
            Op::Certify   => (&[Report], Cert),
            Op::Circulate => (&[ValidSpec], Circ),
            Op::Quantize  => (&[Report, Circ, Policy], Evidence),
            Op::Classify  => (&[Cert, Evidence], Verdict),
            Op::Attribute => (&[Report, DomainGraph, Verdict], Finding),
            Op::Receipt   => (&[Finding, Cert], Receipt),
            Op::Join      => (&[Receipt, Receipt], Receipt),
        }
    }
    pub fn owner(self) -> Owner {
        match self { Op::Validate | Op::Analyze | Op::Extend => Owner::Core, _ => Owner::Products }
    }
    pub const ALL: [Op; 11] = [Op::Encode, Op::Validate, Op::Analyze, Op::Extend, Op::Certify, Op::Circulate,
                               Op::Quantize, Op::Classify, Op::Attribute, Op::Receipt, Op::Join];
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Owner { Core, Products }

/// Where a law lives in the tree.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Scope {
    Root,                         // the whole composite
    Stage(Op),                    // inside one stage
    Edge(Op, Op),                 // the interface from producer to consumer (glue law)
    Product(Prod, Op),            // a product's refinement of a stage
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Prod { Memory, Recon, Consensus, Studio }

/// Touchstone oracle tier: how the law is enforced. Strongest first.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Tier {
    Formal,       // Bend proof (+ --verdict) AND Haskell properties under the same ID; Rust conforms (bridge B)
    Types,        // enforced by Rust types (newtype / enum / no constructor): unrepresentable otherwise
    CoreCited,    // already tested in sheaf-core; pinned by test name + G0 manifest; not re-proved here
    Differential, // two independent implementations must agree (core binary vs typed Rust oracle); disagreement = HOLE
    Property,     // Haskell QuickCheck + Rust proptest, same ID
}

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
    /// For CoreCited: sheaf-core test names. For Formal: Bend law name. Free-form otherwise.
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
    // ============================================================ L0 GLOBAL (the root composite)
    law!("G0", Global, Root, Types, Proposed,
        "Core invariance: crates/**, Cargo.lock, Cargo.toml and rust-toolchain.toml are byte-identical to the pinned manifest; no product law is ever discharged by editing the core.",
        ev["products/gates/core.sha256"], from["I.analyze.pin"], parent None),
    law!("G1", Global, Root, Property, Proposed,
        "Claim soundness: a receipt other than BLOCKED is issued only when the core run is certified (exit 0, spectral and Hodge converged, cut consistent, not near-zero-ambiguous, report well-formed).",
        ev[], from["S.certify.wf", "A.classify.uncertified", "S.receipt.map", "I.analyze.certify"], parent None),
    law!("G2", Global, Root, Types, Proposed,
        "Traceability: every number in a Finding is a Report field verbatim or a Nat produced by Quantize from one; no f64 arithmetic happens outside the core and Circulate.",
        ev[], from["S.attribute.projection", "S.quantize.only-crossing"], parent None),
    law!("G3", Global, Root, CoreCited, Proposed,
        "Determinism: equal input bytes give equal receipt bytes (excluding the timestamp), independent of thread count and of the order in which components or keys are processed.",
        ev["m1_h0_bitwise_independent_of_thread_count", "v1_triangle_harmonic_share_is_one_third_and_report_is_deterministic"],
        from["S.join.monoid", "S.analyze.determinism"], parent None),
    law!("G4", Global, Root, Formal, Proposed,
        "Gauge invariance: the loop verdict is unchanged by adding any exact cochain (delta p) to the observed edge data w. Consequence: a verdict may not depend on harmonic SHARE (Spec v1 L4-L7 violate this; experiment E2a: real arbitrage at share 2e-5).",
        ev[], from["A.circulate.exact-zero", "A.circulate.linear", "A.classify.circ-only"], parent None),
    law!("G5", Global, Root, Property, Proposed,
        "Clean-twin completeness: on exactly consistent data (x a global section, or w = delta p) every product reports consistent / HOLDS, provided the policy tolerance clears the numerical floor.",
        ev["v1_x_only_is_exact_and_rotation_loop_kills_h0"], from["G4", "S.quantize.floor", "S.analyze.x-only-exact"], parent None),
    law!("G6", Global, Root, Differential, Proposed,
        "Compositionality (parallel associativity): the receipt of a multi-component spec equals the join of the receipts of its components analysed separately; findings partition by component.",
        ev["v1_multi_component_h0_matches_dense_count", "v1_two_scale_components_use_local_zero_cut"],
        from["S.analyze.component-local", "S.join.monoid", "S.join.cover"], parent None),
    law!("G7", Global, Root, Types, Proposed,
        "Tier honesty: 'proved' is attached only to Formal-tier law receipts; every numeric field is worded 'certified residual'; approximate layers never write exact fields (ADR-008, ADR-013).",
        ev[], from["S.receipt.tier"], parent None),

    // ============================================================ L1 STAGES (local) — products layer
    law!("S.encode.valid", Stage, Stage(Encode), Formal, Proposed,
        "Encode is valid by construction: for every well-formed domain graph, the emitted spec passes Validate (endpoints exist, u != v, map shapes de x dim, identity only when dims agree, at least one of x, w).",
        ev["encode_valid"], from["A.encode.endpoints", "A.encode.shapes"], parent Some("G1")),
    law!("S.encode.index", Stage, Stage(Encode), Formal, Proposed,
        "Encode fixes a bijection between domain relations and spec edges that preserves order; Attribute uses its inverse.",
        ev["encode_index_bij"], from[], parent Some("G2")),
    law!("S.encode.mode", Stage, Stage(Encode), Formal, Proposed,
        "Mode honesty: the spec carries w iff every relation has an observation; a partial observation set is rejected, never silently mixed.",
        ev["encode_mode"], from[], parent Some("G5")),
    law!("S.circulate.exact", Stage, Stage(Circulate), Differential, Proposed,
        "Circulate computes, for each fundamental cycle of a spanning forest, the transported sum of w around it; for scalar maps with exact rational rates the computation is exact.",
        ev[], from["A.circulate.exact-zero", "A.circulate.linear", "A.circulate.basis-size"], parent Some("G4")),
    law!("S.quantize.only-crossing", Stage, Stage(Quantize), Types, Proposed,
        "Quantize is the only place an f64 becomes a Nat; Classify accepts Ticks<Tol> only (newtype), never f64.",
        ev[], from[], parent Some("G2")),
    law!("S.quantize.threshold", Stage, Stage(Quantize), Property, Proposed,
        "Threshold faithfulness: ticks(v) >= ticks(tol) iff v >= tol, and ticks is monotone; values are measured relative to the domain tolerance (reports span ~41 decades, absolute fixed point fails).",
        ev[], from[], parent Some("G4")),
    law!("S.quantize.floor", Stage, Stage(Quantize), Property, Proposed,
        "Fail-before-freeze for policies: a tolerance below kappa * max(certified residuals, cert_coclosed) is refused (BLOCKED policy), so numerical noise can never cross the threshold.",
        ev[], from[], parent Some("G5")),
    law!("S.quantize.nonfinite", Stage, Stage(Quantize), Property, Proposed,
        "NaN, inf and null report values quantise to BLOCKED, never to 0.",
        ev[], from[], parent Some("G1")),
    law!("S.certify.wf", Stage, Stage(Certify), Property, Proposed,
        "Report well-formedness: h1 = dim_C1 + h0 - dim_C0 (graphs); exact_share + harmonic_share = 1 within cert; unique = (h0_relative == 0); any violation is BLOCKED(report).",
        ev["rank_nullity_and_euler"], from[], parent Some("G1")),
    law!("S.classify.total", Stage, Stage(Classify), Formal, Proposed,
        "Classify is a total, terminating function to a closed verdict enum.",
        ev["classify_total"], from["A.classify.uncertified", "A.classify.undetermined", "A.classify.xonly", "A.classify.mono", "A.classify.circ-only", "A.classify.loop-iff", "A.classify.disc-iff"], parent Some("G1")),
    law!("S.attribute.projection", Stage, Stage(Attribute), Types, Proposed,
        "Attribute only projects: Finding numbers are built through FromReport<T> constructors; there is no public constructor from f64.",
        ev[], from["A.attribute.loop-blame", "A.attribute.fixable", "A.attribute.fill", "A.attribute.x-only"], parent Some("G2")),
    law!("S.receipt.map", Stage, Stage(Receipt), Formal, Proposed,
        "Verdict map: consistent -> HOLDS; fixable or loop -> BROKEN with the finding; inconclusive or underdetermined -> BLOCKED with one question; two oracles disagreeing -> HOLE. Never averaged.",
        ev["receipt_map_hole", "receipt_map_uncertified"], from[], parent Some("G1")),
    law!("S.receipt.binding", Stage, Stage(Receipt), Property, Proposed,
        "A receipt is bound to sha256 of the core manifest, spec bytes, report bytes and interface.lock; changing any one changes the receipt hash.",
        ev[], from[], parent Some("G3")),
    law!("S.receipt.tier", Stage, Stage(Receipt), Types, Proposed,
        "Receipt.tier is an enum; the 'proved' label is a method only on Tier::Formal.",
        ev[], from[], parent Some("G7")),
    law!("S.join.monoid", Stage, Stage(Join), Formal, Proposed,
        "Join of receipts is a bounded join-semilattice on HOLDS < BROKEN < BLOCKED < HOLE: associative, commutative, idempotent, HOLDS is the identity.",
        ev["join_assoc", "join_comm", "join_idem", "join_unit_l", "join_unit_r", "join_hole_absorbs"], from[], parent Some("G6")),
    law!("S.join.cover", Stage, Stage(Join), Formal, Proposed,
        "A system receipt is a join over a declared cover; a missing component contributes HOLE (a child can never claim system-level HOLDS).",
        ev["join_cover_missing_is_hole"], from[], parent Some("G6")),

    // ============================================================ L1 STAGES — core-owned (cited, not re-proved)
    law!("S.analyze.complex", Stage, Stage(Analyze), CoreCited, Proposed,
        "delta delta = 0; a twisted map on a filled face is rejected as not a complex.",
        ev["d_squared_is_zero", "twisted_map_on_filled_face_is_rejected_as_not_a_complex"], from[], parent Some("G1")),
    law!("S.analyze.hodge", Stage, Stage(Analyze), CoreCited, Proposed,
        "Hodge split: exact and harmonic parts are orthogonal and sum to w.",
        ev["hodge_components_orthogonal_and_sum"], from[], parent Some("G4")),
    law!("S.analyze.x-only-exact", Stage, Stage(Analyze), CoreCited, Proposed,
        "With node data only, w := delta x is exact: harmonic share is zero.",
        ev["v1_x_only_is_exact_and_rotation_loop_kills_h0"], from[], parent Some("G5")),
    law!("S.analyze.holonomy", Stage, Stage(Analyze), CoreCited, Proposed,
        "Holonomy around an unfilled loop removes global sections; K4 keeps its loop space.",
        ev["holonomy_kills_h0_and_k4_keeps_loop_space"], from[], parent Some("G4")),
    law!("S.analyze.component-local", Stage, Stage(Analyze), CoreCited, Proposed,
        "H0 cut and spectra are per connected component (L0 block diagonal).",
        ev["v1_multi_component_h0_matches_dense_count", "v1_two_scale_components_use_local_zero_cut"], from[], parent Some("G6")),
    law!("S.analyze.determinism", Stage, Stage(Analyze), CoreCited, Proposed,
        "H0 is bitwise independent of thread count; the report is deterministic.",
        ev["m1_h0_bitwise_independent_of_thread_count"], from[], parent Some("G3")),
    law!("S.validate.mirror", Stage, Stage(Validate), Differential, Proposed,
        "The product-side validator accepts exactly what the core parser accepts, on generated and mutated specs.",
        ev["v1_rejects_bad_specs", "check_edges_refuses_each_bad_edge", "adversarial_overflowing_entries_refused_not_panicking", "adversarial_tiny_spec_huge_dims_refused_before_alloc"], from[], parent Some("G1")),
    law!("S.extend.pins", Stage, Stage(Extend), Differential, Proposed,
        "Pinned boundary values come back bit-identical in the extension output.",
        ev["pinning_kills_sections_and_detects_unpinned_parts"], from[], parent Some("G2")),
    law!("S.extend.unique", Stage, Stage(Extend), CoreCited, Proposed,
        "unique iff dim H0(X,B;F) = 0; an unpinned component is reported not unique.",
        ev["unpinned_component_is_reported_not_unique_and_bad_boundary_rejected", "empty_boundary_equals_absolute_h0_on_the_catalogue"], from[], parent Some("G1")),
    law!("S.extend.minimal", Stage, Stage(Extend), CoreCited, Proposed,
        "The extension minimises energy with the boundary fixed; on a path with constant sheaf it is linear interpolation.",
        ev["extension_minimises_energy", "constant_sheaf_on_path_is_linear_interpolation", "path_extension_is_linear_and_unique"], from[], parent Some("G5")),

    // ============================================================ L1 EDGES (glue / interface laws)
    law!("I.encode.validate", Stage, Edge(Encode, Validate), Differential, Proposed,
        "Every encoder output is accepted by the core binary (exit 0) on the product's fixture set and on generated domain graphs.",
        ev[], from["S.encode.valid", "S.validate.mirror"], parent Some("G1")),
    law!("I.encode.attribute", Stage, Edge(Encode, Attribute), Formal, Proposed,
        "Round trip: Attribute's edge map composed with Encode's index is the identity on relations.",
        ev["index_roundtrip"], from["S.encode.index"], parent Some("G2")),
    law!("I.analyze.pin", Stage, Edge(Analyze, Certify), Types, Proposed,
        "The report binding is generated from report.rs at the pinned commit; a missing or renamed field fails to parse and yields BLOCKED(report-schema).",
        ev[], from[], parent Some("G0")),
    law!("I.analyze.certify", Stage, Edge(Analyze, Certify), Types, Proposed,
        "Certify reads only the certificate fields (converged, cut_consistent, near_zero_ambiguous, residuals, unique, h0_relative) and copies them unchanged.",
        ev[], from["I.analyze.pin"], parent Some("G1")),
    law!("I.circulate.analyze", Stage, Edge(Circulate, Quantize), Differential, Proposed,
        "Two-oracle agreement: all fundamental circulations below tolerance iff the core's harmonic energy is below tolerance squared times the cycle count bound; disagreement is a HOLE, never averaged (F5).",
        ev[], from["S.circulate.exact", "S.analyze.hodge"], parent Some("G4")),
    law!("I.quantize.classify", Stage, Edge(Quantize, Classify), Types, Proposed,
        "Evidence and tolerance share one unit: Ticks<Tol> with the same Tol; mixing units does not type-check.",
        ev[], from[], parent Some("G2")),
    law!("I.classify.receipt", Stage, Edge(Attribute, Receipt), Formal, Proposed,
        "The receipt verdict is a function of the classify verdict and the oracle-agreement bit only.",
        ev["receipt_map_hole"], from["S.receipt.map"], parent Some("G1")),

    // ============================================================ L2 ATOMS
    law!("A.encode.endpoints", Atom, Stage(Encode), Formal, Proposed,
        "Every encoded edge has u < n, v < n and u != v.", ev["encode_ends_ok"], from[], parent Some("S.encode.valid")),
    law!("A.encode.shapes", Atom, Stage(Encode), Formal, Proposed,
        "Every encoded map has de rows and dim(endpoint) columns; omitted maps only when dims agree.", ev["encode_shapes_ok"], from[], parent Some("S.encode.valid")),
    law!("A.circulate.exact-zero", Atom, Stage(Circulate), Formal, Proposed,
        "Telescoping: the circulation of delta p around any closed walk is zero (stated over Nat as: sum of heads = sum of tails).",
        ev["circ_exact_zero"], from[], parent Some("S.circulate.exact")),
    law!("A.circulate.linear", Atom, Stage(Circulate), Formal, Proposed,
        "circ(w1 + w2) = circ(w1) + circ(w2) on every cycle.", ev["circ_linear"], from[], parent Some("S.circulate.exact")),
    law!("A.circulate.basis-size", Atom, Stage(Circulate), Differential, Proposed,
        "The fundamental cycle count equals |E| - |V| + #components, which equals the core's h1 for an injective constant sheaf.",
        ev[], from[], parent Some("S.circulate.exact")),
    law!("A.classify.uncertified", Atom, Stage(Classify), Formal, Proposed,
        "Not certified -> inconclusive, whatever the evidence.", ev["sound_uncertified"], from[], parent Some("S.classify.total")),
    law!("A.classify.undetermined", Atom, Stage(Classify), Formal, Proposed,
        "Certified, boundary given, not unique -> underdetermined.", ev["undetermined"], from[], parent Some("S.classify.total")),
    law!("A.classify.xonly", Atom, Stage(Classify), Formal, Proposed,
        "x-only evidence never yields a loop conflict.", ev["xonly_no_loop"], from[], parent Some("S.classify.total")),
    law!("A.classify.mono", Atom, Stage(Classify), Formal, Proposed,
        "Raising the loop tolerance never turns consistent into loop conflict.", ev["mono_tolerance"], from[], parent Some("S.classify.total")),
    law!("A.classify.circ-only", Atom, Stage(Classify), Formal, Proposed,
        "The loop decision reads the maximum circulation tick only; two evidences with equal max circulation get the same loop decision.",
        ev["circ_only"], from[], parent Some("S.classify.total")),
    law!("A.classify.loop-iff", Atom, Stage(Classify), Formal, Proposed,
        "Certified and bounded: the verdict is LoopConflict exactly when the input is not x-only and the max circulation clears its tolerance.", ev["loop_iff_circ"], from[], parent Some("S.classify.total")),
    law!("A.classify.disc-iff", Atom, Stage(Classify), Formal, Proposed,
        "Certified and bounded: the verdict is FixableNodes exactly when discord clears its tolerance and no loop fires. Added after mutation testing: the swapped discord comparison survived every earlier law.", ev["fixable_iff_disc"], from[], parent Some("S.classify.total")),
    law!("A.attribute.loop-blame", Atom, Stage(Attribute), Types, Proposed,
        "A loop finding names a cycle (two or more edges, closed walk), never a single edge as the cause.", ev[], from[], parent Some("S.attribute.projection")),
    law!("A.attribute.fixable", Atom, Stage(Attribute), Formal, Proposed,
        "Every node named fixable has incident-discord ticks >= the discord tolerance.", ev["fixable_above_tol"], from[], parent Some("S.attribute.projection")),
    law!("A.attribute.fill", Atom, Stage(Attribute), Formal, Proposed,
        "A fill is marked determined only when the extension is unique; write-back requires determined.", ev["fill_determined"], from[], parent Some("S.attribute.projection")),
    law!("A.attribute.x-only", Atom, Stage(Attribute), Types, Proposed,
        "For x-only input the loops field is NotApplicable, never an empty list (an empty list would read as a negative finding).", ev[], from[], parent Some("S.attribute.projection")),

    // ============================================================ L3 PRODUCT REFINEMENTS (leaves)
    law!("P.memory.flag-iff-differs", Product, Product(Memory, Encode), Formal, Proposed,
        "Status in quarter units with gain 4: discord is (level_u - level_v)^2 exactly; a supports-edge flags iff the statuses differ (Spec v1 L8).", ev["flag_iff_status_differs"], from[], parent Some("A.encode.shapes")),
    law!("P.memory.claim-rule", Product, Product(Memory, Encode), Formal, Proposed,
        "claim_of edge flags iff assertiveness exceeds the source's status (Spec v1 L9); equal status never flags (L10).", ev["claim_flag_rule", "equal_status_never_flags"], from[], parent Some("A.encode.shapes")),
    law!("P.memory.contradiction-parity", Product, Product(Memory, Circulate), Formal, Proposed,
        "A cycle of contradicts (-1) maps admits a nonzero consistent assignment iff it has an even number of contradicts edges (E1b: odd cycle -> H0 = 0).", ev["contradiction_parity"], from[], parent Some("A.circulate.exact-zero")),
    law!("P.memory.writeback", Product, Product(Memory, Attribute), Formal, Proposed,
        "Status write-back happens only for fills marked determined.", ev["fill_determined"], from[], parent Some("A.attribute.fill")),
    law!("P.recon.rate-loop", Product, Product(Recon, Circulate), Formal, Proposed,
        "Rates as exact rationals (num, den): a conversion loop closes iff the product of numerators equals the product of denominators.", ev["rate_loop_closes"], from[], parent Some("A.circulate.exact-zero")),
    law!("P.recon.unit-slip-is-exact", Product, Product(Recon, Classify), Formal, Proposed,
        "A single system's unit slip is a node error (exact), never a loop conflict: it shifts x, so circulations are unchanged.", ev["slip_is_exact"], from[], parent Some("A.classify.circ-only")),
    law!("P.recon.minor-units", Product, Product(Recon, Quantize), Property, Proposed,
        "Tolerances are in the currency's minor unit (cents); quantisation is relative to that unit.", ev[], from[], parent Some("S.quantize.threshold")),
    law!("P.recon.per-key-join", Product, Product(Recon, Join), Differential, Proposed,
        "A batch of K keys analysed as one spec gives the join of the K per-key receipts (S1: 5000 keys, 0.49 s).", ev[], from[], parent Some("S.join.monoid")),
    law!("P.consensus.pairwise-for-loops", Product, Product(Consensus, Encode), Formal, Proposed,
        "Loop findings require the pairwise-claim encoding (w); absolute beliefs (x) can only yield outlier findings.", ev["encode_mode"], from[], parent Some("S.encode.mode")),
    law!("P.consensus.frame-holonomy", Product, Product(Consensus, Analyze), CoreCited, Proposed,
        "With rotation frames, consensus exists in some frame iff the loop holonomy fixes a vector (E3b/c).", ev["holonomy_kills_h0_and_k4_keeps_loop_space"], from[], parent Some("S.analyze.holonomy")),
    law!("P.studio.traceable", Product, Product(Studio, Attribute), Property, Proposed,
        "Every number rendered in the Studio maps to a report field path; the page performs no arithmetic.", ev[], from[], parent Some("S.attribute.projection")),
    ]
}

/// Holes found by checking each global law against its derivation (the spec-level consistency pass).
pub fn holes() -> Vec<(&'static str, &'static str)> {
    vec![
        ("H1 G5 needs S.quantize.floor", "Neither Spec v1 nor SPEC.md bound the tolerance from below; on clean data the core returns harmonic share ~1e-33, so a tolerance of 0 makes clean data a loop conflict. The floor law closes it."),
        ("H2 G4 contradicts Spec v1 L4-L7", "Spec v1's proven verdict uses material SHARES (>5%). Adding an exact cochain changes shares, so the verdict is not gauge invariant; E2a shows a real arbitrage at share 2e-5. L4-L7 must be replaced by A.classify.circ-only, and their Bend proofs redone."),
        ("H3 G2 had no crossing stage", "Spec v1 computes per-mille shares in the Python tool and hands them to Bend: the float-to-Nat crossing had no law. Quantize is now an explicit stage with threshold, floor and non-finite laws."),
        ("H4 G6 was only cited", "Component locality was covered by a dense-count core test only; batching keys into one spec (P2) relies on it, so it now needs a differential law (batch vs split)."),
        ("H5 loop localisation", "Harmonic energy spreads over 51 of 52 edges (S2); localisation must come from circulations (Circulate), cross-checked against the core (I.circulate.analyze), with disagreement surfacing as a HOLE."),
    ]
}

pub fn depth(laws: &[Law], id: &str) -> usize {
    let mut d = 0; let mut cur = id;
    while let Some(p) = laws.iter().find(|l| l.id == cur).and_then(|l| l.parent) { d += 1; cur = p; if d > 8 { break; } }
    d
}
