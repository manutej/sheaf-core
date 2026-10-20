"""Pydantic v2 models matching sheaf.schema.json (v1.18). Semantic validators enforce SPEC MUSTs the JSON Schema cannot express.
Strict mode: no str->number / bool->int coercion; JSON integers ARE valid numbers (strict float accepts int), as in the
schema ("number") and Rust (f64). Field patterns are compiled by pydantic-core's rust-regex engine (pinned explicitly with
regex_engine="rust-regex" on Node): without the (?m) flag, `$` matches ONLY at the absolute end of the string, never before a
trailing newline (unlike Python `re`, where `$` also matches just before one final newline). v1.11 (R10 api-typing): this is machine-checked, not asserted:
tests/test_spec_public.py::test_R_ID_01_pydantic_patterns_reject_trailing_newline runs the engine directly on every id/parent
tail case (with the non-claim length guard out of the way), and mutation rows claim_id_parent_plus_newline,
claim_id_suffix_trailing_newline, parent_trailing_newline_on_topic, claim_parent_trailing_newline pin all three mirrors."""
from typing import Literal, Optional, Annotated
from pydantic import BaseModel, ConfigDict, Field, model_validator, ValidationError
from pydantic_core import PydanticCustomError
import math, re
EdgeType = Literal["about", "related", "logged", "claim_of", "same_slug"]
Status = Literal["verified", "partly verified", "search-verified", "claimed", "unverified", "refuted", "not found"]

class Node(BaseModel):
    model_config = ConfigDict(extra="allow", strict=True, regex_engine="rust-regex")   # v1.11: `$` = absolute end (see module doc)
    id: str = Field(pattern=r"^[0-9a-f]{32}(c[0-9]+)?$")
    kind: Literal["source", "topic", "claim", "log"]
    slug: Optional[str] = None; title: Optional[str] = None; status: Optional[str] = None
    url: Optional[str] = None; parent: Optional[str] = Field(default=None, pattern=r"^[0-9a-f]{32}$"); fetched: Optional[str] = None
    @model_validator(mode="after")
    def provenance(self):
        if self.kind != "claim" and len(self.id) != 32: raise PydanticCustomError("R-ID-01", f"R-ID-01: c<k> suffix only on claims ({self.id})")  # M:id_suffix_claims_only
        if self.kind == "source":
            if not self.slug: raise PydanticCustomError("R-PROV-04", f"R-PROV-04: source {self.id} lacks slug")  # M:src_slug
            if not self.url or not self.url.startswith(("http://", "https://")): raise PydanticCustomError("R-PROV-04", f"R-PROV-04: source {self.id} lacks http(s) url")  # M:src_url
            if self.status not in Status.__args__: raise PydanticCustomError("R-STATUS-05", f"R-STATUS-05: {self.slug} status {self.status!r}")  # M:src_status
        if self.kind == "claim" and not self.parent: raise PydanticCustomError("R-CLAIM-01", f"R-CLAIM-01: claim {self.id} has no parent")  # early fast-fail for a MISSING parent only (untagged); NOT equivalent to M:claim_parent_is_source, which also rejects a topic/log/unknown parent and must stay
        return self

class Edge(BaseModel):
    model_config = ConfigDict(extra="forbid", strict=True)
    u: int = Field(ge=0); v: int = Field(ge=0); type: EdgeType

class SheafSpec(BaseModel):
    model_config = ConfigDict(extra="forbid", strict=True)
    d: int = Field(ge=1, le=4096); nodes: list[Node] = Field(min_length=1); x: list[list[float]]; edges: list[Edge]
    maps: dict[EdgeType, Annotated[list[list[list[float]]], Field(min_length=2, max_length=2)]]
    @model_validator(mode="after")
    def invariants(self):
        n, d = len(self.nodes), self.d
        if len(self.x) != n: raise PydanticCustomError("R-DIM-02", "R-DIM-02: x must have n rows")  # M:x_rows
        if any(len(r) != d for r in self.x): raise PydanticCustomError("R-DIM-02", "R-DIM-02: x rows must have d entries")  # M:x_cols
        if not all(math.isfinite(v) for r in self.x for v in r): raise PydanticCustomError("R-FINITE-01", "R-FINITE-01: x")  # M:x_finite
        for t, (Fu, Fv) in self.maps.items():
            for F in (Fu, Fv):
                if len(F) != d: raise PydanticCustomError("R-DIM-02", f"R-DIM-02: map {t} must have {d} rows")  # M:map_rows
                if any(len(r) != d for r in F): raise PydanticCustomError("R-DIM-02", f"R-DIM-02: map {t} rows must have {d} entries")  # M:map_cols
                if not all(math.isfinite(v) for r in F for v in r): raise PydanticCustomError("R-FINITE-01", f"R-FINITE-01: map {t}")  # M:map_finite
        ids = [nd.id for nd in self.nodes]
        if len(set(ids)) != n: raise PydanticCustomError("R-ID-01", "R-ID-01: duplicate node ids")  # M:id_unique
        for j, e in enumerate(self.edges):
            if not e.u < e.v: raise PydanticCustomError("R-ORIENT-01", f"R-ORIENT-01: edge {j} needs u<v")  # M:orient_u_lt_v
            if not e.v < n: raise PydanticCustomError("R-ORIENT-01", f"R-ORIENT-01: edge {j} needs v<n")  # M:orient_v_lt_n
            if e.type not in self.maps: raise PydanticCustomError("R-TYPE-03", f"R-TYPE-03: edge {j} type {e.type} has no map")  # M:type_has_map
        I = [[1.0 if a == b else 0.0 for b in range(d)] for a in range(d)]
        for t, (Fu, Fv) in self.maps.items():
            if Fu == I and Fv == I: raise PydanticCustomError("R-NONID-04", f"R-NONID-04: map pair for {t} is identity/identity")  # M:nonid
        idx = {nd.id: i for i, nd in enumerate(self.nodes)}
        cof = {}
        for e in self.edges:
            if e.type == "claim_of": cof.setdefault(e.v, []).append(e.u)
        for i, nd in enumerate(self.nodes):
            if nd.kind != "claim": continue
            p = idx.get(nd.parent)
            if p is None or self.nodes[p].kind != "source": raise PydanticCustomError("R-CLAIM-01", f"R-CLAIM-01: claim {nd.id} parent is not a source node")  # M:claim_parent_is_source
            if cof.get(i) != [p]: raise PydanticCustomError("R-CLAIM-01", f"R-CLAIM-01: claim {nd.id} needs exactly one claim_of edge from its parent (got {cof.get(i)})")  # M:claim_exactly_one_from_parent
            if not nd.id.startswith(nd.parent): raise PydanticCustomError("R-ID-01", f"R-ID-01: claim id {nd.id} must be its parent id + c<k>")  # M:claim_id_prefix
        for v, us in cof.items():
            if self.nodes[v].kind != "claim": raise PydanticCustomError("R-CLAIM-01", f"R-CLAIM-01: claim_of edge into non-claim node {v}")  # M:claim_of_only_into_claims
        return self


# ---- R-SCHEMA-01 error contract: every rejection message starts with the full requirement ID ----
FIELD_REQ = [(r"^nodes\.\d+\.id$", "R-ID-01"), (r"^nodes\.\d+\.kind$", "R-SCHEMA-01"), (r"^nodes\.\d+\.(url|slug)$", "R-PROV-04"),
             (r"^nodes\.\d+\.status$", "R-STATUS-05"), (r"^nodes\.\d+\.parent$", "R-CLAIM-01"), (r"^edges\.\d+\.(u|v)$", "R-ORIENT-01"),
             (r"^edges\.\d+\.type$", "R-TYPE-03"), (r"^maps\.[^.]+\.\[key\]$", "R-TYPE-03"), (r"^maps\.", "R-DIM-02"), (r"^(d|x)(\.|$)", "R-DIM-02")]
class SpecError(ValueError):
    """Raised by validate_spec; .messages is a list of strings, each beginning with a requirement ID."""
    def __init__(self, messages): self.messages = messages; super().__init__("\n".join(messages))
_SRC_ONLY = re.compile(r"^nodes\.(\d+)\.(url|slug|status)$")
def _rid(err, obj=None):
    if err["type"].startswith("R-"): return err["type"]
    # v1.15 (R14 api-typing): url / slug / status are R-PROV-04 / R-STATUS-05 fields of SOURCE nodes only; the same type
    # error on a topic / claim / log node is a plain structural R-SCHEMA-01 error (like title / fetched)
    m = _SRC_ONLY.match(".".join(str(x) for x in err["loc"]))
    if m and obj is not None:
        try: kind = obj["nodes"][int(m.group(1))].get("kind")
        except Exception: kind = None
        if kind != "source": return "R-SCHEMA-01"
    # v1.7: components are sanitised before joining, so a map key containing '.' (or an empty key) cannot shift the
    # location pattern (an unknown key such as 'cites.x' must still map to R-TYPE-03, not R-DIM-02)
    loc = ".".join(str(x).replace(".", "\u2024") or "\u2205" for x in err["loc"])
    for pat, rid in FIELD_REQ:
        if re.search(pat, loc): return rid
    return "R-SCHEMA-01"
def validate_spec(obj) -> "SheafSpec":
    """Canonical semantic validator (R-SCHEMA-01). Returns the model or raises SpecError(['R-XX-NN: ...', ...])."""
    try: return SheafSpec.model_validate(obj)
    except ValidationError as e:
        msgs = []
        for err in e.errors():
            rid = _rid(err, obj); msg = err["msg"]; loc = ".".join(str(x) for x in err["loc"])
            msgs.append(msg if msg.startswith(rid + ":") else f"{rid}: {loc}: {msg}")
        raise SpecError(msgs) from None

def schema_validator(schema):
    """R-SCHEMA-01 structural check used by every mirror test (v1.6): JSON Schema Draft 2020-12 with ONE declared dialect
    rule — `integer` means a JSON integer literal (Python int, not bool), so 2.0 is rejected exactly as pydantic strict int
    and serde usize reject it. (Plain Draft 2020-12 treats 2.0 as an integer; that alone is not the contract's verdict.)"""
    import jsonschema
    tc = jsonschema.Draft202012Validator.TYPE_CHECKER.redefine("integer", lambda _c, v: isinstance(v, int) and not isinstance(v, bool))
    return jsonschema.validators.extend(jsonschema.Draft202012Validator, type_checker=tc)(schema)
