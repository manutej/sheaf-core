"""[sheaf-core copy, paths retargeted to fixtures/m1/; logic unchanged from frontier contract/tests/gen_mutations.py]
Owner/Eval tool (SPEC v1.18; not run by the gate): regenerates tests/mutations.json, the SHARED rejection table consumed by
pytest (JSON Schema + models.validate_spec) and by cargo test (Rust validate()). Indices are resolved once against
the frozen, hash-locked spec.json, so the table is deterministic. Each case: ops (set/remove/append at a JSON
pointer), expected models result (requirement ID or "ok"), schema verdict ("accept"/"reject"; schema checks
structure only), and Rust result (requirement ID, "ok", or "serde" when serde typing rejects before validate())."""
import json, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[2]; S = json.loads((ROOT / "fixtures/m1/spec.json").read_text())  # sheaf-core: synthetic spec
N = S["nodes"]; E = S["edges"]; d = S["d"]; n = len(N)
src = [i for i, x in enumerate(N) if x["kind"] == "source"]; top = [i for i, x in enumerate(N) if x["kind"] == "topic"]
clm = [i for i, x in enumerate(N) if x["kind"] == "claim"]; logn = [i for i, x in enumerate(N) if x["kind"] == "log"]
jc = next(j for j, e in enumerate(E) if e["type"] == "claim_of"); ec = E[jc]; c0, p0 = ec["v"], ec["u"]
other_src = next(i for i in src if i != p0 and i < c0)
jn = next(j for j, e in enumerate(E) if e["type"] != "claim_of")
I = [[1.0 if a == b else 0.0 for b in range(d)] for a in range(d)]
T = []
def case(name, ops, models, schema, rust): T.append(dict(name=name, ops=ops, models=models, schema=schema, rust=rust))
S_ = lambda p, v: {"op": "set", "path": p, "value": v}; R_ = lambda p: {"op": "remove", "path": p}; A_ = lambda p, v: {"op": "append", "path": p, "value": v}
hexid = N[src[0]]["id"]
# --- valid documents (acceptance side: all three mirrors must accept) ---
case("ok_unchanged", [], "ok", "accept", "ok")
case("ok_int_in_x", [S_("/x/0/0", 1)], "ok", "accept", "ok")
case("ok_int_in_map", [S_("/maps/about/0/0/0", 0)], "ok", "accept", "ok")
# v1.15 (R14 api-typing): url / slug / status type errors on NON-source nodes are structural R-SCHEMA-01 (not R-PROV-04 / R-STATUS-05)
case("topic_slug_number", [S_(f"/nodes/{top[0]}/slug", 5)], "R-SCHEMA-01", "reject", "serde")
case("claim_url_number", [S_(f"/nodes/{clm[0]}/url", 5)], "R-SCHEMA-01", "reject", "serde")
case("log_status_number", [S_(f"/nodes/{logn[0]}/status", 5)], "R-SCHEMA-01", "reject", "serde")
case("ok_null_parent_on_topic", [S_(f"/nodes/{top[0]}/parent", None)], "ok", "accept", "ok")
case("ok_topic_optional_fields_absent", [R_(f"/nodes/{top[0]}/slug"), R_(f"/nodes/{top[0]}/title")], "ok", "accept", "ok")
case("ok_source_fetched_null_title_absent", [S_(f"/nodes/{src[0]}/fetched", None), R_(f"/nodes/{src[0]}/title")], "ok", "accept", "ok")
case("ok_claim_title_slug_absent", [R_(f"/nodes/{c0}/title"), R_(f"/nodes/{c0}/slug")], "ok", "accept", "ok")
case("ok_unknown_node_field", [S_(f"/nodes/{top[0]}/note", "x")], "ok", "accept", "ok")
# v1.7 (R6 test-depth): allow-side of the status enum (each legal value), hex parent on a non-claim, bool map entry
for st_ in ["verified", "partly verified", "search-verified", "claimed", "unverified", "not found", "refuted"]:
    case("ok_status_" + st_.replace(" ", "_").replace("-", "_"), [S_(f"/nodes/{src[0]}/status", st_)], "ok", "accept", "ok")
case("ok_hex_parent_on_topic", [S_(f"/nodes/{top[0]}/parent", N[src[0]]["id"])], "ok", "accept", "ok")
case("claim_id_bare_hex", [S_(f"/nodes/{clm[0]}/id", "0" * 32)], "R-ID-01", "reject", "R-ID-01")          # v1.8: schema requires the c<k> suffix on claims
# v1.10 (R9 test-depth): claim-id near misses built from the claim's OWN parent id (prefix intact), so only the id grammar can reject them
case("claim_id_31hex_c1", [S_(f"/nodes/{c0}/id", N[p0]["id"][:31] + "c1")], "R-ID-01", "reject", "R-ID-01")
case("claim_id_c_no_digits", [S_(f"/nodes/{c0}/id", N[p0]["id"] + "c")], "R-ID-01", "reject", "R-ID-01")
case("claim_id_space_in_suffix", [S_(f"/nodes/{c0}/id", N[p0]["id"] + "c 1")], "R-ID-01", "reject", "R-ID-01")
case("ok_duplicate_source_slug", [S_(f"/nodes/{src[1]}/slug", N[src[0]]["slug"])], "ok", "accept", "ok")   # v1.8: R-PROV-04 allow side
case("ok_http_url", [S_(f"/nodes/{src[0]}/url", "http://example.org/paper")], "ok", "accept", "ok")          # v1.8: http allowed
for t in S["maps"]:
    case(f"ok_one_sided_identity_{t}", [S_(f"/maps/{t}/0", I)], "ok", "accept", "ok")
case("d_4096_in_range_but_x_shape", [S_("/d", 4096)], "R-DIM-02", "accept", "R-DIM-02")
# --- R-SCHEMA-01 structure ---
case("extra_root", [S_("/extra", 1)], "R-SCHEMA-01", "reject", "serde")
case("edge_extra_field", [S_(f"/edges/{jn}/w", 1)], "R-SCHEMA-01", "reject", "serde")
case("bad_kind", [S_(f"/nodes/{top[0]}/kind", "note")], "R-SCHEMA-01", "reject", "serde")
case("title_number", [S_(f"/nodes/{top[0]}/title", 5)], "R-SCHEMA-01", "reject", "serde")
case("fetched_number", [S_(f"/nodes/{src[0]}/fetched", 200)], "R-SCHEMA-01", "reject", "serde")
case("nodes_empty", [S_("/nodes", []), S_("/x", []), S_("/edges", [])], "R-SCHEMA-01", "reject", "R-SCHEMA-01")
# --- R-ID-01 ---
case("id_not_hex", [S_("/nodes/0/id", "XYZ")], "R-ID-01", "reject", "R-ID-01")
case("id_uppercase", [S_(f"/nodes/{src[0]}/id", hexid.upper())], "R-ID-01", "reject", "R-ID-01")
case("id_trailing_newline", [S_(f"/nodes/{top[0]}/id", N[top[0]]["id"] + "\n")], "R-ID-01", "reject", "R-ID-01")
case("id_leading_space", [S_(f"/nodes/{top[0]}/id", " " + N[top[0]]["id"][1:])], "R-ID-01", "reject", "R-ID-01")
case("id_duplicate", [S_("/nodes/1/id", N[0]["id"])], "R-ID-01", "accept", "R-ID-01")
case("id_suffix_on_topic", [S_(f"/nodes/{top[0]}/id", N[top[0]]["id"] + "c1")], "R-ID-01", "reject", "R-ID-01")
case("id_suffix_on_source", [S_(f"/nodes/{src[4]}/id", N[src[4]]["id"] + "c0")], "R-ID-01", "reject", "R-ID-01")
case("claim_id_not_parent_prefix", [S_(f"/nodes/{c0}/id", N[other_src]["id"] + "c99")], "R-ID-01", "accept", "R-ID-01")
# --- R-DIM-02 ---
case("d_zero", [S_("/d", 0)], "R-DIM-02", "reject", "R-DIM-02")
case("d_4097", [S_("/d", 4097)], "R-DIM-02", "reject", "R-DIM-02")
case("d_zero_consistent_shapes", [S_("/d", 0), S_("/x", [[] for _ in range(n)]), S_("/maps", {t: [[], []] for t in S["maps"]})], "R-DIM-02", "reject", "R-DIM-02")
case("d_4097_consistent_shapes", [S_("/d", 4097), S_("/nodes", [N[top[0]]]), S_("/x", [[0.0] * 4097]), S_("/edges", []), S_("/maps", {})], "R-DIM-02", "reject", "R-DIM-02")
case("ok_d_4096_consistent_shapes", [S_("/d", 4096), S_("/nodes", [N[top[0]]]), S_("/x", [[0.5] * 4096]), S_("/edges", []), S_("/maps", {})], "ok", "accept", "ok")
case("d_string", [S_("/d", str(d))], "R-DIM-02", "reject", "serde")
case("d_bool", [S_("/d", True)], "R-DIM-02", "reject", "serde")
case("x_row_long", [A_("/x/0", 0.0)], "R-DIM-02", "accept", "R-DIM-02")
case("x_row_short", [R_("/x/5/0")], "R-DIM-02", "accept", "R-DIM-02")
case("x_missing_row", [R_(f"/x/{n-1}")], "R-DIM-02", "accept", "R-DIM-02")
case("x_string", [S_("/x/0/0", "0.5")], "R-DIM-02", "reject", "serde")
case("x_bool", [S_("/x/0/0", True)], "R-DIM-02", "reject", "serde")
case("map_missing_row", [R_("/maps/related/0/3")], "R-DIM-02", "accept", "R-DIM-02")
case("map_short_row", [R_("/maps/logged/1/2/0")], "R-DIM-02", "accept", "R-DIM-02")
case("map_pair_of_three", [A_("/maps/about", I)], "R-DIM-02", "reject", "serde")
case("map_entry_string", [S_("/maps/about/0/0/0", "1")], "R-DIM-02", "reject", "serde")
case("map_entry_bool", [S_("/maps/about/0/0/0", True)], "R-DIM-02", "reject", "serde")
# --- R-TYPE-03 ---
case("edge_type_unknown", [S_(f"/edges/{jn}/type", "cites")], "R-TYPE-03", "reject", "serde")
case("map_key_unknown", [S_("/maps/cites", S["maps"]["about"])], "R-TYPE-03", "reject", "serde")
case("map_key_unknown_dotted", [S_("/maps/cites.x", S["maps"]["about"])], "R-TYPE-03", "reject", "serde")   # v1.7 (R6 api-typing)
case("map_key_unknown_empty", [S_("/maps/", S["maps"]["about"])], "R-TYPE-03", "reject", "serde")
for t in sorted({e["type"] for e in E}):
    case(f"edge_type_without_map_{t}", [R_(f"/maps/{t}")], "R-TYPE-03", "accept", "R-TYPE-03")
# --- R-NONID-04: every edge type ---
for t in S["maps"]:
    case(f"identity_pair_{t}", [S_(f"/maps/{t}", [I, I])], "R-NONID-04", "accept", "R-NONID-04")
# --- R-PROV-04 / R-STATUS-05 ---
case("src_url_null", [S_(f"/nodes/{src[0]}/url", None)], "R-PROV-04", "reject", "R-PROV-04")
case("src_url_ftp", [S_(f"/nodes/{src[1]}/url", "ftp://x")], "R-PROV-04", "reject", "R-PROV-04")
case("src_url_missing", [R_(f"/nodes/{src[2]}/url")], "R-PROV-04", "reject", "R-PROV-04")
case("src_slug_empty", [S_(f"/nodes/{src[0]}/slug", "")], "R-PROV-04", "reject", "R-PROV-04")
case("src_slug_missing", [R_(f"/nodes/{src[3]}/slug")], "R-PROV-04", "reject", "R-PROV-04")
case("src_status_bad", [S_(f"/nodes/{src[0]}/status", "probably")], "R-STATUS-05", "reject", "R-STATUS-05")
case("src_status_null", [S_(f"/nodes/{src[1]}/status", None)], "R-STATUS-05", "reject", "R-STATUS-05")
case("src_status_case", [S_(f"/nodes/{src[2]}/status", "Verified")], "R-STATUS-05", "reject", "R-STATUS-05")
# --- R-ORIENT-01: every conjunct ---
case("edge_u_negative", [S_(f"/edges/{jn}/u", -1)], "R-ORIENT-01", "reject", "serde")
case("edge_u_gt_v", [S_(f"/edges/{jn}/u", E[jn]["v"]), S_(f"/edges/{jn}/v", E[jn]["u"])], "R-ORIENT-01", "accept", "R-ORIENT-01")
case("edge_u_eq_v", [S_(f"/edges/{jn}/v", E[jn]["u"])], "R-ORIENT-01", "accept", "R-ORIENT-01")
case("edge_v_eq_n", [S_(f"/edges/{jn}/v", n)], "R-ORIENT-01", "accept", "R-ORIENT-01")
case("edge_u_ge_n", [S_(f"/edges/{jn}/u", n), S_(f"/edges/{jn}/v", n + 1)], "R-ORIENT-01", "accept", "R-ORIENT-01")
case("edge_u_string", [S_(f"/edges/{jn}/u", str(E[jn]["u"]))], "R-ORIENT-01", "reject", "serde")
case("edge_u_bool", [S_(f"/edges/{jn}/u", False)], "R-ORIENT-01", "reject", "serde")
case("edge_u_integral_float", [S_(f"/edges/{jn}/u", float(E[jn]["u"]))], "R-ORIENT-01", "reject", "serde")   # v1.6: models.schema_validator integer = JSON integer literal
case("d_integral_float", [S_("/d", float(d))], "R-DIM-02", "reject", "serde")
case("edge_u_float", [S_(f"/edges/{jn}/u", float(E[jn]["u"]) + 0.5)], "R-ORIENT-01", "reject", "serde")
# --- R-CLAIM-01: every conjunct ---
case("claim_parent_null", [S_(f"/nodes/{c0}/parent", None)], "R-CLAIM-01", "reject", "R-CLAIM-01")
case("claim_parent_missing", [R_(f"/nodes/{c0}/parent")], "R-CLAIM-01", "reject", "R-CLAIM-01")
case("claim_parent_is_topic", [S_(f"/nodes/{c0}/parent", N[top[0]]["id"])], "R-CLAIM-01", "accept", "R-CLAIM-01")
case("claim_parent_is_log", [S_(f"/nodes/{c0}/parent", N[logn[0]]["id"])], "R-CLAIM-01", "accept", "R-CLAIM-01")
case("claim_parent_topic_with_matching_edge_and_id", [S_(f"/nodes/{c0}/parent", N[top[0]]["id"]), S_(f"/nodes/{c0}/id", N[top[0]]["id"] + "c0"), S_(f"/edges/{jc}/u", top[0])], "R-CLAIM-01", "accept", "R-CLAIM-01")
case("claim_parent_unknown_id", [S_(f"/nodes/{c0}/parent", "0" * 32)], "R-CLAIM-01", "accept", "R-CLAIM-01")
case("claim_parent_other_source", [S_(f"/nodes/{c0}/parent", N[other_src]["id"])], "R-CLAIM-01", "accept", "R-CLAIM-01")
case("parent_not_hex_on_topic", [S_(f"/nodes/{top[0]}/parent", "zz")], "R-CLAIM-01", "reject", "R-CLAIM-01")
case("parent_suffixed", [S_(f"/nodes/{c0}/parent", N[p0]["id"] + "c0")], "R-CLAIM-01", "reject", "R-CLAIM-01")
# v1.11 (R10 api-typing): trailing-newline tails where ONLY the id/parent grammar can reject (claim ids skip the non-claim
# len==32 guard and keep the parent prefix; a non-claim parent is never rechecked after its Field pattern)
case("claim_id_parent_plus_newline", [S_(f"/nodes/{c0}/id", N[p0]["id"] + "\n")], "R-ID-01", "reject", "R-ID-01")
case("claim_id_suffix_trailing_newline", [S_(f"/nodes/{c0}/id", N[c0]["id"] + "\n")], "R-ID-01", "reject", "R-ID-01")
case("parent_trailing_newline_on_topic", [S_(f"/nodes/{top[0]}/parent", hexid + "\n")], "R-CLAIM-01", "reject", "R-CLAIM-01")
case("claim_parent_trailing_newline", [S_(f"/nodes/{c0}/parent", N[p0]["id"] + "\n")], "R-CLAIM-01", "reject", "R-CLAIM-01")
case("claim_of_dropped", [R_(f"/edges/{jc}")], "R-CLAIM-01", "accept", "R-CLAIM-01")
case("claim_of_duplicate_same_parent", [A_("/edges", dict(ec))], "R-CLAIM-01", "accept", "R-CLAIM-01")
case("claim_of_second_from_other_source", [A_("/edges", {"u": other_src, "v": c0, "type": "claim_of"})], "R-CLAIM-01", "accept", "R-CLAIM-01")
case("claim_of_from_wrong_source_only", [S_(f"/edges/{jc}/u", other_src)], "R-CLAIM-01", "accept", "R-CLAIM-01")
case("claim_of_into_source", [A_("/edges", {"u": src[0], "v": src[1], "type": "claim_of"})], "R-CLAIM-01", "accept", "R-CLAIM-01")
case("claim_of_into_topic", [A_("/edges", {"u": src[0], "v": top[0], "type": "claim_of"})], "R-CLAIM-01", "accept", "R-CLAIM-01")
case("claim_of_into_log", [A_("/edges", {"u": src[0], "v": logn[0], "type": "claim_of"})], "R-CLAIM-01", "accept", "R-CLAIM-01")
case("claim_of_claim_to_claim", [A_("/edges", {"u": clm[0], "v": clm[1], "type": "claim_of"})], "R-CLAIM-01", "accept", "R-CLAIM-01")
(ROOT / "fixtures/m1/mutations.json").write_text(json.dumps(T, indent=0) + "\n"); print(len(T), "cases")
