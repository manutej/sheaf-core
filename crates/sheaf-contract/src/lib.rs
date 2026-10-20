//! Rust mirror of contract/sheaf.schema.json v1.18 (serde) + semantic validation + per-edge discord (R-RUST-20).
//! Every semantic check in `validate()` sits on its own line tagged `// M:<name>`; `../mutants_public.sh` deletes each
//! tagged line in a scratch copy and requires `cargo test` to FAIL (every branch has a killing test).
use serde::Deserialize;
use std::collections::HashMap;
#[derive(Deserialize, Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[serde(rename_all = "snake_case")]
pub enum EdgeType { About, Related, Logged, ClaimOf, SameSlug }
#[derive(Deserialize, Debug, Clone, Copy, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum Kind { Source, Topic, Claim, Log }
/// Unknown node fields are ignored (schema/pydantic: extra allowed); known optional fields are typed (string or null).
#[derive(Deserialize, Debug)]
pub struct Node { pub id: String, pub kind: Kind,
                  // v1.9: explicit #[serde(default)] (serde_derive already maps an absent Option to None; stated, not implied)
                  #[serde(default)] pub slug: Option<String>, #[serde(default)] pub title: Option<String>, #[serde(default)] pub status: Option<String>,
                  #[serde(default)] pub url: Option<String>, #[serde(default)] pub parent: Option<String>, #[serde(default)] pub fetched: Option<String> }
#[derive(Deserialize, Debug)]
#[serde(deny_unknown_fields)]
pub struct Edge { pub u: usize, pub v: usize, #[serde(rename = "type")] pub ty: EdgeType }
#[derive(Deserialize, Debug)]
#[serde(deny_unknown_fields)]
pub struct SheafSpec { pub d: usize, pub nodes: Vec<Node>, pub x: Vec<Vec<f64>>, pub edges: Vec<Edge>, pub maps: HashMap<EdgeType, (Vec<Vec<f64>>, Vec<Vec<f64>>)> }
fn hex32(b: &[u8]) -> bool { b.len() >= 32 && b[..32].iter().all(|c| c.is_ascii_digit() || (b'a'..=b'f').contains(c)) }
/// R-ID-01: 32 lowercase hex (Notion page id) + optional `c<digits>` claim suffix. Same predicate as the schema/pydantic
/// pattern `^[0-9a-f]{32}(c[0-9]+)?$` with end-of-text anchoring (no trailing newline, no whitespace).
pub fn valid_id(id: &str) -> bool {
    let b = id.as_bytes();
    if !hex32(b) { return false; }
    let rest = &b[32..];
    rest.is_empty() || (rest.len() >= 2 && rest[0] == b'c' && rest[1..].iter().all(|c| c.is_ascii_digit()))
}
/// Parent pattern `^[0-9a-f]{32}$` (a source page id, never suffixed).
pub fn valid_parent(p: &str) -> bool { p.len() == 32 && hex32(p.as_bytes()) }
const STATUSES: [&str; 7] = ["verified", "partly verified", "search-verified", "claimed", "unverified", "refuted", "not found"];
fn is_identity(f: &[Vec<f64>]) -> bool { f.iter().enumerate().all(|(a, r)| r.iter().enumerate().all(|(b, v)| *v == if a == b { 1.0 } else { 0.0 })) }
fn err<T>(rid: &str) -> Result<T, String> { Err(rid.to_string()) }
impl SheafSpec {
    /// Returns the first violated requirement ID, if any (models.validate_spec reports the same IDs).
    pub fn validate(&self) -> Result<(), String> {
        let (n, d) = (self.nodes.len(), self.d);
        if n == 0 { return err("R-SCHEMA-01"); }                                                               // M:nodes_nonempty
        if d < 1 || d > 4096 { return err("R-DIM-02"); }                                                       // M:d_range
        if self.x.len() != n { return err("R-DIM-02"); }                                                       // M:x_rows
        if self.x.iter().any(|r| r.len() != d) { return err("R-DIM-02"); }                                     // M:x_cols
        if self.x.iter().flatten().any(|v| !v.is_finite()) { return err("R-FINITE-01"); }                      // M:x_finite
        for (fu, fv) in self.maps.values() { for f in [fu, fv] {
            if f.len() != d { return err("R-DIM-02"); }                                                        // M:map_rows
            if f.iter().any(|r| r.len() != d) { return err("R-DIM-02"); }                                      // M:map_cols
            if f.iter().flatten().any(|v| !v.is_finite()) { return err("R-FINITE-01"); }                       // M:map_finite
        } }
        let mut ids = std::collections::HashSet::new();
        for nd in &self.nodes {
            if !ids.insert(nd.id.as_str()) { return err("R-ID-01"); }                                          // M:id_unique
            if !valid_id(&nd.id) { return err("R-ID-01"); }                                                    // M:id_pattern
            if nd.kind != Kind::Claim && nd.id.len() != 32 { return err("R-ID-01"); }                          // M:id_suffix_claims_only
            if nd.parent.as_deref().map_or(false, |p| !valid_parent(p)) { return err("R-CLAIM-01"); }          // M:parent_pattern
            if nd.kind == Kind::Source {
                if nd.slug.as_deref().unwrap_or("").is_empty() { return err("R-PROV-04"); }                    // M:src_slug
                if !nd.url.as_deref().map_or(false, |u| u.starts_with("http://") || u.starts_with("https://")) { return err("R-PROV-04"); } // M:src_url
                if !STATUSES.contains(&nd.status.as_deref().unwrap_or("")) { return err("R-STATUS-05"); }      // M:src_status
            }
            if nd.kind == Kind::Claim && nd.parent.is_none() { return err("R-CLAIM-01"); }                     // early fast-fail for a MISSING parent only (untagged); NOT equivalent to M:claim_parent_is_source (topic/log/unknown parents), which must stay
        }
        for e in &self.edges {
            if !(e.u < e.v) { return err("R-ORIENT-01"); }                                                     // M:orient_u_lt_v
            if !(e.v < n) { return err("R-ORIENT-01"); }                                                       // M:orient_v_lt_n
            if !self.maps.contains_key(&e.ty) { return err("R-TYPE-03"); }                                     // M:type_has_map
        }
        if self.maps.values().any(|(fu, fv)| is_identity(fu) && is_identity(fv)) { return err("R-NONID-04"); } // M:nonid
        let idx: HashMap<&str, usize> = self.nodes.iter().enumerate().map(|(i, nd)| (nd.id.as_str(), i)).collect();
        let mut cof: HashMap<usize, Vec<usize>> = HashMap::new();
        for e in self.edges.iter().filter(|e| e.ty == EdgeType::ClaimOf) { cof.entry(e.v).or_default().push(e.u); }
        for (i, nd) in self.nodes.iter().enumerate() {
            if nd.kind != Kind::Claim { continue; }
            let p = idx.get(nd.parent.as_deref().unwrap_or("")).copied();
            if !p.map_or(false, |p| self.nodes[p].kind == Kind::Source) { return err("R-CLAIM-01"); }          // M:claim_parent_is_source
            if cof.get(&i).map(|v| v.as_slice()) != Some(&[p.unwrap()][..]) { return err("R-CLAIM-01"); }       // M:claim_exactly_one_from_parent
            if !nd.id.starts_with(nd.parent.as_deref().unwrap()) { return err("R-ID-01"); }                    // M:claim_id_prefix
        }
        if cof.keys().any(|&v| self.nodes[v].kind != Kind::Claim) { return err("R-CLAIM-01"); }                // M:claim_of_only_into_claims
        Ok(())
    }
    /// discord(e) = || F_v x_v - F_u x_u ||^2
    pub fn discord(&self) -> Vec<f64> {
        self.edges.iter().map(|e| { let (fu, fv) = &self.maps[&e.ty];
            (0..self.d).map(|a| { let r: f64 = (0..self.d).map(|b| fv[a][b] * self.x[e.v][b] - fu[a][b] * self.x[e.u][b]).sum(); r * r }).sum() }).collect()
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
    fn raw() -> String { std::fs::read_to_string(std::env::var("SPEC_JSON").unwrap_or("../../spec.json".into())).unwrap() }
    fn load() -> SheafSpec { serde_json::from_str(&raw()).unwrap() }
    #[test] fn r_schema_01_parses_and_validates() { assert_eq!(load().validate(), Ok(())); }
    /// v1.9: every Option field of Node carries an explicit #[serde(default)] (serde_derive would already map an ABSENT
    /// `Option<T>` to None); this test is the executable evidence: absent AND explicit null both parse to None, so
    /// the builder's output (topics/logs without status/url/parent, sources without parent) parses exactly as the schema allows.
    #[test] fn r_schema_01_optional_fields_absent_or_null() {
        let n: Node = serde_json::from_str(r#"{"id":"0123456789abcdef0123456789abcdef","kind":"topic"}"#).unwrap();
        assert!(n.slug.is_none() && n.title.is_none() && n.status.is_none() && n.url.is_none() && n.parent.is_none() && n.fetched.is_none());
        let n: Node = serde_json::from_str(r#"{"id":"0123456789abcdef0123456789abcdef","kind":"log","parent":null,"status":null}"#).unwrap();
        assert!(n.parent.is_none() && n.status.is_none());
        assert!(serde_json::from_str::<Node>(r#"{"id":"0123456789abcdef0123456789abcdef","kind":"topic","title":5}"#).is_err());
        let s = load(); assert!(s.nodes.iter().any(|n| n.kind == Kind::Topic && n.status.is_none() && n.url.is_none() && n.parent.is_none()));
        assert!(s.nodes.iter().any(|n| n.kind == Kind::Source && n.parent.is_none()));
        println!("RUST-EVIDENCE r_schema_01_optional_fields_absent_or_null: {{id,kind}}-only topic parsed, slug/title/status/url/parent/fetched all None; explicit null -> None; title:5 rejected");
    }
    /// v1.9 (R8 security): a JSON integer beyond f64 range is rejected at parse time (serde_json 'number out of range'),
    /// matching models (R-DIM-02 'not a valid number'), sheafcore.from_spec and validate.py (R-DIM-02).
    #[test] fn r_dim_02_int_beyond_f64_rejected() {
        let mut v: Value = serde_json::from_str(&raw()).unwrap();
        v["x"][0][0] = Value::from(0.123456789012345); let base = serde_json::to_string(&v).unwrap();   // unique sentinel = x[0][0]
        assert_eq!(base.matches("0.123456789012345").count(), 1);
        let big = format!("1{}", "0".repeat(309));
        let doc = base.replacen("0.123456789012345", &big, 1); assert_ne!(doc, base);
        assert!(serde_json::from_str::<SheafSpec>(&base).is_ok());
        assert!(serde_json::from_str::<SheafSpec>(&doc).is_err());
        println!("RUST-EVIDENCE r_dim_02_int_beyond_f64_rejected: x[0][0] = 10^309 -> serde error (no Inf reaches validate())");
    }
    /// v1.15 (R14 api-typing): serde_json 1.0.140 parses an integer token that overflows u64 / i64 as f64
    /// (no 'number out of range' below the f64 limit), so FINITE integers such as 2^64, -2^64, 10^20 are accepted exactly
    /// like schema / models / validate.py / sheafcore (they read float(int)); only a token beyond f64 range is a serde error.
    #[test] fn r_dim_02_large_finite_ints_accepted() {
        let mut v: Value = serde_json::from_str(&raw()).unwrap();
        v["x"][0][0] = Value::from(0.123456789012345); let base = serde_json::to_string(&v).unwrap();
        for (tok, want) in [("18446744073709551616", 18446744073709551616.0_f64), ("-18446744073709551616", -18446744073709551616.0_f64),
                            ("100000000000000000000", 1e20_f64), ("18446744073709551615", 18446744073709551615.0_f64), ("-9223372036854775809", -9223372036854775809.0_f64)] {
            let doc = base.replacen("0.123456789012345", tok, 1); assert_ne!(doc, base);
            let s: SheafSpec = serde_json::from_str(&doc).unwrap_or_else(|e| panic!("{tok}: {e}"));
            assert_eq!(s.x[0][0], want, "{tok}"); assert!(s.validate().is_ok(), "{tok}");
        }
        println!("RUST-EVIDENCE r_dim_02_large_finite_ints_accepted: x[0][0] in {{2^64, -2^64, 10^20, u64::MAX, i64::MIN-1}} -> parsed as finite f64, validate() ok");
    }
    #[test] fn r_rust_20_discord_matches_python() {
        let s = load(); let py: Vec<f64> = serde_json::from_str(&std::fs::read_to_string(std::env::var("PY_DISCORD").expect("PY_DISCORD")).unwrap()).unwrap();
        let rs = s.discord(); assert_eq!(rs.len(), py.len());
        for (a, b) in rs.iter().zip(&py) { assert!((a - b).abs() <= 1e-9, "{a} vs {b} (abs tol 1e-9)"); }
    }
    /// Parity as a property: conform.sh draws K fresh random x (seed printed), Python scores each; Rust must agree.
    #[test] fn r_rust_20_discord_parity_fresh_x() {
        let cases: Vec<Value> = serde_json::from_str(&std::fs::read_to_string(std::env::var("PY_DISCORD_RAND").expect("PY_DISCORD_RAND")).unwrap()).unwrap();
        assert!(cases.len() == 6, "R-RUST-20: exactly 6 fresh-x parity cases per run");
        for (k, c) in cases.iter().enumerate() {
            let mut s = load(); s.x = serde_json::from_value(c["x"].clone()).unwrap(); assert_eq!(s.validate(), Ok(()));
            let py: Vec<f64> = serde_json::from_value(c["discord"].clone()).unwrap(); let rs = s.discord(); assert_eq!(rs.len(), py.len());
            for (j, (a, b)) in rs.iter().zip(&py).enumerate() { assert!((a - b).abs() <= 1e-9, "case {k} edge {j}: {a} vs {b}"); }
        }
    }
    #[test] fn r_rust_20_discord_known_answer() {
        // independent hand-computed KAT: d=2, one related edge, F_u = [[1,2],[0,1]], F_v = [[0,1],[1,0]]
        let mut s = load(); s.d = 2; s.edges = vec![Edge { u: 0, v: 1, ty: EdgeType::Related }];
        s.maps = HashMap::from([(EdgeType::Related, (vec![vec![1.0, 2.0], vec![0.0, 1.0]], vec![vec![0.0, 1.0], vec![1.0, 0.0]]))]);
        s.x = vec![vec![1.0, -1.0], vec![3.0, 0.5]];
        // F_v x_v = (0.5, 3); F_u x_u = (-1, -1); residual (1.5, 4) -> 2.25 + 16 = 18.25
        assert_eq!(s.discord(), vec![18.25]);
    }
    #[test] fn r_finite_01_rejects_nonfinite() {
        for (in_x, val) in [(true, f64::NAN), (true, f64::INFINITY), (true, f64::NEG_INFINITY), (false, f64::NAN), (false, f64::INFINITY), (false, f64::NEG_INFINITY)] {
            for side in 0..2 {
                let mut s = load(); assert_eq!(s.validate(), Ok(()));
                if in_x { s.x[3][side] = val; } else { let k = EdgeType::Related; let m = s.maps.get_mut(&k).unwrap(); if side == 0 { m.0[1][2] = val } else { m.1[2][1] = val } }
                assert_eq!(s.validate(), Err("R-FINITE-01".to_string()), "in_x={in_x} val={val} side={side}");
            }
        }
    }
    #[test] fn r_id_01_pattern_kat() {
        assert!(valid_id("0123456789abcdef0123456789abcdef")); assert!(valid_id("0123456789abcdef0123456789abcdefc12"));
        for bad in ["", "0123456789ABCDEF0123456789abcdef", "0123456789abcdef0123456789abcde", "0123456789abcdef0123456789abcdefc",
                    "0123456789abcdef0123456789abcdefx1", "0123456789abcdef0123456789abcdef c1", "0123456789abcdef0123456789abcdef\n",
                    "0123456789abcdef0123456789abcdefc1\n", " 0123456789abcdef0123456789abcdef"] { assert!(!valid_id(bad), "{bad:?}"); }
        assert!(valid_parent("0123456789abcdef0123456789abcdef"));
        for bad in ["0123456789abcdef0123456789abcdefc1", "0123456789abcdef0123456789abcdef\n", "0123456789abcdef0123456789abcdeg"] { assert!(!valid_parent(bad), "{bad:?}"); }
    }
    fn apply(v: &mut Value, op: &Value) {
        let path = op["path"].as_str().unwrap(); let (parent, key) = path.rsplit_once('/').unwrap();
        let tgt = if parent.is_empty() { &mut *v } else { v.pointer_mut(parent).unwrap_or_else(|| panic!("bad path {path}")) };
        match op["op"].as_str().unwrap() {
            "set" => match tgt { Value::Object(m) => { m.insert(key.to_string(), op["value"].clone()); }, Value::Array(a) => { a[key.parse::<usize>().unwrap()] = op["value"].clone(); }, _ => panic!("set {path}") },
            "remove" => match tgt { Value::Object(m) => { m.remove(key); }, Value::Array(a) => { a.remove(key.parse::<usize>().unwrap()); }, _ => panic!("remove {path}") },
            "append" => { v.pointer_mut(path).unwrap().as_array_mut().unwrap().push(op["value"].clone()); },
            o => panic!("unknown op {o}"),
        }
    }
    /// Shared rejection table (tests/mutations.json), also run by pytest against schema + models.validate_spec.
    #[test] fn r_model_mutation_table() {
        let base: Value = serde_json::from_str(&raw()).unwrap();
        let table: Vec<Value> = serde_json::from_str(&std::fs::read_to_string(std::env::var("MUTATIONS_JSON").unwrap_or("../tests/mutations.json".into())).unwrap()).unwrap();
        assert!(table.len() >= 80);  // SPEC R-SCHEMA-01: >= 80 cases
        let mut ok_rows: Vec<String> = Vec::new();
        for case in &table {
            let mut v = base.clone(); for op in case["ops"].as_array().unwrap() { apply(&mut v, op); }
            let name = case["name"].as_str().unwrap(); let want = case["rust"].as_str().unwrap();
            match serde_json::from_value::<SheafSpec>(v) {
                Err(e) => assert_eq!(want, "serde", "{name}: serde rejected ({e}) but table expects {want}"),
                Ok(s) => { let got = s.validate().err().unwrap_or("ok".into()); assert_eq!(got, want, "{name}"); if want == "ok" { ok_rows.push(name.to_string()); } }
            }
        }
        println!("RUST-EVIDENCE r_model_mutation_table: {} rows agree; ok rows parsed by serde AND validated: {}", table.len(), ok_rows.join(","));
    }
}
