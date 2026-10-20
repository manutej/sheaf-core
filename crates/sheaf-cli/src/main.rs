//! Imperative shell (ADR-005): read spec, call the pure core, interpret its effects.
//! usage: sheaf analyze <spec.json> [--out report.json] [--k K] [--null-reps R] [--seed S] [--top N] [--method shift-invert|polynomial] [--boundary FILE]
use sheaf_cli::v1;
use sheaf_cohom::spectral::Method;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 3 || a[1] != "analyze" { eprintln!("usage: sheaf analyze <spec.json> [--out FILE] [--k K] [--null-reps R] [--seed S] [--top N] [--method shift-invert|polynomial] [--boundary FILE]"); std::process::exit(2); }
    let text = std::fs::read_to_string(&a[2]).unwrap_or_else(|e| { eprintln!("P_INPUT_01 read {}: {e}", a[2]); std::process::exit(2) });
    let mut spec = v1::parse(&text).unwrap_or_else(|e| { eprintln!("{e}"); std::process::exit(2) });
    let mut out = None; let mut i = 3;
    while i < a.len() {
        let val = a.get(i + 1).cloned().unwrap_or_default();
        let bad = |f: &str| -> ! { eprintln!("P_INPUT_01 flag {f}: bad or missing value {val:?}"); std::process::exit(2) };
        match a[i].as_str() {
            "--out" => out = Some(val), "--k" => spec.opts.k = val.parse().unwrap_or_else(|_| bad("--k")), "--null-reps" => spec.opts.null_reps = val.parse().unwrap_or_else(|_| bad("--null-reps")),
            "--seed" => spec.opts.seed = val.parse().unwrap_or_else(|_| bad("--seed")), "--top" => spec.opts.top = val.parse().unwrap_or_else(|_| bad("--top")),
            "--method" => spec.opts.method = match val.as_str() { "polynomial" => Method::Polynomial, "shift-invert" => Method::ShiftInvert, _ => bad("--method") },
            "--boundary" => {
                let t = std::fs::read_to_string(&val).unwrap_or_else(|e| { eprintln!("P_INPUT_01 read {val}: {e}"); std::process::exit(2) });
                spec.boundary = Some(v1::parse_boundary(&t, &spec).unwrap_or_else(|e| { eprintln!("{e}"); std::process::exit(2) }));
            }
            f => { eprintln!("P_INPUT_01 unknown flag {f}"); std::process::exit(2) }
        }
        i += 2;
    }
    let (_report, effects) = v1::analyze(&spec, out);
    for e in effects {
        match e {
            v1::Effect::EmitReport { path: Some(p), contents } => { if let Err(e) = std::fs::write(&p, contents) { eprintln!("P_OUTPUT_01 write {p}: {e}"); std::process::exit(1) } eprintln!("wrote {p}"); }
            v1::Effect::EmitReport { path: None, contents } => print!("{contents}"),
        }
    }
}
