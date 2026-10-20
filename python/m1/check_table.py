"""Reference-only: run the synthetic shared mutation table through the Python mirrors (pydantic models + JSON Schema),
like frontier contract/tests/test_spec_public.py does on the frozen spec. usage: python python/m1/check_table.py"""
import copy, json, pathlib, sys
ROOT = pathlib.Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "python"))
import models
S = json.loads((ROOT / "fixtures/m1/spec.json").read_text()); MUT = json.loads((ROOT / "fixtures/m1/mutations.json").read_text())
V = models.schema_validator(json.loads((ROOT / "python/sheaf.schema.json").read_text()))
def apply(doc, op):
    def ptr(p):
        cur = doc
        for k in p.split("/")[1:]: cur = cur[int(k)] if isinstance(cur, list) else cur[k]
        return cur
    par, key = op["path"].rsplit("/", 1); tgt = ptr(par) if par else doc
    if op["op"] == "set":
        if isinstance(tgt, list): tgt[int(key)] = copy.deepcopy(op["value"])
        else: tgt[key] = copy.deepcopy(op["value"])
    elif op["op"] == "remove": tgt.pop(int(key) if isinstance(tgt, list) else key)
    elif op["op"] == "append": ptr(op["path"]).append(copy.deepcopy(op["value"]))
bad = []
for c in MUT:
    T = copy.deepcopy(S)
    for op in c["ops"]: apply(T, op)
    try: models.validate_spec(T); got = "ok"
    except models.SpecError as e: got = sorted({m.split(":")[0] for m in e.messages})
    except Exception as e: got = [type(e).__name__]
    okm = got == "ok" if c["models"] == "ok" else (got != "ok" and c["models"] in got)
    sch = "accept" if V.is_valid(T) else "reject"
    if not okm or sch != c["schema"]: bad.append((c["name"], c["models"], got, c["schema"], sch))
print(f"{len(MUT)} rows; mismatches: {len(bad)}"); [print(b) for b in bad]; sys.exit(1 if bad else 0)
