//! Cursor law for `sheaf_m1_coboundary`.
//! A block is refused when `cursor.checked_add(block)` is `None` or lands past the
//! map buffer. That is `map_or(true, |e| e > buf.len())`, now written `is_none_or`.
use sheaf_cell::m1::{self, M1Edge};
use sheaf_cell::num::sparse::Csr;
use sheaf_cli::ffi::sheaf_m1_coboundary;

const E_INPUT: i64 = -2;

/// Match form of `checked_add(cursor, add).map_or(true, |end| end > len)`.
fn map_or_past(cursor: usize, add: usize, len: usize) -> bool {
    match cursor.checked_add(add) {
        None => true,
        Some(end) => end > len,
    }
}

fn entry(i: usize) -> f64 { [0.0, 1.0, -1.0, 0.5][i % 4] }

/// Admission plus `m1::coboundary`, with the cursor test in match form.
fn oracle(dv: &[usize], u: &[i64], v: &[i64], de: &[usize], fu: &[f64], fv: &[f64]) -> Result<Csr, i64> {
    let nv = dv.len();
    let (mut a, mut b) = (0usize, 0usize);
    let mut edges = Vec::with_capacity(u.len());
    for j in 0..u.len() {
        let uj = usize::try_from(u[j]).map_err(|_| E_INPUT)?;
        let vj = usize::try_from(v[j]).map_err(|_| E_INPUT)?;
        if uj >= nv || vj >= nv { return Err(E_INPUT); }
        let la = de[j].checked_mul(dv[uj]).ok_or(E_INPUT)?;
        let lb = de[j].checked_mul(dv[vj]).ok_or(E_INPUT)?;
        if map_or_past(a, la, fu.len()) || map_or_past(b, lb, fv.len()) { return Err(E_INPUT); }
        edges.push(M1Edge { u: uj, v: vj, de: de[j], fu: &fu[a..a + la], fv: &fv[b..b + lb] });
        a += la;
        b += lb;
    }
    if a != fu.len() || b != fv.len() { return Err(E_INPUT); }
    m1::check_edges(dv, &edges).map_err(|_| E_INPUT)?;
    dv.iter().try_fold(0usize, |s, &k| s.checked_add(k)).ok_or(E_INPUT)?;
    de.iter().try_fold(0usize, |s, &k| s.checked_add(k)).ok_or(E_INPUT)?;
    Ok(m1::coboundary(dv, &edges))
}

fn ffi_call(dv: &[i64], u: &[i64], v: &[i64], de: &[i64], fu: &[f64], fv: &[f64]) -> (i64, Vec<i64>, Vec<i64>, Vec<f64>) {
    let (mut ptr, mut idx, mut val) = (vec![0i64; 64], vec![0i64; 64], vec![0f64; 64]);
    let code = unsafe {
        sheaf_m1_coboundary(
            dv.len() as i64, dv.as_ptr(), u.len() as i64, u.as_ptr(), v.as_ptr(), de.as_ptr(),
            fu.as_ptr(), fu.len() as i64, fv.as_ptr(), fv.len() as i64, 64,
            ptr.as_mut_ptr(), idx.as_mut_ptr(), val.as_mut_ptr(),
        )
    };
    (code, ptr, idx, val)
}

fn assert_same(dv: &[i64], u: &[i64], v: &[i64], de: &[i64], fu: &[f64], fv: &[f64]) {
    let dv_us: Vec<usize> = dv.iter().map(|&x| x as usize).collect();
    let de_us: Vec<usize> = de.iter().map(|&x| x as usize).collect();
    let (code, ptr, idx, val) = ffi_call(dv, u, v, de, fu, fv);
    match oracle(&dv_us, u, v, &de_us, fu, fv) {
        Err(err) => assert_eq!(code, err, "dv={dv:?} u={u:?} v={v:?} de={de:?} fu={} fv={}", fu.len(), fv.len()),
        Ok(csr) => {
            assert_eq!(code, csr.idx.len() as i64, "dv={dv:?} u={u:?} v={v:?} de={de:?}");
            assert_eq!(&ptr[..csr.ptr.len()], &csr.ptr.iter().map(|&x| x as i64).collect::<Vec<_>>());
            assert_eq!(&idx[..csr.idx.len()], &csr.idx.iter().map(|&x| x as i64).collect::<Vec<_>>());
            assert_eq!(
                val[..csr.val.len()].iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                csr.val.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
            );
        }
    }
}

/// Next edge in order `(de, v, u)`. `false` once every triple is exhausted.
fn bump_edges(u: &mut [i64], v: &mut [i64], de: &mut [i64], nv: i64) -> bool {
    for i in 0..u.len() {
        de[i] += 1;
        if de[i] <= 2 { return true; }
        de[i] = 0;
        v[i] += 1;
        if v[i] <= nv { return true; }
        v[i] = 0;
        u[i] += 1;
        if u[i] <= nv { return true; }
        u[i] = 0;
    }
    false
}

/// ∀ cursor, addend, length in `{0, 1, 2, 4, usize::MAX/2, usize::MAX-1, usize::MAX}`:
/// `map_or_past` equals `checked_add(cursor, add).is_none_or(|e| e > len)`, including overflow.
/// ∀ `nv, ne ≤ 2`, ∀ stalk dims in `0..=2`, ∀ endpoints in `0..=nv`, ∀ `de` in `0..=2`,
/// ∀ buffer lengths in `0..=4`, with `fu[i] = entry(i)` and `fv[i] = entry(i + 1)`
/// from `{0, 1, -1, 1/2}`: the FFI returns `-2` exactly when the oracle refuses, and on
/// acceptance writes the same `ptr`, `idx`, and value bits as `m1::coboundary`.
#[test]
fn m1_coboundary_matches_map_or_cursor() {
    let bounds = [0usize, 1, 2, 4, usize::MAX / 2, usize::MAX - 1, usize::MAX];
    for &cursor in &bounds {
        for &add in &bounds {
            for &len in &bounds {
                assert_eq!(map_or_past(cursor, add, len), cursor.checked_add(add).is_none_or(|e| e > len));
            }
        }
    }
    let fu_all: Vec<f64> = (0..4).map(entry).collect();
    let fv_all: Vec<f64> = (0..4).map(|i| entry(i + 1)).collect();
    let mut seen = 0usize;
    for nv in 0..=2 {
        let mut dv = vec![0i64; nv];
        loop {
            for ne in 0..=2 {
                let (mut u, mut v, mut de) = (vec![0i64; ne], vec![0i64; ne], vec![0i64; ne]);
                loop {
                    for fu_len in 0..=4 {
                        for fv_len in 0..=4 {
                            assert_same(&dv, &u, &v, &de, &fu_all[..fu_len], &fv_all[..fv_len]);
                            seen += 1;
                        }
                    }
                    if !bump_edges(&mut u, &mut v, &mut de, nv as i64) { break; }
                }
            }
            let mut i = 0;
            while i < dv.len() {
                dv[i] += 1;
                if dv[i] <= 2 { break; }
                dv[i] = 0;
                i += 1;
            }
            if i == dv.len() { break; }
        }
    }
    // nv=0: 1×(1+3+9); nv=1: 3×157; nv=2: 9×757; times 5×5 buffer lengths.
    assert_eq!(seen, (13 + 3 * 157 + 9 * 757) * 25);
}
