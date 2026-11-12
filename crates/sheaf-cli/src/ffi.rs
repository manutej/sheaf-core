//! C ABI over the m1 math (ADR-016), loaded in-process by the Rust-backed m1 builder shim
//! (`python/m1/sheafimpl_rs.py`, via ctypes: the gate's builder sandbox forbids fork/exec, so a subprocess is not an
//! option). Every entry point:
//! * takes CSR as (nrows, ncols, nnz, indptr: i64[nrows+1], indices: i64[nnz], data: f64[nnz]) and dense operators as
//!   C-order f64 buffers; scipy/numpy arrays are passed as contiguous int64 / float64 copies by the shim;
//! * re-validates the CSR structure (monotone indptr, in-range indices, overflow-checked sizes) and returns -2 on
//!   malformed input; buffer LENGTHS cannot be checked across a C ABI, so the caller must pass buffers of exactly the
//!   stated sizes (the shim's `_csr_args` verifies indptr = nrows+1 and indptr[-1] = len(indices) = len(data) first);
//! * returns -1 if Rust panicked (caught here; nothing unwinds across the boundary), -3 for an R-PERF-18 cap
//!   (the shim checks every cap first with the reference's own messages, so -3 is defence in depth), -4 numeric failure.
//!
//! Admission (R-DIM-02 / R-FINITE-01 / R-PERF-18 with the exact SPEC messages) stays in the shim, unchanged.
#![allow(clippy::missing_safety_doc, clippy::too_many_arguments)]
use sheaf_cell::m1 as cell;
use sheaf_cell::num::sparse::Csr;
use sheaf_cell::num::tol::H0Policy;
use sheaf_cohom::m1 as co;
use std::panic::{catch_unwind, AssertUnwindSafe};

const E_PANIC: i64 = -1; const E_INPUT: i64 = -2; const E_PERF: i64 = -3; const E_NUM: i64 = -4;

fn guard(f: impl FnOnce() -> Result<i64, i64>) -> i64 {
    match catch_unwind(AssertUnwindSafe(f)) { Ok(Ok(v)) => v, Ok(Err(e)) => e, Err(_) => E_PANIC }
}
fn co_err(e: co::M1H0Error) -> i64 { match e { co::M1H0Error::Perf(_) => E_PERF, co::M1H0Error::Policy(_) => E_INPUT, co::M1H0Error::Numeric(_) => E_NUM } }
unsafe fn sl<'a, T>(p: *const T, n: usize) -> Result<&'a [T], i64> {
    if n == 0 { return Ok(&[]); }
    if p.is_null() { return Err(E_INPUT); }
    Ok(std::slice::from_raw_parts(p, n))
}
unsafe fn sl_mut<'a, T>(p: *mut T, n: usize) -> Result<&'a mut [T], i64> {
    if n == 0 { return Ok(&mut []); }
    if p.is_null() { return Err(E_INPUT); }
    Ok(std::slice::from_raw_parts_mut(p, n))
}
fn us(v: i64) -> Result<usize, i64> { usize::try_from(v).map_err(|_| E_INPUT) }

unsafe fn csr(m: i64, n: i64, nnz: i64, ptr: *const i64, idx: *const i64, val: *const f64) -> Result<Csr, i64> {
    let (m, n, nnz) = (us(m)?, us(n)?, us(nnz)?);
    let p = sl(ptr, m.checked_add(1).ok_or(E_INPUT)?)?; let ix = sl(idx, nnz)?; let v = sl(val, nnz)?;
    if p[0] != 0 || p[m] as usize != nnz || p.windows(2).any(|w| w[1] < w[0]) { return Err(E_INPUT); }
    if ix.iter().any(|&j| j < 0 || j as usize >= n) { return Err(E_INPUT); }
    Ok(Csr { nrows: m, ncols: n, ptr: p.iter().map(|&x| x as usize).collect(), idx: ix.iter().map(|&x| x as usize).collect(), val: v.to_vec() })
}
fn write_csr(a: &Csr, cap: usize, op: &mut [i64], oi: &mut [i64], ov: &mut [f64]) -> Result<i64, i64> {
    let nnz = a.idx.len();
    if nnz > cap || op.len() != a.nrows + 1 { return Err(E_INPUT); }
    for (d, &s) in op.iter_mut().zip(&a.ptr) { *d = s as i64; }
    for (d, &s) in oi.iter_mut().zip(&a.idx) { *d = s as i64; }
    ov[..nnz].copy_from_slice(&a.val);
    Ok(nnz as i64)
}

/// δ⁰ of an m1 sheaf with heterogeneous stalks. `fu`/`fv`: every edge's row-major map concatenated in edge order
/// (edge j contributes de[j]*dv[u_j] / de[j]*dv[v_j] entries). Output CSR (rows = Σde, cols = Σdv) into caller buffers
/// out_ptr[Σde+1], out_idx[cap], out_val[cap]; cap ≥ number of nonzero map entries. Returns nnz.
#[no_mangle]
pub unsafe extern "C" fn sheaf_m1_coboundary(nv: i64, dv: *const i64, ne: i64, u: *const i64, v: *const i64, de: *const i64,
    fu: *const f64, fu_len: i64, fv: *const f64, fv_len: i64, cap: i64, out_ptr: *mut i64, out_idx: *mut i64, out_val: *mut f64) -> i64 {
    guard(|| {
        let (nv, ne, cap) = (us(nv)?, us(ne)?, us(cap)?);
        let dv: Vec<usize> = sl(dv, nv)?.iter().map(|&x| us(x)).collect::<Result<_, _>>()?;
        let (u, v) = (sl(u, ne)?, sl(v, ne)?);
        let de: Vec<usize> = sl(de, ne)?.iter().map(|&x| us(x)).collect::<Result<_, _>>()?;
        let (fu, fv) = (sl(fu, us(fu_len)?)?, sl(fv, us(fv_len)?)?);
        let (mut a, mut b) = (0usize, 0usize); let mut edges = Vec::with_capacity(ne);
        for j in 0..ne {
            let (uj, vj) = (us(u[j])?, us(v[j])?);
            if uj >= nv || vj >= nv { return Err(E_INPUT); }
            let la = de[j].checked_mul(dv[uj]).ok_or(E_INPUT)?; let lb = de[j].checked_mul(dv[vj]).ok_or(E_INPUT)?;
            if a.checked_add(la).is_none_or(|e| e > fu.len()) || b.checked_add(lb).is_none_or(|e| e > fv.len()) { return Err(E_INPUT); }
            edges.push(cell::M1Edge { u: uj, v: vj, de: de[j], fu: &fu[a..a + la], fv: &fv[b..b + lb] }); a += la; b += lb;
        }
        if a != fu.len() || b != fv.len() { return Err(E_INPUT); }
        cell::check_edges(&dv, &edges).map_err(|_| E_INPUT)?;
        dv.iter().try_fold(0usize, |s, &k| s.checked_add(k)).ok_or(E_INPUT)?;
        de.iter().try_fold(0usize, |s, &k| s.checked_add(k)).ok_or(E_INPUT)?;
        let d = cell::coboundary(&dv, &edges);
        write_csr(&d, cap, sl_mut(out_ptr, d.nrows + 1)?, sl_mut(out_idx, cap)?, sl_mut(out_val, cap)?)
    })
}

/// Upper bound on nnz(δᵀδ) (exact structural count incl. products that cancel to 0).
#[no_mangle]
pub unsafe extern "C" fn sheaf_m1_laplacian_nnz(m: i64, n: i64, nnz: i64, ptr: *const i64, idx: *const i64, val: *const f64) -> i64 {
    guard(|| { let d = csr(m, n, nnz, ptr, idx, val)?; Ok(cell::laplacian_nnz(&d) as i64) })
}

/// L₀ = δᵀδ into out_ptr[n+1], out_idx[cap], out_val[cap]; exact zeros dropped (scipy SpGEMM), columns sorted. Returns nnz.
#[no_mangle]
pub unsafe extern "C" fn sheaf_m1_laplacian(m: i64, n: i64, nnz: i64, ptr: *const i64, idx: *const i64, val: *const f64,
    cap: i64, out_ptr: *mut i64, out_idx: *mut i64, out_val: *mut f64) -> i64 {
    guard(|| {
        let d = csr(m, n, nnz, ptr, idx, val)?; let cap = us(cap)?;
        let l = cell::drop_zeros(cell::laplacian(&d));
        write_csr(&l, cap, sl_mut(out_ptr, l.nrows + 1)?, sl_mut(out_idx, cap)?, sl_mut(out_val, cap)?)
    })
}

fn de_vec(de: &[i64], rows: usize) -> Result<Vec<usize>, i64> {
    let v: Vec<usize> = de.iter().map(|&x| us(x)).collect::<Result<_, _>>()?;
    if v.iter().try_fold(0usize, |s, &k| s.checked_add(k)) != Some(rows) { return Err(E_INPUT); }
    Ok(v)
}

#[no_mangle]
pub unsafe extern "C" fn sheaf_m1_edge_discord_csr(m: i64, n: i64, nnz: i64, ptr: *const i64, idx: *const i64, val: *const f64,
    x: *const f64, ne: i64, de: *const i64, out: *mut f64) -> i64 {
    guard(|| {
        let d = csr(m, n, nnz, ptr, idx, val)?; let x = sl(x, d.ncols)?; let ne = us(ne)?;
        let de = de_vec(sl(de, ne)?, d.nrows)?;
        sl_mut(out, ne)?.copy_from_slice(&cell::edge_discord(&d, x, &de)); Ok(0)
    })
}

#[no_mangle]
pub unsafe extern "C" fn sheaf_m1_edge_discord_dense(m: i64, n: i64, a: *const f64, x: *const f64, ne: i64, de: *const i64, out: *mut f64) -> i64 {
    guard(|| {
        let (m, n, ne) = (us(m)?, us(n)?, us(ne)?);
        let a = sl(a, m.checked_mul(n).ok_or(E_INPUT)?)?; let x = sl(x, n)?; let de = de_vec(sl(de, ne)?, m)?;
        sl_mut(out, ne)?.copy_from_slice(&cell::edge_discord_dense(a, m, n, x, &de)); Ok(0)
    })
}

/// E(x) = xᵀLx for a square CSR L.
#[no_mangle]
pub unsafe extern "C" fn sheaf_m1_energy_csr(m: i64, n: i64, nnz: i64, ptr: *const i64, idx: *const i64, val: *const f64, x: *const f64, out: *mut f64) -> i64 {
    guard(|| {
        let l = csr(m, n, nnz, ptr, idx, val)?; if l.nrows != l.ncols { return Err(E_INPUT); }
        *sl_mut(out, 1)?.first_mut().unwrap() = cell::energy(&l, sl(x, l.ncols)?); Ok(0)
    })
}

#[no_mangle]
pub unsafe extern "C" fn sheaf_m1_energy_dense(n: i64, a: *const f64, x: *const f64, out: *mut f64) -> i64 {
    guard(|| {
        let n = us(n)?; let a = sl(a, n.checked_mul(n).ok_or(E_INPUT)?)?;
        *sl_mut(out, 1)?.first_mut().unwrap() = cell::energy_dense(a, n, sl(x, n)?); Ok(0)
    })
}

/// Components of δ (reference order) as (ncols, nrows) pairs into out_sizes[2*ncols]. Returns the count.
#[no_mangle]
pub unsafe extern "C" fn sheaf_m1_components(m: i64, n: i64, nnz: i64, ptr: *const i64, idx: *const i64, val: *const f64, out_sizes: *mut i64) -> i64 {
    guard(|| {
        let d = csr(m, n, nnz, ptr, idx, val)?; let comps = co::components(&d); let o = sl_mut(out_sizes, 2 * d.ncols)?;
        for (k, c) in comps.iter().enumerate() { o[2 * k] = c.cols.len() as i64; o[2 * k + 1] = c.rows.len() as i64; }
        Ok(comps.len() as i64)
    })
}

/// Pass 1 of global_sections under `H0Policy::m1_svd(rel)`: per-component ranks into out_ranks[ncols]. Returns the count.
#[no_mangle]
pub unsafe extern "C" fn sheaf_m1_gs_ranks(m: i64, n: i64, nnz: i64, ptr: *const i64, idx: *const i64, val: *const f64, rel: f64, out_ranks: *mut i64) -> i64 {
    guard(|| {
        let d = csr(m, n, nnz, ptr, idx, val)?; let comps = co::components(&d);
        co::check_caps(&comps).map_err(co_err)?;
        let (rk, _) = co::ranks(&d, &comps, H0Policy::m1_svd(rel)).map_err(co_err)?;
        let o = sl_mut(out_ranks, d.ncols)?; for (k, r) in rk.iter().enumerate() { o[k] = *r as i64; }
        Ok(rk.len() as i64)
    })
}

/// Pass 2: write each component's kernel basis into the zero-initialised C-order g[ncols × big_n] (component order,
/// column offset accumulating), given the pass-1 ranks. Returns the number of columns written (must equal big_n).
#[no_mangle]
pub unsafe extern "C" fn sheaf_m1_gs_fill(m: i64, n: i64, nnz: i64, ptr: *const i64, idx: *const i64, val: *const f64,
    ranks: *const i64, ncomp: i64, g: *mut f64, big_n: i64) -> i64 {
    guard(|| {
        let d = csr(m, n, nnz, ptr, idx, val)?; let comps = co::components(&d); let big_n = us(big_n)?;
        let rk = sl(ranks, us(ncomp)?)?; if rk.len() != comps.len() { return Err(E_INPUT); }
        co::check_caps(&comps).map_err(co_err)?;
        let w: usize = comps.iter().zip(rk).map(|(c, &r)| c.cols.len().saturating_sub(r as usize)).sum();
        if w != big_n || rk.iter().zip(&comps).any(|(&r, c)| r < 0 || r as usize > c.cols.len()) { return Err(E_INPUT); }
        let out = sl_mut(g, d.ncols.checked_mul(big_n).ok_or(E_INPUT)?)?; let mut k = 0;
        for (c, &r) in comps.iter().zip(rk) {
            let kb = co::block_kernel(&d, c, r as usize).map_err(co_err)?;     // one dense block alive at a time
            for (a, &row) in c.cols.iter().enumerate() { for b in 0..kb.ncols() { out[row * big_n + k + b] = kb[(a, b)]; } }
            k += kb.ncols();
        }
        Ok(k as i64)
    })
}

/// Component sizes of a square L (undirected nonzero adjacency) into out_sizes[n]. Returns the count.
#[no_mangle]
pub unsafe extern "C" fn sheaf_m1_l_components(n: i64, nnz: i64, ptr: *const i64, idx: *const i64, val: *const f64, out_sizes: *mut i64) -> i64 {
    guard(|| {
        let l = csr(n, n, nnz, ptr, idx, val)?; let g = co::l_components(&l); let o = sl_mut(out_sizes, l.nrows)?;
        for (k, c) in g.iter().enumerate() { o[k] = c.len() as i64; }
        Ok(g.len() as i64)
    })
}

/// dim ker L under `H0Policy::m1_eig(rel)`.
#[no_mangle]
pub unsafe extern "C" fn sheaf_m1_h0_dim(n: i64, nnz: i64, ptr: *const i64, idx: *const i64, val: *const f64, rel: f64) -> i64 {
    guard(|| { let l = csr(n, n, nnz, ptr, idx, val)?; co::h0_dim_from_l(&l, H0Policy::m1_eig(rel)).map(|v| v as i64).map_err(co_err) })
}

/// ABI version, so the shim refuses a stale library.
#[no_mangle]
pub extern "C" fn sheaf_m1_abi_version() -> i64 { 1 }
