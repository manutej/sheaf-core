//! Block law for `m1::edge_discord` and `m1::edge_discord_dense`.
//! Each edge owns a contiguous residual block. The value is the left fold
//! `s = 0; s += v * v` over that block in row order.
use sheaf_cell::m1::{edge_discord, edge_discord_dense};
use sheaf_cell::num::sparse::Csr;

fn fold_squares(block: &[f64]) -> f64 {
    block.iter().fold(0.0, |s, v| s + v * v)
}

fn blocks(residual: &[f64], de: &[usize]) -> Vec<f64> {
    let mut out = Vec::with_capacity(de.len());
    let mut row = 0;
    for &k in de {
        out.push(fold_squares(&residual[row..row + k]));
        row += k;
    }
    assert_eq!(row, residual.len());
    out
}

/// Every composition of `total` into at most `total + 1` parts, zeros allowed.
/// `total == 0` also includes the empty list and one or two explicit zero widths.
fn parts(total: usize) -> Vec<Vec<usize>> {
    if total == 0 {
        return vec![vec![], vec![0], vec![0, 0]];
    }
    let radix = total + 1;
    let mut out = Vec::new();
    for len in 1..=total + 1 {
        let space = radix.pow(u32::try_from(len).expect("len fits u32"));
        for mask in 0..space {
            let mut de = Vec::with_capacity(len);
            let mut m = mask;
            let mut sum = 0usize;
            for _ in 0..len {
                let p = m % radix;
                m /= radix;
                de.push(p);
                sum += p;
            }
            if sum == total {
                out.push(de);
            }
        }
    }
    out
}

fn grid(n: usize, alpha: &[f64]) -> Vec<Vec<f64>> {
    let mut out = vec![Vec::new()];
    for _ in 0..n {
        let mut next = Vec::new();
        for prefix in &out {
            for &a in alpha {
                let mut v = prefix.clone();
                v.push(a);
                next.push(v);
            }
        }
        out = next;
    }
    out
}

fn assert_bits(got: &[f64], want: &[f64]) {
    assert_eq!(
        got.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
        want.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
    );
}

/// ∀ n ≤ 3, ∀ x with entries in `{0, -0, 1, -1, 1/2, 2^27, 1e16}`, ∀ `de` a
/// composition of n (zeros allowed, at most n+1 parts): on the identity, both
/// discord functions equal the left fold of squares of each block of x, bit for bit.
/// ∀ x in `{-2..2}^2`, ∀ composition `de` of 3: on
/// `A = [[2, -1], [3, 0], [-2, 4]]`, both equal that fold of `Ax`.
/// Every product and block sum in the second family is an integer with absolute
/// value ≤ 2187, so the f64 fold is the real sum of squares.
#[test]
fn discord_blocks_match_left_fold_of_squares() {
    let alpha = [0.0, -0.0, 1.0, -1.0, 0.5, 2f64.powi(27), 1e16];
    for n in 0..=3 {
        let id = Csr::from_triplets(n, n, (0..n).map(|i| (i, i, 1.0)).collect());
        let mut dense = vec![0.0; n * n];
        for i in 0..n {
            dense[i * n + i] = 1.0;
        }
        let widths = parts(n);
        for x in grid(n, &alpha) {
            for de in &widths {
                let want = blocks(&x, de);
                assert_bits(&edge_discord(&id, &x, de), &want);
                assert_bits(&edge_discord_dense(&dense, n, n, &x, de), &want);
            }
        }
    }

    let a = [2.0, -1.0, 3.0, 0.0, -2.0, 4.0];
    let csr = Csr::from_triplets(
        3,
        2,
        vec![(0, 0, 2.0), (0, 1, -1.0), (1, 0, 3.0), (2, 0, -2.0), (2, 1, 4.0)],
    );
    let widths = parts(3);
    for x0 in -2..=2 {
        for x1 in -2..=2 {
            let x = [x0 as f64, x1 as f64];
            let residual: Vec<f64> = (0..3).map(|i| a[i * 2] * x[0] + a[i * 2 + 1] * x[1]).collect();
            for de in &widths {
                let want = blocks(&residual, de);
                assert_bits(&edge_discord(&csr, &x, de), &want);
                assert_bits(&edge_discord_dense(&a, 3, 2, &x, de), &want);
            }
        }
    }
}
