//! Nonzero column-component law for `m1::pattern_components`.
//! A stored 0 (including -0.0) is not a coupling. NaN and infinities are.
use proptest::prelude::*;
use sheaf_cell::m1::pattern_components;
use sheaf_cell::num::sparse::Csr;

fn couples(v: f64) -> bool {
    v != 0.0
}

/// Columns are adjacent when one row stores a nonzero in both. BFS starts at the
/// least unseen column, so component ids follow the minimum column index.
/// Rows with no coupling are omitted. Rows and columns inside a group are ascending.
fn nonzero_column_bfs(nrows: usize, ncols: usize, entries: &[(usize, usize, f64)]) -> Vec<(Vec<usize>, Vec<usize>)> {
    let mut row_nz: Vec<Vec<usize>> = vec![Vec::new(); nrows];
    for &(r, c, v) in entries {
        if couples(v) {
            row_nz[r].push(c);
        }
    }
    let mut adj = vec![Vec::new(); ncols];
    for cols in &mut row_nz {
        cols.sort_unstable();
        cols.dedup();
        if let Some(&hub) = cols.first() {
            for &c in &cols[1..] {
                adj[hub].push(c);
                adj[c].push(hub);
            }
        }
    }
    let mut col_comp = vec![usize::MAX; ncols];
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut start = 0;
    while start < ncols {
        if col_comp[start] != usize::MAX {
            start += 1;
            continue;
        }
        let id = groups.len();
        col_comp[start] = id;
        let mut members = vec![start];
        let mut q = vec![start];
        let mut h = 0;
        while h < q.len() {
            let x = q[h];
            h += 1;
            for &y in &adj[x] {
                if col_comp[y] == usize::MAX {
                    col_comp[y] = id;
                    members.push(y);
                    q.push(y);
                }
            }
        }
        members.sort_unstable();
        groups.push(members);
        start += 1;
    }
    let mut out: Vec<(Vec<usize>, Vec<usize>)> = groups.into_iter().map(|cols| (Vec::new(), cols)).collect();
    for (i, cols) in row_nz.iter().enumerate() {
        if let Some(&c) = cols.first() {
            out[col_comp[c]].0.push(i);
        }
    }
    out
}

fn arb_value() -> impl Strategy<Value = f64> {
    prop_oneof![
        Just(0.0),
        Just(-0.0),
        Just(1.0),
        Just(-1.0),
        Just(f64::NAN),
        Just(f64::INFINITY),
        Just(f64::NEG_INFINITY),
    ]
}

/// Every shape in `0..5` × `0..5`, and every triplet list of length `< 12`.
/// Out-of-range endpoints are dropped, so `from_triplets` sees only legal indices.
fn arb_triplets() -> impl Strategy<Value = (usize, usize, Vec<(usize, usize, f64)>)> {
    (0usize..5, 0usize..5, prop::collection::vec((0usize..5, 0usize..5, arb_value()), 0..12)).prop_map(
        |(nrows, ncols, entries)| {
            let entries = entries.into_iter().filter(|(r, c, _)| *r < nrows && *c < ncols).collect();
            (nrows, ncols, entries)
        },
    )
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]

    /// ∀ nrows, ncols < 5, ∀ triplet lists with endpoints in range and values in
    /// {0, -0, ±1, NaN, ±∞}: `pattern_components` equals nonzero-column BFS.
    /// Every column appears once. A row appears iff it stores a coupling.
    #[test]
    fn pattern_components_match_nonzero_column_bfs((nrows, ncols, entries) in arb_triplets()) {
        let want = nonzero_column_bfs(nrows, ncols, &entries);
        let got = pattern_components(&Csr::from_triplets(nrows, ncols, entries.clone()));
        prop_assert_eq!(&got, &want);
        let mut seen = vec![false; ncols];
        for (rows, cols) in &got {
            prop_assert!(rows.windows(2).all(|w| w[0] < w[1]));
            prop_assert!(cols.windows(2).all(|w| w[0] < w[1]));
            for &c in cols {
                prop_assert!(!seen[c]);
                seen[c] = true;
            }
        }
        prop_assert!(seen.iter().all(|&s| s));
        let mut row_hit = vec![false; nrows];
        for &(r, _, v) in &entries {
            if couples(v) {
                row_hit[r] = true;
            }
        }
        let mut seen_row = vec![false; nrows];
        for (rows, _) in &got {
            for &r in rows {
                prop_assert!(row_hit[r]);
                prop_assert!(!seen_row[r]);
                seen_row[r] = true;
            }
        }
        prop_assert_eq!(seen_row, row_hit);
    }
}
