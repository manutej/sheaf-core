//! Block-count law for `h0_relative`'s interior union-find.
//! `n_blocks` is the number of connected components of the interior columns.
//! Two interior columns are adjacent when one row stores a nonzero in both.
//! Explicit zeros are not edges. Boundary columns are not components.
use proptest::prelude::*;
use sheaf_cell::num::sparse::Csr;
use sheaf_cell::num::tol::H0Policy;
use sheaf_cohom::relative::{h0_relative, RelativeError};

/// BFS over interior columns. A new component starts at the least unvisited interior column.
/// Self-pairs and repeated nonzeros in one row do not add components.
fn interior_blocks(nrows: usize, ncols: usize, trips: &[(usize, usize, i8)], boundary: &[usize]) -> usize {
    let mut is_b = vec![false; ncols];
    for &b in boundary {
        is_b[b] = true;
    }
    let mut adj = vec![Vec::new(); ncols];
    let mut by_row = vec![Vec::new(); nrows];
    for &(r, c, v) in trips {
        if v != 0 && !is_b[c] {
            by_row[r].push(c);
        }
    }
    for cols in &by_row {
        if let Some((&first, rest)) = cols.split_first() {
            for &c in rest {
                adj[first].push(c);
                adj[c].push(first);
            }
        }
    }
    let mut seen = vec![false; ncols];
    let mut nc = 0usize;
    for i in 0..ncols {
        if is_b[i] || seen[i] {
            continue;
        }
        nc += 1;
        seen[i] = true;
        let mut q = vec![i];
        let mut h = 0;
        while h < q.len() {
            let x = q[h];
            h += 1;
            for &y in &adj[x] {
                if !seen[y] {
                    seen[y] = true;
                    q.push(y);
                }
            }
        }
    }
    nc
}

type Entry = (usize, usize, i8);
type Case = (usize, usize, Vec<Entry>, Vec<usize>);

/// `nrows, ncols` in `0..6`. Entries are `-1`, `0`, or `1`. Boundary indices may sit one past `ncols`.
fn arb_case() -> impl Strategy<Value = Case> {
    prop_oneof![
        1 => Just((0usize, 0usize, Vec::new(), Vec::new())),
        1 => Just((0usize, 0usize, Vec::new(), vec![0usize])),
        1 => (1usize..5).prop_map(|n| (0usize, n, Vec::new(), Vec::new())),
        1 => (1usize..5).prop_map(|m| (m, 0usize, Vec::new(), Vec::new())),
        8 => (1usize..6, 1usize..6).prop_flat_map(|(nrows, ncols)| {
            (
                Just(nrows),
                Just(ncols),
                prop::collection::vec((0..nrows, 0..ncols, -1i8..=1), 0..10),
                prop::collection::vec(0..ncols + 2, 0..5),
            )
        }),
    ]
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]

    /// ∀ sizes below 6, ∀ triplets in range with values in {-1, 0, 1}, ∀ boundary lists
    /// (indices may be `≥ ncols`):
    /// an out-of-range index yields `BoundaryIndex` of the first such index;
    /// otherwise `n_boundary_dofs` equals the number of distinct boundary indices (a repeat is
    /// the same dof) and `n_blocks` equals the
    /// interior-column component count (nonzero row-sharing; explicit zeros do not connect).
    #[test]
    fn h0_relative_n_blocks_matches_interior_components(
        (nrows, ncols, trips, boundary) in arb_case(),
    ) {
        let t: Vec<(usize, usize, f64)> = trips.iter().map(|&(r, c, v)| (r, c, v as f64)).collect();
        let d0 = Csr::from_triplets(nrows, ncols, t);
        let got = h0_relative(&d0, &boundary, H0Policy::adr012(), ncols.max(1));
        match boundary.iter().copied().find(|&b| b >= ncols) {
            Some(b) => prop_assert_eq!(got.unwrap_err(), RelativeError::BoundaryIndex(b)),
            None => {
                let rel = got.unwrap().value;
                let mut distinct = boundary.clone();
                distinct.sort_unstable();
                distinct.dedup();
                prop_assert_eq!(rel.n_boundary_dofs, distinct.len());
                prop_assert_eq!(rel.n_blocks, interior_blocks(nrows, ncols, &trips, &boundary));
            }
        }
    }
}
