use std::collections::BTreeSet;

use faer::Mat;
use proptest::prelude::*;
use sheaf_cell::complex::{CellComplex, ComplexError};
use sheaf_cell::sheaf::Sheaf;

/// (n, edges, faces) on n in 3..7 vertices; every face has three distinct vertices, edges may repeat or be reversed.
fn arb_simplicial() -> impl Strategy<Value = (usize, Vec<(usize, usize)>, Vec<[usize; 3]>)> {
    (3usize..7).prop_flat_map(|n| {
        let edges = prop::collection::vec((0..n, 0..n).prop_filter("no loop", |(a, b)| a != b), 0..8);
        let faces = prop::collection::vec(prop::sample::subsequence((0..n).collect::<Vec<_>>(), 3), 0..6)
            .prop_map(|fs| fs.into_iter().map(|f| [f[2], f[0], f[1]]).collect());
        (Just(n), edges, faces)
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    /// `simplicial` always yields a complex that passes `validate` (indices in range, boundary of boundary = 0), with
    /// n_cells[1] = #distinct undirected edges (given ∪ implied by faces), n_cells[2] = #distinct vertex triples,
    /// and exactly 2 incidences per edge and 3 per triangle.
    #[test]
    fn simplicial_is_valid_with_exact_cell_counts((n, edges, faces) in arb_simplicial()) {
        let cx = CellComplex::simplicial(n, &edges, &faces);
        prop_assert_eq!(cx.validate(), Ok(()));
        let mut es: BTreeSet<(usize, usize)> = edges.iter().map(|&(a, b)| (a.min(b), a.max(b))).collect();
        let ts: BTreeSet<[usize; 3]> = faces.iter().map(|f| { let mut t = *f; t.sort(); t }).collect();
        for t in &ts { es.extend([(t[0], t[1]), (t[0], t[2]), (t[1], t[2])]); }
        prop_assert_eq!(cx.n_cells[0], n);
        prop_assert_eq!(cx.n_cells[1], es.len());
        prop_assert_eq!(cx.inc[0].len(), 2 * es.len());
        prop_assert_eq!(cx.dim(), if ts.is_empty() { 1 } else { 2 });
        if !ts.is_empty() {
            prop_assert_eq!(cx.n_cells[2], ts.len());
            prop_assert_eq!(cx.inc[1].len(), 3 * ts.len());
        }
    }

    /// A functorial sheaf on a simplicial 2-complex — F(c) = R^d, F(σ ≤ τ) = (s_τ / s_σ) I with a power-of-two
    /// scale s_c per cell — has δ¹ δ⁰ = 0 EXACTLY (each term of the product is ± a power of two, so f64 is exact).
    #[test]
    fn coboundary_squares_to_zero_on_a_functorial_sheaf(
        (n, edges, faces) in arb_simplicial(), d in 1usize..3, exps in prop::collection::vec(-4i32..=4, 64),
    ) {
        let complex = CellComplex::simplicial(n, &edges, &faces);
        prop_assume!(complex.dim() == 2);
        let base = [0, complex.n_cells[0], complex.n_cells[0] + complex.n_cells[1]];
        let s = |k: usize, i: usize| 2f64.powi(exps[(base[k] + i) % exps.len()]);
        let mut sheaf = Sheaf::constant(complex, d);
        for k in 0..2 {
            for (j, inc) in sheaf.complex.inc[k].iter().enumerate() {
                let r = s(k + 1, inc.cofacet) / s(k, inc.facet);
                sheaf.maps[k][j] = Mat::<f64>::identity(d, d) * faer::Scale(r);
            }
        }
        prop_assert_eq!(sheaf.validate(), Ok(()));
        let (d0, d1) = (sheaf.coboundary_dense(0), sheaf.coboundary_dense(1));
        prop_assert_eq!((d1.ncols(), d0.nrows()), (sheaf.cochain_dim(1), sheaf.cochain_dim(1)));
        for i in 0..d1.nrows() { for j in 0..d0.ncols() {
            let v: f64 = (0..d0.nrows()).map(|m| d1[(i, m)] * d0[(m, j)]).sum();
            prop_assert_eq!(v, 0.0, "(d1 d0)[{}, {}] = {}", i, j, v);
        }}
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    /// Flipping ONE incidence sign of a valid simplicial complex: `validate` fails iff some (k-1, k+1) pair now sums
    /// to ±2, and reports the lexicographically least such (facet, cocofacet). Flip [e:t] (e = [a, b], a < b) =>
    /// witness (a, t). Flip [v:e] => Ok iff no triangle contains e, else witness (v, least triangle containing e).
    #[test]
    fn validate_reports_the_least_witness_of_a_flipped_sign((n, edges, faces) in arb_simplicial(), k in 0usize..2, pick in any::<prop::sample::Index>()) {
        let mut cx = CellComplex::simplicial(n, &edges, &faces);
        prop_assume!(cx.inc.get(k).is_some_and(|l| !l.is_empty()));
        let i = pick.index(cx.inc[k].len());
        cx.inc[k][i].sign = -cx.inc[k][i].sign;
        let (facet, cofacet) = (cx.inc[k][i].facet, cx.inc[k][i].cofacet);
        let expected = if k == 1 {
            Err(ComplexError::BoundaryOfBoundaryNonzero { dim: 0, facet: cx.verts[1][facet][0], cocofacet: cofacet })
        } else {
            let tri = cx.inc.get(1).and_then(|l| l.iter().filter(|b| b.facet == cofacet).map(|b| b.cofacet).min());
            tri.map_or(Ok(()), |t| Err(ComplexError::BoundaryOfBoundaryNonzero { dim: 0, facet, cocofacet: t }))
        };
        prop_assert_eq!(cx.validate(), expected);
    }
}
