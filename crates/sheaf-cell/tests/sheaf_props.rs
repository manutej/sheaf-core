use faer::Mat;
use proptest::prelude::*;
use sheaf_cell::sheaf::{Sheaf, SheafError};

/// Edge j = (u, v, d_e, ku, kv) on n vertices, u != v. F_u is d_e x (dims[u] + [ku == 0]); F_v is
/// (d_e + [kv == 1]) x (dims[v] + [kv == 0]); any other kind is the well-shaped map.
type Spec = (usize, usize, usize, u8, u8);

fn arb_graph() -> impl Strategy<Value = (Vec<usize>, Vec<Spec>)> {
    (2usize..6).prop_flat_map(|n| {
        let e = (0..n, 1..n, 0usize..3, 0u8..6, 0u8..8).prop_map(move |(u, k, de, ku, kv)| (u, (u + k) % n, de, ku, kv));
        (prop::collection::vec(0usize..3, n), prop::collection::vec(e, 0..8))
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]
    /// `Sheaf::from_graph(dims, E)` with valid, loop-free endpoints is Ok iff every edge has F_u: R^{dims[u]} -> R^{d_e}
    /// and F_v: R^{dims[v]} -> R^{d_e} (d_e = rows of F_u). Otherwise it is `MapShape` naming the FIRST bad map in the
    /// order F_u(e_0), F_v(e_0), F_u(e_1), ... (incidence 2j / 2j+1) with the expected and actual shape. On Ok,
    /// C^0 = sum dims, C^1 = sum d_e, offsets are their prefix sums, and δ⁰ has shape C^1 x C^0.
    #[test]
    fn from_graph_accepts_exactly_the_well_shaped_maps((dims, es) in arb_graph()) {
        let mats: Vec<_> = es.iter().map(|&(u, v, de, ku, kv)| (
            u, v,
            Mat::<f64>::zeros(de, dims[u] + usize::from(ku == 0)),
            Mat::<f64>::zeros(de + usize::from(kv == 1), dims[v] + usize::from(kv == 0)),
        )).collect();
        let first_bad = mats.iter().enumerate().find_map(|(j, (u, v, fu, fv))| {
            let de = fu.nrows();
            if fu.ncols() != dims[*u] { return Some((2 * j, (de, dims[*u]), (fu.nrows(), fu.ncols()))); }
            if (fv.nrows(), fv.ncols()) != (de, dims[*v]) { return Some((2 * j + 1, (de, dims[*v]), (fv.nrows(), fv.ncols()))); }
            None
        });
        match (Sheaf::from_graph(&dims, &mats), first_bad) {
            (Ok(s), None) => {
                let des: Vec<usize> = es.iter().map(|e| e.2).collect();
                let prefix = |xs: &[usize]| std::iter::once(0).chain(xs.iter().scan(0, |a, &x| { *a += x; Some(*a) })).collect::<Vec<_>>();
                prop_assert_eq!(s.offsets(0), prefix(&dims));
                prop_assert_eq!(s.offsets(1), prefix(&des));
                prop_assert_eq!(s.cochain_dim(0), dims.iter().sum::<usize>());
                prop_assert_eq!(s.cochain_dim(1), des.iter().sum::<usize>());
                prop_assert_eq!(s.cochain_dim(2), 0);
                let d0 = s.coboundary_dense(0);
                prop_assert_eq!((d0.nrows(), d0.ncols()), (s.cochain_dim(1), s.cochain_dim(0)));
            }
            (Err(e), Some((incidence, expected, got))) => {
                prop_assert_eq!(e, SheafError::MapShape { dim: 0, incidence, expected, got });
            }
            (r, want) => prop_assert!(false, "from_graph = {:?}, expected first bad map {:?}", r.map(|_| ()), want),
        }
    }
}
