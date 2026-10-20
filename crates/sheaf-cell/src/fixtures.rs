//! Small complexes with known answers (used by tests, the reference dump and the Python cross-check).
use crate::complex::CellComplex;
use crate::sheaf::Sheaf;
use faer::Mat;

pub fn triangle_graph() -> CellComplex { CellComplex::simplicial(3, &[(0, 1), (1, 2), (0, 2)], &[]) }
pub fn triangle_filled() -> CellComplex { CellComplex::simplicial(3, &[], &[[0, 1, 2]]) }
pub fn cycle_graph(n: usize) -> CellComplex {
    let e: Vec<(usize, usize)> = (0..n).map(|i| (i, (i + 1) % n)).collect();
    CellComplex::simplicial(n, &e, &[])
}
pub fn k4_graph() -> CellComplex {
    CellComplex::simplicial(4, &[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)], &[])
}
/// boundary of the tetrahedron = S^2: betti (1,0,1)
pub fn sphere() -> CellComplex { CellComplex::simplicial(4, &[], &[[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]) }
/// 7-vertex (Moebius-Csaszar) torus: betti (1,2,1)
pub fn torus7() -> CellComplex {
    let mut f = Vec::new();
    for i in 0..7 { f.push([i, (i + 1) % 7, (i + 3) % 7]); f.push([i, (i + 2) % 7, (i + 3) % 7]); }
    CellComplex::simplicial(7, &[], &f)
}
/// 6-vertex RP^2: betti over R = (1,0,0); over F_2 = (1,1,1)
pub fn rp2() -> CellComplex {
    let f = [[1,2,3],[1,2,4],[1,3,5],[1,4,6],[1,5,6],[2,3,6],[2,4,5],[2,5,6],[3,4,5],[3,4,6]];
    let f: Vec<[usize; 3]> = f.iter().map(|t| [t[0] - 1, t[1] - 1, t[2] - 1]).collect();
    CellComplex::simplicial(6, &[], &f)
}
pub fn rotation(theta: f64) -> Mat<f64> {
    let (c, s) = (theta.cos(), theta.sin());
    Mat::from_fn(2, 2, |i, j| [[c, -s], [s, c]][i][j])
}
/// R^2 sheaf on `cx` (graph part) that is constant except: the map from vertex `v` into edge `e` is `r`.
/// Holonomy around any cycle through e = r^{+-1}. If r has no eigenvalue 1, no global section passes.
pub fn with_twist(cx: CellComplex, e: usize, v: usize, r: Mat<f64>) -> Sheaf {
    let mut s = Sheaf::constant(cx, r.nrows());
    let j = s.complex.inc[0].iter().position(|i| i.cofacet == e && i.facet == v).expect("incidence");
    s.maps[0][j] = r;
    s
}

pub struct Named { pub name: &'static str, pub sheaf: Sheaf, pub expect_r: Vec<usize> }

pub fn catalogue() -> Vec<Named> {
    let pi3 = std::f64::consts::PI / 3.0;
    vec![
        Named { name: "triangle_graph_const1", sheaf: Sheaf::constant(triangle_graph(), 1), expect_r: vec![1, 1] },
        Named { name: "triangle_filled_const1", sheaf: Sheaf::constant(triangle_filled(), 1), expect_r: vec![1, 0, 0] },
        Named { name: "cycle5_const2", sheaf: Sheaf::constant(cycle_graph(5), 2), expect_r: vec![2, 2] },
        Named { name: "sphere_const1", sheaf: Sheaf::constant(sphere(), 1), expect_r: vec![1, 0, 1] },
        Named { name: "torus7_const1", sheaf: Sheaf::constant(torus7(), 1), expect_r: vec![1, 2, 1] },
        Named { name: "rp2_const1", sheaf: Sheaf::constant(rp2(), 1), expect_r: vec![1, 0, 0] },
        // holonomy = rotation by 60 deg: kills H0 (2 -> 0); Euler char 0 forces H1 2 -> 0 too
        Named { name: "triangle_rot60", sheaf: with_twist(triangle_graph(), 0, 1, rotation(pi3)), expect_r: vec![0, 0] },
        // reflection holonomy (eigs +1,-1): one direction survives
        Named { name: "triangle_reflect", sheaf: with_twist(triangle_graph(), 0, 1, Mat::from_fn(2, 2, |i, j| if i == j { if i == 0 { 1.0 } else { -1.0 } } else { 0.0 })), expect_r: vec![1, 1] },
        // K4 graph, chi = 8 - 12 = -4: rotation kills H0 and H1 stays 4 (nonzero loop-level space)
        Named { name: "k4_rot60", sheaf: with_twist(k4_graph(), 0, 1, rotation(pi3)), expect_r: vec![0, 4] },
    ]
}
