//! Regular cell complexes, stored as signed incidences between consecutive dimensions.
//!
//! ASSUMPTIONS (A1-A3, see DESIGN in research/rust-cohomology-design.md):
//! A1 finite, regular cell complex of dimension <= 2 in the spike (code is k-generic);
//! A2 every cell carries a fixed orientation; the incidence number [s:t] in {-1,+1}
//!    is supplied by the builder (simplicial builder derives it from sorted vertex order);
//! A3 boundary of boundary is zero: sum_t [s:t][t:r] = 0 for every s<t<r (checked).

#[derive(Debug, Clone, PartialEq)]
pub struct Incidence {
    /// index of the (k)-cell (the face, "sigma")
    pub facet: usize,
    /// index of the (k+1)-cell (the coface, "tau")
    pub cofacet: usize,
    /// incidence number [sigma : tau]
    pub sign: i8,
}

#[derive(Debug, Clone)]
pub struct CellComplex {
    /// n_cells[k] = number of k-cells
    pub n_cells: Vec<usize>,
    /// inc[k] = incidences between k-cells and (k+1)-cells
    pub inc: Vec<Vec<Incidence>>,
    /// optional vertex labels of each cell (simplicial builder fills these)
    pub verts: Vec<Vec<Vec<usize>>>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ComplexError {
    BadIndex { dim: usize, index: usize },
    BoundaryOfBoundaryNonzero { dim: usize, facet: usize, cocofacet: usize },
    SelfLoop { edge: usize, vertex: usize },
    /// `inc` has `inc_levels` incidence levels but `n_cells` only has `n_levels` entries; a valid complex needs
    /// `inc.len() < n_cells.len()` (inc[k] links k-cells to (k+1)-cells).
    LevelMismatch { inc_levels: usize, n_levels: usize },
}

impl CellComplex {
    pub fn dim(&self) -> usize { self.n_cells.len().saturating_sub(1) }

    /// Graph (1-complex) that KEEPS the given edge order, multi-edges and orientation: edge j = (u, v) gets incidences
    /// [u : e_j] = -1 and [v : e_j] = +1 even when u > v, so (delta x)_j = F_v x_v - F_u x_u exactly as in m1/v1
    /// (SPEC R-COB-06, section 3k). `simplicial` instead sorts, dedups and orients by vertex order. Self-loops are rejected.
    pub fn graph(n_vertices: usize, edges: &[(usize, usize)]) -> Result<Self, ComplexError> {
        let mut inc0 = Vec::with_capacity(2 * edges.len());
        for (j, &(u, v)) in edges.iter().enumerate() {
            if u >= n_vertices { return Err(ComplexError::BadIndex { dim: 0, index: u }); }
            if v >= n_vertices { return Err(ComplexError::BadIndex { dim: 0, index: v }); }
            if u == v { return Err(ComplexError::SelfLoop { edge: j, vertex: u }); }
            inc0.push(Incidence { facet: u, cofacet: j, sign: -1 });
            inc0.push(Incidence { facet: v, cofacet: j, sign: 1 });
        }
        let verts = vec![(0..n_vertices).map(|v| vec![v]).collect(), edges.iter().map(|&(u, v)| vec![u, v]).collect()];
        Ok(CellComplex { n_cells: vec![n_vertices, edges.len()], inc: vec![inc0], verts })
    }

    /// Simplicial 2-complex (or graph if `faces` is empty) from vertex count, edges and triangles.
    /// Edges implied by triangles are added automatically. Orientation = increasing vertex order.
    pub fn simplicial(n_vertices: usize, edges: &[(usize, usize)], faces: &[[usize; 3]]) -> Self {
        let mut es: Vec<(usize, usize)> = edges.iter().map(|&(a, b)| (a.min(b), a.max(b))).collect();
        for f in faces {
            let mut t = *f; t.sort();
            es.push((t[0], t[1])); es.push((t[0], t[2])); es.push((t[1], t[2]));
        }
        es.sort(); es.dedup();
        let eidx = |a: usize, b: usize| es.binary_search(&(a, b)).expect("edge present");
        let mut inc0 = Vec::new();
        for (j, &(a, b)) in es.iter().enumerate() {
            // boundary(e=[a,b]) = b - a
            inc0.push(Incidence { facet: a, cofacet: j, sign: -1 });
            inc0.push(Incidence { facet: b, cofacet: j, sign: 1 });
        }
        let mut fs: Vec<[usize; 3]> = faces.iter().map(|f| { let mut t = *f; t.sort(); t }).collect();
        fs.sort(); fs.dedup();
        let mut inc1 = Vec::new();
        for (j, t) in fs.iter().enumerate() {
            // boundary([a,b,c]) = [b,c] - [a,c] + [a,b]
            inc1.push(Incidence { facet: eidx(t[1], t[2]), cofacet: j, sign: 1 });
            inc1.push(Incidence { facet: eidx(t[0], t[2]), cofacet: j, sign: -1 });
            inc1.push(Incidence { facet: eidx(t[0], t[1]), cofacet: j, sign: 1 });
        }
        let mut n_cells = vec![n_vertices, es.len()];
        let mut inc = vec![inc0];
        let mut verts = vec![(0..n_vertices).map(|v| vec![v]).collect(), es.iter().map(|&(a, b)| vec![a, b]).collect()];
        if !fs.is_empty() {
            n_cells.push(fs.len());
            inc.push(inc1);
            verts.push(fs.iter().map(|t| t.to_vec()).collect());
        }
        CellComplex { n_cells, inc, verts }
    }

    /// Check indices and boundary-of-boundary = 0 (A3). Pure; returns typed error.
    pub fn validate(&self) -> Result<(), ComplexError> {
        // Fields are pub, so a hand-built complex may carry more incidence levels than cell counts: reject it instead
        // of panicking on self.n_cells[k + 1] below.
        if self.inc.len() >= self.n_cells.len() && !self.inc.is_empty() {
            return Err(ComplexError::LevelMismatch { inc_levels: self.inc.len(), n_levels: self.n_cells.len() });
        }
        for (k, list) in self.inc.iter().enumerate() {
            for i in list {
                if i.facet >= self.n_cells[k] { return Err(ComplexError::BadIndex { dim: k, index: i.facet }); }
                if i.cofacet >= self.n_cells[k + 1] { return Err(ComplexError::BadIndex { dim: k + 1, index: i.cofacet }); }
            }
        }
        for k in 0..self.inc.len().saturating_sub(1) {
            use std::collections::HashMap;
            let mut acc: HashMap<(usize, usize), i32> = HashMap::new();
            for a in &self.inc[k] {
                for b in self.inc[k + 1].iter().filter(|b| b.facet == a.cofacet) {
                    *acc.entry((a.facet, b.cofacet)).or_default() += (a.sign as i32) * (b.sign as i32);
                }
            }
            if let Some((&(s, r), _)) = acc.iter().find(|(_, &v)| v != 0) {
                return Err(ComplexError::BoundaryOfBoundaryNonzero { dim: k, facet: s, cocofacet: r });
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn inc(facet: usize, cofacet: usize, sign: i8) -> Incidence { Incidence { facet, cofacet, sign } }

    /// Regression (sheaf-validate-index-001 observation): more incidence levels than cell counts used to panic.
    #[test]
    fn validate_rejects_more_inc_levels_than_cell_counts() {
        // one level of cell counts, one level of incidences: inc[0] would need n_cells[1]
        let c = CellComplex { n_cells: vec![2], inc: vec![vec![inc(0, 0, -1), inc(1, 0, 1)]], verts: vec![] };
        assert_eq!(c.validate(), Err(ComplexError::LevelMismatch { inc_levels: 1, n_levels: 1 }));
        // no cell counts at all
        let c = CellComplex { n_cells: vec![], inc: vec![vec![]], verts: vec![] };
        assert_eq!(c.validate(), Err(ComplexError::LevelMismatch { inc_levels: 1, n_levels: 0 }));
        // graph shape with an extra (empty) incidence level
        let c = CellComplex { n_cells: vec![2, 1], inc: vec![vec![inc(0, 0, -1), inc(1, 0, 1)], vec![]], verts: vec![] };
        assert_eq!(c.validate(), Err(ComplexError::LevelMismatch { inc_levels: 2, n_levels: 2 }));
    }

    #[test]
    fn validate_still_accepts_consistent_shapes() {
        let c = CellComplex { n_cells: vec![2, 1], inc: vec![vec![inc(0, 0, -1), inc(1, 0, 1)]], verts: vec![] };
        assert_eq!(c.validate(), Ok(()));
        let c = CellComplex { n_cells: vec![3], inc: vec![], verts: vec![] };
        assert_eq!(c.validate(), Ok(()));
        let c = CellComplex { n_cells: vec![], inc: vec![], verts: vec![] };
        assert_eq!(c.validate(), Ok(()));
        assert_eq!(CellComplex::simplicial(3, &[], &[[0, 1, 2]]).validate(), Ok(()));
        // index errors keep their order and payload
        let c = CellComplex { n_cells: vec![2, 1], inc: vec![vec![inc(0, 1, -1)]], verts: vec![] };
        assert_eq!(c.validate(), Err(ComplexError::BadIndex { dim: 1, index: 1 }));
    }
}
