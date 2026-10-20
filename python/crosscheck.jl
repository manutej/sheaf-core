# Cross-check graph fixtures against AlgebraicJulia CellularSheaves.jl (graph-only: H0 via nullspace_ldlt).
# H1 on a graph = dim C1 - rank d0 = dim C1 - (dim C0 - H0).  Rust expectations hard-coded from fixtures.rs.
using CellularSheaves, Graphs, LinearAlgebra
const CS = CellularSheaves
rot(t) = [cos(t) -sin(t); sin(t) cos(t)]
function build(n, edges, d, twist)   # twist = (edge (a,b), matrix) applied to the map from b into that edge
    s = CS.EuclideanSheaf{Float64}(fill(d, n))
    for (a, b) in edges
        r2 = ((a, b) == twist[1]) ? twist[2] : Matrix{Float64}(I, d, d)
        CS.add_sheaf_edge!(s, a, b, Matrix{Float64}(I, d, d), r2)
    end
    s
end
Id2 = Matrix{Float64}(I, 2, 2)
tri = [(1,2),(1,3),(2,3)]; k4 = [(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]; c5 = [(i, i % 5 + 1) for i in 1:5]
cases = [("triangle_graph_const1", build(3, tri, 1, ((0,0), zeros(1,1))), [1,1]),
         ("cycle5_const2", build(5, c5, 2, ((0,0), Id2)), [2,2]),
         ("triangle_rot60", build(3, tri, 2, ((1,2), rot(pi/3))), [0,0]),
         ("triangle_reflect", build(3, tri, 2, ((1,2), [1.0 0; 0 -1])), [1,1]),
         ("k4_rot60", build(4, k4, 2, ((1,2), rot(pi/3))), [0,4])]
allok = true
for (name, s, exp) in cases
    L = Matrix(CS.sheaf_laplacian_matrix(s))
    c0 = size(L, 1); c1 = sum(CS.edge_stalk_dimensions(s))
    h0_ldlt = size(CS.nullspace_ldlt(s), 2)
    h0_svd = c0 - rank(L)
    b = [h0_ldlt, c1 - (c0 - h0_ldlt)]
    ok = b == exp && h0_ldlt == h0_svd
    global allok &= ok
    println(rpad(name, 24), " CellularSheaves.jl H0(ldlt)=", h0_ldlt, " H0(svd)=", h0_svd, " -> ", b, " rust=", exp, ok ? " OK" : " MISMATCH")
end
println(allok ? "ALL_OK" : "FAIL")
