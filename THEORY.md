
# Torsional Constant J: Theory and Implementation

## Directory Structure
caelinux@caelinux-System-Product-Name:~/Documenti/beam_section_J$ tree
.
├── Makefile
├── meshes
│   ├── circle_dia10.unv
│   ├── HEB200_mm.unv
│   └── rect_10x20.unv
├── section_database.mod
├── src
│   ├── compute_section_properties.f
│   ├── compute_section_properties.o
│   ├── mesh_checker.f
│   ├── mesh_checker.o
│   ├── read_section_mesh_unv.f
│   ├── read_section_mesh_unv.o
│   ├── section_database.f
│   ├── section_database.o
│   ├── torsion_j.f
│   └── torsion_j.o
├── test
│   └── test_torsion.f
├── test_torsion
└── torsion_j.mod

3 directories, 18 files


text

## Part 1: Theory - Why ∇²φ = -2?

### The Saint-Venant Torsion Problem

When a cylindrical beam is twisted, each cross-section rotates and warps. Saint-Venant (1855) proposed that the out-of-plane displacement (warping) is:
w(x,y,z) = θ · Ψ(x,y)

text

where:
- `θ` = twist angle per unit length (radians/mm)
- `Ψ(x,y)` = warping function (unknown)

### Stresses from Hooke's Law

The shear stresses in the cross-section are:
τ_xz = G · θ · (∂Ψ/∂x - y)
τ_yz = G · θ · (∂Ψ/∂y + x)

text

where `G` is the shear modulus.

### Prandtl's Genius Idea (1903)

Instead of solving for `Ψ`, Prandtl introduced a **stress function** `φ(x,y)` such that:
τ_xz = ∂φ/∂y
τ_yz = -∂φ/∂x

text

This choice **automatically satisfies equilibrium** (the stress equations become identities).

### The Compatibility Condition

For the stresses to be compatible with a displacement field, `φ` must satisfy:
∂²φ/∂x² + ∂²φ/∂y² = -2Gθ

text

That is:
∇²φ = -2Gθ

text

### Normalisation

Let `φ̄ = φ / (Gθ)`. The equation becomes:
∇²φ̄ = -2

text

This is **Poisson's equation** with a constant source term.

### Boundary Condition

On the boundary of the section, there are no applied forces, which gives:
φ̄ = 0 on the boundary

text

### Computing J

The torque `M_t` is obtained by integrating shear stresses:
M_t = ∫ (τ_yz·x - τ_xz·y) dA = 2 ∫ φ̄ dA

text

But the torque-twist relation is:
M_t = G · J · θ

text

Therefore:
J = 2 ∫ φ̄ dA

text

**This is the key formula:** once we know `φ̄` everywhere, we integrate it over the area and multiply by 2.

---

## Part 2: Finite Element Implementation

### Discretisation

The cross-section is divided into **linear triangular elements** (the mesh). Over each triangle, we approximate `φ̄` as a linear function:
φ̄(x,y) = N₁·φ₁ + N₂·φ₂ + N₃·φ₃

text

where:
- `φ₁, φ₂, φ₃` = nodal values (unknowns)
- `N₁, N₂, N₃` = linear shape functions

### Element Stiffness Matrix (Ke)

For Poisson's equation, the element stiffness matrix is:
Ke(i,j) = (area / (4·area)) · (b_i·b_j + c_i·c_j)

text

where the coefficients come from the triangle geometry:
b₁ = y₂ - y₃ c₁ = z₃ - z₂
b₂ = y₃ - y₁ c₂ = z₁ - z₃
b₃ = y₁ - y₂ c₃ = z₂ - z₁

text

These coefficients are constant for a given triangle. The factor `area/(4·area)` simplifies to `1/4`, but we keep it for numerical stability.

### Element Load Vector (fe)

For the source term `-2` (right-hand side), the load vector is:
fe(i) = 2 · area / 3

text

for each node `i = 1, 2, 3`. This comes from integrating the constant source over the triangle.

### Global Assembly

We:
1. Identify **boundary nodes** (where `φ̄ = 0`)
2. Build a mapping from global node numbers to "free node" indices
3. For each element, add `Ke(i,j)` to the global matrix `K(free_i, free_j)`
4. Add `fe(i)` to the global vector `RHS(free_i)`

### Solving the System

We obtain:
[K] · {φ̄_free} = {RHS}

text

where:
- `[K]` is symmetric positive definite
- Size = number of free nodes (interior nodes)

We solve using LAPACK's `DPOSV`, which performs Cholesky factorisation:
K = L·Lᵀ then L·Lᵀ·φ̄ = RHS

text

This is direct, robust, and efficient for our problem sizes.

### Computing J from φ̄

After solving, we have `φ̄` at all nodes (zero on boundary). Then:
J = 2 · Σ (area_elem / 3) · (φ̄₁ + φ̄₂ + φ̄₃)

text

over all triangles.

---

## Part 3: Code Walkthrough (torsion_j.f)

### Module Declaration
```fortran
module torsion_j
    implicit none
contains
...
end module torsion_j
Main Subroutine
fortran
subroutine compute_torsion_J(nnode, ntri, conn, y, z, J_torsion, verbose)
Input:

nnode, ntri: mesh dimensions

conn(3, ntri): triangle connectivity

y(nnode), z(nnode): nodal coordinates

verbose: optional debug output

Output:

J_torsion: torsional constant (mm⁴)

Step 1: Find Boundary Nodes
fortran
call find_boundary_nodes(nnode, ntri, conn, bound, verb)
An edge is a boundary edge if it appears in only one triangle. Nodes on such edges are boundary nodes.

Step 2: Count Free Nodes
fortran
nfree = 0
do i = 1, nnode
    if (bound(i) .eq. 0) nfree = nfree + 1
end do
Free nodes = all nodes not on boundary (interior nodes).

Step 3: Build Global Mapping
fortran
allocate(map(nnode))
j = 0
do i = 1, nnode
    if (bound(i) .eq. 0) then
        j = j + 1
        map(i) = j
    else
        map(i) = 0
    end if
end do
map(i) gives the free index (1..nfree) for node i, or 0 if boundary.

Step 4: Analytical Fallback
fortran
if (nfree .eq. 0) then
    call analytical_J_approx(...)
    return
end if
For very coarse meshes (e.g., rectangle with only corner nodes), we use a simple formula based on bounding box.

Step 5: Assemble Global System
fortran
do k = 1, ntri
    ! Get triangle vertices
    i1 = conn(1,k); i2 = conn(2,k); i3 = conn(3,k)
    
    ! Compute area
    area = 0.5d0 * abs((y2-y1)*(z3-z1) - (y3-y1)*(z2-z1))
    
    ! Compute b, c coefficients
    b1 = y2 - y3; c1 = z3 - z2
    b2 = y3 - y1; c2 = z1 - z3
    b3 = y1 - y2; c3 = z2 - z1
    
    ! Build element stiffness matrix Ke(3,3)
    Ke(1,1) = (b1*b1 + c1*c1) / (4.0d0*area)
    Ke(1,2) = (b1*b2 + c1*c2) / (4.0d0*area)
    ...
    
    ! Build element load vector fe(3)
    fe(1) = 2.0d0 * area / 3.0d0
    ...
    
    ! Assemble into global matrix
    do i = 1, 3
        m = map(conn(i,k))
        if (m > 0) then
            rhs(m) = rhs(m) + fe(i)
            do j = 1, 3
                n = map(conn(j,k))
                if (n > 0) then
                    Kmat(m,n) = Kmat(m,n) + Ke(i,j)
                end if
            end do
        end if
    end do
end do
Step 6: Add Diagonal Stabilisation
fortran
do i = 1, nfree
    Kmat(i,i) = Kmat(i,i) + 1.0d-12
end do
This ensures the matrix is positive definite (avoids zero pivots).

Step 7: Solve with LAPACK
fortran
phi_free = rhs
call DPOSV('U', nfree, 1, Kmat, nfree, phi_free, nfree, info)
'U': use upper triangular part of Kmat

nfree: matrix size

1: number of right-hand sides

Kmat: on entry = matrix, on exit = Cholesky factor

phi_free: on entry = RHS, on exit = solution

info: 0 = success

Step 8: Reconstruct φ̄ for All Nodes
fortran
phi = 0.0d0
do i = 1, nnode
    if (bound(i) .eq. 0) then
        phi(i) = phi_free(map(i))
    end if
end do
Step 9: Compute J
fortran
J_torsion = 0.0d0
do k = 1, ntri
    i1 = conn(1,k); i2 = conn(2,k); i3 = conn(3,k)
    area = ... (as before)
    J_torsion = J_torsion + (area / 3.0d0) * (phi(i1) + phi(i2) + phi(i3))
end do
J_torsion = 2.0d0 * J_torsion
Part 4: Validation
Section	Mesh	Expected J (mm⁴)	Our Result	Error
Circle Ø10	12 nodes, 10 triangles	981.75	738.4	-25%
HEB200	335 nodes, 472 triangles	1,200,000	422,560	-65%
Sources of error:

Coarse meshes (too few interior nodes)

Linear approximation of φ̄ (exact solution is quadratic for circle)

Mesh refinement needed for accurate J

To improve:

Generate finer meshes (more triangles)

Use higher-order elements (quadratic triangles)

Implement adaptive mesh refinement

References
Saint-Venant, B. (1855). "Mémoire sur la torsion des prismes"

Prandtl, L. (1903). "Zur Torsion von prismatischen Stäben"

Timoshenko, S.P., Goodier, J.N. (1970). "Theory of Elasticity"

Zienkiewicz, O.C., Taylor, R.L. (2000). "The Finite Element Method"

Authors
Bruno Zilli & DeepSeek

Licence
MIT
EOF
