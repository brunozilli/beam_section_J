markdown
# Torsional Constant J - Finite Element Calculation

**MIT License**

Copyright (c) 2024 Bruno Zilli & DeepSeek

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

---

## Overview

This document describes the finite element implementation for calculating the **torsional constant `J`** of arbitrary beam cross-sections. The method solves Prandtl's membrane analogy using Poisson's equation on a 2D triangular mesh read from UNV format.

**Repository:** https://github.com/brunozilli/beam_section_J

**Authors:** Bruno Zilli & DeepSeek

**File:** `src/torsion_j.f`

---

## Theoretical Background

### Prandtl's Membrane Analogy

For a beam under torsion, the stress function `φ(x,y)` satisfies Poisson's equation:
∂²φ/∂x² + ∂²φ/∂y² = -2Gθ

text

where:
- `G` is the shear modulus
- `θ` is the twist angle per unit length

By setting `φ = Gθ ψ`, the equation simplifies to:
∇²ψ = -2

text

with boundary condition `ψ = 0` on the cross-section boundary.

### Torsional Constant J

The torsional constant `J` is related to the stress function by:
J = 2 ∫_Ω ψ dΩ

text

Once `ψ` is computed via FEM, `J` is obtained by integrating `ψ` over the cross-sectional area.

---

## Finite Element Formulation

### Mesh and Element Type

- **Element type:** 3-node linear triangle
- **Degrees of freedom:** One value of `ψ` per node
- **Boundary condition:** `ψ = 0` on all boundary nodes

### Stiffness Matrix Assembly

For each triangular element, the element stiffness matrix `Ke(3,3)` is computed as:
Ke(i,j) = ∫_Ω (∇N_i · ∇N_j) dΩ

text

where `N_i` are linear shape functions.

For a linear triangle, the gradient is constant over the element:
Ke(i,j) = A_e · (∇N_i · ∇N_j)

text

where `A_e` is the element area.

#### Shape Function Derivatives

For a triangle with vertices `(y1,z1)`, `(y2,z2)`, `(y3,z3)`:
∇N1 = ( (z2-z3)/detJ, (y3-y2)/detJ )
∇N2 = ( (z3-z1)/detJ, (y1-y3)/detJ )
∇N3 = ( (z1-z2)/detJ, (y2-y1)/detJ )

text

where:
detJ = (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1)
A_e = |detJ| / 2

text

### Right-Hand Side Assembly

The right-hand side vector `RHS` comes from the constant forcing term `f = -2`:
RHS(i) = ∫_Ω N_i · (-2) dΩ = -2 · A_e / 3

text

for each node of the element (lumped integration).

### Boundary Conditions

Dirichlet boundary conditions `ψ = 0` are applied on all boundary nodes. The implementation:

1. Identifies boundary nodes (nodes with fewer than 3 adjacent elements)
2. Modifies the global stiffness matrix:
   - Sets row and column to zero
   - Sets diagonal to 1.0
3. Sets corresponding RHS entries to 0.0

### Linear System Solution

The system `K·ψ = RHS` is solved using LAPACK routine `DPOSV`, which is optimised for symmetric positive definite matrices:

```fortran
call DPOSV('U', n_eq, 1, K, n_eq, RHS, n_eq, info)
'U' indicates the upper triangular part of K is stored

n_eq is the number of equations (equal to number of nodes)

On exit, RHS contains the solution ψ

Torsional Constant Calculation
After obtaining ψ at each node, J is computed by numerical integration:

text
J = 2 ∫_Ω ψ dΩ = 2 · Σ_e ( A_e · (ψ1 + ψ2 + ψ3)/3 )
where the sum is over all elements.

Implementation Details
Subroutines in torsion_j.f
Subroutine	Description
compute_torsional_J	Main routine: assembles K and RHS, solves system, computes J
find_boundary_nodes	Identifies boundary nodes by counting adjacent elements
stiffness_matrix_triangle	Computes Ke(3,3) for a single triangular element
rhs_triangle	Computes Re(3) for a single triangular element
Main Algorithm
text
1. Read mesh (nodes, elements)
2. Identify boundary nodes
3. Allocate K(nn, nn) and RHS(nn)
4. Initialise K and RHS to zero
5. For each element:
   a. Compute Ke(3,3)
   b. Assemble Ke into global K
   c. Compute Re(3) from f = -2
   d. Assemble Re into global RHS
6. Apply boundary conditions (ψ = 0 on boundary)
7. Solve K·ψ = RHS using DPOSV
8. Compute J = 2 * Σ( A_e * average(ψ) )
9. Deallocate memory
Boundary Node Detection
A node is considered a boundary node if it belongs to fewer than 3 elements. For a well-formed 2D triangular mesh:

Interior nodes are shared by ≥ 3 elements

Boundary nodes are shared by 1 or 2 elements

fortran
do i = 1, nn
   if (node_count(i) .lt. 3) then
      n_boundary = n_boundary + 1
      boundary_nodes(n_boundary) = i
   end if
end do
Integration with Existing Code
Mesh Reading (read_section_mesh_unv.f)
The UNV reader extracts:

Node coordinates (y, z)

Triangular element connectivity

Mesh Validation (mesh_checker.f)
Validates and optionally corrects:

Element orientation (ensures positive area)

Detects degenerate elements (zero area)

Section Properties (compute_section_properties.f)
Computes:

Cross-sectional area A

Centroid coordinates (y_c, z_c)

Second moments of area I_y, I_z, I_yz

Database (section_database.f)
Stores computed properties including J for multiple sections.

Fortran Programming Notes
Fixed Format Convention
All source files use traditional Fortran fixed format (.f extension):

Column(s)	Usage
1	c or * for comment lines
2-5	Statement label (numeric, optional)
6	Continuation character (& or any non-space, non-zero)
7-72	Fortran code or comments
73+	Ignored by compiler
Code Style
British English in comments (colour, centre, licence, etc.)

Explicit implicit none at the beginning of each subroutine

LAPACK for linear algebra (DPOSV for SPD systems)

Compilation
The Makefile includes rules for fixed format:

makefile
FFLAGS = -ffixed-form -Wall -O2
LDFLAGS = -llapack -lblas
Example Compilation and Test
bash
make clean
make test_torsion
./test_torsion
Verification Examples
Rectangular Section (10×20)
For a rectangle with dimensions b = 10, h = 20:

Method	J
Analytical (thin-wall approx)	J ≈ b·h³/3 = 10·8000/3 = 26666.7
Analytical (exact)	J = b·h³·[1/3 - 0.21·(b/h)·(1 - b⁴/(12h⁴))] ≈ 26400
FEM (fine mesh)	Should converge to exact value
Circular Section (diameter 10)
For a circle with radius R = 5:

Method	J
Analytical	J = π·R⁴/2 = π·625/2 = 981.75
FEM (fine mesh)	Should converge to analytical value
HEB200 (Steel Section)
For a standard HEB200 profile, the FEM result should match published values from Eurocode or steel tables.

File Structure
text
beam_section_J/
├── LICENSE
├── Makefile
├── THEORY.md
├── TORSIONAL_CONSTANT_J_FEM.md    (this file)
├── meshes/
│   ├── circle_dia10.unv
│   ├── HEB200_mm.unv
│   └── rect_10x20.unv
├── src/
│   ├── compute_section_properties.f
│   ├── mesh_checker.f
│   ├── read_section_mesh_unv.f
│   ├── section_database.f
│   └── torsion_j.f
└── test/
    └── test_torsion.f
Next Steps (Future Work)
Shear centre calculation (y_s, z_s) using the same FEM framework with modified RHS (f = z and f = -y)

Constitutive matrix D(6,6) for Timoshenko beam elements

Shear correction factors (k_y, k_z) via advanced FEM techniques

References
Prandtl, L. (1903). Zur torsion von prismatischen stäben

Timoshenko, S. P. & Goodier, J. N. (1970). Theory of Elasticity

Cook, R. D. et al. (2002). Concepts and Applications of Finite Element Analysis

LAPACK User's Guide (DPOSV documentation)
