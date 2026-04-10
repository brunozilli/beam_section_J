markdown
# Shear Centre and Constitutive Matrix D(6,6) Implementation

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

This document summarises the implementation of shear centre calculation and the 6×6 constitutive matrix `D` for Timoshenko beam sections. The work extends an existing Fortran 90/95 tool (fixed format `.f`) that reads UNV meshes and computes cross-sectional properties including the torsional constant `J` via finite element method.

**Repository:** https://github.com/brunozilli/beam_section_J

**Authors:** Bruno Zilli & DeepSeek

---

## Completed Work

### 1. Shear Centre Calculation (`src/shear_center.f`)

The shear centre `(y_s, z_s)` is computed by solving two Poisson problems on the same finite element mesh used for the torsional constant `J`.

#### Theoretical Background

For a beam cross-section subjected to shear forces `V_y` and `V_z`, the shear centre is the point where an applied shear force produces no twisting moment. The location is found by solving:

| Case | Condition | Right-hand side `f` | Result |
|------|-----------|---------------------|--------|
| 1 | `V_y = 1, V_z = 0` | `f = z` | `M_torque1 = y_s` |
| 2 | `V_y = 0, V_z = 1` | `f = -y` | `M_torque2 = -z_s` |

The governing equation is Poisson's equation:
∇²φ = f

text

with Dirichlet boundary conditions `φ = 0` on the boundary of the cross-section.

#### Finite Element Formulation

**Stiffness matrix `K`:** Same as for the torsion problem:
K(i,j) = ∫_Ω (∇N_i · ∇N_j) dΩ

text

where `N_i` are linear shape functions for triangular elements.

**Right-hand side vector `RHS`:** For each case:
RHS(i) = ∫_Ω N_i · f dΩ

text

**Solution:** The linear system `K·φ = RHS` is solved using LAPACK routine `DPOSV` (for symmetric positive definite matrices).

**Torque calculation:**
M_torque = ∫_Ω φ dΩ

text

**Shear centre coordinates:**
y_s = M_torque1
z_s = -M_torque2

text

#### Implementation Details

- **Element type:** 3-node linear triangle
- **Integration:** One-point (centroid) integration for RHS; exact integration for stiffness matrix (constant gradient over element)
- **Boundary condition identification:** Nodes with fewer than 3 adjacent elements are considered boundary nodes
- **Solver:** LAPACK `DPOSV` (Cholesky factorisation for SPD matrices)

#### Subroutines in `shear_center.f`

| Subroutine | Description |
|------------|-------------|
| `compute_shear_center` | Main routine: assembles K, solves for both cases, computes y_s, z_s |
| `find_boundary_nodes` | Identifies boundary nodes by counting adjacent elements |

#### Input/Output

**Input:**
- `nn` - number of nodes
- `ne` - number of elements
- `nodes(nn,2)` - nodal coordinates (y, z)
- `elements(ne,3)` - element connectivity

**Output:**
- `y_s` - shear centre y-coordinate
- `z_s` - shear centre z-coordinate

---

### 2. Constitutive Matrix D(6,6) (`src/build_D_full.f`)

The 6×6 constitutive matrix for a Timoshenko beam element relates generalised strains to generalised forces:
[F] = [D] · [ε]

text

where:
- `F = [N, V_y, V_z, M_x, M_z, M_y]^T` (axial force, shear forces, torque, bending moments)
- `ε = [ε_x, γ_y, γ_z, φ_x, κ_z, κ_y]^T` (axial strain, shear strains, twist, curvatures)

#### Matrix Definition

The matrix `D` is symmetric and has the following non-zero entries:

| D(i,j) | Expression | Description |
|--------|------------|-------------|
| D(1,1) | `E·A` | Axial stiffness |
| D(2,2) | `k_y·G·A` | Shear stiffness in y-direction |
| D(3,3) | `k_z·G·A` | Shear stiffness in z-direction |
| D(4,4) | `G·J` | Torsional stiffness |
| D(5,5) | `E·I_z` | Bending stiffness about z-axis |
| D(6,6) | `E·I_y` | Bending stiffness about y-axis |
| D(1,5) | `-E·A·e_z` | Axial-bending coupling (z-offset) |
| D(1,6) | `E·A·e_y` | Axial-bending coupling (y-offset) |
| D(5,6) | `E·I_yz` | Bending-bending coupling (product of inertia) |

where:
e_y = y_s - y_c
e_z = z_s - z_c

text
are the offsets between shear centre `(y_s, z_s)` and centroid `(y_c, z_c)`.

#### Packing for CalculiX (`elcon` vector)

For use with CalculiX generalised beam elements, the matrix is packed into a 12-element vector:

| elcon index | D(i,j) | Description |
|-------------|--------|-------------|
| 1 | D(1,1) | Axial stiffness |
| 2 | D(2,2) | Shear y stiffness |
| 3 | D(3,3) | Shear z stiffness |
| 4 | D(4,4) | Torsional stiffness |
| 5 | D(5,5) | Bending Iz stiffness |
| 6 | D(6,6) | Bending Iy stiffness |
| 7 | D(1,5) | Axial-bending Iz coupling |
| 8 | D(1,6) | Axial-bending Iy coupling |
| 9 | D(5,6) | Bending-bending coupling (Iyz) |
| 10 | D(1,2) | (Typically zero) |
| 11 | D(1,3) | (Typically zero) |
| 12 | D(2,3) | (Typically zero) |

#### Subroutines in `build_D_full.f`

| Subroutine | Description |
|------------|-------------|
| `build_D_matrix` | Constructs the full 6×6 D matrix from section properties |
| `pack_elcon` | Packs D(6,6) into elcon(12) for CalculiX |
| `zero_matrix` | Utility to initialise a matrix to zero |

---

### 3. Database Extension (`src/section_database.f`)

The section database has been extended to include the following fields:

| Field | Type | Description |
|-------|------|-------------|
| `y_s`, `z_s` | `double precision` | Shear centre coordinates |
| `k_y`, `k_z` | `double precision` | Shear correction factors (default = 1.0) |
| `D(6,6)` | `double precision` | Full constitutive matrix |
| `elcon(12)` | `double precision` | Packed vector for CalculiX |

---

## Fortran Programming Notes

### Fixed Format Convention

All source files use traditional Fortran fixed format (`.f` extension):

| Column(s) | Usage |
|-----------|-------|
| 1 | `c` or `*` for comment lines |
| 2-5 | Statement label (numeric, optional) |
| 6 | Continuation character (`&` or any non-space, non-zero) |
| 7-72 | Fortran code or comments |
| 73+ | Ignored by compiler |

### Code Style

- **British English** in comments (colour, centre, licence, etc.)
- **Uppercase** for keywords (optional but traditional)
- **Explicit `implicit none`** at the beginning of each subroutine
- **LAPACK** used for linear algebra (`DPOSV` for symmetric positive definite systems)

### Compilation

The `Makefile` includes rules for fixed format:

```makefile
FFLAGS = -ffixed-form -Wall -O2
LDFLAGS = -llapack -lblas
Example Compilation and Test
bash
make clean
make test_shear_center
./test_shear_center
Expected Verification (HEB200 Section)
For a doubly symmetric section like HEB200:

y_s ≈ y_c (shear centre coincides with centroid)

z_s ≈ z_c

I_yz ≈ 0 (product of inertia is zero)

Matrix D is block-diagonal (no coupling terms)

File Structure
text
beam_section_J/
├── LICENSE
├── Makefile
├── THEORY.md
├── SHEAR_CENTER_IMPLEMENTATION.md   (this file)
├── meshes/
│   ├── circle_dia10.unv
│   ├── HEB200_mm.unv
│   └── rect_10x20.unv
├── src/
│   ├── compute_section_properties.f
│   ├── mesh_checker.f
│   ├── read_section_mesh_unv.f
│   ├── section_database.f
│   ├── torsion_j.f
│   ├── shear_center.f                (new)
│   └── build_D_full.f                (new)
└── test/
    ├── test_torsion.f
    └── test_shear_center.f           (to be created)
Next Steps (Optional)
Implement test_shear_center.f to validate the implementation on HEB200 mesh

Add advanced shear correction factor calculation (k_y, k_z) via FEM

Integrate shear centre calculation into main program

Add support for additional mesh formats
