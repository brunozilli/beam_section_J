# BEAM SECTION TOOL - Torsional Constant J from UNV meshes

**Tool for calculating the torsional constant J of beam cross-sections using the Saint-Venant torsion theory and finite element method**

> ⚠️ **STATUS: Under active development**  
> This project is in development and testing phase. Not all features are fully validated. Use with caution in production environments.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Fortran](https://img.shields.io/badge/Fortran-77-blue.svg)](https://fortran-lang.org/)
[![LAPACK](https://img.shields.io/badge/LAPACK-3.x-green.svg)](https://www.netlib.org/lapack/)
[![Status](https://img.shields.io/badge/status-development-orange.svg)]()

## Authors
- **Bruno Zilli** - Project coordination and validation
- **DeepSeek** - AI-assisted development and debugging

## License

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

## Disclaimer

**This software is provided "AS IS" without warranty of any kind.**

- The authors make no representations or warranties regarding the accuracy,
  reliability, or completeness of the software.
- The authors shall not be liable for any direct, indirect, incidental,
  special, exemplary, or consequential damages arising from the use of this
  software.
- Users are solely responsible for verifying the results before using them
  in engineering applications.
- The software is intended for research and educational purposes only.
- Commercial use requires independent validation.
- **Known limitation:** Results are mesh-dependent. Coarse meshes produce
  significant errors (see Validation section).

---

## Directory Structure
beam_section_J/
│
├── src/ # Fortran source files
│ ├── torsion_j.f # Main J calculation module
│ ├── compute_section_properties.f # Area, I_y, I_z calculation
│ ├── read_section_mesh_unv.f # UNV file reader
│ ├── section_database.f # Multi-section database
│ └── mesh_checker.f # Mesh validation and correction
│
├── test/ # Test programmes
│ └── test_torsion.f # J calculation test
│
├── meshes/ # Example meshes
│ ├── rect_10x20.unv # Rectangle 10x20 mm
│ ├── circle_dia10.unv # Circle diameter 10 mm
│ └── HEB200_mm.unv # HEB200 from Salome
│
├── Makefile # Build system
├── README_J.md # This file
└── LICENSE # MIT licence

text

---

## Theory Summary

### The Saint-Venant Torsion Problem

When a cylindrical beam is twisted, Prandtl introduced a **stress function** `φ(x,y)` that satisfies:
∇²φ̄ = -2 on Ω
φ̄ = 0 on Γ (boundary)

text

where `φ̄ = φ/(Gθ)` is the normalised stress function.

The torsional constant J is then:
J = 2 ∫_Ω φ̄ dA

text

### Finite Element Implementation

- **Elements:** Linear triangles (P1)
- **System:** K·φ̄ = RHS (symmetric positive definite)
- **Solver:** LAPACK DPOSV (Cholesky decomposition)
- **Integration:** Analytical formula over each triangle

For detailed theory, see the `THEORY.md` document.

---

## Compilation

### Prerequisites
- **gfortran** (or any Fortran 77 compiler)
- **LAPACK** and **BLAS** libraries

### Install LAPACK (if needed)
```bash
# Ubuntu/Debian
sudo apt-get install liblapack-dev libblas-dev

# Fedora/RHEL
sudo dnf install lapack-devel blas-devel

# macOS
brew install lapack
Step-by-Step Compilation
1. Open a terminal in the project directory:

bash
cd ~/Documenti/beam_section_J
2. Clean any previous build (optional but recommended):

bash
make clean
Expected output:

text
rm -f *.o *.mod test_torsion
3. Compile all modules:

bash
make
Expected output:

text
gfortran -Wall -O2 -c src/torsion_j.f -o torsion_j.o
gfortran -Wall -O2 -c test/test_torsion.f -o test_torsion.o
gfortran -Wall -O2 -o test_torsion test_torsion.o torsion_j.o \
         read_section_mesh_unv.o mesh_checker.o compute_section_properties.o
4. Verify the executable was created:

bash
ls -la test_torsion
Expected output:

text
-rwxr-xr-x 1 user user 35000 test_torsion
5. Run the test:

bash
./test_torsion
Usage
Basic Test (all sections)
bash
./test_torsion
Expected Output
text
=== Torsional Constant J (FEM) ===
  Nodes:       335
  Triangles:   472
  Free nodes:  141
  Boundary nodes: 194
  J = 422560.12 mm⁴
Validation
Test Results
SectionMeshExpected J (mm⁴)Our ResultErrorNotes
Circle Ø1012 nodes, 10 triangles981.75738.4-25%Coarse mesh
HEB200335 nodes, 472 triangles~1,200,000422,560-65%Mesh too coarse
Sources of Error
Coarse meshes - Too few interior nodes to capture quadratic solution

Linear approximation - Exact solution for circle is quadratic, linear elements underestimate

Mesh quality - Irregular triangles affect accuracy

How to Improve
Generate finer meshes (more triangles)

Use higher-order elements (quadratic triangles)

Implement adaptive mesh refinement

Refine mesh near boundaries where gradients are high

Analytical Approximations
For very coarse meshes with no interior nodes, the code falls back to analytical approximations:

ShapeFormula
CircleJ = π·r⁴/2
Rectangle (thin)J ≈ b·h³/3
SquareJ ≈ 0.1406·a⁴
Code Structure (torsion_j.f)
Main Subroutine
fortran
subroutine compute_torsion_J(nnode, ntri, conn, y, z, J_torsion, verbose)
Steps
find_boundary_nodes() - Identify boundary edges

map() - Build free node mapping

analytical_J_approx() - Fallback for coarse meshes

assemble() - Build K matrix and RHS vector

DPOSV() - Solve linear system

integrate() - Compute J = 2∫φ dA

Troubleshooting
ProblemSolution
"No triangles found"Check UNV file has section 2412 with type 41
DPOSV fails (info > 0)Matrix singular; check boundary conditions
J = 0No free nodes; refine mesh
Large errorMesh too coarse; use more triangles
Segmentation faultIncrease array dimensions in source
"make: command not found"Install build-essential: sudo apt install build-essential
"gfortran: command not found"Install gfortran: sudo apt install gfortran
Performance Notes
NodesFree NodesMatrix SizeMemorySolve Time
335141141×141~0.16 MB<0.01 sec
1000500500×500~2 MB~0.05 sec
500030003000×3000~72 MB~1 sec
Known Issues
Coarse mesh accuracy - For accurate J, need fine mesh (1000+ triangles)

Circle approximation - Linear elements cannot represent quadratic solution exactly

HEB200 error - 65% error with current mesh (needs refinement)

No mesh refinement - Adaptive refinement not yet implemented

References
Saint-Venant, B. (1855). "Mémoire sur la torsion des prismes"

Prandtl, L. (1903). "Zur Torsion von prismatischen Stäben"

Timoshenko, S.P., Goodier, J.N. (1970). Theory of Elasticity

Zienkiewicz, O.C., Taylor, R.L. (2000). The Finite Element Method

Version History
VersionDateChanges
0.12024Initial release: J calculation with DPOSV
0.22024Added analytical fallback for coarse meshes
0.32024Validation and documentation
Contributing
Contributions are welcome! Areas needing help:

Finer test meshes

Higher-order element implementation

Adaptive mesh refinement

Validation against analytical solutions

Please:

Fork the repository

Create a feature branch

Submit a pull request

Acknowledgements
Salome Platform for mesh generation

Code_Aster for reference UNV format

LAPACK team for numerical libraries

Happy engineering! 🚀

Last updated: 2024
