# P2-P1 Finite Element Solver for Stokes Equations

A rigorous MATLAB implementation of a P2-P1 (Taylor-Hood) Finite Element Method (FEM) solver designed to simulate incompressible, steady, viscous fluid flows at low Reynolds numbers. The solver handles the full pipeline from mesh generation via `gmsh` to the assembly of global matrices and numerical validation.

## 📄 Project Documentation
* **[Read the Full Scientific Report (PDF)](stokes_report.pdf)** *(Note: The report and code comments are written in French, but the variable naming and matrix structures follow standard FEM mathematical conventions).*

## Mathematical & Technical Framework
* **Governing Equations:** Solves the simplified Navier-Stokes equations neglecting the convection term: `-νΔu + ∇p = f` and `div(u) = 0`.
* **Mixed Finite Elements:** Uses stable P2-P1 Taylor-Hood elements (quadratic interpolation for velocity, linear for pressure) to ensure numerical stability (LBB inf-sup condition).
* **Numerical Core:** Custom numerical integration using a 6-point Gauss-Legendre quadrature (order 4) for the mass matrix and a 3-point quadrature (order 2) for the stiffness matrix on reference triangles.
* **Boundary Conditions:** Implementation of Dirichlet and Neumann boundary conditions using a pseudo-elimination method to modify the global algebraic system without destroying matrix symmetry.

## Validation & Test Cases
* **Analytical Validation:** Convergence analysis using exact analytical solutions (e.g., `u(x,y) = 3*cos(π*x)*cos(2*π*y)`) to compute L² and H¹ norms of the error.
* **Complex Geometry (Forward Facing Step):** Simulation of flow over a geometric step with an imposed parabolic Poiseuille inlet velocity profile. Includes a parametric study on the influence of dynamic viscosity (ν) on the pressure field.

## Requirements
* MATLAB (R2022a or later recommended)
* GMSH API or standalone installation
