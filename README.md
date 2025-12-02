2D Finite Element Analysis (Plane Stress)

A MATLAB implementation of a 2D Linear Static Finite Element solver using 4-node Isoparametric Quadrilateral (Q4) elements. This project was developed for the ME 5310 Finite Element Method course.

Features

Element Type: 4-Node Isoparametric Quadrilateral (Q4)

Physics: Plane Stress Elasticity

Integration: 2x2 Gaussian Quadrature

Post-Processing: Stress calculation at Gauss Points and Von Mises Stress evaluation

Visualization: Deformed mesh superposition

Validation

The code has been verified against two test cases:

1. Single Element Patch Test

Verifies constant stress states and Poisson locking effects under clamped boundary conditions.

2. Plate with Central Hole ($K_t$ Analysis)

Models a quarter-symmetry plate with a circular hole ($d/w = 0.1$).

Theoretical Max Stress: 3022 MPa ($K_t \approx 2.72$)

FEA Computed Stress: 2952 MPa

Error: ~2.3%

Usage

Clone the repository.

Open MATLAB and navigate to the examples/ folder.

Run Run_TestCase2_PlateHole.m.

Mathematical Formulation

The stiffness matrix is derived using:


$$k^e = \int_{-1}^{1} \int_{-1}^{1} B^T D B \cdot t \cdot \det(J) \, d\xi d\eta$$

Acknowledgements

Author: [Your Name]

Course: ME 5310, Fall 2025

Assisted by Gemini for code debugging and documentation structure.
