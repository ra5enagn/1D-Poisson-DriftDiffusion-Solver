# 1D Poisson–Drift-Diffusion Solver (PN Diode) – MATLAB

## Overview
This repository contains a MATLAB implementation of a one-dimensional, self-consistent semiconductor device solver for a PN diode under applied bias. The code solves the electrostatic potential using Poisson’s equation and evaluates carrier concentrations and transport quantities consistent with a drift–diffusion formulation. A bias sweep is supported to study non-equilibrium behavior.

## Key Capabilities
- 1D Poisson solution for electrostatic potential and electric field
- Electron and hole carrier concentration evaluation from the potential
- Bias sweep support (forward bias range configurable in the script)
- Extraction of transport-related quantities such as current density and quasi-Fermi level trends (as implemented in the script)
- Uses a lightweight LU-based linear solver for the Poisson/linearized system

## Files
- `DD_solver.m`  
  Main script that defines constants, device/doping setup, mesh, bias sweep, solves the coupled equations, and generates outputs/plots.
- `LU_Decomposition.m`  
  Helper routine used as the linear solver (commonly used for tridiagonal or banded systems, depending on implementation).

## Requirements
- MATLAB (base installation)
- No additional MATLAB toolboxes are required based on the current implementation

## How to Run
1. Clone the repository or download the files.
2. Open MATLAB and set the working directory to the repository folder.
3. Run:
   ```matlab
   DD_solver
