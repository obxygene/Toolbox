# RecursiveGreenFunction Module

This module provides MATLAB scripts and tools for recursive Green's function calculations, which are widely used in quantum transport and condensed matter physics. The recursive Green's function method is a powerful technique for computing electronic properties, transmission, and surface Green's functions in low-dimensional systems.

## Contents

- `FiniteTemperatureTransmission.m`  
  Calculate transmission at finite temperature using Green's functions. (Not complete now)
- `GreenFunction_Transmission.m`  
  Compute transmission coefficients via Green's function formalism.
- `RecursiveGreenFunction_1L.m`  
  **Main program**. Main part of Recursive Green's function.
- `SurfaceGreenFunction.m`
  Calculate the surface green function for a given onsite and hopping hamiltonian(*This version uses iteration number as convergence criterion*)
- `SurfaceGreenFunction_Broadening.m` 
  Calculate the Broadening of Surface Green's function
- `SurfaceGreenFunction_SelfEnergy.m`
  Calculate the Broadening of Surface Green's function
- `SurfaceGreenFunction_V2.m`  
  Same as `SurfaceGreenFunction.m`. (*This version uses relative error as convergence criterion*)
- `tridiag_block_inv.m`  
  Inversion of block tridiagonal matrices, useful for recursive algorithms. For general cases, this will be slower than `RecursiveGreenFunction_1L.m` 
- `Example/`  
  Example scripts demonstrating usage of the above tools.

## Features

- Efficient calculation of Green's functions for large or quasi-1D systems.
- Transmission and conductance calculations for quantum transport.
- Surface Green's function and self-energy evaluation for open systems.

## References
- Katharina Laubscher, Clara S. Weber, Maximilian Hünenberger, Herbert Schoeller, Dante M.
Kennes, Daniel Loss, and Jelena Klinovaja. RKKY interaction in one-dimensional flat-band lat
tices. Physical Review B, 108(15):155429, October 2023.

- Wei Ren, Cai-Zhuang Wang, Kai-Ming Ho, and C. T. Chan. Transport in a metallic nanotube at
finite temperature. Physical Review B, 79(16):161404, April 2009.

- Juntao Song, Haiwen Liu, Jie Liu, Yu-Xian Li, Robert Joynt, Qing-feng Sun, and X. C. Xie. Quan-
tum interference in topological insulator Josephson junctions. Physical Review B, 93(19):195302, May 2016

...To be added