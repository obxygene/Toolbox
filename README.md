
# Toolbox project

This is a **personal collection** of MATLAB scripts and tools for quantum transport, condensed matter physics, and various numerical simulations. The repository is modular, with each folder focusing on a specific class of physical models or computational techniques.

Some codes are adapted or referenced from online sources; where possible, original URLs are cited in the code comments.


## Directory Overview

- **RecursiveGreenFunction/**
	- Tools for recursive Green's function calculations, including surface Green's functions and transmission.
- **Kernel_Polynomial_Method/**
	- Implementation of the Kernel Polynomial Method (KPM) for density of states, conductivity, and spectral functions.
	- Includes Chebyshev expansion, kernel corrections, and test examples for graphene and other lattices.
	- **Note:** Currently, only the DOS (density of states) calculation yields reliable results. The conductivity calculation does not show quantized behavior and should be used with caution.
- **Lattice/**
	- Scripts for constructing and analyzing lattice models in both momentum and real space.
	- Subfolders:
		- `Momentum_Space/`: Graphene, Haldane, SSH, and other models in k-space.
		- `Real_Space/`: Real-space Hamiltonians, hopping matrices, and visualization tools.
- **Numerical Simulations/**
	- General numerical experiments, including disorder, magnetism, quantum transport, and more.
	- Subfolders for Anderson localization, Drude weight, exact diagonalization, magnetism, quantum transport, superconductivity, topology, transfer matrix, and unclassified simulations.



## Features

- Efficient calculation of density of states using KPM.
- Quantum transport calculations via lattice Green's functions in real frequency.



## Requirements

- **MATLAB R2018b or later** is required. Most scripts are written for the Live Script environment (R2025a or newer recommended for best compatibility).
- Some scripts may contain extra formatting characters for figures.

> **Note:** Starting in R2025a, the Live Editor supports a new plain text Live Code file format (`.m`) for live scripts, as an alternative to the default binary format (`.mlx`). This format supports all Live Editor features and is compatible with version control.

For more information, see: [MATLAB Plain Text File Format for Live Scripts](https://ww2.mathworks.cn/help/matlab/matlab_prog/plain-text-file-format-for-live-scripts.html)