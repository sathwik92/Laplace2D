# Solving Laplace's Equation in two-dimensions (2D) using the Finite Element Method

An Example code to solve the scalar Laplace's problem in two-dimensions (2D) using finite element method (FEM) implemented within an MPI-based Parallel Computing Environment. Here we have applied Dirichlet boundary conditions with 3 sides of
the rectangle as zero and one side as a sine function. 

We perform discretization within unstructured 2D triangles. Mesh files are included in the mesh/ directory. 

This program can handle both linear and quintic Hermite polynomials as the finite element basis functions. 

The Gauss-Dunavant quadrature method is used to compute the numerical integration. 

The FEM matrix can be solved here using several parallel solvers listed below. 

## Available Sparse Linear Solvers				
Richardson,	Chebyshev,	Conjugate Gradient,
BiConjugate Gradient,	Generalized Minimal,	Residual,	
Flexible Generalized,	Minimal Residual,	Deflated Generalized,	
Minimal Residual,	Generalized Conjugate,	Residual,	
Squared,	Transpose-Free,	Quasi-Minimal Residual,	Conjugate Residual,	
Least Squares Method	

# Getting Started:

## Prerequisites:

Petsc
Petsc: Quick Installation Instructions in Ubuntu

git clone -b release https://gitlab.com/petsc/petsc.git petsc
git pull

cd petsc/

./configure --with-cc=mpicc --with-cxx=mpicxx --download-mpich --download-fblaslapack
make PETSC_DIR=/mnt/c/Users/sathw/Documents/petsc PETSC_ARCH=arch-linux-c-debug all
make PETSC_DIR=/mnt/c/Users/sathw/Documents/petsc PETSC_ARCH=arch-linux-c-debug check
make all check

## Installation:

cd src

make

## Example Input File:
FEMstruct2d.inp

This file can be found in the input directory.

## Output:
Output files are written in the output directory. 
Solutions are written in:
output.out

FEM matrix parameters are written in:
input_params.out

## Example:
./Laplace2d FEMstruct2d.inp

# Author
Sathwik Bharadwaj

## License
This project is licensed under the MIT License.

## References
Balay, S. et al. (2015) PETSc Web page. Available at: http://www.mcs.anl.gov/petsc. 
