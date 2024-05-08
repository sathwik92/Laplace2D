//=================================================================
//  File petsc_diagonalizer.cpp
//  Created:  March 8, 2017: Sathwik Bharadwaj
//  Modified: 
//=================================================================
//
//
// ****************************************************************
// This files diagonalizes the matrix, output the eigen values 
// and calls solution.
// condition number is calculated 
// ****************************************************************
//
//
//=================================================================
#include <petscksp.h>
#include "global_params.h"
PetscErrorCode petsc_diagonalizer(global_matrices & gmat , data & dat)
{
	Vec            x,tx;  /* approx solution, exact solution */
	KSP            ksp;     /* linear solver context */
	PetscInt       its; /* Number of iternations */
	PetscErrorCode ierr; /* Petsc Error Code */


	//Create the solution vector x

	ierr = MatCreateVecs(gmat.A,NULL,&x);CHKERRQ(ierr);
	ierr = VecSetFromOptions(x);CHKERRQ(ierr);


	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Create the linear solver and set various options
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	/*
	   Create linear solver context
	   */
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

	/*
	   Set operators. Here the matrix that defines the linear system
	   also serves as the preconditioning matrix.
	   */
	ierr = KSPSetOperators(ksp,gmat.A,gmat.A);CHKERRQ(ierr);

	/*
	   Set linear solver defaults for this problem (optional).
	   - By extracting the KSP and PC contexts from the KSP context,
	   we can then directly call any KSP and PC routines to set
	   various options.
	   - The following two statements are optional; all of these
	   parameters could alternatively be specified at runtime via
	   KSPSetFromOptions().  All of these defaults can be
	   overridden at runtime, as indicated below.
	   */
	ierr = KSPSetTolerances(ksp, PETSC_DEFAULT,  PETSC_DEFAULT, PETSC_DEFAULT , PETSC_DEFAULT);CHKERRQ(ierr);

	/*
	   KSP Solver type
	   */ 

	char solvetype[50];
	strcpy(solvetype, dat.solvertype.c_str());
	solvetype[sizeof(solvetype) - 1] = 0;
	ierr = KSPSetType(ksp, solvetype);

	/*
	   Set runtime options, e.g.,
	   -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
	   These options will override those specified above as long as
	   KSPSetFromOptions() is called _after_ any other customization
	   routines.
	   */
	ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);


	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Solve the linear system
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	ierr = KSPSolve(ksp, gmat.rhs, x);CHKERRQ(ierr);

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Check solution and clean up
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	/*
	   Check the error
	   */
	ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);

	/*
	   Print convergence information.  PetscPrintf() produces a single
	   print statement from all processes that share a communicator.
	   An alternative is PetscFPrintf(), which prints to a file.
	   */
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Size of the global matrix and vector:\t%d\n", dat.nglobal);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Total number of iterations taken by the linear solver:\t%d\n",its);CHKERRQ(ierr);
	ierr = KSPGetSolution(ksp, &x);
	/*
	   Copy the parallel solution vector to a new vector in the zeroth process 
	   */	
	VecScatter     ctx;
	ierr = VecScatterCreateToAll(x,&ctx, &tx);CHKERRQ(ierr);     
	ierr = VecScatterBegin( ctx, x, tx, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
	ierr = VecScatterEnd( ctx, x, tx, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
	ierr = VecScatterDestroy(&ctx);
	//Interpolation of the solution
	solution(gmat, dat, tx); 
	/*
	   Free work space.  All PETSc objects should be destroyed when they
	   are no longer needed.
	   */
	ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
	ierr = MatDestroy(&gmat.A);CHKERRQ(ierr);
	ierr = VecDestroy(&gmat.rhs);CHKERRQ(ierr);
	ierr = VecDestroy(&x);CHKERRQ(ierr);
	return ierr;
}
