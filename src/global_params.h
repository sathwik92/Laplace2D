//================================================================
//    File global.h
//    Created: Feb 28, 2017: Sathwik Bharadwaj 
//    Define global physical constants and variables
//================================================================
//
//
// ***************************************************************
// This file defines the global structure. 
// It also defines data and global_matrices data structures. 
// It includes miscellaneous other files.
// ***************************************************************
//
//
/* if not already defined, define the variable GLOBAL_H*/
#ifndef GLOBAL_H
#define GLOBAL_H

/*#undef max
#undef min*/
#include "petsc.h"
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <vector>
#include "constants.h"
#include <cstring>   // Needed for strcpy
#include <sstream>   // Needed for string conversions

//=====================================================
//*****************************************************
//=====================================================

//==============================================================
//  The class "data" contains input data. 
//  Usually in the codes, the class is labeled 
//  "dat" or "indat"
//==============================================================
//==============================================================
class data //class efining the structure data
{ 
	public:

		int ndebug;
		//==============================================  

		//==============================================
		// FEM related and geometry related parameters:
		// =============================================

		//This first set has to do with the parameters having to do with "size"
		int nglobal;   //=ngnodes*ndof*nband
		int ngnodes;   // global number of nodes
		int nelem;     // total number of elements
		int nbelem;
		int nbnode;
		int node_elem; // nodes per element
		int ndof;      // number of deg of freedom at each node
		int nele_mat;  // size of element dmat =node_elem*ndof
		int ndeg;      // Degree of interp. polynomial=(node_elem*ndof -1)
		int nband;
		int nboundary; //number of boundarires. 
		double xmin, xmax, ymin, ymax; //minimum and maximum values along x and y axis


		//==============================================
		//These are the arrays used to store the mesh data
		//The sizes of the arrays are dependent on some of the above variables
		//==============================================
		double **node;
		int **elem;
		int *material; //material of each element
		int **belem;
		double **bnode;
		int **bmaterial;
		double **bvalue; //values assigned to the nodes in mesh generator
		int *btob; //type of boundary condition
		int *bnum; // bounndary number
		
		//=====================================================
		// Gauss Quadrature parameters
		//=====================================================
		int ngaus; 
		double *xigaus; 
		double *etagaus;
		double *wgaus;
		//=====================================================


		//=====================================================
		// Parameters for output
		//=====================================================

		int ndzy; //number of output points along ZY axis
		int ndzx; //number of output points along ZX axis
		double dzx; //distance between interpolation points along ZX axis
		double dzy; //distance between interpolation points along Zy axis

		//=====================================================
		//We create a subdirectory into which we 
		//put all the output files:
		//=====================================================

		std::string output_files_path;

		//=====================================================
		// path to the output files folder:
		//=====================================================

		FILE* log1; // for reading in input band parameters

		PetscMPIInt rank, size;//petsc parameters

		std::string 	  solvertype; //solver type used in the calculation
		double  rtol, atol, divtol;//right, absolute and relative tolerence
		int     maxit;//maximum number of allowed iterations

};

//=====================================================
class global_matrices //defining the global_matrices struct
{
	public:
		Mat A; //Left Hand Matrix used to store Laplace data
		Vec rhs; //Right hand side vector        
};

//=====================================================
//*****************************************************
//=====================================================

// Note: the prototypes use the definitions of the various 
// structures made in the above: So we have to include the 
// prototypes below HERE ***:

#include "prototypes.h"

#endif


//=====================================================
//*****************************************************
//=====================================================

