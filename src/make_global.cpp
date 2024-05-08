//=================================================================
//  File: make_global.cpp
//  Created:  March 17, 2017: Sathwik Bharadwaj
//  Modified: 
//  Last change: 
//
//=================================================================
/* This file is used to defined the element matrices and then construct the
 * global matrices from it
 * the element matrices are defined according to the action principle and 
 * the Lagrangian */
//Function in this file:
//	PetscErrorCode make_global(global_matrices& gmat, data& dat);
//=================================================================

//header file with all external libraries and user defined classes.
#include "global_params.h"
#include "global_params.h"

PetscErrorCode make_global(global_matrices& gmat, data& dat)
{
	/* in the case of the Laplace equation we end up with 
	 * a matrix in the LHS and a vector on the RHS */
	//Declare local variables
	//========================

	int iel, ie, je, igaus; //declaring indices used for loops
	double phipx[dat.ndeg], phipy[dat.ndeg]; //Shape functions
	double B[dat.ndeg][18]; //Coefficient matrix for the shape functions
	generateShapeCoeffMatrix(B, dat.ndeg);
	double T[6][6];
	double xc[3];
	double yc[3];
	memset(xc, 0.0, sizeof(xc));//clea
	memset(yc, 0.0, sizeof(yc));
	// wpp and wff are temp variables used to store
	// the gauss weight*xjacobian multiplied by phi_primes
	// or by phi_s

	double wpp; //in Laplace there is no use for wff
	PetscErrorCode ierr; //Petsc error code
	int start_ind, end_ind, n_loc_elem;//start and end index of the elementsin

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &dat.rank); CHKERRQ(ierr);//rank of the current process/
	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &dat.size); CHKERRQ(ierr);//rank of the current process/
									 //Creating the matrix with appropriate dimensions
	ierr = MatCreate(PETSC_COMM_WORLD,&gmat.A);CHKERRQ(ierr);
	ierr = MatSetSizes(gmat.A,PETSC_DECIDE,PETSC_DECIDE,dat.nglobal,dat.nglobal);CHKERRQ(ierr);
	ierr = MatSetFromOptions(gmat.A);CHKERRQ(ierr);
	ierr = MatSetType(gmat.A, MATMPIAIJ); CHKERRQ(ierr);
	ierr = MatSetUp(gmat.A); CHKERRQ(ierr);
	//rhs vector is created with same matrix parallel layout
	ierr = MatCreateVecs(gmat.A,NULL,&gmat.rhs);CHKERRQ(ierr);
	ierr = VecSetFromOptions(gmat.rhs);CHKERRQ(ierr);

	if(dat.rank < dat.nelem%dat.size)
	{
		start_ind= (dat.nelem/dat.size+ 1)*dat.rank;
		end_ind = (dat.nelem/dat.size+ 1)*(dat.rank+1)-1;
	}
	else
	{
		start_ind= dat.nelem/dat.size*dat.rank +dat.nelem%dat.size;
		end_ind=  dat.nelem/dat.size*(dat.rank+1)+ dat.nelem%dat.size-1;
	}
	n_loc_elem = end_ind-start_ind+1;


	//==========================
	//Begin loop over elements;
	//==========================
	for(iel=start_ind; iel<end_ind+1; iel++)
	{

		//start of element calculations
		//==============================

		int n1 = dat.elem[iel][0]; //these gives us the nodes of the element
		int n2 = dat.elem[iel][1]; //the indexing is of course in a counterclockwise way
		int n3 = dat.elem[iel][2];
		xc[0] = dat.node[n1][0]; //xc[0] - x coordinate of first node
		yc[0] = dat.node[n1][1]; //yc[0] - y coordinate of first node
		xc[1] = dat.node[n2][0]; //xc[1] - x coordinate of 2nd node
		yc[1] = dat.node[n2][1];
		xc[2] = dat.node[n3][0];
		yc[2] = dat.node[n3][1];
		double jacobian = 2*abs(xc[0]*(yc[1]-yc[2])+xc[1]*(yc[2]-yc[0])+xc[2]*(yc[0]-yc[1]));
		//the jacobian here corresponds to twice the area of the triangle - it is
		//easier if you just take the area of rectangle but...
		generateTransformationMatrix(xc, yc, T);
		//start of Gauss loop
		//=====================

		for(igaus=0;igaus<dat.ngaus;igaus++)
		{
			double xi = dat.xigaus[igaus];
			double eta = dat.etagaus[igaus];
			double wt = dat.wgaus[igaus];

			//Caluculate the shape and deriv functions
			//========================================

			deriv1(phipx, phipy, xi , eta, T, B, dat.ndeg);
			for(int a=0;a<dat.node_elem;a++) //loop over each entry along element matrix row
			{
				for(int b=0;b<dat.node_elem;b++) //loop over each entry along element matrix column
				{
					for(int c=0;c<dat.ndof;c++) //loop over degrees of freedom along row - necessary for Hermite
					{
						for(int d=0;d<dat.ndof;d++) //loop over deg of freedom along column - necessary for Hermite
						{
							int f_i = dat.elem[iel][a]; //f_i corresponds to the "row node number"

							int f_j = dat.elem[iel][b]; //f_j corresponds to the "column node number"

							ie = a*dat.ndof+c;
							je = b*dat.ndof+d;
							wpp = 0.5*wt*jacobian*(phipx[ie]*phipx[je]+phipy[ie]*phipy[je]); 
							//	cout << "x derivative value: " << phipx[je] << "y derivative: "<<phipy[je] <<endl;
							//      wpp = 0.5*wt*jacobian*(phipx[ie]*phipx[je]+phipy[ie]*phipy[je]); 
							// the jacobian here arises from the integration which is done in the
							// "local" element and therefore requires a jacobian(scale) factor!!
							// wff = wt*jacobian*(phi[ie]*phi[je]+phiy[ie]*phiy[je]);
							ierr = MatSetValue(gmat.A,(f_i*dat.ndof)+c,(f_j*dat.ndof)+d,wpp,ADD_VALUES); CHKERRQ(ierr);

							//    Aelement[ie][je]+=wpp;
						} //end of je loop - loop over columns
					} //end of ie loop - loop over rows
				}
			}
		} //end of igaus loop - loop over gauss points
	}//end of iel loop
	ierr = MatAssemblyBegin(gmat.A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(gmat.A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	return ierr;
}// End of make_global subroutine

//==========================================================
//==========================================================

