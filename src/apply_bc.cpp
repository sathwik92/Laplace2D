//=================================================================
//  File apply_bc.h
//  Created:  March 17, 2017: Sathwik Bharadwaj
//  Modified: 
//  Last change: 
//
//=================================================================
// In this function we apply the boundary conditions.
/* We know the boundary nodes coordinates from reading the mesh output file
 * bnode.dat
 * We then first use the locelem function to locate in which element the node
 * is and then we compare the 3 nodes of the element with our bnodes
 * coordinates. This allows us to obtain the node number of the "bnode"
 * In the bnode.dat we have the function value, derivative values...etc given
 * and therefore each of these values can be assigned a global index which can
 * be related to the node number [(node_number*ndof*nband)+dof)]. 
 * Given the global index we then proceed to the "half-benediction process"
 */
//Function in this file:
//	PetscErrorCode apply_bc(global_matrices & gmat , data & dat);
//=================================================================
#include "global_params.h"

PetscErrorCode apply_bc(global_matrices& gmat, data& dat)
{
	int i, j, iel, ig;
	int n1, n2, n3, in;
	PetscErrorCode ierr; //Petsc error code
	int id[1];
	//Here we obtain the node number of the boundary node
	for(i=0;i<dat.nbnode;i++){
		locelem(dat.bnode[i][0], dat.bnode[i][1], dat, iel); //locating which element bnode is found in
		n1 = dat.elem[iel][0];
		n2 = dat.elem[iel][1]; 
		n3 = dat.elem[iel][2];
		if((dat.node[n1][0]==dat.bnode[i][0]) && (dat.node[n1][1]==dat.bnode[i][1]))
			in = n1;
		else{ 
			if((dat.node[n2][0]==dat.bnode[i][0]) && (dat.node[n2][1]==dat.bnode[i][1]))
				in = n2;
			else
				in = n3;
		} //conditional operators are used to determining the correct node number by explicitly comparing
		  //the coordinates 

		  // The "half-benediction" process consist of assigning the value of the
		  // global node on the rhs vector and deleting the appropriate row on the
		  // lhs(matrix) and assigning a value of 1 to the diagonal entry

		  // However because of the additional degrees of freedom used in Hermite
		  // interpolation, we seperate our "half-benediction" process into two
		  // possible case

		for(j=0;j<dat.ndof*dat.nband;j++){ //in the case of linear it still works because ndot would be 1 and j=0(no loop)
			if(dat.bvalue[i][j]!=99)
			{
				ig = (in*dat.ndof*dat.nband)+j; 
				id[0] = ig;
				ierr = MatZeroRows(gmat.A, 1, id, 1.0, NULL,  NULL);CHKERRQ(ierr);
				ierr = VecSetValue(gmat.rhs, ig , dat.bvalue[i][j], INSERT_VALUES); 	   	   
			}
		} //end of j loop
	} //end of i loop
	ierr = VecAssemblyBegin(gmat.rhs);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(gmat.rhs);CHKERRQ(ierr);
	return ierr;

	//===========================================
} //end apply_bc subroutine

//===========================================
//===========================================
