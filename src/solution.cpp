//=================================================================
//  File solution.cpp
//    Created:  March 10, 2017: Sathwik bharadwaj
//    Modified: 
//    Last change: 
//  
//=================================================================
//  
//  
// ****************************************************************
// This file interpolates the ouput vector and writes the solution on output.out
// ****************************************************************
// Function in this file:
// 	 PetscErrorCode solution(global_matrices & gmat, data & dat, Vec psi);
//  
// =================================================================
//header file with all external libraries and user defined classes.
#include "global_params.h"

PetscErrorCode solution(global_matrices & gmat ,data & dat, Vec psi)                 
{

	double phi[dat.ndeg];
	double B[dat.ndeg][18];
	generateShapeCoeffMatrix(B, dat.ndeg);
	double T[6][6];
	double xc[3], yc[3];



	//========================================================================
	// Interpolate the nodal solutions obtained to reconstruct the actual
	//   wavefunctions.
	//========================================================================

	double zx, zy, xi, eta, det;
	int izx, izy, iel, i;
	int n1, n2, n3; //used for the three nodes index in next loop
	char buffer[200]; //used for output file path
	PetscScalar fz;
	PetscErrorCode ierr;
	int ind[dat.nele_mat], sind[dat.nele_mat];
	PetscScalar elemvec[dat.nele_mat];
	double rfz, imfz, afz;
	char outputpath[50];
	strncpy(outputpath, dat.output_files_path.c_str(), sizeof(outputpath));
	outputpath[sizeof(outputpath) - 1] = 0;
	sprintf(buffer,"%soutput.out", outputpath);
	std::ofstream fout(buffer);
	fout<<"#	x	y	value.real	value.imag	value.abs"<<std::endl;

	fout.setf(std::ios::scientific);
	fout.setf(std::ios::showpos);

	// we have the number of interpolation point along each axis
	// i.e ndzx and ndzy!
	for(izy=0;izy<dat.ndzy+1;izy++)
	{
		for(izx=0;izx<dat.ndzx+1;izx++)
		{
			zx = dat.node[0][0] + double(izx)*dat.dzx; //dat.node[0][0] corresponds to xmin
			zy = dat.node[0][1] + double(izy)*dat.dzy; //dat.node[0][1] corresponds to ymin

			locelem(zx, zy, dat, iel); //locate the element in which the interpolation points is
			n1 = dat.elem[iel][0]; //this corresponds to the first node
			n2 = dat.elem[iel][1]; //2nd node
			n3 = dat.elem[iel][2]; //3rd node
			xc[0] = dat.node[n1][0]; 
			yc[0] = dat.node[n1][1];
			xc[1] = dat.node[n2][0];
			yc[1] = dat.node[n2][1]; 
			xc[2] = dat.node[n3][0];
			yc[2] = dat.node[n3][1];
			generateTransformationMatrix(xc, yc, T);
			det = (xc[1]-xc[0])*(yc[2]-yc[0])-(xc[2]-xc[0])*(yc[1]-yc[0]); //determinant used in inversion

			xi = ((yc[2]-yc[0])*(zx-xc[0])+(xc[0]-xc[2])*(zy-yc[0]))/det;
			eta =((yc[0]-yc[1])*(zx-xc[0])+(xc[1]-xc[0])*(zy-yc[0]))/det;

			//we have the matrix equation for x and y in terms of eta and xi
			//the determinant is used because we take the inverse to solve for eta and xi

			shape(phi, xi, eta, T, B, dat.ndeg); //creating shape function for specific interpolation point

			rfz = 0.0;
			imfz= 0.0;
			fz  = 0.0;
			i =0;
			for(int ie=0;ie<dat.node_elem; ie++) //first loop over number of nodes in 1 element
			{
				for(int id=0;id<dat.ndof;id++) //loop over degrees of freedom - necessary for Hermite
				{
					ind[i] =(dat.elem[iel][ie]*dat.ndof+id); 
					sind[i] =ie*dat.ndof+id;
					i = i+1;
				}
			}
			ierr = VecGetValues(psi, dat.nele_mat , ind , elemvec);CHKERRQ(ierr);

			for(int ie=0;ie<dat.nele_mat; ie++){ //first loop over number of nodes in 1 element
				fz += phi[sind[ie]]*elemvec[ie]; 
			} //end of ie loop
			rfz = PetscRealPart(fz);
			imfz = 0;
			afz = sqrt(rfz*rfz+ imfz*imfz);

			fout<<zx<<"\t"<<zy<<"\t"<<rfz<<"\t"<<imfz<<"\t"<<afz<<std::endl; //printing out solutions
		} //end of izx loop
		fout<<std::endl;
	} //end of izy loop
	fout.close();

	return ierr;
}



