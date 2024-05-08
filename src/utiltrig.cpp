//=================================================================
//  File utiltrig.h
//  Created:  March 17, 2017: Sathwik Bharadwaj
//  Modified: 
//  Last change: 
//
//=================================================================
// ****************************************************************
// This function contains the finite element shape functions, and their derivatives. 
// ****************************************************************
// Functions within this file:
//	void shape(double phi[], double xi_val, double eta_val, double Transform[][6], double shapeCoeff[][18], int& ndeg1);
//	void deriv1(double derivx[], double derivy[], double xi_val, double eta_val, double Transform[][6], double shapeCoeff[][18], int& ndeg1);
//	void deriv2(double derivxx[], double derivxy[], double derivyy[], double xi_val, double eta_val, double Transform[][6], double shapeCoeff[][18]);
//	void generateTransformationMatrix(double xp[3], double yp[3], double T[][6]);
//	void generateShapeCoeffMatrix(double B[][18], int& ndeg1);
//	void locelem(double x, double y, data &dat, int &iel);
//	void trigauss(int& ngaus, double *x, double *y, double *wt);
//
//=========================================================

#include <stdio.h>
#include "global_params.h"

//==========================================================
// Locate element
//==========================================================

//structured version
/* In the case of the structured grid
 * we know already what is the x_min, x_max and y_min, y_max of
 * each triangular element
 * It will correspond each time to the x-coordinate and y-coordinate
 * of the same nodes everytime
 * even though the indexing is different for "even" and "odd" triangular
 * elements
 * for any element x_min = node[0][0] i.e the x coordinate of the first node
 * y_min = node[0][1] i.e the y coordinate of the first node
 * x_max = node[1][0] i.e the x_coordinate of the 2nd node
 * y_max = node[2][1] i.e the y_coordinate of the 3rd node */

/* Now we know that we have a total number of elements nelem
 * loop over all elements so that x_min<x<x_max and y_min<y<y_max
 * however normally there would be two triangles satisfying the conditions
 * so we only deal with one type of triangle
 * from the element structure, for each rectangular element, the 2 triangles
 * would satisfy the conditions -> consider only one type of triangle
 * either "odd" or "even" triangles and then later determine which triangle
 * is the point in */
//the +1 is added only because of loop and < sign
//for unstructured nelem could be odd but still because i is int type nelem/2
//would return the smallest integer eg n<nelem/2<m -> nelem/2=n
//                                          n3
//                                          /|
//                                         / |
//                                        /__|
//                                       n1  n2
//the above triangle shows the type of triangle we are limiting ourselves to  

//the reason for dividing into two in loop is because we want to 
//consider only the "even" triangular elements

//now we know that the point is either in the i*2 or (i*2)+1 element
//to check in which triangle the point is in, we use the equation of the
//line seperating the two triangle
//this has the same form for all the triangle because of the structure of
//the grid
//equation of line: y-y_min = ((y_max-y_min)/(x_max-x_min))*(x-x_min)
/* points to consider - what if the point is on one of the edges??
 * should be solved now - hopefully 
 * need to look into it again
 * worry about rounding errors, where to put the fabs(float... thing?? */

//=================================================================================

void locelem(double x, double y , data &dat, int &iel)
{
	/* need to do a version for structured and another for unstructured */

	//structured version
	/* In the case of the structured grid
	 * we know already what is the x_min, x_max and y_min, y_max of
	 * each triangular element
	 * It will correspond each time to the x-coordinate and y-coordinate
	 * of the same nodes everytime
	 * even though the indexing is different for "even" and "odd" triangular
	 * elements
	 * for any element x_min = node[0][0] i.e the x coordinate of the first node
	 * y_min = node[0][1] i.e the y coordinate of the first node
	 * x_max = node[1][0] i.e the x_coordinate of the 2nd node
	 * y_max = node[2][1] i.e the y_coordinate of the 3rd node */

	/* Now we know that we have a total number of elements nelem
	 * loop over all elements so that x_min<x<x_max and y_min<y<y_max
	 * however normally there would be two triangles satisfying the conditions
	 * so we only deal with one type of triangle
	 * from the element structure, for each rectangular element, the 2 triangles
	 * would satisfy the conditions -> consider only one type of triangle
	 * either "odd" or "even" triangles and then later determine which triangle
	 * is the point in */


	int i, n1, n2, n3;
	/* 

	//the +1 is added only because of loop and < sign
	//for unstructured nelem could be odd but still because i is int type nelem/2
	//would return the smallest integer eg n<nelem/2<m -> nelem/2=n
	//                                          n3
	//                                          /|
	//                                         / |
	//                                        /__|
	//                                       n1  n2
	//the above triangle shows the type of triangle we are limiting ourselves to  

*/

	// Unstructured algorithm
	double a1, a2, a3, a4;
	double tol = 1e-10;
	for(i = 0; i < dat.nelem; i++)
	{
		n1 = dat.elem[i][0];
		n2 = dat.elem[i][1];
		n3 = dat.elem[i][2];

		a1 = dat.node[n1][0] * (dat.node[n2][1] - dat.node[n3][1])
			+ dat.node[n2][0] * (dat.node[n3][1] - dat.node[n1][1])
			+ dat.node[n3][0] * (dat.node[n1][1] - dat.node[n2][1]);
		a2 = x * (dat.node[n2][1] - dat.node[n3][1])
			+ dat.node[n2][0] * (dat.node[n3][1] - y)
			+ dat.node[n3][0] * (y - dat.node[n2][1]);
		a3 = dat.node[n1][0] * (y - dat.node[n3][1])
			+ x * (dat.node[n3][1] - dat.node[n1][1])
			+ dat.node[n3][0] * (dat.node[n1][1] - y);
		a4 = dat.node[n1][0] * (dat.node[n2][1] - y)
			+ dat.node[n2][0] * (y - dat.node[n1][1])
			+ x * (dat.node[n1][1] - dat.node[n2][1]);

		if((a1 > -tol) && (a2 > -tol) && (a3 > -tol) && (a4 > -tol))
		{
			iel = i;
			return;
		}
		if((a1 < tol) && (a2 < tol) && (a3 < tol) && (a4 < tol))
		{
			iel = i;
			return;
		}
	}
	std::cout << "PROBLEM: " << x << " " << y << " coordinates are out of range." << ".\n";
} //end of locelem

//=======================================================================================
//the reason for dividing into two in loop is because we want to 
//consider only the "even" triangular elements

//now we know that the point is either in the i*2 or (i*2)+1 element
//to check in which triangle the point is in, we use the equation of the
//line seperating the two triangle
//this has the same form for all the triangle because of the structure of
//the grid
//equation of line: y-y_min = ((y_max-y_min)/(x_max-x_min))*(x-x_min)
/* points to consider - what if the point is on one of the edges??
 * should be solved now - hopefully 
 * need to look into it again
 * worry about rounding errors, where to put the fabs(float... thing?? */
//=======================================================================================

//=======================================================================================
//shape functions and derivative shape functions
//======================================================================================

void shape(double phi[], double xi_val, double eta_val, double Transform[][6], double shapeCoeff[][18], int& ndeg1)
{
	double x = xi_val;
	double y = eta_val;
	double point[ndeg1]; //18 corresponding to the largest possible size of point
	int i,j;
	switch(ndeg1){
		case 3:
			//point corresponds to the basis functions of the linear shape functions
			point[0] = 1;
			point[1] = x;
			point[2] = y;

			//the linear shape functions are then given in terms of the basis functions
			//where the coefficient associated to each basis function for every shape
			//function is obtained from the shapeCoeff
			for(i=0;i< 3; i++){
				phi[i] = 0.0;
				for(j=0; j<3;j++){
					phi[i] += point[j]*shapeCoeff[j][i];
				}
			}
			break;

		case 18:
			double x = xi_val;
			double y = eta_val;
			//works the same way as for the linear case
			point[0]  = 1;
			point[1]  = x;
			point[2]  = y;
			point[3]  = x*x;
			point[4]  = x*y;
			point[5]  = y*y;
			point[6]  = x*x*x;
			point[7]  = x*x*y;
			point[8]  = y*y*x;
			point[9]  = y*y*y;
			point[10] = x*x*x*x;
			point[11] = x*x*x*y;
			point[12] = x*x*y*y;
			point[13] = x*y*y*y;
			point[14] = y*y*y*y;
			point[15] = x*x*x*x*x-5*x*x*x*y*y;
			point[16] = x*x*y*y*y-x*x*x*y*y;
			point[17] = y*y*y*y*y-5*x*x*x*y*y;

			double phi2[18];
			//This is in essence the multiplication of each basis function to its
			//corresponding coefficient in the shape functions
			//phi2 are the "unscaled" basis functions - pg 101 Dhatt and Touzot
			for(int i=0; i<18 ;i++){
				phi2[i]=0.0;
				for(int j=0; j<18 ; j++){
					phi2[i]+=point[j]*shapeCoeff[j][i];
				}
			}

			//This is where we carry out the transformation (but instead of doing 
			//u_xi = T*u as in D&T, we did u_xi = u*T for the same T as in D&T 
			for(int i=0; i<3 ;i++){
				phi[6*i] = phi2[6*i];
				phi[6*i+1] = phi2[6*i+1]*Transform[1][1]+phi2[6*i+2]*Transform[2][1];
				phi[6*i+2] = phi2[6*i+1]*Transform[1][2]+phi2[6*i+2]*Transform[2][2];
				phi[6*i+3] = phi2[6*i+3]*Transform[3][3]+phi2[6*i+4]*Transform[4][3]+phi2[6*i+5]*Transform[5][3];
				phi[6*i+4] = phi2[6*i+3]*Transform[3][4]+phi2[6*i+4]*Transform[4][4]+phi2[6*i+5]*Transform[5][4];
				phi[6*i+5] = phi2[6*i+3]*Transform[3][5]+phi2[6*i+4]*Transform[4][5]+phi2[6*i+5]*Transform[5][5];
			}

			break; 
	}
}

//The deriv1 and deriv2 shape functions are constructed in a similar way
//except for some "jacobian factors" arising from the taking the derivatives

void deriv1(double derivx[], double derivy[], double xi_val, double eta_val, double Transform[][6], double shapeCoeff[][18], int& ndeg1)
{
	double x = xi_val;
	double y = eta_val;
	int i,j;
	double pointdxi[ndeg1];
	double pointdeta[ndeg1];
	memset(pointdxi, 0.0,sizeof(pointdxi));
	memset(pointdeta, 0.0,sizeof(pointdeta));
	double detadx = Transform[1][2]/(Transform[2][1]*Transform[1][2]-Transform[2][2]*Transform[1][1]);
	double dxidx  = -Transform[2][2]/(Transform[2][1]*Transform[1][2]-Transform[2][2]*Transform[1][1]);

	double detady = -Transform[1][1]/(Transform[2][1]*Transform[1][2]-Transform[2][2]*Transform[1][1]);
	double dxidy  = Transform[2][1]/(Transform[2][1]*Transform[1][2]-Transform[2][2]*Transform[1][1]);
	double temp1[ndeg1];
	double temp2[ndeg1];
	switch(ndeg1){
		case 3:
			pointdxi[0] = 0;
			pointdxi[1] = 1;
			pointdxi[2] = 0;

			pointdeta[0] = 0;
			pointdeta[1] = 0;
			pointdeta[2] = 1;
			for(i=0;i<3 ; i++){
				temp1[i] =0.0;
				temp2[i] =0.0;
				for(j=0;j<3 ; j++){
					temp1[i] += pointdxi[j]*shapeCoeff[j][i];
					temp2[i] += pointdeta[j]*shapeCoeff[j][i];
				}
			}
			for(i=0;i<3 ; i++){
				derivx[i] = temp1[i]*dxidx+ temp2[i]*detadx;
				derivy[i] = temp1[i]*dxidy+ temp2[i]*detady;
			}
			break;

		case 18:
			pointdxi[0]  = 0;
			pointdxi[1]  = 1;
			pointdxi[2]  = 0;
			pointdxi[3]  = 2*x;
			pointdxi[4]  = 1*y;
			pointdxi[5]  = 0;
			pointdxi[6]  = 3*x*x;
			pointdxi[7]  = 2*x*y;
			pointdxi[8]  = y*y*1;
			pointdxi[9]  = 0;
			pointdxi[10] = 4*x*x*x;
			pointdxi[11] = 3*x*x*y;
			pointdxi[12] = 2*x*y*y;
			pointdxi[13] = 1*y*y*y;
			pointdxi[14] = 0;      
			pointdxi[15] = 5*x*x*x*x-5*3*x*x*y*y;
			pointdxi[16] = 2*x*y*y*y-3*x*x*y*y;
			pointdxi[17] = -5*3*x*x*y*y;         

			pointdeta[0]  = 0;
			pointdeta[1]  = 0;
			pointdeta[2]  = 1;
			pointdeta[3]  = 0;
			pointdeta[4]  = x*1;
			pointdeta[5]  = 2*y;
			pointdeta[6]  = 0;
			pointdeta[7]  = x*x*1;
			pointdeta[8]  = 2*y*x;
			pointdeta[9]  = 3*y*y;
			pointdeta[10] = 0;
			pointdeta[11] = x*x*x*1;
			pointdeta[12] = x*x*y*2;
			pointdeta[13] = x*y*y*3;
			pointdeta[14] = y*y*y*4;
			pointdeta[15] = -5*x*x*x*2*y;
			pointdeta[16] = x*x*3*y*y-x*x*x*2*y;
			pointdeta[17] = 5*y*y*y*y-5*x*x*x*2*y;

			double phipxi[18];
			double phipeta[18];
			memset(phipxi, 0.0, sizeof(phipxi));
			memset(phipeta, 0.0, sizeof(phipxi));

			for( i=0; i<18 ;i++){ 
				for( j=0; j<18 ; j++){
					phipxi[i]+=pointdxi[j]*shapeCoeff[j][i];
					phipeta[i]+=pointdeta[j]*shapeCoeff[j][i];
				}
			}

			for( i=0;i<3;i++){
				temp1[6*i] = phipxi[6*i];
				temp1[6*i+1] = phipxi[6*i+1]*Transform[1][1]+phipxi[6*i+2]*Transform[2][1];
				temp1[6*i+2] = phipxi[6*i+1]*Transform[1][2]+phipxi[6*i+2]*Transform[2][2];
				temp1[6*i+3] = phipxi[6*i+3]*Transform[3][3]+phipxi[6*i+4]*Transform[4][3]+phipxi[6*i+5]*Transform[5][3];
				temp1[6*i+4] = phipxi[6*i+3]*Transform[3][4]+phipxi[6*i+4]*Transform[4][4]+phipxi[6*i+5]*Transform[5][4];
				temp1[6*i+5] = phipxi[6*i+3]*Transform[3][5]+phipxi[6*i+4]*Transform[4][5]+phipxi[6*i+5]*Transform[5][5];

				temp2[6*i] = phipeta[6*i];
				temp2[6*i+1] = phipeta[6*i+1]*Transform[1][1]+phipeta[6*i+2]*Transform[2][1];
				temp2[6*i+2] = phipeta[6*i+1]*Transform[1][2]+phipeta[6*i+2]*Transform[2][2];
				temp2[6*i+3] = phipeta[6*i+3]*Transform[3][3]+phipeta[6*i+4]*Transform[4][3]+phipeta[6*i+5]*Transform[5][3];
				temp2[6*i+4] = phipeta[6*i+3]*Transform[3][4]+phipeta[6*i+4]*Transform[4][4]+phipeta[6*i+5]*Transform[5][4];
				temp2[6*i+5] = phipeta[6*i+3]*Transform[3][5]+phipeta[6*i+4]*Transform[4][5]+phipeta[6*i+5]*Transform[5][5];
			}

			for(i=0;i<18 ; i++){
				derivx[i] = temp1[i]*dxidx + temp2[i]*detadx;//again we take care of the jacobian factors
				derivy[i] = temp1[i]*dxidy + temp2[i]*detady;
			}
			break;
	}
}

void deriv2(double derivxx[], double derivxy[], double derivyy[], double xi_val, double eta_val, double Transform[][6], double shapeCoeff[][18])
{
	double x = xi_val;
	double y = eta_val;

	double pointddxi[18];
	pointddxi[0]  = 0;
	pointddxi[1]  = 0;
	pointddxi[2]  = 0;
	pointddxi[3]  = 2*1;
	pointddxi[4]  = 0;
	pointddxi[5]  = 0;
	pointddxi[6]  = 3*2*x;
	pointddxi[7]  = 2*1*y;
	pointddxi[8]  = 0;
	pointddxi[9]  = 0;
	pointddxi[10] = 4*3*x*x;
	pointddxi[11] = 3*2*x*y;
	pointddxi[12] = 2*1*y*y;
	pointddxi[13] = 0;
	pointddxi[14] = 0;
	pointddxi[15] = 5*4*x*x*x-5*3*2*x*y*y;
	pointddxi[16] = 2*1*y*y*y-3*2*x*y*y;
	pointddxi[17] = -5*3*2*x*y*y;

	double pointdxideta[18];
	memset(pointdxideta, 0.0, sizeof(pointdxideta));
	pointdxideta[0]  = 0;
	pointdxideta[1]  = 0;
	pointdxideta[2]  = 0;
	pointdxideta[3]  = 0;
	pointdxideta[4]  = 1*1;
	pointdxideta[5]  = 0;
	pointdxideta[6]  = 0;
	pointdxideta[7]  = 2*x*1;
	pointdxideta[8]  = 2*y*1;
	pointdxideta[9]  = 0;
	pointdxideta[10] = 0;
	pointdxideta[11] = 3*x*x*1;
	pointdxideta[12] = 2*x*y*2;
	pointdxideta[13] = 1*y*y*3;
	pointdxideta[14] = 0;
	pointdxideta[15] = -5*3*x*x*2*y;
	pointdxideta[16] = 2*x*3*y*y-3*x*x*2*y;
	pointdxideta[17] = -5*3*x*x*2*y;

	double pointddeta[18];
	memset(pointddeta, 0.0, sizeof(pointddeta));
	pointddeta[0]  = 0;
	pointddeta[1]  = 0;
	pointddeta[2]  = 0;
	pointddeta[3]  = 0;
	pointddeta[4]  = 0;
	pointddeta[5]  = 2*1;
	pointddeta[6]  = 0;
	pointddeta[7]  = 0;
	pointddeta[8]  = 2*1*x;
	pointddeta[9]  = 3*2*y;
	pointddeta[10] = 0;
	pointddeta[11] = 0;
	pointddeta[12] = x*x*1*2;
	pointddeta[13] = x*2*y*3;
	pointddeta[14] = 3*y*y*4;
	pointddeta[15] = -5*x*x*x*2*1;
	pointddeta[16] = x*x*3*2*y-x*x*x*2*1;
	pointddeta[17] = 5*4*y*y*y-5*x*x*x*2*1;

	double detadx = Transform[1][2]/(Transform[2][1]*Transform[1][2]-Transform[2][2]*Transform[1][1]);
	double dxidx  = -Transform[2][2]/(Transform[2][1]*Transform[1][2]-Transform[2][2]*Transform[1][1]);

	double detady = -Transform[1][1]/(Transform[2][1]*Transform[1][2]-Transform[2][2]*Transform[1][1]);
	double dxidy  = Transform[2][1]/(Transform[2][1]*Transform[1][2]-Transform[2][2]*Transform[1][1]);

	double phip2xi[18];
	double phip2eta[18];
	double phipxieta[18];
	double temp1[18];
	double temp2[18];
	double temp3[18];
	memset(phip2xi, 0.0, sizeof(phip2xi));
	memset(phip2eta, 0.0, sizeof(phip2xi));
	memset(phipxieta, 0.0, sizeof(phip2xi));
	memset(temp1, 0.0, sizeof(phip2xi));
	memset(temp2, 0.0, sizeof(phip2xi));
	memset(temp3, 0.0, sizeof(phip2xi));
	for(int i=0; i<18 ;i++){
		for(int j=0; j<18 ; j++){
			phip2xi[i]+=pointddxi[j]*shapeCoeff[j][i];
			phip2eta[i]+=pointddeta[j]*shapeCoeff[j][i];
			phipxieta[i]+=pointdxideta[j]*shapeCoeff[j][i];
		}
	}

	for(int i=0;i<3;i++){
		temp1[6*i] = phip2xi[6*i];
		temp1[6*i+1] = phip2xi[6*i+1]*Transform[1][1]+phip2xi[6*i+2]*Transform[2][1];
		temp1[6*i+2] = phip2xi[6*i+1]*Transform[1][2]+phip2xi[6*i+2]*Transform[2][2];
		temp1[6*i+3] = phip2xi[6*i+3]*Transform[3][3]+phip2xi[6*i+4]*Transform[4][3]+phip2xi[6*i+5]*Transform[5][3];
		temp1[6*i+4] = phip2xi[6*i+3]*Transform[3][4]+phip2xi[6*i+4]*Transform[4][4]+phip2xi[6*i+5]*Transform[5][4];
		temp1[6*i+5] = phip2xi[6*i+3]*Transform[3][5]+phip2xi[6*i+4]*Transform[4][5]+phip2xi[6*i+5]*Transform[5][5];

		temp2[6*i] = phipxieta[6*i];
		temp2[6*i+1] = phipxieta[6*i+1]*Transform[1][1]+phipxieta[6*i+2]*Transform[2][1];
		temp2[6*i+2] = phipxieta[6*i+1]*Transform[1][2]+phipxieta[6*i+2]*Transform[2][2];
		temp2[6*i+3] = phipxieta[6*i+3]*Transform[3][3]+phipxieta[6*i+4]*Transform[4][3]+phipxieta[6*i+5]*Transform[5][3];
		temp2[6*i+4] = phipxieta[6*i+3]*Transform[3][4]+phipxieta[6*i+4]*Transform[4][4]+phipxieta[6*i+5]*Transform[5][4];
		temp2[6*i+5] = phipxieta[6*i+3]*Transform[3][5]+phipxieta[6*i+4]*Transform[4][5]+phipxieta[6*i+5]*Transform[5][5];

		temp3[6*i] = phip2eta[6*i];
		temp3[6*i+1] = phip2eta[6*i+1]*Transform[1][1]+phip2eta[6*i+2]*Transform[2][1];
		temp3[6*i+2] = phip2eta[6*i+1]*Transform[1][2]+phip2eta[6*i+2]*Transform[2][2];
		temp3[6*i+3] = phip2eta[6*i+3]*Transform[3][3]+phip2eta[6*i+4]*Transform[4][3]+phip2eta[6*i+5]*Transform[5][3];
		temp3[6*i+4] = phip2eta[6*i+3]*Transform[3][4]+phip2eta[6*i+4]*Transform[4][4]+phip2eta[6*i+5]*Transform[5][4];
		temp3[6*i+5] = phip2eta[6*i+3]*Transform[3][5]+phip2eta[6*i+4]*Transform[4][5]+phip2eta[6*i+5]*Transform[5][5];
	}
	for(int i=0; i < 18; i++){
		derivxx[i] = phip2xi[i]*dxidx*dxidx + phipxieta[i]*2*dxidx*detadx + phip2eta[i]*detadx*detadx; //same as for deriv1, we have to take care of the jacobian factors

		derivxy[i] = phip2xi[i]*dxidx*dxidy + phipxieta[i]*(dxidx*detady + dxidy*detadx) + phip2eta[i]*detadx*detady;

		derivyy[i] = phip2xi[i]*dxidy*dxidy + phipxieta[i]*2*dxidy*detady + phip2eta[i]*detady*detady;
	}
}

// NOTE: This function should be simplified
void generateShapeCoeffMatrix(double B[][18],int& ndeg1)
{
	//This functions gives us directly the shape coefficients depending on whether
	//we are using hermite or linear
	memset(B, 0.0, sizeof(B[0][0])*ndeg1*18);
	switch(ndeg1){
		case 3:
			B[0][0] = 1 , B[0][1] = 0, B[0][2] = 0;
			B[1][0] = -1, B[1][1] = 1, B[1][2] = 0;
			B[2][0] = -1, B[2][1] = 0, B[2][2] = 1;
			break;

		case 18:
			B[0][0] = 1;
			B[1][1] = 1;
			B[2][2] = 1;
			B[3][3] = 0.5;
			B[4][4] = 1;
			B[5][5] = 0.5;
			B[6][0] = -10,     B[6][1] = -6,     B[6][3] = -1.5,    B[6][6] = 10,     B[6][7] = -4,     B[6][9] = 0.5;
			B[7][2] = -3,      B[7][4] = -2,     B[7][8] = 3,	     B[7][10] = -1;
			B[8][1] = -3,      B[8][4] = -2,     B[8][5] = 0.0,     B[8][13] = 3,     B[8][16] = -1;
			B[9][0] = -10,     B[9][2] = -6,     B[9][5] = -1.5,    B[9][12] = 10,    B[9][14] = -4,    B[9][17] = 0.5;
			B[10][0] = 15,     B[10][1] = 8,     B[10][3] = 1.5,    B[10][6] = -15,   B[10][7] = 7,	 B[10][9] = -1;
			B[11][2] = 2,      B[11][4] = 1,     B[11][8] = -2,     B[11][10] = 1;
			B[12][0] = -30,    B[12][1] = -6,    B[12][2] = -6,     B[12][3] = -1.5,  B[12][4] = 2,	 B[12][5] = -1.5,
				B[12][6] = 15,     B[12][7] = -7.5,  B[12][8] = -1.5,   B[12][9] = 1.25,  B[12][10] = 0.5,
				B[12][11] = 0.25,  B[12][12] = 15,   B[12][13] = -1.5,  B[12][14] = -7.5, B[12][15] = 0.25,
				B[12][16] = 0.5,   B[12][17] = 1.25;
			B[13][1] = 2,      B[13][4] = 1,     B[13][13] = -2,    B[13][16] = 1;
			B[14][0] = 15,     B[14][2] = 8,     B[14][5] = 1.5,    B[14][12] = -15,  B[14][14] = 7,	 B[14][17] = -1;
			B[15][0] = -6,     B[15][1] = -3,    B[15][3] = -0.5,   B[15][6] = 6,     B[15][7] = -3,	 B[15][9] = 0.5;
			B[16][0] = 30,     B[16][1] = 6,     B[16][2] = 9,	     B[16][3] = 1,     B[16][5] = 1.5,	 B[16][6] = -15,
				B[16][7] = 7.5,    B[16][8] = -1.5,  B[16][9] = -1.25,  B[16][10] = 0.5,  B[16][11] = -0.25,
				B[16][12] = -15,   B[16][13] = 1.5,  B[16][14] = 7.5,   B[16][15] = 0.25, B[16][16] = -0.5,
				B[16][17] = -1.25;
			B[17][0] = -6,     B[17][2] = -3,    B[17][5] = -0.5,   B[17][12] = 6,    B[17][14] = -3,	 B[17][17] = 0.5;
			break;
	}
}

// NOTE: it might be better to have this function return the transformed 
// matrix instead of the transformation matrix, so that the matrix 
// multiplication is only done once per element.
void generateTransformationMatrix(double xp[3], double yp[3], double T[][6])
{
	//This is the Transformation Matrix as described D&T
	memset(T, 0.0, sizeof(T[0][0])*6*6);
	T[0][0] = 1;

	T[1][1] = xp[1]-xp[0]; T[1][2] = yp[1]-yp[0];
	T[2][1] = xp[2]-xp[0]; T[2][2] = yp[2]-yp[0];

	T[3][3] = (xp[1]-xp[0])*(xp[1]-xp[0]); 
	T[3][4] = 2*(xp[1]-xp[0])*(yp[1]-yp[0]);
	T[3][5] = (yp[1]-yp[0])*(yp[1]-yp[0]);

	T[4][3] = (xp[1]-xp[0])*(xp[2]-xp[0]); 
	T[4][4] = (yp[1]-yp[0])*(xp[2]-xp[0]) + (xp[1]-xp[0])*(yp[2]-yp[0]);
	T[4][5] = (yp[1]-yp[0])*(yp[2]-yp[0]);

	T[5][3] = (xp[2]-xp[0])*(xp[2]-xp[0]); 
	T[5][4] = 2*(xp[2]-xp[0])*(yp[2]-yp[0]);
	T[5][5] = (yp[2]-yp[0])*(yp[2]-yp[0]); 
}   


//=========================================================
// Sets gauss quadrature weights and coords 
// used for integrations.
//=========================================================


void trigauss(int& ng, double* xi, double* eta, double* w)
{

	/*=========================================================
	  Gauss-Dunavant quadrature is exact for polynomials in
	  xi, eta, with a maximum power of xi^m*eta^n, m+n =maxpower.
	  This requires the number of gauss-points to correspond to
	  ng. This relation is as in the following table:

	  maxpower = n + m,  where f(xi,eta) ~=~ [xi^n * eta^m]
	  ================         ============================

	  for maxpower=   1,  2,  3,  4,  5,  6,  7,  8,  9, 10
	  use ng=         1,  3,  4,  6,  7, 12, 13, 16, 19, 25

	  for maxpower=  11, 12, 13, 14, 15, 16, || 17, 18, 19, 20
	  use ng=        27, 33, 37, 42, 48, 52, || 61, 70, 73, 79
	  ---------------------------------------------------------

	  Gauss points/wts for triangles:

Reference:
D. A. Dunavant, Int. J. Num. Methods. Eng.
{\bf 21}, 1129 (1985).

I(Global_coord_triangle)

= Int(x,y,triangle)

= |Jacobian|*Int_0^1 d_eta Int_0^{1-eta} d_xi

= 2*area*(1/2) = Sum_i  w_i = area

where: |Jacobian| = 2* area

Putting xi = (1+u)/2,  and
eta = (1-xi)(1+v)/2 = (1-u)(1-v)/4

gives:

I(Global_coord_triangle)

= (A/4)*Int_{u=-1}^{u=1} (1-u)
	 * Int_{v=-1}^{v=1} f(xi,eta)*du*dv

	 = (A/4)* Sum_i wi (1-u_i) *Sum_j wj f(xi_i, eta_i)

	 This then reduces to a regular double integral given
	 in terms of a double sum over gauss pts/wts


	 In the following, p = 2 means that if f(xi,eta) has
	 powers of xi and eta upto 2 the integration is exact.

	 ---------------------------------------------------------
	 */

switch(ng)
{
	case 1://p=1
		w[0]     =  1.000000000000000e0;
		xi[0]    =  0.333333333333333e0;
		eta[0]   =  0.333333333333333e0;
		break;

	case 3://p=2
		w[0]     =  0.333333333333333e0;
		w[1]     =  0.333333333333333e0;
		w[2]     =  0.333333333333333e0;

		xi[0]    =  0.666666666666667e0;
		xi[1]    =  0.166666666666667e0;
		xi[2]    =  0.166666666666667e0;

		eta[0]   =  0.166666666666667e0;
		eta[1]   =  0.666666666666667e0;
		eta[2]   =  0.166666666666667e0;

		break;

	case 4: //p =3
		w[0]     = -0.562500000000000e0;

		xi[0]    =  0.333333333333333e0;
		eta[0]   =  0.333333333333333e0;

		w[1]     =  0.520833333333333e0;
		w[2]     =  0.520833333333333e0;
		w[3]     =  0.520833333333333e0;

		xi[1]    =  0.600000000000000e0;
		xi[2]    =  0.200000000000000e0;
		xi[3]    =  0.200000000000000e0;

		eta[1]   =  0.200000000000000e0;
		eta[2]   =  0.600000000000000e0;
		eta[3]   =  0.200000000000000e0;

		break;

	case  6://p=4
		w[0]     =  0.223381589678011e0;
		w[1]     =  0.223381589678011e0;
		w[2]     =  0.223381589678011e0;

		xi[0]    =  0.108103018168070e0;
		xi[1]    =  0.445948490915965e0;
		xi[2]    =  0.445948490915965e0;

		eta[0]   =  0.445948490915965e0;
		eta[1]   =  0.445948490915965e0;
		eta[2]   =  0.108103018168070e0;

		w[3]     =  0.109951743655322e0;
		w[4]     =  0.109951743655322e0;
		w[5]     =  0.109951743655322e0;

		xi[3]    =  0.816847572980459e0;
		xi[4]    =  0.091576213509771e0;
		xi[5]    =  0.091576213509771e0;

		eta[3]   =  0.091576213509771e0;
		eta[4]   =  0.816847572980459e0;
		eta[5]   =  0.091576213509771e0;

		break; 

	case      7://p=5
		w[0]     =  0.225000000000000e0;

		xi[0]    =  0.333333333333333e0;
		eta[0]   =  0.333333333333333e0;

		w[1]     =  0.132394152788506e0;
		w[2]     =  0.132394152788506e0;
		w[3]     =  0.132394152788506e0;

		xi[1]    =  0.059715871789770e0;
		xi[2]    =  0.470142064105115e0;
		xi[3]    =  0.470142064105115e0;

		eta[1]   =  0.470142064105115e0;
		eta[2]   =  0.059715871789770e0;
		eta[3]   =  0.470142064105115e0;

		w[4]     =  0.125939180544827e0;
		w[5]     =  0.125939180544827e0;
		w[6]     =  0.125939180544827e0;

		xi[4]    =  0.797426985353087e0;
		xi[5]    =  0.101286507323456e0;
		xi[6]    =  0.101286507323456e0;


		eta[4]   =  0.101286507323456e0;
		eta[5]   =  0.797426985353087e0;
		eta[6]   =  0.101286507323456e0;

		break;

	case      12: //p =6

		w[0]     =   0.116786275726379e0;
		w[1]     =   0.116786275726379e0;
		w[2]     =   0.116786275726379e0;

		xi[0]    =   0.501426509658179e0;
		xi[1]    =   0.249286745170910e0;
		xi[2]    =   0.249286745170910e0;

		eta[0]   =   0.249286745170910e0;
		eta[1]   =   0.501426509658179e0;
		eta[2]   =   0.249286745170910e0;

		w[3]     =   0.050844906370207e0;
		w[4]     =   0.050844906370207e0;
		w[5]     =   0.050844906370207e0;

		xi[3]    =   0.873821971016996e0;
		xi[4]    =   0.063089014491502e0;
		xi[5]    =   0.063089014491502e0;

		eta[3]   =   0.063089014491502e0;
		eta[4]   =   0.873821971016996e0;
		eta[5]   =   0.063089014491502e0;

		w[6]     =   0.082851075618374e0;
		w[7]     =   0.082851075618374e0;
		w[8]     =   0.082851075618374e0;
		w[9]     =   0.082851075618374e0;
		w[10]    =   0.082851075618374e0;
		w[11]    =   0.082851075618374e0;

		xi[6]    =   0.053145049844817e0;
		xi[7]    =   0.053145049844817e0;
		xi[8]    =   0.310352451033784e0;
		xi[9]    =   0.310352451033784e0;
		xi[10]   =   0.636502499121399e0;
		xi[11]   =   0.636502499121399e0;

		eta[6]   =   0.310352451033784e0;
		eta[7]   =   0.636502499121399e0;
		eta[8]   =   0.053145049844817e0;
		eta[9]   =   0.636502499121399e0;
		eta[10]  =   0.053145049844817e0;
		eta[11]  =   0.310352451033784e0;

		break;

	case 13:
		//p = 7
		w[0]     =  -0.149570044467682e0;
		xi[0]    =   0.333333333333333e0;
		eta[0]   =   0.333333333333333e0;

		w[1]     =   0.175615257433208e0;
		w[2]     =   0.175615257433208e0;
		w[3]     =   0.175615257433208e0;

		xi[1]    =   0.479308067841920e0;
		xi[2]    =   0.260345966079040e0;
		xi[3]    =   0.260345966079040e0;

		eta[1]   =   0.260345966079040e0;
		eta[2]   =   0.479308067841920e0;
		eta[3]   =   0.260345966079040e0;

		w[4]     =   0.053347235608838e0;
		w[5]     =   0.053347235608838e0;
		w[6]     =   0.053347235608838e0;

		xi[4]    =   0.869739794195568e0;
		xi[5]    =   0.065130102902216e0;
		xi[6]    =   0.065130102902216e0;

		eta[4]   =   0.065130102902216e0;
		eta[5]   =   0.869739794195568e0;
		eta[6]   =   0.065130102902216e0;

		w[7]     =   0.077113760890257e0;
		w[8]     =   0.077113760890257e0;
		w[9]     =   0.077113760890257e0;
		w[10]    =   0.077113760890257e0;
		w[11]    =   0.077113760890257e0;
		w[12]    =   0.077113760890257e0;

		xi[7]    =   0.048690315425316e0;
		xi[8]    =   0.048690315425316e0;
		xi[9]    =   0.638444188569810e0;
		xi[10]   =   0.638444188569810e0;
		xi[11]   =   0.312865496004874e0;
		xi[12]   =   0.312865496004874e0;

		eta[7]   =   0.312865496004874e0;
		eta[8]   =   0.638444188569810e0;
		eta[9]   =   0.048690315425316e0;
		eta[10]  =   0.312865496004874e0;
		eta[11]  =   0.048690315425316e0;
		eta[12]  =   0.638444188569810e0;

		break;
	case     16: //p =8 
		w[0]     =   0.144315607677787e0;
		xi[0]    =   0.333333333333333e0;
		eta[0]   =   0.333333333333333e0;

		w[1]     =   0.095091634267285e0;
		w[2]     =   0.095091634267285e0;
		w[3]     =   0.095091634267285e0;

		xi[1]    =   0.081414823414554e0;
		xi[2]    =   0.459292588292723e0;
		xi[3]    =   0.459292588292723e0;

		eta[1]   =   0.459292588292723e0;
		eta[2]   =   0.081414823414554e0;
		eta[3]   =   0.459292588292723e0;

		w[4]     =   0.103217370534718e0;
		w[5]     =   0.103217370534718e0;
		w[6]     =   0.103217370534718e0;

		xi[4]    =   0.170569307751760e0;
		xi[5]    =   0.658861384496480e0;
		xi[6]    =   0.170569307751760e0;

		eta[4]   =   0.658861384496480e0;
		eta[5]   =   0.170569307751760e0;
		eta[6]   =   0.170569307751760e0;

		w[7]     =   0.032458497623198e0;
		w[8]     =   0.032458497623198e0;
		w[9]     =   0.032458497623198e0;

		xi[7]    =   0.898905543365938e0;
		xi[8]    =   0.050547228317031e0;
		xi[9]    =   0.050547228317031e0;

		eta[7]   =   0.050547228317031e0;
		eta[8]   =   0.898905543365938e0;
		eta[9]   =   0.050547228317031e0;

		w[10]    =   0.027230314174435e0;
		w[11]    =   0.027230314174435e0;
		w[12]    =   0.027230314174435e0;
		w[13]    =   0.027230314174435e0;
		w[14]    =   0.027230314174435e0;
		w[15]    =   0.027230314174435e0;

		xi[10]   =   0.263112829634638e0;
		xi[11]   =   0.728492392955404e0;
		xi[12]   =   0.008394777409958e0;
		xi[13]   =   0.728492392955404e0;
		xi[14]   =   0.263112829634638e0;
		xi[15]   =   0.008394777409958e0;

		eta[10]  =   0.008394777409958e0;
		eta[11]  =   0.008394777409958e0;
		eta[12]  =   0.263112829634638e0;
		eta[13]  =   0.263112829634638e0;
		eta[14]  =   0.728492392955404e0;
		eta[15]  =   0.728492392955404e0;

		break;

	case     19://p = 9

		w[0]     =   0.097135796282799e0;
		xi[0]    =   0.333333333333333e0;
		eta[0]   =   0.333333333333333e0;

		w[1]     =   0.031334700227139e0;
		w[2]     =   0.031334700227139e0;
		w[3]     =   0.031334700227139e0;

		xi[1]    =   0.020634961602525e0;
		xi[2]    =   0.489682519198738e0;
		xi[3]    =   0.489682519198738e0;

		eta[1]   =   0.489682519198738e0;
		eta[2]   =   0.020634961602525e0;
		eta[3]   =   0.489682519198738e0;

		w[4]     =   0.077827541004774e0;
		w[5]     =   0.077827541004774e0;
		w[6]     =   0.077827541004774e0;

		xi[4]    =   0.125820817014127e0;
		xi[5]    =   0.437089591492937e0;
		xi[6]    =   0.437089591492937e0;

		eta[4]   =   0.437089591492937e0;
		eta[5]   =   0.125820817014127e0;
		eta[6]   =   0.437089591492937e0;

		w[7]     =   0.079647738927210e0;
		w[8]     =   0.079647738927210e0;
		w[9]     =   0.079647738927210e0;

		xi[7]    =   0.623592928761935e0;
		xi[8]    =   0.188203535619033e0;
		xi[9]    =   0.188203535619033e0;

		eta[7]   =   0.188203535619033e0;
		eta[8]   =   0.623592928761935e0;
		eta[9]   =   0.188203535619033e0;

		w[10]    =   0.025577675658698e0;
		w[11]    =   0.025577675658698e0;
		w[12]    =   0.025577675658698e0;

		xi[10]   =   0.910540973211095e0;
		xi[11]   =   0.044729513394453e0;
		xi[12]   =   0.044729513394453e0;

		eta[10]  =   0.044729513394453e0;
		eta[11]  =   0.910540973211095e0;
		eta[12]  =   0.044729513394453e0;

		w[13]    =   0.043283539377289e0;
		w[14]    =   0.043283539377289e0;
		w[15]    =   0.043283539377289e0;
		w[16]    =   0.043283539377289e0;
		w[17]    =   0.043283539377289e0;
		w[18]    =   0.043283539377289e0;

		xi[13]   =   0.036838412054736e0;
		xi[14]   =   0.036838412054736e0;
		xi[15]   =   0.221962989160766e0;
		xi[16]   =   0.221962989160766e0;
		xi[17]   =   0.741198598784498e0;
		xi[18]   =   0.741198598784498e0;

		eta[13]  =   0.221962989160766e0;
		eta[14]  =   0.741198598784498e0;
		eta[15]  =   0.036838412054736e0;
		eta[16]  =   0.741198598784498e0;
		eta[17]  =   0.221962989160766e0;
		eta[18]  =   0.036838412054736e0;

		break;

	case 25://p=10

		w[0]     =   0.090817990382754e0;
		xi[0]    =   0.333333333333333e0;
		eta[0]   =   0.333333333333333e0;

		w[1]     =   0.036725957756467e0;
		w[2]     =   0.036725957756467e0;
		w[3]     =   0.036725957756467e0;

		xi[1]    =   0.028844733232685e0;
		xi[2]    =   0.485577633383657e0;
		xi[3]    =   0.485577633383657e0;

		eta[1]   =   0.485577633383657e0;
		eta[2]   =   0.028844733232685e0;
		eta[3]   =   0.485577633383657e0;

		w[4]     =   0.045321059435528e0;
		w[5]     =   0.045321059435528e0;
		w[6]     =   0.045321059435528e0;

		xi[4]    =   0.781036849029926e0;
		xi[5]    =   0.109481575485037e0;
		xi[6]    =   0.109481575485037e0;

		eta[4]   =   0.109481575485037e0;
		eta[5]   =   0.781036849029926e0;
		eta[6]   =   0.109481575485037e0;

		w[7]     =   0.072757916845420e0;
		w[8]     =   0.072757916845420e0;
		w[9]     =   0.072757916845420e0;
		w[10]    =   0.072757916845420e0;
		w[11]    =   0.072757916845420e0;
		w[12]    =   0.072757916845420e0;

		xi[7]    =   0.141707219414880e0;
		xi[8]    =   0.141707219414880e0;
		xi[9]    =   0.550352941820999e0;
		xi[10]   =   0.550352941820999e0;
		xi[11]   =   0.307939838764121e0;
		xi[12]   =   0.307939838764121e0;

		eta[7]   =   0.307939838764121e0;
		eta[8]   =   0.550352941820999e0;
		eta[9]   =   0.141707219414880e0;
		eta[10]  =   0.307939838764121e0;
		eta[11]  =   0.550352941820999e0;
		eta[12]  =   0.141707219414880e0;

		w[13]    =   0.028327242531057e0;
		w[14]    =   0.028327242531057e0;
		w[15]    =   0.028327242531057e0;
		w[16]    =   0.028327242531057e0;
		w[17]    =   0.028327242531057e0;
		w[18]    =   0.028327242531057e0;

		xi[13]   =  0.025003534762686e0;
		xi[14]   =  0.025003534762686e0;
		xi[15]   =  0.246672560639903e0;
		xi[16]   =  0.246672560639903e0;
		xi[17]   =  0.728323904597411e0;
		xi[18]   =  0.728323904597411e0;

		eta[13]  =   0.246672560639903e0;
		eta[14]  =   0.728323904597411e0;
		eta[15]  =   0.025003534762686e0;
		eta[16]  =   0.728323904597411e0;
		eta[17]  =   0.025003534762686e0;
		eta[18]  =   0.246672560639903e0;

		w[19]    =   0.009421666963733e0;
		w[20]    =   0.009421666963733e0;
		w[21]    =   0.009421666963733e0;
		w[22]    =   0.009421666963733e0;
		w[23]    =   0.009421666963733e0;
		w[24]    =   0.009421666963733e0;

		xi[19]   =   0.009540815400299e0;
		xi[20]   =   0.009540815400299e0;
		xi[21]   =   0.923655933587500e0;
		xi[22]   =   0.923655933587500e0;
		xi[23]   =   0.066803251012200e0;
		xi[24]   =   0.066803251012200e0;

		eta[19]  =   0.066803251012200e0;
		eta[20]  =   0.923655933587500e0;
		eta[21]  =   0.009540815400299e0;
		eta[22]  =   0.066803251012200e0;
		eta[23]  =   0.009540815400299e0;
		eta[24]  =   0.923655933587500e0;

		break;
	case  27: //p =11

		w[0]     =   0.000927006328961e0;
		w[1]     =   0.000927006328961e0;
		w[2]     =   0.000927006328961e0;

		xi[0]    =  -0.069222096541517e0;
		xi[1]    =   0.534611048270758e0;
		xi[2]    =   0.534611048270758e0;

		eta[0]   =   0.534611048270758e0;
		eta[1]   =  -0.069222096541517e0;
		eta[2]   =   0.534611048270758e0;

		w[3]     =   0.077149534914813e0;
		w[4]     =   0.077149534914813e0;
		w[5]     =   0.077149534914813e0;

		xi[3]    =   0.202061394068290e0;
		xi[4]    =   0.398969302965855e0;
		xi[5]    =   0.398969302965855e0;

		eta[3]   =   0.398969302965855e0;
		eta[4]   =   0.202061394068290e0;
		eta[5]   =   0.398969302965855e0;

		w[6]     =   0.059322977380774e0;
		w[7]     =   0.059322977380774e0;
		w[8]     =   0.059322977380774e0;

		xi[6]    =   0.593380199137435e0;
		xi[7]    =   0.203309900431282e0;
		xi[8]    =   0.203309900431282e0;

		eta[6]   =   0.203309900431282e0;
		eta[7]   =   0.593380199137435e0;
		eta[8]   =   0.203309900431282e0;

		w[9]     =   0.036184540503418e0;
		w[10]    =   0.036184540503418e0;
		w[11]    =   0.036184540503418e0;

		xi[9]    =   0.761298175434837e0;
		xi[10]   =   0.119350912282581e0;
		xi[11]   =   0.119350912282581e0;

		eta[9]   =   0.119350912282581e0;
		eta[10]  =   0.761298175434837e0;
		eta[11]  =   0.119350912282581e0;

		w[12]    =   0.013659731002678e0;
		w[13]    =   0.013659731002678e0;
		w[14]    =   0.013659731002678e0;

		xi[12]   =   0.935270103777448e0;
		xi[13]   =   0.032364948111276e0;
		xi[14]   =   0.032364948111276e0;

		eta[12]  =   0.032364948111276e0;
		eta[13]  =   0.935270103777448e0;
		eta[14]  =   0.032364948111276e0;

		w[15]    =   0.052337111962204e0;
		w[16]    =   0.052337111962204e0;
		w[17]    =   0.052337111962204e0;
		w[18]    =   0.052337111962204e0;
		w[19]    =   0.052337111962204e0;
		w[20]    =   0.052337111962204e0;

		xi[15]   =   0.050178138310495e0;
		xi[16]   =   0.050178138310495e0;
		xi[17]   =   0.593201213428213e0;
		xi[18]   =   0.593201213428213e0;
		xi[19]   =   0.356620648261293e0;
		xi[20]   =   0.356620648261293e0;

		eta[15]  =   0.356620648261293e0;
		eta[16]  =   0.593201213428213e0;
		eta[17]  =   0.050178138310495e0;
		eta[18]  =   0.356620648261293e0;
		eta[19]  =   0.593201213428213e0;
		eta[20]  =   0.050178138310495e0;

		w[21]    =   0.020707659639141e0;
		w[22]    =   0.020707659639141e0;
		w[23]    =   0.020707659639141e0;
		w[24]    =   0.020707659639141e0;
		w[25]    =   0.020707659639141e0;
		w[26]    =   0.020707659639141e0;

		xi[21]   =   0.021022016536166e0;
		xi[22]   =   0.021022016536166e0;
		xi[23]   =   0.807489003159792e0;
		xi[24]   =   0.807489003159792e0;
		xi[25]   =   0.171488980304042e0;
		xi[26]   =   0.171488980304042e0;

		eta[21]  =   0.171488980304042e0;
		eta[22]  =   0.807489003159792e0;
		eta[23]  =   0.021022016536166e0;
		eta[24]  =   0.171488980304042e0;
		eta[25]  =   0.021022016536166e0;
		eta[26]  =   0.807489003159792e0;

		break;
	case 33: // p = 12

		w[0]     =   0.025731066440455e0;
		w[1]     =   0.025731066440455e0;
		w[2]     =   0.025731066440455e0;

		xi[0]    =   0.023565220452390e0;
		xi[1]    =   0.488217389773805e0;
		xi[2]    =   0.488217389773805e0;

		eta[0]   =   0.488217389773805e0;
		eta[1]   =   0.023565220452390e0;
		eta[2]   =   0.488217389773805e0;

		w[3]     =   0.043692544538038e0;
		w[4]     =   0.043692544538038e0;
		w[5]     =   0.043692544538038e0;

		xi[3]    =   0.120551215411079e0;
		xi[4]    =   0.439724392294460e0;
		xi[5]    =   0.439724392294460e0;

		eta[3]   =   0.439724392294460e0;
		eta[4]   =   0.120551215411079e0;
		eta[5]   =   0.439724392294460e0;

		w[6]     =   0.062858224217885e0;
		w[7]     =   0.062858224217885e0;
		w[8]     =   0.062858224217885e0;

		xi[6]    =   0.457579229975768e0;
		xi[7]    =   0.271210385012116e0;
		xi[8]    =   0.271210385012116e0;

		eta[6]   =   0.271210385012116e0;
		eta[7]   =   0.457579229975768e0;
		eta[8]   =   0.271210385012116e0;

		w[9]     =   0.034796112930709e0;
		w[10]    =   0.034796112930709e0;
		w[11]    =   0.034796112930709e0;

		xi[9]    =   0.744847708916828e0;
		xi[10]   =   0.127576145541586e0;
		xi[11]   =   0.127576145541586e0;

		eta[9]   =   0.127576145541586e0;
		eta[10]  =   0.744847708916828e0;
		eta[11]  =   0.127576145541586e0;

		w[12]    =   0.006166261051559e0;
		w[13]    =   0.006166261051559e0;
		w[14]    =   0.006166261051559e0;

		xi[12]   =   0.957365299093579e0;
		xi[13]   =   0.021317350453210e0;
		xi[14]   =   0.021317350453210e0;

		eta[12]  =   0.021317350453210e0;
		eta[13]  =   0.957365299093579e0;
		eta[14]  =   0.021317350453210e0;

		w[15]    =   0.040371557766381e0;
		w[16]    =   0.040371557766381e0;
		w[17]    =   0.040371557766381e0;
		w[18]    =   0.040371557766381e0;
		w[19]    =   0.040371557766381e0;
		w[20]    =   0.040371557766381e0;

		xi[15]   =   0.115343494534698e0;
		xi[16]   =   0.115343494534698e0;
		xi[17]   =   0.608943235779788e0;
		xi[18]   =   0.608943235779788e0;
		xi[19]   =   0.275713269685514e0;
		xi[20]   =   0.275713269685514e0;

		eta[15]  =   0.608943235779788e0;
		eta[16]  =   0.275713269685514e0;
		eta[17]  =   0.115343494534698e0;
		eta[18]  =   0.275713269685514e0;
		eta[19]  =   0.115343494534698e0;
		eta[20]  =   0.608943235779788e0;

		w[21]    =   0.022356773202303e0;
		w[22]    =   0.022356773202303e0;
		w[23]    =   0.022356773202303e0;
		w[24]    =   0.022356773202303e0;
		w[25]    =   0.022356773202303e0;
		w[26]    =   0.022356773202303e0;

		xi[21]   =   0.022838332222257e0;
		xi[22]   =   0.022838332222257e0;
		xi[23]   =   0.695836086787803e0;
		xi[24]   =   0.695836086787803e0;
		xi[25]   =   0.281325580989940e0;
		xi[26]   =   0.281325580989940e0;

		eta[21]  =   0.281325580989940e0;
		eta[22]  =   0.695836086787803e0;
		eta[23]  =   0.022838332222257e0;
		eta[24]  =   0.281325580989940e0;
		eta[25]  =   0.022838332222257e0;
		eta[26]  =   0.695836086787803e0;

		w[27]    =   0.017316231108659e0;
		w[28]    =   0.017316231108659e0;
		w[29]    =   0.017316231108659e0;
		w[30]    =   0.017316231108659e0;
		w[31]    =   0.017316231108659e0;
		w[32]    =   0.017316231108659e0;

		xi[27]   =   0.025734050548330e0;
		xi[28]   =   0.025734050548330e0;
		xi[29]   =   0.858014033544073e0;
		xi[30]   =   0.858014033544073e0;
		xi[31]   =   0.116251915907597e0;
		xi[32]   =   0.116251915907597e0;

		eta[27]  =   0.116251915907597e0;
		eta[28]  =   0.858014033544073e0;
		eta[29]  =   0.025734050548330e0;
		eta[30]  =   0.116251915907597e0;
		eta[31]  =   0.858014033544073e0;
		eta[32]  =   0.025734050548330e0;

		break;

	case 37: // p = 13

		w[0]     =   0.052520923400802e0;
		xi[0]    =   0.333333333333333e0;
		eta[0]   =   0.333333333333333e0;

		w[1]     =   0.011280145209330e0;
		w[2]     =   0.011280145209330e0;
		w[3]     =   0.011280145209330e0;

		xi[1]    =   0.009903630120591e0;
		xi[2]    =   0.495048184939705e0;
		xi[3]    =   0.495048184939705e0;

		eta[1]   =   0.495048184939705e0;
		eta[2]   =   0.009903630120591e0;
		eta[3]   =   0.495048184939705e0;

		w[4]     =   0.031423518362454e0;
		w[5]     =   0.031423518362454e0;
		w[6]     =   0.031423518362454e0;

		xi[4]    =   0.062566729780852e0;
		xi[5]    =   0.468716635109574e0;
		xi[6]    =   0.468716635109574e0;

		eta[4]   =   0.468716635109574e0;
		eta[5]   =   0.062566729780852e0;
		eta[6]   =   0.468716635109574e0;

		w[7]     =   0.047072502504194e0;
		w[8]     =   0.047072502504194e0;
		w[9]     =   0.047072502504194e0;

		xi[7]    =   0.170957326397447e0;
		xi[8]    =   0.414521336801277e0;
		xi[9]    =   0.414521336801277e0;

		eta[7]   =   0.414521336801277e0;
		eta[8]   =   0.170957326397447e0;
		eta[9]   =   0.414521336801277e0;

		w[10]    =   0.047363586536355e0;
		w[11]    =   0.047363586536355e0;
		w[12]    =   0.047363586536355e0;

		xi[10]   =   0.541200855914337e0;
		xi[11]   =   0.229399572042831e0;
		xi[12]   =   0.229399572042831e0;

		eta[10]  =   0.229399572042831e0;
		eta[11]  =   0.541200855914337e0;
		eta[12]  =   0.229399572042831e0;

		w[13]    =   0.031167529045794e0;
		w[14]    =   0.031167529045794e0;
		w[15]    =   0.031167529045794e0;

		xi[13]   =   0.771151009607340e0;
		xi[14]   =   0.114424495196330e0;
		xi[15]   =   0.114424495196330e0;

		eta[13]  =   0.114424495196330e0;
		eta[14]  =   0.771151009607340e0;
		eta[15]  =   0.114424495196330e0;

		w[16]    =   0.007975771465074e0;
		w[17]    =   0.007975771465074e0;
		w[18]    =   0.007975771465074e0;

		xi[16]   =   0.950377217273082e0;
		xi[17]   =   0.024811391363459e0;
		xi[18]   =   0.024811391363459e0;

		eta[16]  =   0.024811391363459e0;
		eta[17]  =   0.950377217273082e0;
		eta[18]  =   0.024811391363459e0;

		w[19]    =   0.036848402728732e0;
		w[20]    =   0.036848402728732e0;
		w[21]    =   0.036848402728732e0;
		w[22]    =   0.036848402728732e0;
		w[23]    =   0.036848402728732e0;
		w[24]    =   0.036848402728732e0;

		xi[19]   =   0.094853828379579e0;
		xi[20]   =   0.094853828379579e0;
		xi[21]   =   0.636351174561660e0;
		xi[22]   =   0.636351174561660e0;
		xi[23]   =   0.268794997058761e0;
		xi[24]   =   0.268794997058761e0;

		eta[19]  =   0.268794997058761e0;
		eta[20]  =   0.636351174561660e0;
		eta[21]  =   0.094853828379579e0;
		eta[22]  =   0.268794997058761e0;
		eta[23]  =   0.636351174561660e0;
		eta[24]  =   0.094853828379579e0;

		w[25]    =   0.017401463303822e0;
		w[26]    =   0.017401463303822e0;
		w[27]    =   0.017401463303822e0;
		w[28]    =   0.017401463303822e0;
		w[29]    =   0.017401463303822e0;
		w[30]    =   0.017401463303822e0;

		xi[25]   =   0.018100773278807e0;
		xi[26]   =   0.018100773278807e0;
		xi[27]   =   0.690169159986905e0;
		xi[28]   =   0.690169159986905e0;
		xi[29]   =   0.291730066734288e0;
		xi[30]   =   0.291730066734288e0;

		eta[25]  =   0.291730066734288e0;
		eta[26]  =   0.690169159986905e0;
		eta[27]  =   0.018100773278807e0;
		eta[28]  =   0.291730066734288e0;
		eta[29]  =   0.690169159986905e0;
		eta[30]  =   0.018100773278807e0;

		w[31]    =   0.015521786839045e0;
		w[32]    =   0.015521786839045e0;
		w[33]    =   0.015521786839045e0;
		w[34]    =   0.015521786839045e0;
		w[35]    =   0.015521786839045e0;
		w[36]    =   0.015521786839045e0;

		xi[31]   =   0.022233076674090e0;
		xi[32]   =   0.022233076674090e0;
		xi[33]   =   0.851409537834241e0;
		xi[34]   =   0.851409537834241e0;
		xi[35]   =   0.126357385491669e0;
		xi[36]   =   0.126357385491669e0;

		eta[31]  =   0.126357385491669e0;
		eta[32]  =   0.851409537834241e0;
		eta[33]  =   0.022233076674090e0;
		eta[34]  =   0.126357385491669e0;
		eta[35]  =   0.851409537834241e0;
		eta[36]  =   0.022233076674090e0;

		break;

	case 42: //p =14
		w[0]     =   0.021883581369429e0;
		w[1]     =   0.021883581369429e0;
		w[2]     =   0.021883581369429e0;

		xi[0]    =   0.022072179275643e0;
		xi[1]    =   0.488963910362179e0;
		xi[2]    =   0.488963910362179e0;

		eta[0]   =   0.488963910362179e0;
		eta[1]   =   0.022072179275643e0;
		eta[2]   =   0.488963910362179e0;

		w[3]     =    0.032788353544125e0;
		w[4]     =    0.032788353544125e0;
		w[5]     =    0.032788353544125e0;

		xi[3]    =    0.164710561319092e0;
		xi[4]    =    0.417644719340454e0;
		xi[5]    =    0.417644719340454e0;

		eta[3]   =    0.417644719340454e0;
		eta[4]   =    0.164710561319092e0;
		eta[5]   =    0.417644719340454e0;

		w[6]     =   0.051774104507292e0;
		w[7]     =   0.051774104507292e0;
		w[8]     =   0.051774104507292e0;

		xi[6]    =   0.453044943382323e0;
		xi[7]    =   0.273477528308839e0;
		xi[8]    =   0.273477528308839e0;

		eta[6]   =   0.273477528308839e0;
		eta[7]   =   0.453044943382323e0;
		eta[8]   =   0.273477528308839e0;

		w[9]    =   0.042162588736993e0;
		w[10]    =   0.042162588736993e0;
		w[11]    =   0.042162588736993e0;

		xi[9]   =   0.645588935174913e0;
		xi[10]   =   0.177205532412543e0;
		xi[11]   =   0.177205532412543e0;

		eta[9]  =   0.177205532412543e0;
		eta[10]  =   0.645588935174913e0;
		eta[11]  =   0.177205532412543e0;

		w[12]    =   0.014433699669777e0;
		w[13]    =   0.014433699669777e0;
		w[14]    =   0.014433699669777e0;

		xi[12]   =   0.876400233818255e0;
		xi[13]   =   0.061799883090873e0;
		xi[14]   =   0.061799883090873e0;

		eta[12]  =   0.061799883090873e0;
		eta[13]  =   0.876400233818255e0;
		eta[14]  =   0.061799883090873e0;

		w[15]    =   0.004923403602400e0;
		w[16]    =   0.004923403602400e0;
		w[17]    =   0.004923403602400e0;

		xi[15]   =   0.961218077502598e0;
		xi[16]   =   0.019390961248701e0;
		xi[17]   =   0.019390961248701e0;

		eta[15]  =   0.019390961248701e0;
		eta[16]  =   0.961218077502598e0;
		eta[17]  =   0.019390961248701e0;

		w[18]    =   0.024665753212564e0;
		w[19]    =   0.024665753212564e0;
		w[20]    =   0.024665753212564e0;
		w[21]    =   0.024665753212564e0;
		w[22]    =   0.024665753212564e0;
		w[23]    =   0.024665753212564e0;

		xi[18]   =   0.057124757403648e0;
		xi[19]   =   0.057124757403648e0;
		xi[20]   =   0.770608554774996e0;
		xi[21]   =   0.770608554774996e0;
		xi[22]   =   0.172266687821356e0;
		xi[23]   =   0.172266687821356e0;

		eta[18]  =   0.172266687821356e0;
		eta[19]  =   0.770608554774996e0;
		eta[20]  =   0.057124757403648e0;
		eta[21]  =   0.172266687821356e0;
		eta[22]  =   0.057124757403648e0;
		eta[23]  =   0.770608554774996e0;

		w[24]    =   0.038571510787061e0;
		w[25]    =   0.038571510787061e0;
		w[26]    =   0.038571510787061e0;
		w[27]    =   0.038571510787061e0;
		w[28]    =   0.038571510787061e0;
		w[29]    =   0.038571510787061e0;

		xi[24]   =   0.092916249356972e0;
		xi[25]   =   0.092916249356972e0;
		xi[26]   =   0.570222290846683e0;
		xi[27]   =   0.570222290846683e0;
		xi[28]   =   0.336861459796345e0;
		xi[29]   =   0.336861459796345e0;

		eta[24]  =   0.336861459796345e0;
		eta[25]  =   0.570222290846683e0;
		eta[26]  =   0.092916249356972e0;
		eta[27]  =   0.336861459796345e0;
		eta[28]  =   0.092916249356972e0;
		eta[29]  =   0.570222290846683e0;

		w[30]    =   0.014436308113534e0;
		w[31]    =   0.014436308113534e0;
		w[32]    =   0.014436308113534e0;
		w[33]    =   0.014436308113534e0;
		w[34]    =   0.014436308113534e0;
		w[35]    =   0.014436308113534e0;

		xi[30]   =   0.014646950055654e0;
		xi[31]   =   0.014646950055654e0;
		xi[32]   =   0.686980167808088e0;
		xi[33]   =   0.686980167808088e0;
		xi[34]   =   0.298372882136258e0;
		xi[35]   =   0.298372882136258e0;

		eta[30]  =   0.298372882136258e0;
		eta[31]  =   0.686980167808088e0;
		eta[32]  =   0.014646950055654e0;
		eta[33]  =   0.298372882136258e0;
		eta[34]  =   0.014646950055654e0;
		eta[35]  =   0.686980167808088e0;

		w[36]    =   0.005010228838501e0;
		w[37]    =   0.005010228838501e0;
		w[38]    =   0.005010228838501e0;
		w[39]    =   0.005010228838501e0;
		w[40]    =   0.005010228838501e0;
		w[41]    =   0.005010228838501e0;

		xi[36]   =   0.001268330932872e0;
		xi[37]   =   0.001268330932872e0;
		xi[38]   =   0.118974497696957e0;
		xi[39]   =   0.118974497696957e0;
		xi[40]   =   0.879757171370171e0;
		xi[41]   =   0.879757171370171e0;

		eta[36]  =   0.118974497696957e0;
		eta[37]  =   0.879757171370171e0;
		eta[38]  =   0.001268330932872e0;
		eta[39]  =   0.879757171370171e0;
		eta[40]  =   0.001268330932872e0;
		eta[41]  =   0.118974497696957e0;

		break;

	case 48://p =15

		w[0]     =   0.001916875642849e0;
		w[1]     =   0.001916875642849e0;
		w[2]     =   0.001916875642849e0;

		xi[0]    =   -0.013945833716486e0;
		xi[1]    =    0.506972916858243e0;
		xi[2]    =    0.506972916858243e0;

		eta[0]   =    0.506972916858243e0;
		eta[1]   =   -0.013945833716486e0;
		eta[2]   =    0.506972916858243e0;

		w[3]     =    0.044249027271145e0;
		w[4]     =    0.044249027271145e0;
		w[5]     =    0.044249027271145e0;

		xi[3]    =    0.137187291433955e0;
		xi[4]    =    0.431406354283023e0;
		xi[5]    =    0.431406354283023e0;

		eta[3]   =    0.431406354283023e0;
		eta[4]   =    0.137187291433955e0;
		eta[5]   =    0.431406354283023e0;

		w[6]     =   0.051186548718852e0;
		w[7]     =   0.051186548718852e0;
		w[8]     =   0.051186548718852e0;

		xi[6]    =   0.444612710305711e0;
		xi[7]    =   0.277693644847144e0;
		xi[8]    =   0.277693644847144e0;

		eta[6]   =   0.277693644847144e0;
		eta[7]   =   0.444612710305711e0;
		eta[8]   =   0.277693644847144e0;

		w[9]    =   0.023687735870688e0;
		w[10]    =   0.023687735870688e0;
		w[11]    =   0.023687735870688e0;

		xi[9]   =   0.747070217917492e0;
		xi[10]   =   0.126464891041254e0;
		xi[11]   =   0.126464891041254e0;

		eta[9]  =   0.126464891041254e0;
		eta[10]  =   0.747070217917492e0;
		eta[11]  =   0.126464891041254e0;

		w[12]    =   0.013289775690021e0;
		w[13]    =   0.013289775690021e0;
		w[14]    =   0.013289775690021e0;

		xi[12]   =   0.858383228050628e0;
		xi[13]   =   0.070808385974686e0;
		xi[14]   =   0.070808385974686e0;

		eta[12]  =   0.070808385974686e0;
		eta[13]  =   0.858383228050628e0;
		eta[14]  =   0.070808385974686e0;

		w[15]    =   0.004748916608192e0;
		w[16]    =   0.004748916608192e0;
		w[17]    =   0.004748916608192e0;

		xi[15]   =   0.962069659517853e0;
		xi[16]   =   0.018965170241073e0;
		xi[17]   =   0.018965170241073e0;

		eta[15]  =   0.018965170241073e0;
		eta[16]  =   0.962069659517853e0;
		eta[17]  =   0.018965170241073e0;

		w[18]    =   0.038550072599593e0;
		w[19]    =   0.038550072599593e0;
		w[20]    =   0.038550072599593e0;
		w[21]    =   0.038550072599593e0;
		w[22]    =   0.038550072599593e0;
		w[23]    =   0.038550072599593e0;

		xi[18]   =   0.133734161966621e0;
		xi[19]   =   0.133734161966621e0;
		xi[20]   =   0.604954466893291e0;
		xi[21]   =   0.604954466893291e0;
		xi[22]   =   0.261311371140087e0;
		xi[23]   =   0.261311371140087e0;

		eta[18]  =   0.261311371140087e0;
		eta[19]  =   0.604954466893291e0;
		eta[20]  =   0.133734161966621e0;
		eta[21]  =   0.261311371140087e0;
		eta[22]  =   0.133734161966621e0;
		eta[23]  =   0.604954466893291e0;

		w[24]    =   0.027215814320624e0;
		w[25]    =   0.027215814320624e0;
		w[26]    =   0.027215814320624e0;
		w[27]    =   0.027215814320624e0;
		w[28]    =   0.027215814320624e0;
		w[29]    =   0.027215814320624e0;

		xi[24]   =   0.036366677396917e0;
		xi[25]   =   0.036366677396917e0;
		xi[26]   =   0.575586555512814e0;
		xi[27]   =   0.575586555512814e0;
		xi[28]   =   0.388046767090269e0;
		xi[29]   =   0.388046767090269e0;


		eta[24]  =   0.388046767090269e0;
		eta[25]  =   0.575586555512814e0;
		eta[26]  =   0.036366677396917e0;
		eta[27]  =   0.388046767090269e0;
		eta[28]  =   0.036366677396917e0;
		eta[29]  =   0.575586555512814e0;

		w[30]    =   0.002182077366797e0;
		w[31]    =   0.002182077366797e0;
		w[32]    =   0.002182077366797e0;
		w[33]    =   0.002182077366797e0;
		w[34]    =   0.002182077366797e0;
		w[35]    =   0.002182077366797e0;

		xi[30]   =  -0.010174883126571e0;
		xi[31]   =  -0.010174883126571e0;
		xi[32]   =   0.724462663076655e0;
		xi[33]   =   0.724462663076655e0;
		xi[34]   =   0.285712220049916e0;
		xi[35]   =   0.285712220049916e0;

		eta[30]  =   0.285712220049916e0;
		eta[31]  =   0.724462663076655e0;
		eta[32]  =  -0.010174883126571e0;
		eta[33]  =   0.285712220049916e0;
		eta[34]  =  -0.010174883126571e0;
		eta[35]  =   0.724462663076655e0;

		w[36]    =   0.021505319847731e0;
		w[37]    =   0.021505319847731e0;
		w[38]    =   0.021505319847731e0;
		w[39]    =   0.021505319847731e0;
		w[40]    =   0.021505319847731e0;
		w[41]    =   0.021505319847731e0;

		xi[36]   =   0.036843869875878e0;
		xi[37]   =   0.036843869875878e0;
		xi[38]   =   0.747556466051838e0;
		xi[39]   =   0.747556466051838e0;
		xi[40]   =   0.215599664072284e0;
		xi[41]   =   0.215599664072284e0;

		eta[36]  =   0.215599664072284e0;
		eta[37]  =   0.747556466051838e0;
		eta[38]  =   0.036843869875878e0;
		eta[39]  =   0.215599664072284e0;
		eta[40]  =   0.036843869875878e0;
		eta[41]  =   0.747556466051838e0;

		w[42]    =   0.007673942631049e0;
		w[43]    =   0.007673942631049e0;
		w[44]    =   0.007673942631049e0;
		w[45]    =   0.007673942631049e0;
		w[46]    =   0.007673942631049e0;
		w[47]    =   0.007673942631049e0;

		xi[42]   =   0.012459809331199e0;
		xi[43]   =   0.012459809331199e0;
		xi[44]   =   0.883964574092416e0;
		xi[45]   =   0.883964574092416e0;
		xi[46]   =   0.103575616576386e0;
		xi[47]   =   0.103575616576386e0;

		eta[42]  =   0.103575616576386e0;
		eta[43]  =   0.883964574092416e0;
		eta[44]  =   0.012459809331199e0;
		eta[45]  =   0.103575616576386e0;
		eta[46]  =   0.012459809331199e0;
		eta[47]  =   0.883964574092416e0;

		break;

	case 52: //p =16

		w[0]     =   0.046875697427642e0;
		xi[0]    =   0.333333333333333e0;
		eta[0]   =   0.333333333333333e0;

		w[1]     =   0.006405878578585e0;
		w[2]     =   0.006405878578585e0;
		w[3]     =   0.006405878578585e0;

		xi[1]    =   0.005238916103123e0;
		xi[2]    =   0.497380541948438e0;
		xi[3]    =   0.497380541948438e0;

		eta[1]   =   0.497380541948438e0;
		eta[2]   =   0.005238916103123e0;
		eta[3]   =   0.497380541948438e0;

		w[4]     =   0.041710296739387e0;
		w[5]     =   0.041710296739387e0;
		w[6]     =   0.041710296739387e0;

		xi[4]    =   0.173061122901295e0;
		xi[5]    =   0.413469438549352e0;
		xi[6]    =   0.413469438549352e0;

		eta[4]   =   0.413469438549352e0;
		eta[5]   =   0.173061122901295e0;
		eta[6]   =   0.413469438549352e0;

		w[7]     =   0.026891484250064e0;
		w[8]     =   0.026891484250064e0;
		w[9]    =   0.026891484250064e0;

		xi[7]    =   0.059082801866017e0;
		xi[8]    =   0.470458599066991e0;
		xi[9]   =   0.470458599066991e0;

		eta[7]   =   0.470458599066991e0;
		eta[8]   =   0.059082801866017e0;
		eta[9]  =   0.470458599066991e0;

		w[10]    =   0.042132522761650e0;
		w[11]    =   0.042132522761650e0;
		w[12]    =   0.042132522761650e0;

		xi[10]   =   0.518892500060958e0;
		xi[11]   =   0.240553749969521e0;
		xi[12]   =   0.240553749969521e0;

		eta[10]  =   0.240553749969521e0;
		eta[11]  =   0.518892500060958e0;
		eta[12]  =   0.240553749969521e0;

		w[13]    =   0.030000266842773e0;
		w[14]    =   0.030000266842773e0;
		w[15]    =   0.030000266842773e0;

		xi[13]   =   0.704068411554854e0;
		xi[14]   =   0.147965794222573e0;
		xi[15]   =   0.147965794222573e0;

		eta[13]  =   0.147965794222573e0;
		eta[14]  =   0.704068411554854e0;
		eta[15]  =   0.147965794222573e0;

		w[16]    =   0.014200098925024e0;
		w[17]    =   0.014200098925024e0;
		w[18]    =   0.014200098925024e0;

		xi[16]   =   0.849069624685052e0;
		xi[17]   =   0.075465187657474e0;
		xi[18]   =   0.075465187657474e0;

		eta[16]  =   0.075465187657474e0;
		eta[17]  =   0.849069624685052e0;
		eta[18]  =   0.075465187657474e0;

		w[19]    =   0.003582462351273e0;
		w[20]    =   0.003582462351273e0;
		w[21]    =   0.003582462351273e0;

		xi[19]   =   0.966807194753950e0;
		xi[20]   =   0.016596402623025e0;
		xi[21]   =   0.016596402623025e0;

		eta[19]  =   0.016596402623025e0;
		eta[20]  =   0.966807194753950e0;
		eta[21]  =   0.016596402623025e0;

		w[22]    =   0.032773147460627e0;
		w[23]    =   0.032773147460627e0;
		w[24]    =   0.032773147460627e0;
		w[25]    =   0.032773147460627e0;
		w[26]    =   0.032773147460627e0;
		w[27]    =   0.032773147460627e0;

		xi[22]   =   0.103575692245252e0;
		xi[23]   =   0.103575692245252e0;
		xi[24]   =   0.296555596579887e0;
		xi[25]   =   0.296555596579887e0;
		xi[26]   =   0.599868711174861e0;
		xi[27]   =   0.599868711174861e0;

		eta[22]  =   0.296555596579887e0;
		eta[23]  =   0.599868711174861e0;
		eta[24]  =   0.103575692245252e0;
		eta[25]  =   0.599868711174861e0;
		eta[26]  =   0.103575692245252e0;
		eta[27]  =   0.296555596579887e0;

		w[28]    =   0.015298306248441e0;
		w[29]    =   0.015298306248441e0;
		w[30]    =   0.015298306248441e0;
		w[31]    =   0.015298306248441e0;
		w[32]    =   0.015298306248441e0;
		w[33]    =   0.015298306248441e0;

		xi[28]   =   0.020083411655416e0;
		xi[29]   =   0.020083411655416e0;
		xi[30]   =   0.642193524941505e0;
		xi[31]   =   0.642193524941505e0;
		xi[32]   =   0.337723063403079e0;
		xi[33]   =   0.337723063403079e0;

		eta[28]  =   0.337723063403079e0;
		eta[29]  =   0.642193524941505e0;
		eta[30]  =   0.020083411655416e0;
		eta[31]  =   0.337723063403079e0;
		eta[32]  =   0.020083411655416e0;
		eta[33]  =   0.642193524941505e0;

		w[34]    =   0.002386244192839e0;
		w[35]    =   0.002386244192839e0;
		w[36]    =   0.002386244192839e0;
		w[37]    =   0.002386244192839e0;
		w[38]    =   0.002386244192839e0;
		w[39]    =   0.002386244192839e0;

		xi[34]   =  -0.004341002614139e0;
		xi[35]   =  -0.004341002614139e0;
		xi[36]   =   0.799592720971327e0;
		xi[37]   =   0.799592720971327e0;
		xi[38]   =   0.204748281642812e0;
		xi[39]   =   0.204748281642812e0;

		eta[34]  =   0.204748281642812e0;
		eta[35]  =   0.799592720971327e0;
		eta[36]  =  -0.004341002614139e0;
		eta[37]  =   0.204748281642812e0;
		eta[38]  =   0.799592720971327e0;
		eta[39]  =  -0.004341002614139e0;

		w[40]    =   0.019084792755899e0;
		w[41]    =   0.019084792755899e0;
		w[42]    =   0.019084792755899e0;
		w[43]    =   0.019084792755899e0;
		w[44]    =   0.019084792755899e0;
		w[45]    =   0.019084792755899e0;

		xi[40]   =   0.041941786468010e0;
		xi[41]   =   0.041941786468010e0;
		xi[42]   =   0.768699721401368e0;
		xi[43]   =   0.768699721401368e0;
		xi[44]   =   0.189358492130623e0;
		xi[45]   =   0.189358492130623e0;

		eta[40]  =   0.189358492130623e0;
		eta[41]  =   0.768699721401368e0;
		eta[42]  =   0.041941786468010e0;
		eta[43]  =   0.189358492130623e0;
		eta[44]  =   0.041941786468010e0;
		eta[45]  =   0.768699721401368e0;

		w[46]    =   0.006850054546542e0;
		w[47]    =   0.006850054546542e0;
		w[48]    =   0.006850054546542e0;
		w[49]    =   0.006850054546542e0;
		w[50]    =   0.006850054546542e0;
		w[51]    =   0.006850054546542e0;

		xi[46]   =   0.014317320230681e0;
		xi[47]   =   0.014317320230681e0;
		xi[48]   =   0.900399064086661e0;
		xi[49]   =   0.900399064086661e0;
		xi[50]   =   0.085283615682657e0;
		xi[51]   =   0.085283615682657e0;

		eta[46]  =   0.085283615682657e0;
		eta[47]  =   0.900399064086661e0;
		eta[48]  =   0.014317320230681e0;
		eta[49]  =   0.085283615682657e0;
		eta[50]  =   0.900399064086661e0;
		eta[51]  =   0.014317320230681e0;
		break;

	case 61://p =17

		w[0]    =   0.033437199290803;
		xi[0]   =   0.333333333333333;
		eta[0]  =   0.333333333333333;

		w[1]    =   0.005093415440507;
		xi[1]   =   0.005658918886452;
		eta[1]  =   0.497170540556774;

		w[2] =   0.005093415440507;
		xi[2] =   0.497170540556774;
		eta[2] =   0.005658918886452;

		w[3] =   0.005093415440507;
		xi[3] =   0.497170540556774;
		eta[3] =   0.497170540556774;

		w[4]    =  0.014670864527638;
		xi[4]   =  0.035647354750751;
		eta[4]  =  0.482176322624625;

		w[5]    =  0.014670864527638;
		xi[5]   =  0.482176322624625;
		eta[5]  =  0.035647354750751;

		w[6]    =  0.014670864527638;
		xi[6]   =  0.482176322624625;
		eta[6]  =  0.482176322624625;

		w[7]    =  0.024350878353672;
		xi[7]   =  0.099520061958437;
		eta[7]  =  0.450239969020782;

		w[8]    =  0.024350878353672;
		xi[8]   =  0.450239969020782;
		eta[8]  =  0.099520061958437;

		w[9]    =  0.024350878353672;
		xi[9]   =  0.450239969020782;
		eta[9]  =  0.450239969020782;

		w[10]    =  0.031107550868969;
		xi[10]   =  0.199467521245206;
		eta[10]  =  0.400266239377397;

		w[11]    =  0.031107550868969;
		xi[11]   =  0.400266239377397;
		eta[11]  =  0.199467521245206;

		w[12]    =  0.031107550868969;
		xi[12]   =  0.400266239377397;
		eta[12]  =  0.400266239377397;


		w[13]    =  0.031257111218620 ;
		xi[13]   =  0.495717464058095;
		eta[13]  =  0.252141267970953  ;

		w[14]    =  0.031257111218620 ;
		xi[14]   =  0.252141267970953 ;
		eta[14]  =  0.495717464058095;

		w[15]    =  0.031257111218620 ;
		xi[15]   =  0.252141267970953 ;
		eta[15]  =  0.252141267970953 ;

		w[16]    =  0.024815654339665;
		xi[16]   =  0.675905990683077;
		eta[16]  =  0.162047004658461;

		w[17]    =  0.024815654339665;
		xi[17]   =  0.162047004658461;
		eta[17]  =  0.675905990683077;

		w[18]    =  0.024815654339665;
		xi[18]   =  0.162047004658461;
		eta[18]  =  0.162047004658461;


		w[19]    =  0.014056073070557;
		xi[19]   =  0.848248235478508;
		eta[19]  =  0.075875882260746;

		w[20]    =  0.014056073070557;
		xi[20]   =  0.075875882260746;
		eta[20]  =  0.848248235478508;

		w[21]    =  0.014056073070557;
		xi[21]   =  0.075875882260746;
		eta[21]  =  0.075875882260746;

		w[22]    =  0.003194676173779;
		xi[22]   =  0.968690546064356 ;
		eta[22]  =  0.015654726967822;

		w[23]    =  0.003194676173779;
		xi[23]   =  0.015654726967822;
		eta[23]  =  0.968690546064356;

		w[24]    =  0.003194676173779;
		xi[24]   =  0.015654726967822;
		eta[24]  =  0.015654726967822;

		w[25]    =  0.008119655318993;
		xi[25]   =  0.010186928826919;
		eta[25]  =  0.334319867363658;

		w[26]    =  0.008119655318993;
		xi[26]   =  0.010186928826919;
		eta[26]  =  0.655493203809423;

		w[27]    =  0.008119655318993;
		xi[27]   =  0.334319867363658;
		eta[27]  =  0.010186928826919;

		w[28]    =  0.008119655318993;
		xi[28]   =  0.334319867363658;
		eta[28]  =  0.655493203809423;

		w[29]    =  0.008119655318993;
		xi[29]   =  0.655493203809423;
		eta[29]  =  0.010186928826919;

		w[30]    =  0.008119655318993;
		xi[30]   =  0.655493203809423;
		eta[30]  =  0.334319867363658;

		w[31]    =  0.026805742283163 ;
		xi[31]   =  0.135440871671036 ;
		eta[31]  =  0.292221537796944;

		w[32]    =  0.026805742283163;
		xi[32]   =  0.135440871671036;
		eta[32]  =  0.572337590532020;

		w[33]    =  0.026805742283163;
		xi[33]   =  0.292221537796944;
		eta[33]  =  0.135440871671036;

		w[34]    =  0.026805742283163;
		xi[34]   =  0.292221537796944;
		eta[34]  =  0.572337590532020;

		w[35]    =  0.026805742283163;
		xi[35]   =  0.572337590532020;
		eta[35]  =  0.135440871671036;

		w[36]    =  0.026805742283163;
		xi[36]   =  0.572337590532020 ;
		eta[36]  =  0.292221537796944;

		w[37]    =  0.018459993210822;
		xi[37]   =  0.054423924290583;
		eta[37]  =  0.319574885423190;

		w[38]    =  0.018459993210822;
		xi[38]   =  0.054423924290583;
		eta[38]  =  0.626001190286228;

		w[39]    =  0.018459993210822;
		xi[39]   =  0.319574885423190;
		eta[39]  =  0.054423924290583;

		w[40]    =  0.018459993210822;
		xi[40]   =  0.319574885423190;
		eta[40]  =  0.626001190286228;

		w[41]    =  0.018459993210822;
		xi[41]   =  0.626001190286228;
		eta[41]  =  0.054423924290583;

		w[42]    =  0.018459993210822;
		xi[42]   =  0.626001190286228;
		eta[42]  =  0.319574885423190;

		w[43]    =  0.008476868534328;
		xi[43]   =  0.012868560833637;
		eta[43]  =  0.190704224192292;

		w[44]    =  0.008476868534328;
		xi[44]   =  0.012868560833637;
		eta[44]  =  0.796427214974071;

		w[45]    =  0.008476868534328;
		xi[45]   =  0.190704224192292;
		eta[45]  =  0.012868560833637;

		w[46]    =  0.008476868534328;
		xi[46]   =  0.190704224192292;
		eta[46]  =  0.796427214974071;

		w[47]    =  0.008476868534328;
		xi[47]   =  0.796427214974071;
		eta[47]  =  0.012868560833637;

		w[48]    =  0.008476868534328;
		xi[48]   =  0.796427214974071;
		eta[48]  =  0.190704224192292;

		w[49]    =  0.018292796770025;
		xi[49]   =  0.067165782413524;
		eta[49]  =  0.180483211648746;

		w[50]    =  0.018292796770025;
		xi[50]   =  0.067165782413524;
		eta[50]  =  0.752351005937729;

		w[51]    =  0.018292796770025;
		xi[51]   =  0.180483211648746;
		eta[51]  =  0.067165782413524;

		w[52]    =  0.018292796770025;
		xi[52]   =  0.180483211648746;
		eta[52]  =  0.752351005937729;

		w[53]    =  0.018292796770025;
		xi[53]   =  0.752351005937729;
		eta[53]  =  0.067165782413524;

		w[54]    =  0.018292796770025;
		xi[54]   =  0.752351005937729;
		eta[54]  =  0.180483211648746;

		w[55]    =  0.006665632004165;
		xi[55]   =  0.014663182224828;
		eta[55]  =  0.080711313679564;

		w[56]    =  0.006665632004165;
		xi[56]   =  0.014663182224828;
		eta[56]  =  0.904625504095608;

		w[57]    =  0.006665632004165;
		xi[57]   =  0.080711313679564;
		eta[57]  =  0.014663182224828;

		w[58]    =  0.006665632004165;
		xi[58]   =  0.080711313679564;
		eta[58]  =  0.904625504095608;

		w[59]    =  0.006665632004165;
		xi[59]   =  0.904625504095608;
		eta[59]  =  0.014663182224828;

		w[60]    =  0.006665632004165;
		xi[60]   =  0.904625504095608;
		eta[60]  =  0.080711313679564;

		break;

	case 70: // p =18

		w[0]    =  0.030809939937647;
		xi[0]   =  0.333333333333333  ;
		eta[0]  =  0.333333333333333  ;


		w[1]    =  0.009072436679404;
		xi[1]   =  0.013310382738157;
		eta[1]  =  0.493344808630921   ;

		w[2]    =  0.009072436679404;
		xi[2]   =  0.493344808630921  ;
		eta[2]  =  0.013310382738157;

		w[3]    =  0.009072436679404;
		xi[3]   =  0.493344808630921  ;
		eta[3]  =  0.493344808630921  ;

		w[4]    =  0.018761316939594;
		xi[4]   =  0.061578811516086;
		eta[4]  =  0.469210594241957;

		w[5]    =  0.018761316939594;
		xi[5]   =  0.469210594241957;
		eta[5]  =  0.061578811516086;

		w[6]    =  0.018761316939594;
		xi[6]   =  0.469210594241957;
		eta[6]  =  0.469210594241957;

		w[7]    =  0.019441097985477;
		xi[7]   =  0.127437208225989;
		eta[7]  =  0.436281395887006;

		w[8]    =  0.019441097985477;
		xi[8]   =  0.436281395887006;
		eta[8]  =  0.127437208225989;

		w[9]    =  0.019441097985477;
		xi[9]   =  0.436281395887006;
		eta[9]  =  0.436281395887006;

		w[10]    =  0.027753948610810;
		xi[10]   =  0.210307658653168;
		eta[10]  =  0.394846170673416;

		w[11]    =  0.027753948610810;
		xi[11]   =  0.394846170673416;
		eta[11]  =  0.210307658653168;

		w[12]    =  0.027753948610810;
		xi[12]   =  0.394846170673416;
		eta[12]  =  0.394846170673416;

		w[13]    =  0.032256225351457;
		xi[13]   =  0.500410862393686;
		eta[13]  =  0.249794568803157  ;

		w[14]    =  0.032256225351457;
		xi[14]   =  0.249794568803157  ;
		eta[14]  =  0.500410862393686;

		w[15]    =  0.032256225351457;
		xi[15]   =  0.249794568803157  ;
		eta[15]  =  0.249794568803157  ;



		w[16]    =  0.025074032616922;
		xi[16]   =  0.677135612512315;
		eta[16]  =  0.161432193743843;

		w[17]    =  0.025074032616922;
		xi[17]   =  0.161432193743843;
		eta[17]  =  0.677135612512315;

		w[18]    =  0.025074032616922;
		xi[18]   =  0.161432193743843;
		eta[18]  =  0.161432193743843;



		w[19]    =  0.015271927971832;
		xi[19]   =  0.846803545029257;
		eta[19]  =  0.076598227485371   ;

		w[20]    =  0.015271927971832;
		xi[20]   =  0.076598227485371  ;
		eta[20]  =  0.846803545029257;

		w[21]    =  0.015271927971832;
		xi[21]   =  0.076598227485371  ;
		eta[21]  =  0.076598227485371  ;

		w[22]    =  0.006793922022963;
		xi[22]   =  0.951495121293100;
		eta[22]  =  0.024252439353450;

		w[23]    =  0.006793922022963;
		xi[23]   =  0.024252439353450;
		eta[23]  =  0.951495121293100;

		w[24]    =  0.006793922022963;
		xi[24]   =  0.024252439353450;
		eta[24]  =  0.024252439353450;

		w[25]    =  -0.002223098729920;
		xi[25]   =  0.913707265566071;
		eta[25]  =  0.043146367216965;

		w[26]    =  -0.002223098729920;
		xi[26]   =  0.043146367216965;
		eta[26]  =  0.913707265566071;

		w[27]    =  -0.002223098729920;
		xi[27]   =  0.043146367216965;
		eta[27]  =  0.043146367216965;

		w[28]    =  0.006331914076406;
		xi[28]   =  0.008430536202420;
		eta[28]  =  0.358911494940944;

		w[29]    =  0.006331914076406;
		xi[29]   =  0.008430536202420;
		eta[29]  =  0.632657968856636;

		w[30]    =  0.006331914076406;
		xi[30]   =  0.358911494940944;
		eta[30]  =  0.008430536202420;

		w[31]    =  0.006331914076406;
		xi[31]   =  0.358911494940944;
		eta[31]  =  0.632657968856636;

		w[32]    =  0.006331914076406;
		xi[32]   =  0.632657968856636;
		eta[32]  =  0.008430536202420;

		w[33]    =  0.006331914076406;
		xi[33]   =  0.632657968856636;
		eta[33]  =  0.358911494940944;

		w[34]    =  0.027257538049138;
		xi[34]   =  0.131186551737188;
		eta[34]  =  0.294402476751957;

		w[35]    =  0.027257538049138;
		xi[35]   =  0.131186551737188;
		eta[35]  =  0.574410971510855;

		w[36]    =  0.027257538049138;
		xi[36]   =  0.294402476751957;
		eta[36]  =  0.131186551737188;

		w[37]    =  0.027257538049138;
		xi[37]   =  0.294402476751957;
		eta[37]  =  0.574410971510855;

		w[38]    =  0.027257538049138;
		xi[38]   =  0.574410971510855;
		eta[38]  =  0.131186551737188;

		w[39]    =  0.027257538049138;
		xi[39]   =  0.574410971510855;
		eta[39]  =  0.294402476751957;


		w[40]    =  0.017676785649465;
		xi[40]   =  0.050203151565675;
		eta[40]  =  0.325017801641814;

		w[41]    =  0.017676785649465;
		xi[41]   =  0.050203151565675;
		eta[41]  =  0.624779046792512;

		w[42]    =  0.017676785649465;
		xi[42]   =  0.325017801641814;
		eta[42]  =  0.050203151565675;

		w[43]    =  0.017676785649465;
		xi[43]   =  0.325017801641814;
		eta[43]  =  0.624779046792512;

		w[44]    =  0.017676785649465;
		xi[44]   =  0.624779046792512;
		eta[44]  =  0.050203151565675;

		w[45]    =  0.017676785649465;
		xi[45]   =  0.624779046792512;
		eta[45]  =  0.325017801641814;



		w[46]    =  0.018379484638070;
		xi[46]   =  0.066329263810916;
		eta[46]  =  0.184737559666046;

		w[47]    =  0.018379484638070;
		xi[47]   =  0.066329263810916;
		eta[47]  =  0.748933176523037;

		w[48]    =  0.018379484638070;
		xi[48]   =  0.184737559666046;
		eta[48]  =  0.066329263810916;

		w[49]    =  0.018379484638070;
		xi[49]   =  0.184737559666046;
		eta[49]  =  0.748933176523037;

		w[50]    =  0.018379484638070;
		xi[50]   =  0.748933176523037;
		eta[50]  =  0.066329263810916;

		w[51]    =  0.018379484638070;
		xi[51]   =  0.748933176523037;
		eta[51]  =  0.184737559666046;

		w[52]    =  0.008104732808192;
		xi[52]   =  0.011996194566236;
		eta[52]  =  0.218796800013321  ;

		w[53]    =  0.008104732808192;
		xi[53]   =  0.011996194566236;
		eta[53]  =  0.769207005420443;

		w[54]    =  0.008104732808192;
		xi[54]   =  0.218796800013321  ;
		eta[54]  =  0.011996194566236;

		w[55]    =  0.008104732808192;
		xi[55]   =  0.218796800013321  ;
		eta[55]  =  0.769207005420443;

		w[56]    =  0.008104732808192;
		xi[56]   =  0.769207005420443;
		eta[56]  =  0.011996194566236;

		w[57]    =  0.008104732808192;
		xi[57]   =  0.769207005420443;
		eta[57]  =  0.218796800013321  ;

		w[58]    =  0.007634129070725;
		xi[58]   =  0.014858100590125;
		eta[58]  =  0.101179597136408;

		w[59]    =  0.007634129070725;
		xi[59]   =  0.014858100590125;
		eta[59]  =  0.883962302273467;

		w[60]    =  0.007634129070725;
		xi[60]   =  0.101179597136408;
		eta[60]  =  0.014858100590125;

		w[61]    =  0.007634129070725;
		xi[61]   =  0.101179597136408;
		eta[61]  =  0.883962302273467;

		w[62]    =  0.007634129070725;
		xi[62]   =  0.883962302273467;
		eta[62]  =  0.014858100590125;

		w[63]    =  0.007634129070725;
		xi[63]   =  0.883962302273467;
		eta[63]  =  0.101179597136408;

		w[64]    =  0.000046187660794;
		xi[64]   =  -0.035222015287949;
		eta[64]  =  0.020874755282586;

		w[65]    =  0.000046187660794;
		xi[65]   =  -0.035222015287949;
		eta[65]  =  1.014347260005363;

		w[66]    =  0.000046187660794;
		xi[66]   =  0.020874755282586;
		eta[66]  =  -0.035222015287949;

		w[67]    =  0.000046187660794;
		xi[67]   =  0.020874755282586;
		eta[67]  =  1.014347260005363;

		w[68]    =  0.000046187660794;
		xi[68]   =  1.014347260005363;
		eta[68]  =  -0.035222015287949;

		w[69]    =  0.000046187660794;
		xi[69]   =  1.014347260005363;
		eta[69]  =  0.020874755282586;

		break;

	case 73: // p=19

		w[0]    =  0.032906331388919;
		xi[0]   =  0.333333333333333;
		eta[0]  =  0.333333333333333;

		w[1]    =  0.010330731891272;
		xi[1]   =  0.020780025853987;
		eta[1]  =  0.489609987073006;

		w[2]    =  0.010330731891272;
		xi[2]   =  0.489609987073006;
		eta[2]  =  0.020780025853987;

		w[3]    =  0.010330731891272;
		xi[3]   =  0.489609987073006;
		eta[3]  =  0.489609987073006;



		w[4]    =  0.022387247263016;
		xi[4]   =  0.090926214604215;
		eta[4]  =  0.454536892697893;

		w[5]    =  0.022387247263016;
		xi[5]   =  0.454536892697893;
		eta[5]  =  0.090926214604215;

		w[6]    =  0.022387247263016;
		xi[6]   =  0.454536892697893;
		eta[6]  =  0.454536892697893;

		w[7]    =  0.030266125869468;
		xi[7]   =  0.197166638701138;
		eta[7]  =  0.401416680649431;

		w[8]    =  0.030266125869468;
		xi[8]   =  0.401416680649431;
		eta[8]  =  0.197166638701138;

		w[9]    =  0.030266125869468;
		xi[9]   =  0.401416680649431;
		eta[9]  =  0.401416680649431;

		w[10]    =  0.030490967802198;
		xi[10]   =  0.488896691193805;
		eta[10]  =  0.255551654403098;

		w[11]    =  0.030490967802198;
		xi[11]   =  0.255551654403098;
		eta[11]  =  0.488896691193805;

		w[12]    =  0.030490967802198;
		xi[12]   =  0.255551654403098;
		eta[12]  =  0.255551654403098;


		w[13]    =  0.024159212741641;
		xi[13]   =  0.645844115695741;
		eta[13]  =  0.177077942152130;

		w[14]    =  0.024159212741641;
		xi[14]   =  0.177077942152130;
		eta[14]  =  0.645844115695741;

		w[15]    =  0.024159212741641;
		xi[15]   =  0.177077942152130;
		eta[15]  =  0.177077942152130;

		w[16]    =  0.016050803586801;
		xi[16]   =  0.779877893544096;
		eta[16]  =  0.110061053227952;

		w[17]    =  0.016050803586801;
		xi[17]   =  0.110061053227952;
		eta[17]  =  0.779877893544096;

		w[18]    =  0.016050803586801;
		xi[18]   =  0.110061053227952;
		eta[18]  =  0.110061053227952;


		w[19]    =  0.008084580261784;
		xi[19]   =  0.888942751496321;
		eta[19]  =  0.055528624251840;

		w[20]    =  0.008084580261784;
		xi[20]   =  0.055528624251840;
		eta[20]  =  0.888942751496321;

		w[21]    =  0.008084580261748;
		xi[21]   =  0.055528624251840;
		eta[21]  =  0.055528624251840;


		w[22]    =  0.002079362027485;
		xi[22]   =  0.974756272445543;
		eta[22]  =  0.012621863777229;

		w[23]    =  0.002079362027485;
		xi[23]   =  0.012621863777229;
		eta[23]  =  0.974756272445543;

		w[24]    =  0.002079362027485;
		xi[24]   =  0.012621863777229;
		eta[24]  =  0.012621863777229;


		w[25]    =  0.003884876904981;
		xi[25]   =  0.003611417848412;
		eta[25]  =  0.395754787356943;

		w[26]    =  0.003884876904981;
		xi[26]   =  0.003611417848412;
		eta[26]  =  0.600633794794645;

		w[27]    =  0.003884876904981;
		xi[27]   =  0.395754787356943;
		eta[27]  =  0.003611417848412;

		w[28]    =  0.003884876904981;
		xi[28]   =  0.395754787356943;
		eta[28]  =  0.600633794794645;

		w[29]    =  0.003884876904981;
		xi[29]   =  0.600633794794645;
		eta[29]  =  0.003611417848412;

		w[30]    =  0.003884876904981;
		xi[30]   =  0.600633794794645;
		eta[30]  =  0.395754787356943;

		w[31]    =  0.025574160612022;
		xi[31]   =  0.134466754530780;
		eta[31]  =  0.307929983880436;

		w[32]    =  0.025574160612022;
		xi[32]   =  0.134466754530780;
		eta[32]  =  0.557603261588784;

		w[33]    =  0.025574160612022;
		xi[33]   =  0.307929983880436;
		eta[33]  =  0.134466754530780;

		w[34]    =  0.025574160612022;
		xi[34]   =  0.307929983880436;
		eta[34]  =  0.557603261588784;

		w[35]    =  0.025574160612022;
		xi[35]   =  0.557603261588784;
		eta[35]  =  0.134466754530780;

		w[36]    =  0.025574160612022;
		xi[36]   =  0.557603261588784;
		eta[36]  =  0.307929983880436;


		w[37]    =  0.008880903573338;
		xi[37]   =  0.014446025776115;
		eta[37]  =  0.264566948406520;

		w[38]    =  0.008880903573338;
		xi[38]   =  0.014446025776115;
		eta[38]  =  0.720987025817365;

		w[39]    =  0.008880903573338;
		xi[39]   =  0.264566948406520;
		eta[39]  =  0.014446025776115;

		w[40]    =  0.008880903573338;
		xi[40]   =  0.264566948406520;
		eta[40]  =  0.720987025817365;

		w[41]    =  0.008880903573338;
		xi[41]   =  0.720987025817365;
		eta[41]  =  0.014446025776115;

		w[42]    =  0.008880903573338;
		xi[42]   =  0.720987025817365;
		eta[42]  =  0.264566948406520;

		w[43]    =  0.016124546761731;
		xi[43]   =  0.046933578838178;
		eta[43]  =  0.358539352205951;

		w[44]    =  0.016124546761731;
		xi[44]   =  0.046933578838178;
		eta[44]  =  0.594527068955871;

		w[45]    =  0.016124546761731;
		xi[45]   =  0.358539352205951;
		eta[45]  =  0.046933578838178;

		w[46]    =  0.016124546761731;
		xi[46]   =  0.358539352205951;
		eta[46]  =  0.594527068955871;

		w[47]    =  0.016124546761731;
		xi[47]   =  0.594527068955871;
		eta[47]  =  0.046933578838178;

		w[48]    =  0.016124546761731;
		xi[48]   =  0.594527068955871;
		eta[48]  =  0.358539352205951;

		w[49]    =  0.002491941817491;
		xi[49]   =  0.002861120350567;
		eta[49]  =  0.157807405968595;

		w[50]    =  0.002491941817491;
		xi[50]   =  0.002861120350567;
		eta[50]  =  0.839331473680839;

		w[51]    =  0.002491941817491;
		xi[51]   =  0.157807405968595;
		eta[51]  =  0.002861120350567;

		w[52]    =  0.002491941817491;
		xi[52]   =  0.157807405968595;
		eta[52]  =  0.839331473680839;

		w[53]    =  0.002491941817491;
		xi[53]   =  0.839331473680839;
		eta[53]  =  0.002861120350567;

		w[54]    =  0.002491941817491;
		xi[54]   =  0.839331473680839;
		eta[54]  =  0.157807405968595;



		w[55]    =  0.018242840118951;
		xi[55]   =  0.223861424097916;
		eta[55]  =  0.075050596975911;

		w[56]    =  0.018242840118951;
		xi[56]   =  0.223861424097916;
		eta[56]  =  0.701087978926173;

		w[57]    =  0.018242840118951;
		xi[57]   =  0.075050596975911;
		eta[57]  =  0.223861424097916;

		w[58]    =  0.018242840118951;
		xi[58]   =  0.075050596975911;
		eta[58]  =  0.701087978926173;

		w[59]    =  0.018242840118951;
		xi[59]   =  0.701087978926173;
		eta[59]  =  0.223861424097916;

		w[60]    =  0.018242840118951;
		xi[60]   =  0.701087978926173;
		eta[60]  =  0.075050596975911;


		w[61]    =  0.010258563736199;
		xi[61]   =  0.034647074816760;
		eta[61]  =  0.142421601113383;

		w[62]    =  0.010258563736199;
		xi[62]   =  0.034647074816760;
		eta[62]  =  0.822931324069857;

		w[63]    =  0.010258563736199;
		xi[63]   =  0.142421601113383;
		eta[63]  =  0.034647074816760;

		w[64]    =  0.010258563736199;
		xi[64]   =  0.142421601113383;
		eta[64]  =  0.822931324069857;

		w[65]    =  0.010258563736199;
		xi[65]   =  0.822931324069857;
		eta[65]  =  0.034647074816760;

		w[66]    =  0.010258563736199;
		xi[66]   =  0.822931324069857;
		eta[66]  =  0.142421601113383;


		w[67]    =  0.003799928855302;
		xi[67]   =  0.010161119296278;
		eta[67]  =  0.065494628082938;

		w[68]    =  0.003799928855302;
		xi[68]   =  0.010161119296278;
		eta[68]  =  0.924344252620784;

		w[69]    =  0.003799928855302;
		xi[69]   =  0.065494628082938;
		eta[69]  =  0.010161119296278;

		w[70]    =  0.003799928855302;
		xi[70]   =  0.065494628082938 ;
		eta[70]  =  0.924344252620784;

		w[71]    =  0.003799928855302;
		xi[71]   =  0.924344252620784;
		eta[71]  =  0.010161119296278;

		w[72]    =  0.003799928855302;
		xi[72]   =  0.924344252620784;
		eta[72]  =  0.065494628082938;

		break;
	case 79:// p =20

		w[0]    =  0.033057055541624;
		xi[0]   =  0.333333333333333;
		eta[0]  =  0.333333333333333;

		w[1]    =  0.000867019185663;
		xi[1]   =  -0.001900928704400;
		eta[1]  =  0.500950464352200;

		w[2]    =  0.000867019185663;
		xi[2]   =  0.500950464352200;
		eta[2]  =  -0.001900928704400;

		w[3]    =  0.000867019185663;
		xi[3]   =  0.500950464352200;
		eta[3]  =  0.500950464352200;

		w[4]    =  0.011660052716448;
		xi[4]   =  0.023574084130543;
		eta[4]  =  0.488212957934729;

		w[5]    =  0.011660052716448;
		xi[5]   =  0.488212957934729;
		eta[5]  =  0.023574084130543;

		w[6]    =  0.011660052716448;
		xi[6]   =  0.488212957934729;
		eta[6]  =  0.488212957934729;

		w[7]    =  0.022876936356421;
		xi[7]   =  0.089726636099435;
		eta[7]  =  0.455136681950283;

		w[8]    =  0.022876936356421;
		xi[8]   =  0.455136681950283;
		eta[8]  =  0.089726636099435;

		w[9]    =  0.022876936356421;
		xi[9]   =  0.455136681950283;
		eta[9]  =  0.455136681950283;


		w[10]    =  0.030448982673938;
		xi[10]   =  0.196007481363421;
		eta[10]  =  0.401996259318289;

		w[11]    =  0.030448982673938;
		xi[11]   =  0.401996259318289;
		eta[11]  =  0.196007481363421;

		w[12]    =  0.030448982673938;
		xi[12]   =  0.401996259318289;
		eta[12]  =  0.401996259318289;

		w[13]    =  0.030624891725355;
		xi[13]   =  0.488214180481157;
		eta[13]  =  0.255892909759421;

		w[14]    =  0.030624891725355;
		xi[14]   =  0.255892909759421;
		eta[14]  =  0.488214180481157;

		w[15]    =  0.030624891725355;
		xi[15]   =  0.255892909759421;
		eta[15]  =  0.255892909759421;

		w[16]    =  0.024368057676800;
		xi[16]   =  0.647023488009788;
		eta[16]  =  0.176488255995106;

		w[17]    =  0.024368057676800;
		xi[17]   =  0.176488255995106;
		eta[17]  =  0.647023488009788;

		w[18]    =  0.024368057676800;
		xi[18]   =  0.176488255995106;
		eta[18]  =  0.176488255995106;


		w[19]    =  0.015997432032024;
		xi[19]   =  0.791658289326483;
		eta[19]  =  0.104170855336758;

		w[20]    =  0.015997432032024;
		xi[20]   =  0.104170855336758;
		eta[20]  =  0.791658289326483;

		w[21]    =  0.015997432032024;
		xi[21]   =  0.104170855336758;
		eta[21]  =  0.104170855336758;

		w[22]    =  0.007698301815602;
		xi[22]   =  0.893862072318140;
		eta[22]  =  0.053068963840930;

		w[23]    =  0.007698301815602;
		xi[23]   =  0.053068963840930;
		eta[23]  =  0.893862072318140;

		w[24]    =  0.007698301815602;
		xi[24]   =  0.053068963840930;
		eta[24]  =  0.053068963840930;


		w[25]    =  -0.000632060497488;
		xi[25]   =  0.916762569607942;
		eta[25]  =  0.041618715196029;

		w[26]    =  -0.000632060497488;
		xi[26]   =  0.041618715196029;
		eta[26]  =  0.916762569607942;

		w[27]    =  -0.000632060497488;
		xi[27]   =  0.041618715196029;
		eta[27]  =  0.041618715196029;


		w[28]    =  0.001751134301193;
		xi[28]   =  0.976836157186356;
		eta[28]  =  0.011581921406822;

		w[29]    =  0.001751134301193;
		xi[29]   =  0.011581921406822;
		eta[29]  =  0.976836157186356;

		w[30]    =  0.001751134301193;
		xi[30]   =  0.011581921406822;
		eta[30]  =  0.011581921406822;


		w[31]    =  0.016465839189576;
		xi[31]   =  0.048741583664839;
		eta[31]  =  0.344855770229001;

		w[32]    =  0.016465839189576;
		xi[32]   =  0.048741583664839;
		eta[32]  =  0.606402646106160;

		w[33]    =  0.016465839189576;
		xi[33]   =  0.344855770229001;
		eta[33]  =  0.048741583664839;

		w[34]    =  0.016465839189576;
		xi[34]   =  0.344855770229001;
		eta[34]  =  0.606402646106160;

		w[35]    =  0.016465839189576;
		xi[35]   =  0.606402646106160;
		eta[35]  =  0.048741583664839;

		w[36]    =  0.016465839189576;
		xi[36]   =  0.606402646106160;
		eta[36]  =  0.344855770229001;

		w[37]    =  0.004839033540485;
		xi[37]   =  0.006314115948605;
		eta[37]  =  0.377843269594854;

		w[38]    =  0.004839033540485;
		xi[38]   =  0.006314115948605;
		eta[38]  =  0.615842614456541;

		w[39]    =  0.004839033540485;
		xi[39]   =  0.377843269594854;
		eta[39]  =  0.006314115948605;

		w[40]    =  0.004839033540485;
		xi[40]   =  0.377843269594854;
		eta[40]  =  0.615842614456541;

		w[41]    =  0.004839033540485;
		xi[41]   =  0.615842614456541;
		eta[41]  =  0.006314115948605;

		w[42]    =  0.004839033540485;
		xi[42]   =  0.615842614456541;
		eta[42]  =  0.377843269594854;

		w[43]    =  0.025804906534650;
		xi[43]   =  0.134316520547348;
		eta[43]  =  0.306635479062357;

		w[44]    =  0.025804906534650;
		xi[44]   =  0.134316520547348;
		eta[44]  =  0.559048000390295;

		w[45]    =  0.025804906534650;
		xi[45]   =  0.306635479062357;
		eta[45]  =  0.134316520547348;

		w[46]    =  0.025804906534650;
		xi[46]   =  0.306635479062357;
		eta[46]  =  0.559048000390295;

		w[47]    =  0.025804906534650;
		xi[47]   =  0.559048000390295;
		eta[47]  =  0.134316520547348;

		w[48]    =  0.025804906534650;
		xi[48]   =  0.559048000390295;
		eta[48]  =  0.306635479062357;


		w[49]    =  0.008471091054441;
		xi[49]   =  0.013973893962392;
		eta[49]  =  0.249419362774742;

		w[50]    =  0.008471091054441;
		xi[50]   =  0.013973893962392;
		eta[50]  =  0.736606743262866;

		w[51]    =  0.008471091054441;
		xi[51]   =  0.249419362774742;
		eta[51]  =  0.013973893962392;

		w[52]    =  0.008471091054441;
		xi[52]   =  0.249419362774742;
		eta[52]  =  0.736606743262866;

		w[53]    =  0.008471091054441;
		xi[53]   =  0.736606743262866;
		eta[53]  =  0.013973893962392;

		w[54]    =  0.008471091054441;
		xi[54]   =  0.736606743262866;
		eta[54]  =  0.249419362774742;

		w[55]    =  0.018354914106280;
		xi[55]   =  0.075549132909764;
		eta[55]  =  0.212775724802802;

		w[56]    =  0.018354914106280;
		xi[56]   =  0.075549132909764;
		eta[56]  =  0.711675142287434;

		w[57]    =  0.018354914106280;
		xi[57]   =  0.212775724802802;
		eta[57]  =  0.075549132909764;

		w[58]    =  0.018354914106280;
		xi[58]   =  0.212775724802802;
		eta[58]  =  0.711675142287434;

		w[59]    =  0.018354914106280;
		xi[59]   =  0.711675142287434;
		eta[59]  =  0.075549132909764;

		w[60]    =  0.018354914106280;
		xi[60]   =  0.711675142287434;
		eta[60]  =  0.212775724802802;

		w[61]    =  0.000704404677908;
		xi[61]   =  -0.008368153208227;
		eta[61]  =  0.146965436053239;

		w[62]    =  0.000704404677908;
		xi[62]   =  -0.008368153208227;
		eta[62]  =  0.861402717154987;

		w[63]    =  0.000704404677908;
		xi[63]   =  0.146965436053239;
		eta[63]  =  -0.008368153208227;

		w[64]    =  0.000704404677908;
		xi[64]   =  0.146965436053239;
		eta[64]  =  0.861402717154987;

		w[65]    =  0.000704404677908;
		xi[65]   =  0.861402717154987;
		eta[65]  =  -0.008368153208227;

		w[66]    =  0.000704404677908;
		xi[66]   =  0.861402717154987;
		eta[66]  =  0.146965436053239;

		w[67]    =  0.010112684927462;
		xi[67]   =  0.026686063258714;
		eta[67]  =  0.137726978828923;

		w[68]    =  0.010112684927462;
		xi[68]   =  0.026686063258714;
		eta[68]  =  0.835586957912363;

		w[69]    =  0.010112684927462;
		xi[69]   =  0.137726978828923;
		eta[69]  =  0.026686063258714;

		w[70]    =  0.010112684927462;
		xi[70]   =  0.137726978828923;
		eta[70]  =  0.835586957912363;

		w[71]    =  0.010112684927462;
		xi[71]   =  0.835586957912363;
		eta[71]  =  0.026686063258714;

		w[72]    =  0.010112684927462;
		xi[72]   =  0.835586957912363;
		eta[72]  =  0.137726978828923;

		w[73]    =  0.003573909385950;
		xi[73]   =  0.010547719294141;
		eta[73]  =  0.059696109149007;

		w[74]    =  0.003573909385950;
		xi[74]   =  0.010547719294141;
		eta[74]  =  0.929756171556853;

		w[75]    =  0.003573909385950;
		xi[75]   =  0.059696109149007;
		eta[75]  =  0.010547719294141;

		w[76]    =  0.003573909385950;
		xi[76]   =  0.059696109149007;
		eta[76]  =  0.929756171556853;


		w[77]    =  0.003573909385950;
		xi[77]   =  0.929756171556853;
		eta[77]  =  0.010547719294141;

		w[78]    =  0.003573909385950;
		xi[78]   =  0.929756171556853;
		eta[78]  =  0.059696109149007;

		break;
}
}
