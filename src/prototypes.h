//=================================================================
//  File Prototypes.cpp
//  Created:  March 17, 2017
//  Author: Sathwik Bharadwaj
//  Modified: 
//  Last change: 
//  This file contains the definitation of all functions used in our codes. 
//  Prototype.h header function wil be included in all files. 
//=================================================================
#ifndef PROTOTYPESGAN_H
#define PROTOTYPESGAN_H
//======================================================

// ==================
//   PROTOTYPES  
// ==================
// In Solve.cpp
// ==================
PetscErrorCode solution(global_matrices & gmat, data & dat, Vec psi);

// ==================
// In apply_BC.cpp
//====================
PetscErrorCode apply_bc(global_matrices & gmat , data & dat);

// ==================
//In make_global.cpp
//==================
PetscErrorCode make_global(global_matrices& gmat, data& dat);

//==================

PetscErrorCode petsc_diagonalizer(global_matrices & gmat , data & dat);

// ==================
// From input_reader.cpp
// ==================

void get_line(FILE* inputfile,char* buffer);
void input_reader(std::ifstream &in, data &dat);
PetscErrorCode mesh_input(data &dat);

// ==================
// From utiltrig.cpp
// ===============

void shape(double phi[], double xi_val, double eta_val, double Transform[][6], double shapeCoeff[][18], int& ndeg1);
void deriv1(double derivx[], double derivy[], double xi_val, double eta_val, double Transform[][6], double shapeCoeff[][18], int& ndeg1);
void deriv2(double derivxx[], double derivxy[], double derivyy[], double xi_val, double eta_val, double Transform[][6], double shapeCoeff[][18]);
void generateTransformationMatrix(double xp[3], double yp[3], double T[][6]);
void generateShapeCoeffMatrix(double B[][18], int& ndeg1);
void locelem(double x, double y, data &dat, int &iel);
void trigauss(int& ngaus, double *x, double *y, double *wt);

// ==================
// From readin.cpp
// ==================

template <typename T> T readin(std::ifstream& file, const char* description, int ndebug);
void ignore_until ( std::ifstream& file, char c );
int determine_current_line_number( std::ifstream& file);

//======================================================
#endif
