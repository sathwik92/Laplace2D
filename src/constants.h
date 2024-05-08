//=================================================================
//  File Constants.h
//=================================================================

//=================================================================
#ifndef CONSTANTS_H 
#define CONSTANTS_H
//=================================================================

/* linux has PI defined in /usr/include/math.h */

#ifndef   PI
#define PI  3.14159265358979323846
#endif

//////////////////////////////////////////////////////
//
// Note: the complex constant ei=(0.0,1.0) is already
//       defined in the dmat library: in constants.h
//
//////////////////////////////////////////////////////

//====================================================
//    Physical constants: 
// m_o*c^2 (eV), hbar_c (ev*cm), xscale:1 Angstrom
//====================================================

#define xm0     0.5109990615E+6
#define hbarc   1.9732705359E-5
#define xscale  1.0E-8
#define vac_epsilon 8.854187817e-12 // in Farad/meter

// LRR's C_scale_energy, del = 2*cscale, 
//       xl0=length scale:1 Angstrom

#define cscale 	3.809984038754390
#define del 	7.619968077508781
#define xl0     1.0e-8
#define Kboltz  (1.0e0/11604.5e0) // Boltzmann's constant

//=====================================================
// B-field parameters:
//
// Square of Landau orbit r^2/B_0, hbar*omega_B
// r2b1 is divided by b0 (in Tesla) to get Landau R_o**2 
//        in units of cm^2 
// hwb is multiplied by b0 (in Tesla) to get hbar*omega_b 
//        in units of eV 
//=====================================================

#define r2b1 	6.58212224061E-12
#define hwb	1.1576764E-4

// Some numerical constants appearing k.P Hamiltonian
//
#define sqi2	0.70710678118654752440
#define sqi3	0.57735026918962576451
#define sqi6	0.40824829046386301637
#define sq2	1.41421356237309504880
#define sq3	1.73205080756887729353

#define sq6   	2.449489742783178098197
#define sq32	1.224744871391589049099
#define sq23	0.816496580927726032732

//==================================================================
//==================================================================
//==================================================================
#endif


