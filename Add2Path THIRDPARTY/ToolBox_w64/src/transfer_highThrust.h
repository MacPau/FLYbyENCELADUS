/*******************************************************************************
*                          transfer_highThrust.h                               *
*                   Functions for high-thrust transfers                        *
*                                                                              *
*                              Space Toolbox                                   *
*******************************************************************************/

#ifndef TRANSFER_HIGHTHRUST_H
#define TRANSFER_HIGHTHRUST_H

/* Include libraries
 * ----------------- */
#ifdef DEBUG
#include <stdio.h>
#endif /* DEBUG */

#include <math.h>
#include "mathUtils.h"

/* Mex compiler switch
 * ------------------- */
#ifdef MEXCOMPILE
#define printf mexPrintf
#define malloc mxMalloc
#define calloc mxCalloc
#define free mxFree
#endif

#ifndef MEXCOMPILE
#include <stdio.h>
#endif /* MEXCOMPILE */

/* Constants
 * ---------- */
#define LAMBERTMR_TOL (1e-14)
#define LAMBERTMR_NITERMAX 2000

/* Prototypes
 * ---------- */

int lambertMR(const double ri[3], const double rf[3], double tof, double mu, int orbittype, int nrev, int ncase, int optionslmr,
	double *a_out, double *p_out, double *e_out, int *error_out, double vi_out[3], double vf_out[3], double *tpar_out, double *theta_out);
/* lambertMR - Multiple revolution Lambert problem solver.
* 
* PROTOTYPE:
*	int lambertMR(const double ri[3], const double rf[3], double tof,
*		double mu, int orbittype, int nrev, int ncase, int optionslmr,
*		double *a_out, double *p_out, double *e_out, int *error_out,
*		double vi_out[3], double vf_out[3], double *tpar_out, double *theta_out)
*
* DESCRIPTION:
*	Lambert's problem solver for all possible transfers:
*		1 - zero-revolution (for all possible types of orbits: circles,
*			ellipses, parabolas and hyperbolas)
*		2 - multirevolution case
*		3 - inversion of the motion
* 
*	1- ZERO-REVOLUTION LAMBERT'S PROBLEM
*
*	For the solution of Lambert's problem with number of revolution = 0 the
*	subroutine by Chris D'Souza is included here.
*	This subroutine is a Lambert Algorithm which given two radius vectors and
*	the time to get from one to the other, it finds the orbit connecting the
*	two. It solves the problem using a new algorithm developed by R. Battin.
*	It solves the Lambert problem for all possible types of orbits (circles, 
*	ellipses, parabolas and hyperbolas).
*	The only singularity is for the case of a transfer angle of 360 degrees,
*	which is a rather obscure case.
*	It computes the velocity vectors corresponding to the given radius
*	vectors except for the case when the transfer angle is 180 degrees
*	in which case the orbit plane is ambiguous (an infinite number of
*	transfer orbits exist).
* 
*	2- MULTIREVOLUTION LAMBERT'S PROBLEM
*
*	For the solution of Lambert's problem with Nrev>0 number of revolution,
*	Battin's formulation has been extended to accomodate N-revolution
*	transfer orbits, by following the paper: "USING BATTIN METHOD TO OBTAIN 
*	MULTIPLE-REVOLUTION LAMBERT'S SOLUTIONS" by Shen and Tsiotras.
*
*	When Nrev>0 the possible orbits are just ellipses.
*	If 0<=Nrev<=Nmax, there are two Nrev-revolution transfer orbits.
*	These two transfer orbits have different semi-major axis and they may be 
*	all combinations of large-e and small-e transfer orbits.
*	The Original Successive Substitution Method by Battin converges to one of
*	the two possible solution with a viable initial guest, however it
*	diverges from the other one. Then a Reversed Successive Substitution is
*	used to converge to the second solution.
*	A procedure is implemented in order to guarantee to provide initial
*	guesses in the convergence region. If Nrew exceeds the maximum number of
*	revolution an ERROR is given:
*	warning('off','lambertMR:SuccessiveSubstitutionDiverged') to take out the
*	warnings or use optionsLMR(1) = 0.
* 
*	3- INVERSION OF THE MOTION
* 
*	Direct or retrograde option can be selected for the transfer
*
*	The algorithm computes the semi-major axis, the parameter (semi-latus 
*	rectum), the eccentricity and the velocity vectors.
* 
*	Note: If ERROR occurs or the 360 or 180 degree transfer case is 
*	encountered. 
* 
*	Issues to be solved:
*		- 180 degrees transfer indetermination
*		- 360 degrees transfer singularity
*       - Nmax number of max revolution for a given TOF:
*			work in progress - Camilla Colombo
*
*	Notes:
*		The semi-major axis, positions, times, & gravitational parameter
*		must be in compatible units.
*
* INPUT:
*	(const double) ri[3]	Pointer to vector of initial position [L].
*   (const double) rf[3]    Pointer to vector of final position [L].
*   (double) tof            Transfer time, time of flight [T]
*   (double) mu             Planetary constant of the planet
*                           (mu = mass * G) [L^3/T^2].
*   (int) orbittype         Logical variable defining whether transfer is
*								0: direct transfer from R1 to R2
*									(counterclockwise);
*                               1: retrograde transfer from R1 to R2
*									(clockwise).
*   (int) nrev              Number of revolutions:
*								0: ZERO-REVOLUTION transfer is calculated.
*                               if nrev > 0 two transfers are possible.
*									Ncase should be defined to select one of
*									the two.
*   (int) ncase             Logical variable defining the small-a or
*                           large-a option in case of nrev > 0:
*								0: small-a option;
*                               1: large-a option.
*   (int) optionslmr        Display options:
*								0: no display;
*                               1: warnings are displayed only when
*									the algorithm does not converge;
*                               2: full warnings displayed.
*      
* OUTPUT:
*	(double *) a_out		Semi-major axis of the transfer orbit [L].
*   (double *) p_out        Semi-latus rectum of the transfer orbit [L].
*   (double *) e_out        Eccentricity of the transfer orbit.
*   (int *) error_out       Error flag:
*								0: No error
*                               1: Error, routine failed to converge
*                               -1: 180 degrees transfer
*                               2: 360 degrees transfer
*                               3: the algorithm doesn't converge because
*									the number of revolutions is bigger
*                                   than nrevmax for that TOF.
*                               4: Routine failed to converge, maximum
*									number of iterations exceeded.
*   (double) vi_out[3]		Pointer to initial velocity vector [L].
*   (double) vf_out[3]      Pointer to final velocity vector [L].
*   (double) *tpar_out      Parabolic flight time between ri and rf [T].
*   (double) *theta_out     Transfer angle [rad].
*  
* NON-STANDARD LIBRARIES:
*	mathUtils.h
*
* REFERENCES:
*	"USING BATTIN METHOD TO OBTAIN MULTIPLE-REVOLUTION LAMBERT'S SOLUTIONS"
*		by Shen and Tsiotras.
*
* ORIGINAL VERSION:
*	Chris D'Souza, 20/01/1989, MATLAB, lambert.m
*
* AUTHOR:
*	Camilla Colombo, 10/11/2006, MATLAB, lambertMR.m
*
* PORTING:
*	Matteo Ceriotti, 29/01/2009, from MATLAB, lambertMR.m
*  
* CHANGELOG:
*	13/11/2006, Camilla Colombo: added ERROR = 3 if Nrev > NrevMAX
*   21/11/2006, Camilla Colombo: added another case of ERROR = 3 (index
*		N3) corresponding to the limit case when small-a solution =
*       large-a solution. No solution is given in this case.
*   06/08/2007, Camilla Colombo: optionsLMR added as an input
*   28/11/2007, Camilla Colombo: minor changes
*   29/01/2009, Matteo Ceriotti:
*		- Introduced variable for maximum number of iterations nitermax.
*       - Corrected final check on maximum number of iterations
*			exceeded, from "==" to ">=" (if N1 >= nitermax || N >= nitermax).
*       - Increased maxumum number of iterations to 2000, not to lose
*			some solutions.
*       - In OSS loop, added check for maximum number of iterations
*			exceeded, which then sets checkNconvOSS = 0.
*       - Changed the way of coumputing X given Y1 in RSS. Now the
*			Newton-Raphson method with initial guess suggested by Shen,
*           Tsiotras is used. This should guarantee convergence without
*           the need of an external zero finder (fsolve).
*       - Changed absolute tolerance into relative tolerance in all
*			loops X0-X.
*           Now the condition is: while "abs(X0-X) >= abs(X)*TOL+TOL".
*       - Added return immediately when any error is detected. - Moved
*           check on 4*TOF*LAMBDA==0 after computing LAMBDA. - Moved check
*           on THETA==0 || THETA==2*PI after computing THETA. - Added
*           error code 4 (number of iterations exceeded). - Removed
*           variable Nwhile, as not strictly needed. - Removed variable
*           PIE=pi.
*	29/01/2009, REVISION: Matteo Ceriotti.
*	29/01/2009, Matteo Ceriotti: Conversion of the M-file into C
*
*	Note: Please if you have got any changes that you would like to be done,
*       do not change the function, please contact the authors.
*
* -------------------------------------------------------------------------
*/

double h_E(const double E, const double y, const double m, const int nrev);
/* Used by lambertMR
*/

double h_E2(const double E, const double y, const double m, const int nrev, double *dh);
/* Used by lambertMR
*/

#endif /* TRANSFER_HIGHTHRUST_H */
