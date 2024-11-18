/******************************************************************************
*                               astroConstants.h                              *
*                  Definition of the astronautical constants                  *
*																			  *
*								 Space Toolbox								  *
******************************************************************************/

#ifndef ASTROCONSTANTS_H
#define ASTROCONSTANTS_H

/* Generic astronomical constants:
 * ------------------------------- */
/* Universal gravity constant (G) [km^3/(kg*s^2)], from Ditan */
#define GRAVITY_CONSTANT 6.67259e-20
/* Astronomical Unit (AU) [km], from Ditan */
#define AU 149597870.7
/* Speed of light in the vacuum [km/s], definition in the SI */
#define V_LIGHT 299792.458
/* Standard free fall (the acceleration due to gravity on the Earth's surface at
 * sea level) [m/s^2], definition in Wertz */
#define FREE_FALL 9.80665
/* Mean distance Earth-Moon [km], definition in Wertz */
#define D_EARTH_MOON 384401.0
/* Obliquity (angle) of the ecliptic at Epoch 2000 [rad], definition in Wertz */
#define OBLIQUITY_ECLIPTIC_2000 (23.43928111/180.0 * 3.141592653589793)
/* Energy flux density of the Sun [W/m2 at 1 AU], definition from SMAD/Wertz */
#define FLUX_DENSITY 1367.0


/* Planetary constants of the solar system bodies (mu = mass * G) [km^3/s^2],
 * from DITAN
 * ------------------------------------------------------------------------- */
#define MU_SUN  (0.19891000000000E+31*GRAVITY_CONSTANT) /* SUN */
#define MU_ME   (0.33020000000000E+24*GRAVITY_CONSTANT) /* MERCURY */
#define MU_V    (0.48685000000000E+25*GRAVITY_CONSTANT) /* VENUS */
#define MU_E    (0.59736990612667E+25*GRAVITY_CONSTANT) /* EARTH */
#define MU_M    (0.64184999247389E+24*GRAVITY_CONSTANT) /* MARS */
#define MU_J    (0.18986000000000E+28*GRAVITY_CONSTANT) /* JUPITER */
#define MU_S    (0.56846000000000E+27*GRAVITY_CONSTANT) /* SATURN */
#define MU_U    (0.86832000000000E+26*GRAVITY_CONSTANT) /* URANUS */
#define MU_N    (0.10243000000000E+27*GRAVITY_CONSTANT) /* NEPTUNE */
#define MU_P    (0.14120000000000E+23*GRAVITY_CONSTANT) /* PLUTON */
#define MU_MOON (0.73476418263373E+23*GRAVITY_CONSTANT) /* MOON */
#define MU_PLANETS {MU_ME, MU_V, MU_E, MU_M, MU_J, MU_S, MU_U, MU_N, MU_P, MU_MOON}

/* Mean radius of the solar system bodies [km], from DITAN
 * ------------------------------------------------------- */
#define R_SUN  700000.0 /* SUN */
#define R_ME   0.24400000000000E+04 /* MERCURY */
#define R_V    0.60518000000000E+04 /* VENUS */
#define R_E    0.63781600000000E+04 /* EARTH */
#define R_M    0.33899200000000E+04 /* MARS */
#define R_J    0.69911000000000E+05 /* JUPITER */
#define R_S    0.58232000000000E+05 /* SATURN */
#define R_U    0.25362000000000E+05 /* URANUS */
#define R_N    0.24624000000000E+05 /* NEPTUNE */
#define R_P    0.11510000000000E+04 /* PLUTON */
#define R_MOON 0.17380000000000E+04 /* MOON */
#define R_PLANETS {R_ME, R_V, R_E, R_M, R_J, R_S, R_U, R_N, R_P, R_MOON}

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
#endif


/* Prototype
 * --------- */
 
int astroConstants(const int n, const int in[], double out[]);
/*  astroConstants - Astronautical constants of the solar system
*
* PROTOTYPE:
*     int astroConstants(const int n, const int in[], double out[])
*
* DESCRIPTION:
*	This function returns a vector of constants corresponding to the input
*	vector query. For single value, it is recommended to use the #define
*	conastants instead.
* 
*	List of identifiers:
*		Generic astronomical constants:
*			1   Universal gravity constant (G) [km^3/(kg*s^2)]
*			2   Astronomical Unit (AU) [km]
*		Sun related:
*			3   Sun mean radius [km]
*			4   Sun planetary constant (mu = mass * G) [km^3/s^2]
*			31  Energy flux density of the Sun [W/m2 at 1 AU]
*		Other:
*			5   Speed of light in the vacuum [km/s]
*			6   Standard free fall (the acceleration due to gravity on the
*				Earth's surface at sea level) [m/s^2]
*			7   Mean distance Earth-Moon [km]
*			8   Obliquity (angle) of the ecliptic at Epoch 2000 [rad]
*		Planetary constants of the planets (mu = mass * G) [km^3/s^2]:
*			11  Me
*			12  V
*			13  E
*			14  Ma
*			15  J
*			16  S
*			17  U
*			18  N
*			19  P
*			20  Moon
*		Mean radius of the planets [km]:
*			21  Me
*			22  V
*			23  E
*			24  Ma
*			25  J
*			26  S
*			27  U
*			28  N
*			29  P
*			30  Moon
*
*	Instructions for upgrading this function:
*		1) It is possible to add new constants. Please DO NOT change the
*			structure of this function, as well as its prototype.
*		2) DO NOT change the identifiers of the constants that have already
*			been defined in this function. If you want to add a new constant,
*			please use an unused identifier.
*		3) DO NOT add constants that can be easily computed starting form
*			other ones (avoid redundancy).
*		4) ALWAYS check that you are modifying the latest version of the
*			function, and make the new version available to everyone as soon
*			as possible.
*		5) When adding new constants, please:
*           a) Add the constant as a #define
*           b) Upgrade the help of the function
*		Thanks!
*			                                 Matteo Ceriotti, Nicolas Croisard
*
* INPUT:
*	(int) n                 Integer of the number of elements of the vector in.
*   (int) in[n]             Pointer to the vector of identifier of the
*                           required constants.
* 
* OUTPUT:
*   (double) out[n]			Pointer to the vector of constants
*                           corresponding to the vector of identifiers.
*                           The vector must be preallocated outside of
*                           length n.
*   (int) astroConstants    Error code: 1 if an element of the query
*                           vector does not correspond to any available
*                           identifier, 0 otherwise.
*
* NON-STANDARD LIBRARIES:
*	(none)
*
* REFERENCES:
*   James R. Wertz and Wiley J. Larson (editors), "Space Mission
*		Analysis and Design", Third Edition, 1999
*   DITAN
*   System International (SI)
*
* AUTHOR:
*	Matteo Ceriotti, 2006, MATLAB, astroConstants.m
*
* PORTING:
*   Nicolas Croisard, 17/09/2008, from MATLAB, astroConstants.m
*
* CHANGELOG:
*
* -------------------------------------------------------------------------
*/

     
#endif /* ASTROCONSTANTS_H */



    

