/******************************************************************************
*                                  dynamics.c                                 *
*                             Functions for dynamics                          *
*                                                                             *
*                                Space Toolbox                                *
******************************************************************************/

#ifndef DYNAMICS_H
#define DYNAMICS_H

/* Include libraries
 * ----------------- */
#include "mathUtils.h"        /* In-house C mathematical utilities */
/*#include "conversion.h"*/       /* In-house reference frame conversion functions */
/*#include "ephemerides.h"*/      /* In-house ephemerides of the solar system */
/*#include "astroConstants.h"*/   /* In-house astronautical constants */
 
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


/* Prototypes
 * ---------- */

int dyn_2BP(const double t, const double s[6], const double mu, double out[6]);
/*  dyn_2BP - Equation of the dynamics of the two body problem.
*
* DESCRIPTION:
*	All units to be consistent with each other.
*
* PROTOTYPE:
*   int d2b_eq(const double t, const double s[6], const double mu,
*		double out[6])
*
* INPUT:
*   (double) t          Double containing the time
*   (double) s[6]       Pointer to a vector of 6 doubles containing the
*                       position and velocity vectors.
*   (double) mu         Planetary constant of the planet. It is
*                       recommended to use the constants defined in
*                       astro_constants.h
*
* OUTPUT:
*   (double) out[6]     Pointer to a vector of 6 doubles containing the
*                       derivative of the state vector.
*                       out should be pre-allocated.
*   (int) dyn_2BP		Error code (always 0)
*
* NON-STANDARD LIBRARIES:
*   mathUtils
*
* REFERENCES:
*   (none)
*
* AUTHOR:
*   Matteo Ceriotti, 2006, MATLAB
*
* PORTING:
*   Nicolas Croisard, 23/10/2008
*
* CHANGELOG:
*
* -------------------------------------------------------------------------
*/

int dyn_2BP_u(const double t, const double s[6], const double u[3], const double mu, double out[6]);
/*  dyn_2BP_u - Equation of the dynamics of the two body problem,
*		with control.
*
* DESCRIPTION:
*	All units to be consistent with each other.
*
* PROTOTYPE:
*   int d2b_eq_u(const double t, const double s[6], const double u[3],
*		const double kp, double out[3])
*
* INPUT:
*	(double) t          Double containing the time
*   (double) s[6]       Pointer to a vector of 6 doubles containing the
*                       position and velocity vectors.
*   (double) u[3]       Pointer to a vector of 3 doubles containing the
*                       control vector.
*   (double) mu         Planetary constant of the planet. It is
*                       recommended to use the constants defined in
*                       astro_constants.h
*
* OUTPUT:
*   (double) out[6]		Pointer to a vector of 6 doubles containing the
*						derivative of the state vector.
*						out should be pre-allocated.
*   (int) dyn_2BP_u		Error code (always 0)
*
* NON-STANDARD LIBRARIES:
*   mathUtils
*
* REFERENCES:
*   (none)
*
* AUTHOR:
*   Camilla Colombo, 2005, MATLAB
*
* PORTING:
*   Nicolas Croisard, 23/10/2008
*
* CHANGELOG:
*
* -------------------------------------------------------------------------
*/

#endif /* DYNAMICS_H */
