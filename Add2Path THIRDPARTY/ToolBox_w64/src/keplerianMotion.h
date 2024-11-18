/******************************************************************************
*                              keplerianMotion.h                         *
*                       Functions for Keplerian motion                   *
*                                                                        *
*                                Space Toolbox                           *
******************************************************************************/

#ifndef KEPLERIANMOTION_H
#define KEPLERIANMOTION_H

/* Include libraries
 * ----------------- */
#include "mathUtils.h"		/* In-house C mathematical utilities */
#include "ephemerides.h"    /* In-house ephemerides of the solar system */
#include "astroConstants.h"	/* In-house astronautical constants */
#include "conversion.h"

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

/* Constants */
#define KEPPRO_TINYNUMBER (1e-10)
#define KEPPRO_HUGENUMBER (1e299)
#define KEPPRO_DEFAULTMAXITER 1000
#define KEPPRO_DEFAULTTOL (1e-6)

/* Inline definitions
 * ------------------ */
static double semiperiod_hohmarg1, semiperiod_hohmarg2;

#define SEMIPERIOD_HOHM(mu,a) (semiperiod_hohmarg1=(mu),semiperiod_hohmarg2=(a),(PI*sqrt(DCUBE(semiperiod_hohmarg2)/semiperiod_hohmarg1)))
/* SEMIPERIOD_HOHM - Semi period of an Hohmann transfer
* 
* PROTOTYPE:
*   semiperiod = SEMIPERIOD_HOHM(mu,a);
*
* INPUT:
*   (double) mu                 Double containing the planetary constant
*                               [L^3/T^2]
*   (double) a                  Double containing the semi-major axis
*                               [L]
*
* OUTPUT:
*   (double) SEMIPERIOD_HOHM    Semi period of the orbit [T]
*
* NON-STANDARD LIBRARIES:
*   mathUtils
*
* REFERENCES:
*   (none)
*
* AUTHOR:
*   Camila Colombo, 02/08/2006, MATLAB
*
* PORTING:
*   Nicolas Croisard, 30/10/2008
*
* CHANGELOG:
*
* -------------------------------------------------------------------------
*/

/* Prototypes
 * ---------- */

int kepPro(double out[6], const double in[6], const double dt, const double mu, const int npar, double tol, int maxiter);
/* kepPro - Analytic Keplerian propagator with universal variables.
*  
* PROTOTYPE:
*	int kepPro(double out[6], const double in[6], double dt, double mu,
*		const int npar, double tol, int maxiter)
*
* DESCRIPTION:
*	Analytic orbit propagator for Keplerian orbits.
*	Uses an algorithm from D. A. Vallado (see references section).
*  
* INPUT:
*	(double) in[6]		Pointer to the input state vector. 6 components
*                       (position and velocity).
*   (double) dt         Propagation time.
*   (double) mu         Planetary constant.
*   (int) npar          Number of following inputs actually passed.
*   (double) tol        Tolerance on the time law. Uses the same number
*                       to check if the parameter alpha is zero.
*                       Optional. If not provided, uses 1e-6.
*   (int) maxiter       Maximum number of itarations for the loop.
*                       Optional. If not provided, uses 1000.
*      
* OUTPUT:
*   (double) out[6]     Pointer to the output state vector. 6 components
*                       (position and velocity), preallocated.
*   (int) kepPro        Error code.
*							0 = no error;
*                           1 = check on (f*gd - fd*g == 1) failed;
*                           2 = number of iterations exceeded before
*								convergence.
*  
* NON-STANDARD LIBRARIES:
*   (none)
*
* REFERENCES:
*	D. A. Vallado, "Fundamentals of Astrodynamics and Applications, Second
*	Edition", Microcosm Press, pp. 101-102.
*
* AUTHOR:
*   Matteo Ceriotti, 13-02-2007, MATLAB, kepPro3.m
*
* PORTING:
*   Matteo Ceriotti, 02/07/2008, from MATLAB, kepPro3.m
*  
* CHANGELOG:
*	21/08/2008, Matteo Ceriotti: Corrected a bug which led to non-sense
*		output when fabs(dt) < TINYNUMBER (or dt == 0, of course).
*   25/06/2009, Matteo Ceriotti: Corrected a bug for initial condition in
*       the parabolic case, at line 64.
*   21/07/2009, Matteo Ceriotti: Added condition:
*		(fabs((s-c)/s) < LAMBERTMR_TOL)
*		at line 193 to prevent the case in which 4.*tof*lambda not exactly
*		zero, but very close to it.
*
* -------------------------------------------------------------------------
*/

int dvInsertion(const int ibody, const double t, const double vf[3], const double rp, const double ecc, double *dv);
/*  dvInsertion - Delta-v to insert in an orbit around a planet.
*
* PROTOTYPE:
*   int dvInsertion(const int ibody, const double t, const double vf[3],
*		const double rp, const double ecc, double *dv)
*
* INPUT:
*   (int) ibody         Integer containing the principal body index.
*                       Valid indeces: 1 to 9 (planets).
*   (double) t          Double containing the time of the insertion
*                       (final time of the trajectory) [d, MJD2000].
*   (double) vf[3]      Pointer to a vector of 3 doubles containing the
*                       absolute cartesian velocity at the planet
*                       arrival [km/s].
*   (double) rp         Double containing the radius of the pericentre
*                       of the orbit to achieve around the planet [km].
*   (double) ecc        Double containing the eccentricity of the orbit
*                       to achieve.
*
* OUTPUT:
*   (double *) dv       Pointer to a double containing the delta-v
*                       needed to insert into the defined orbit [km/s].
*   (int) dvInsertion   Error code: 1 if ibody is not within [1, 9]
*                                   2 if an error occurs in EphSS
*                                   0 otherwise
*
* NON-STANDARD LIBRARIES:
*   mathUtils
*   ephemerides
*   astroConstants
*
* REFERENCES:
*   (none)
*
* AUTHOR:
*   Matteo Ceriotti, 19/05/2007, MATLAB
*
* PORTING:
*   Nicolas Croisard, 30/10/2008
*
* CHANGELOG:
*
* -------------------------------------------------------------------------
*/

int kepEq_t(double f, const double a, const double ecc, const double mu,
        double f0, const double t0, double *t);
/*  kepEq_t - Keplerian equation. It finds the time corresponding to a certain
*		true anomaly.
*
* PROTOTYPE:
*	int kepEq_t(double f, const double a, const double ecc, const double mu,
*		double f0, const double t0, double *t);
*
* INPUT:
*   (double) f      Double containing the true anomaly [rad].
*   (double) a      Double containing the semi-major axis [L].
*   (double) ecc    Double containing the eccentricity.
*   (double) mu     Double containing the planetary constant
*                       (mu = mass * G) [L^3/T^2].
*   (double) f0     Double containing a given true anomaly [rad].
*   (double) t0     Double containing the time corresponding to f0 [T].
*
* OUTPUT:
*   (double *) t    Pointer to a double containing the time
*                   corresponding to f [T] between [-Inf, Inf].
*   (int) kepEq_t   Error code: 1 if ecc>=1
*                               0 otherwise
*
* NON-STANDARD LIBRARIES:
*   mathUtils
*
* REFERENCES:
*   (none)
*
* AUTHOR:
*   Camilla Colombo, Nicolas Croisard, 20/11/2007, MATLAB
*
* PORTING:
*   Nicolas Croisard, 04/12/2008
*
* CHANGELOG:
*
* -------------------------------------------------------------------------
*/

int kepEq_f(const double t, const double a, const double ecc, const double mu,
            const double f0, const double t0, int const imax, double const tol, 
            double *f);
/*  kepEq_f - Keplerian equation. It finds the true anomaly corresponding
*		to a given time.
*
* PROTOTYPE:
*   int kepEq_f(const double t, const double a, const double ecc,
*		const double mu, const double f0, const double t0,
*		int const imax, double const tol, double *f);
*
* DESCRIPTION:
*	It finds the true anomaly corresponding to a time t. Multiple
*   revolutions are considered. Elliptic orbits (e < 1).
*	
*	Future work:
*		Implementation of other methods to be applied if the Newton loop
*		does not converge.
*
* INPUT:
*	(double) t      Double containing the time when the true anomaly is
*                   required [T].
*   (double) a      Double containing the semi-major axis [L].
*   (double) ecc    Double containing the eccentricity.
*   (double) mu     Double containing the planetary constant
*                   (mu = mass * G) [L^3/T^2].
*   (double) f0     Double containing a given true anomaly [rad].
*   (double) t0     Double containing the time corresponding to f0 [T].
*   (int) imax      Integer containing the maximum number of iteration.
*   (double) tol    Double containing the convergence tolerance for the
*                   Newton loop on abs(E-E0) [rad].
*
* OUTPUT:
*   (double *) f    Pointer to a double containing the true anomaly at
*                   time t [rad] in [-Inf, Inf].
*                   Note: Multiple revolutions are considered.
*   (int) kepEq_f   Error code: 1 The Newton loop did not converge
*                                 within the required tolerance;
*                               2 Non elliptic orbit;
*                               0 otherwise.
*
* NON-STANDARD LIBRARIES:
*   mathUtils
*
* REFERENCES:
*   (none)
*
* AUTHOR:
*   Camilla Colombo, 20/11/2007, MATLAB
*
* PORTING:
*   Nicolas Croisard, 04/12/2008
*
* CHANGELOG:
*
* -------------------------------------------------------------------------
*/

int swingby(const double v1[3], const double rp, const double mu, const double gamma, const double n_r[3], double v2[3]);
/*  swingby - outgoing velocity after the swingby of a planet, given
*           incoming velocity and swingby parameters.
*
* PROTOTYPE:
*   int swingby(const double v1[3], const double rp, const double mu, const double gamma, const double n_r[3], double v2[3])
*
* DESCRIPTION:
*	All units to be consistent with each other, angles in radians.
*
* INPUT:
*   (double) v1[3]      Pointer to a vector of 3 doubles containing the
*                       incoming velocity vector, before the swingby,
*                       relative to the planet.
*   (double) rp         Double containing the radius of pericentre of
*                       the hyperbola.
*   (double) mu         Planetary constant of the planet. It is
*                       recommended to use the constants defined in
*                       astro_constants.h
*   (double) gamma      Plane angle [rad]. This angle identifies the
*                       inclination of the hyperbola plane around the
*                       incoming velocity vector, and is the angle
*                       between the vector n_r and the vector normal to
*                       the hyperbola plane.
*   (double) n_r[3]     Pointer to a vector of 3 doubles containing the
*                       reference vector, used as a origin to measure
*                       gamma. In principle, this vector is arbitrary. A
*                       choice can be to use the normal to the plane
*                       containing the incoming velocity and the
*                       heliocentric velocity of the planet:
*                           n_r = dcross(v1,vp)/dnorm(dcross(v1,vp),3);
*                               where vp is the heliocentric velocity of
*                               the planet.
*
* OUTPUT:
*   (double) v2[3]      Pointer to a vector of 3 doubles containing the
*                       outgoing velocity, after the swingby, relative
*                       to the planet.
*                       v2 should be pre-allocated.
*   (int) swingby       Error code: 1 or 2 if an error occurs while
*                                   calling euler_axis_angle internally,
*                                   0 otherwise.
*
* NON-STANDARD LIBRARIES:
*   mathUtils
*   conversions
*
* REFERENCES:
*   Kaplan, "Modern spacecraft dynamics and control", pag. 93
*
* AUTHOR:
*   Matteo Ceriotti, 11/01/2007
*
* PORTING:
*   Nicolas Croisard, 30/10/2008
*
* CHANGELOG:
*
* -------------------------------------------------------------------------
*/

#endif /* KEPLERIANMOTION_H */
