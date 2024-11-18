/*******************************************************************************
*                              conversion.h                                   *
*               Functions for time and frame conversion                       *
*                                                                             *
*                             Space Toolbox                                   *
*******************************************************************************/

#ifndef CONVERSION_H
#define CONVERSION_H

/* Include libraries
 * ----------------- */
#include "mathUtils.h"     /* In-house C mathematical utilities */
#include <math.h>

#ifdef MEXCOMPILE
#include "mex.h"            /* For MATLAB MEX files */
#endif

#ifndef MEXCOMPILE
#include <stdio.h>
#endif

/* Constants
 * --------- */
#define MJD_ZERO_IN_JD 2400000.5    /* Date of modified Julian Day = 0
                                       in Julian day */
#define MJD2000_ZERO_IN_MJD 51544.5 /* Date of modified Julian day 2000 = 0
                                       in modified Julian day */
#define CAR2KEP_ELIMIT_CIR 0.00000001
/* Threshold on eccentricity for considering the orbit to be circular
* Value determined comparing the relative error on state and position
* between using the circular case and the elliptic case. With this elimit
* the relative error on position and velocity is always less then 1e-7.
*/
#define CAR2KEP_ELIMIT_PAR 0.00000001   /* This is set arbitrarly */

#define KEP2CAR_ELIMIT_CIR CAR2KEP_ELIMIT_CIR  /* This is set arbitrarly */
#define KEP2CAR_ELIMIT_PAR CAR2KEP_ELIMIT_PAR  /* This is set arbitrarly */

#define RTH2TNH_TOL_ECC 1e-10 /* Tolerance on eccentricity to identify a parabola */
#define TNH2RTH_TOL_ECC 1e-10 /* Tolerance on eccentricity to identify a parabola */

/* Mex compiler switch
 * ------------------- */
#ifdef MEXCOMPILE
#define printf mexPrintf
#define malloc mxMalloc
#define calloc mxCalloc
#define free mxFree
#endif

/* Prototypes
 * ---------- */


/******************************************************************************
 *                   REFERENCE FRAME CONVERSION FUNCTIONS                    *
 ******************************************************************************/

int kep2car(const double kep[6], const double mu, const double p, double out[6]);
/* kep2car - Converts from keplerian orbital elements to Cartesian coordinates
* 
* PROTOTYPE:
*  int kep2cart(const double kep[6], const double mu, double out[6])
* 
* DESCRIPTION:
*  Converts from keplerian orbital elements to Cartesian coordinates.
*  All units to be consistent each other. Angles in radians.
*  Notes:
*      - Tolerance on eccentricity is used to identify the circular and
*          parabolic cases.
*      - In case of a parabola, kep[0] must be the radius of pericentre,
*          and not the semi-major axis (Inf for a parabola).
*      - In the case of hyperbola, theta must be such that the point is on
*          the physical leg of the hyperbola (the leg around the attracting
*          body).
*
* INPUT:
*  (double) kep[6] Pointer to a vector of 6 doubles containing the
*                  mean keplerian elements at date.
*                      kep[0] = semimajor axis [L]
*                      kep[1] = eccentricity
*                      kep[2] = inclination, rad
*                      kep[3] = right ascension of the ascending
*                               node, rad
*                      kep[4] = argument of perigee, rad
*                      kep[5] = true anomaly, rad
*  (double) p      Semi-latus rectum [L]. It is used only for the case of
*                  parabola (in which case kep[0] is not used).
*  (double) mu     Planetary constant of the central body. [L^3/(M*T^2)]
*    
* OUTPUT:
*  (double) out[6] Pointer to a vector of 6 doubles containing the
*                  position and velocity vector in cartesian
*                  coordinates. out should be pre-allocated. [L, L/T]
*  (int) kep2car   Integer giving the type of orbit:
*                      0 for circular orbit
*                      1 for elliptic orbit
*                      2 for parabola
*                      3 for hyperbola
*
* NON-STANDARD LIBRARIES:
*  mathUtils
*  
* AUTHOR:
*  Massimiliano Vasile, 2002, MATLAB, kep2cart.m
*
* PORTING:
*  Nicolas Croisard, 17/09/2008, from MATLAB, kep2car.m
*
* CHANGELOG:
*   07/10/2010, Matteo Ceriotti:
*       Added the semi-latus rectum to the input to eliminate the
*       assumption kep(1)=rp.
*
* -------------------------------------------------------------------------
*/

int car2kep(const double in[6], const double mu, double kep[6], double *eccAnom, double *meanAnom, double *dt, double *p);
/* car2kep - Convertion from cartesian position and velocity to keplerian
*  elements.
* 
* PROTOTYPE:
*  int cart2kep(const double in[6], const double mu,
*      double kep[6], double *eccAnom, double *meanAnom,
*      double *dt, double *p)
* 
* DESCRIPTION:
*  Convertion from cartesian position and velocity to keplerian
*  elements. All the units have to be consistent, angles in radians.
*
* INPUT:
*  (double) in[6]  Pointer to a vector of 6 doubles containing the
*                  position and velocity vector in cartesian
*                  coordinates. [L, L/T]
*  (double) mu     Planetary constant of the central body. [L^3/(M*T^2)]
*    
* OUTPUT:
*  (double) kep[6]     Pointer to a vector of 6 doubles containing the
*                      keplerian elements at date. kep should be
*                      pre-allocated.
*                          kep[0] = semimajor axis [L]
*                          kep[1] = eccentricity
*                          kep[2] = inclination, rad (0 <= kep[2] <= PI)
*                          kep[3] = right ascension of the ascending
*                              node in rad (0<=kep[3]<PI2)
*                          kep[4] = argument of perigee, rad
*                              (0 <= kep[4] < PI2)
*                          kep[5] = true anomaly, rad
*                              (0 <= kep[5] < PI2)
*  (double *) eccAnom  Pointer to a double containing the eccentric
*                      anomaly, hyperbolic anomaly or parabolic anomaly
*                      (for definitions see Vallado pag. 49).
*  (double *) meanAnom Pointer to a double containing the mean anomaly
*  (double *) dt       Pointer to a double containing the time from the
*                      pericentre passage.
*  (double *) p        Pointer to a double containing the parameter.
*  (int) kep2cart      Integer giving the type of orbit:
*                          0 for circular orbit
*                          1 for elliptic orbit
*                          2 for parabola
*                          3 for hyperbola
*
* REFERENCES:
*  D. A. Vallado, "Fundamentals of Astrodynamics and Applications, Second
*  Edition", Microcosm Press, 2001, pp. 101-102.
*
* NON-STANDARD LIBRARIES:
*  mathUtils
*
* ORIGINAL VERSION:
*  Massimiliano Vasile, 2002, MATLAB, ca2kep.m
*  
* AUTHOR:
*  Matteo Ceriotti, 08/02/2007
*
* PORTING:
*  Nicolas Croisard, 17/09/2008, from MATLAB, car2kep.m
*  
* CHANGELOG:
*  02/12/2008, Matteo Ceriotti: Added correction for om when i == 0.
*  04/12/2008, Matteo Ceriotti: Added correction for om for planar
*      retrograde orbits (i == pi). Changed help for bounds on Om, om,
*      th, now in [0, 2pi).
*  30/09/2009, Matteo Ceriotti: Header and function name in accordance
*      with guidlines.
*  16/10/2009, Matteo Ceriotti ported modifications by Jeannette Heiligers,
*      Camilla Colombo:
*      - Added threshold value for eccentricity in case of circular orbit.
*      - Modified condition for determinimg is om = 2*pi-om. New condition
*          condition works for prograde/retrograde, planar/inclined case.
*      - Modified condition for determinimg is th = 2*pi-om. New condition
*          condition works for prograde/retrograde, planar/inclined case.
*          Previous condition was wrong in the circular case because
*          dot(r,v) jumps between <0 and >0.
*  19/10/2009, Matteo Ceriotti, Camilla Colombo:
*      - Added threshold value for eccentricity in case of circular orbit,
*          value determined with test_car2kep.m.
*      - Modified condition for determinimg is om = 2*pi-om. New condition
*          condition works for prograde/retrograde, planar/inclined case.
*      - Modified condition for determinimg is th = 2*pi-om. New condition
*          condition works for prograde/retrograde, planar/inclined case.
*          Previous condition was wrong in the circular case because
*          dot(r,v) jumps between <0 and >0.
*      
* -------------------------------------------------------------------------
*/

int car2bpl(const double x_car[3], const double U_car[3], const double vp_car[3], double x_bpl[3]);
/* car2bpl - Vector reference frame transformation from Cartesian
*	reference frame to b-plane reference frame.
*
* DESCRIPTION:
*	Cartesian reference frame: {x,y,z} inertial reference frame.
*	The b-plane is the plane perpendicular to the incoming relative velocity
*	of the small body at the planet arrival, containing the planet.
*	b-plane reference frame: {xi,eta,zeta} where
*		eta-axis: as the incoming relative velocity of the small body on its
*			arrival.
*		zeta-axis: in the b-plane, direction opposite to the projection of
*			the heliocentric velocity of the planet on the b-plane.
*		xi-axis: in the b-plane, completes the reference frame.
*
* PROTOTYPE:
*   int car2bpl(const double x_car[3], const double U_car[3],
*		const double vp_car[3], double x_bpl[3])
*
* INPUT:
*   (double) x_car[3]   Pointer to a vector of 3 doubles containing the
*                       coordinates of the vector to be transformed,
*                       expressed in {x,y,z}.
*   (double) U_car[3]   Pointer to a vector of 3 doubles containing the
*                       velocity of the small body relative to the
*                       planet, expressed in {x,y,z}.
*   (double) vp_car[3]  Pointer to a vector of 3 doubles containing the
*                       velocity of the planet, expressed in {x,y,z}.
*
* OUTPUT:
*   (double) x_bpl[3]   Pointer to a vector of 3 doubles containing the
*                       coordinates of the vector x_car expressed in
*                       {xi,eta,zeta}.
*                       x_bpl should be pre-allocated.
*   (int) car_bplT      Error code: 1 if the norm of U_car or vp_car is
*						null, 0 otherwise.
*
* NON-STANDARD LIBRARIES:
*   mathUtils
*
* REFERENCES:
*   (none)
*
* AUTHOR:
*   Camilla Colombo, 04/05/2007, MATLAB, car2bpl.m
*
* PORTING:
*   Nicolas Croisard, 19/09/2008, from MATLAB, car2bpl.m
* 
* CHANGELOG:
*
* -----------------------------------------------------------------------------
*/

int eulerAxisAngle(const double v[3],const double n[3], const double theta, double v1[3]);
/* EulerAxisAngle - Rotates a vector about an axis of a given angle
*	(Eules axis and angle rotation).
*
* DESCRIPTION:
*   Rotates a vector about an axis of a given angle (counterclockwise
%   according to the right-hand rule).
%   Note: If you want to rotate the coordinate system of a given angle
%   (counterclockwise according to the right-hand rule), use the minus in
%   front of the angle (i.e., -theta). (See reference).
*
* PROTOTYPE:
*   int eulerAxisAngle(const double v[3],const double n[3],
*		const double theta, double v1[3])
*
* INPUT:
*   (double) v[3]			Pointer to a vector of 3 doubles containing
*                           the coordinates of the vector to be rotated.
*   (double) n[3]           Pointer to a vector of 3 doubles containing
*                           the coordinates of the axis of rotation.
*   (double) theta          Angle of rotation in radians.
*
* OUTPUT:
*   (double) vi[3]          Pointer to a vector of 3 doubles containing
*	                        the coordinates of the rotated vector.
*   (int) eulerAxisAngle	Error code:
*								1 if the norm of n is null
*								0 otherwise.
*
* NON-STANDARD LIBRARIES:
*   mathUtils
*
* REFERENCES:
*	Schaub and Junkins, "Analytical Mechanics of Space Systems", AIAA
*	Education Series, 2003, pp. 90.
*   see http://mathworld.wolfram.com/RotationMatrix.html for Matrix, vector
*   rotation.
*
* AUTHOR:
*   Matteo Ceriotti, 11/01/2007, MATLAB, eulerAxisAngle.m
*
* PORTING:
*   Nicolas Croisard, 19/09/2008
* 
* CHANGELOG:
*   28/03/2011, Matteo Ceriotti, Camilla COlombo:
*       Now the function rotates a vector, not a reference frame.
*  
* -------------------------------------------------------------------------
*/
     
int rth2car(const double x_rth[3], const double s_car[6], double x_car[3]);
/* rth2car - Vector reference frame transformation from radial-trasversal-h
*	reference frame to Cartesian reference frame.
*
* DESCRIPTION:
*	Cartesian reference frame: {x,y,z} inertial reference frame.
*	rth reference frame: {r,t,h}
*		r-axis: direction of the orbit radius
*		h-axis: direction of angular momentum
*		t-axis: in the orbit plane, completes the reference frame (inward)
*
* PROTOTYPE:
*   int rth2car(const double x_rth[3], const double s_car[6],
*		double x_car[3])
*
* INPUT:
*   (double) x_rth[3]   Pointer to a vector of 3 doubles containing the
*                       coordinates of the vector to be transformed,
*                       expressed in {r,t,h}.
*   (double) s_car[3]   Pointer to a vector of 6 doubles containing the
*                       position and velocity of the orbiting body,
*                       expressed in {x,y,z}
*
* OUTPUT:
*   (double) x_car[3]   Pointer to a vector of 3 doubles containing the
*                       coordinates of the vector x_rth expressed in
*                       {x,y,z}.
*                       x_car should be pre-allocated.
*   (int) rth2car       Error code:
*							1 if the position or velocity of the orbiting body
*								is null,
*                           0 otherwise.
*
* NON-STANDARD LIBRARIES:
*   mathUtils
*
* REFERENCES:
*   (none)
*
* AUTHOR:
*   Camilla Colombo, 03/03/2006, MATLAB, rth2car.m
*
* PORTING:
*   Nicolas Croisard, 19/09/2008
* 
* CHANGELOG:
* -------------------------------------------------------------------------
*/

int car2rth(const double x_car[3], const double s_car[6], double x_rth[3]);
/* car2rth - Vector reference frame transformation from Cartesian
*	reference frame to radial-trasversal-h reference frame.
*
* DESCRIPTION:
*	Cartesian reference frame: {x,y,z} inertial reference frame.
*	rth reference frame: {r,t,h}
*		r-axis: direction of the orbit radius
*		h-axis: direction of angular momentum
*		t-axis: in the orbit plane, completes the reference frame (inward)
*
* PROTOTYPE:
*   int car2rth(const double x_car[3], const double s_car[6],
*		double x_rth[3])
*
* INPUT:
*   (double) x_car[3]   Pointer to a vector of 3 doubles containing the
*                       coordinates of the vector to be transformed,
*                       expressed in {x,y,z}.
*   (double) s_car[3]   Pointer to a vector of 6 doubles containing the
*                       position and velocity of the orbiting body,
*                       expressed in {x,y,z}
*
* OUTPUT:
*   (double) x_rth[3]   Pointer to a vector of 3 doubles containing the
*                       coordinates of the vector x_car expressed in
*                       {r,t,h}.
*                       x_rth should be pre-allocated.
*   (int) car2rth       Error code:
*							1 if the position or velocity of the orbiting
*								body is null,
*                           0 otherwise.
*
* NON-STANDARD LIBRARIES:
*   mathUtils
*
* REFERENCES:
*   (none)
*
* AUTHOR
*   Camilla Colombo, 03/03/2006, MATLAB, car2rth.m
*
* PORTING:
*   Nicolas Croisard, 19/09/2008
* 
* CHANGELOG:
* -------------------------------------------------------------------------
*/

int tnh2car(const double x_tnh[3], const double s_car[6], double x_car[3]);
/* tnh2car - Vector reference frame transformation from tangent-normal-h
*	reference frame to Cartesian reference frame.
*
* DESCRIPTION:
*	Cartesian reference frame: {x,y,z} inertial reference frame.
*	thn reference frame: {t,n,h}
*		t-axis: tangent to the motion
*		h-axis: direction of angular momentum
*		n-axis: inward normal to t, in the orbit plane
*
* PROTOTYPE:
*   int tnh_carT(const double x_tnh[3], const double s_car[6], double x_car[3])
*
* INPUT:
*   (double) x_tnh[3]   Pointer to a vector of 3 doubles containing the
*                       coordinates of the vector to be transformed,
*                       expressed in {t,n,h}.
*   (double) s_car[3]   Pointer to a vector of 6 doubles containing the
*                       position and velocity of the orbiting body,
*                       expressed in {x,y,z}
*
* OUTPUT:
*   (double) x_car[3]   Pointer to a vector of 3 doubles containing the
*                       coordinates of the vector x_tnh expressed in {x,y,z}.
*                       x_car should be pre-allocated.
*   (int) tnh2car       Error code:
*							1 if the position or velocity of the orbiting body
*								is null,
*                           0 otherwise.
*
* NON-STANDARD LIBRARIES:
*   mathUtils
*
* REFERENCES:
*   (none)
*
* AUTHOR:
*   Camilla Colombo, 03/03/2006, MATLAB, tnh2car.m
*
* PORTING:
*   Nicolas Croisard, 20/09/2008
* 
* CHANGELOG:
* -------------------------------------------------------------------------
*/

int car2tnh(const double x_car[3], const double s_car[6], double x_tnh[3]);
/* car2tnh - Vector reference frame transformation from Cartesian
*	reference frame to tangent-normal-h reference frame.
*
* DESCRIPTION:
*	Cartesian reference frame: {x,y,z} inertial reference frame.
*	thn reference frame: {t,n,h}
*		t-axis: tangent to the motion
*		h-axis: direction of angular momentum
*		n-axis: inward normal to t, in the orbit plane
*
* PROTOTYPE:
*   int car2tnh(const double x_car[3], const double s_car[6], double x_tnh[3])
*
* INPUT:
*   (double) x_car[3]   Pointer to a vector of 3 doubles containing the
*                       coordinates of the vector to be transformed,
*                       expressed in {x,y,z}.
*   (double) s_car[3]   Pointer to a vector of 6 doubles containing the
*                       position and velocity of the orbiting body,
*                       expressed in {x,y,z}
*
* OUTPUT:
*   (double) x_tnh[3]   Pointer to a vector of 3 doubles containing the
*                       coordinates of the vector x_car expressed in
*                       {t,n,h}.
*                       x_tnh should be pre-allocated.
*   (int) car2tnh       Error code:
*							1 if the position or velocity of the
*								orbiting body is null,
*                           0 otherwise.
*
* NON-STANDARD LIBRARIES:
*   mathUtils
*
* REFERENCES:
*   (none)
*
* AUTHOR:
*   Camilla Colombo, 03/03/2006, MATLAB, car2tnh.m
*
* PORTING:
*   Nicolas Croisard, 20/09/2008
* 
* CHANGELOG:
* -------------------------------------------------------------------------
*/

int radec2tnh(const double x_radec[3], double x_tnh[3]);
/* radec2tnh - Vector reference frame transformation from Right Ascension and
*	Declination reference frame to tangent-normal-h reference frame.
*
* DESCRIPTION:
*	radec reference frame: {r,alpha,delta} right ascension and declination
*	(spherical equatorial) reference frame.
*		r = modulus of the vector
*		alpha = in-plane right ascention angle, counted from the tangential
*			direction to the projection of the vector on the orbital plane [rad]
*		delta = out-of-plane declination angle from the projection of the
*			vector on the orbital plane up to the vector itself [rad]
*	thn reference frame: {t,n,h}
*		t-axis: tangent to the motion
*		h-axis: direction of angular momentum
*		n-axis: inward normal to t, in the orbit plane
*
* PROTOTYPE:
*   int radec2tnh(const double x_radec[3], double x_tnh[3])
*
* INPUT:
*   (double) x_radec[3] Pointer to a vector of 3 doubles containing the
*                       coordinates of the vector to be transformed,
*                       expressed in {r,alpha,delta}.
*
* OUTPUT:
*   (double) x_tnh[3]   Pointer to a vector of 3 doubles containing the
*                       coordinates of the vector x_radec expressed in
*                       {t,n,h}.
*                       x_tnh should be pre-allocated.
*   (int) radec2tnh     Error code (0)
*
* NON-STANDARD LIBRARIES:
*   mathUtils
*
* REFERENCES:
*   (none)
*
* AUTHOR:
*   Camilla Colombo, 07/12/2007, MATLAB, radec2tnh.m
*
* PORTING:
*   Nicolas Croisard, 20/09/2008
* 
* CHANGELOG
* -------------------------------------------------------------------------
*/

int tnh2radec(const double x_tnh[3], double x_radec[3]);
/* tnh2radec - Vector reference frame transformation from tangent-normal-h
*	reference frame to Right Ascension and DEClination reference frame.
*
* DESCRIPTION:
*	radec reference frame: {r,alpha,delta} right ascension and declination
*	(spherical equatorial) reference frame.
*		r = modulus of the vector
*		alpha = in-plane right ascention angle, counted from the tangential
*			direction to the projection of the vector on the orbital
*           plane [rad]
*		delta = out-of-plane declination angle from the projection of the
*           vector on the orbital plane up to the vector itself [rad]
*	thn reference frame: {t,n,h}
*		t-axis: tangent to the motion
*		h-axis: direction of angular momentum
*		n-axis: inward normal to t, in the orbit plane
*
* PROTOTYPE:
*   tnh_radecT(const double x_tnh[3], double x_radec[3])
*
* INPUT:
*   (double) x_tnh[3] Pointer to a vector of 3 doubles containing the
*                     coordinates of the vector to be transformed,
*                     expressed in {t,n,h}.
*
* OUTPUT:
*   (double) x_radec[3] Pointer to a vector of 3 doubles containing the
*                       coordinates of the vector x_tnh expressed in
*                       {r,alpha,delta}.
*                       x_radec should be pre-allocated.
*   (int) tnh2radec     Error code: always 0, no error possible.
*
* NON-STANDARD LIBRARIES:
*   mathUtils
*
* REFERENCES:
*   (none)
*
* AUTHOR:
*   Camilla Colombo, 23/10/2007, MATLAB, tnh2radec.m
*
* PORTING:
*   Nicolas Croisard, 20/09/2008
* 
* CHANGELOG:
*   02/04/2008, Matteo Ceriotti: Removed error when input is [0,0,0].
* -------------------------------------------------------------------------
*/

int radec2car(const double x_radec[3], const double s_car[6], double x_car[3]);
/* radec2car - Vector reference frame transformation from Right Ascension and
*	DEClination reference frame to Cartesian reference frame.
*
* DESCRIPTION:
*	Cartesian reference frame: {x,y,z} inertial reference frame.
*	radec reference frame: {r,alpha,delta} right ascension and declination
*	(spherical equatorial) reference frame.
*		r = modulus of the vector
*		alpha = in-plane right ascention angle, counted from the tangential
*		    direction to the projection of the vector on the orbital
*           plane [rad]
*		delta = out-of-plane declination angle from the projection of the
*           vector on the orbital plane up to the vector itself [rad]
*
* PROTOTYPE:
*   int radec2car(const double x_radec[3], const double s_car[6],
*		double x_car[3])
*
* INPUT:
*   (double) x_radec[3] Pointer to a vector of 3 doubles containing the
*                       coordinates of the vector to be transformed,
*                       expressed in {r,alpha,delta}.
*   (double) s_car[6]   Pointer to a vector of 6 doubles containing the
*                       position and velocity of the orbiting body,
*                       expressed in {x,y,z}
*
* OUTPUT:
*   (double) x_car[3]   Pointer to a vector of 3 doubles containing the
*                       coordinates of the vector x_radec expressed in
*                       {x,y,z}.
*                       x_car should be pre-allocated.
*   (int) radec2car     Error code:
*							1 if the position or velocity of the orbiting body
*								is null,
*                           0 otherwise.
*
* NON-STANDARD LIBRARIES:
*   mathUtils
*
* REFERENCES:
*   (none)
*
* AUTHOR:
*   Camilla Colombo, 07/12/2007, MATLAB, radec2car.m
*
* PORTING:
*   Nicolas Croisard, 20/09/2008
* 
* CHANGELOG:
* -------------------------------------------------------------------------
*/

int car2radec(const double x_car[3], const double s_car[6], double x_radec[3]);
/* car2radec - Vector reference frame transformation from Cartesian reference
*	frame to Right Ascension and DEClination reference frame.
*
* DESCRIPTION:
*	Cartesian reference frame: {x,y,z} inertial reference frame.
*	radec reference frame: {r,alpha,delta} right ascension and declination
*	(spherical equatorial) reference frame.
*		r = modulus of the vector
*		alpha = in-plane right ascention angle, counted from the tangential
*           direction to the projection of the vector on the orbital
*           plane [rad]
*		delta = out-of-plane declination angle from the projection of the
*           vector on the orbital plane up to the vector itself [rad]
*
* PROTOTYPE:
*   int car2radec(const double x_car[3], const double s_car[6],
*		double x_radec[3])
*
* INPUT:
*   (double) x_car[3]   Pointer to a vector of 3 doubles containing the
*                       coordinates of the vector to be transformed,
*                       expressed in {x,y,z}.
*   (double) s_car[6]   Pointer to a vector of 6 doubles containing the
*                       position and velocity of the orbiting body,
*                       expressed in {x,y,z}
*
* OUTPUT:
*   (double) x_radec[3] Pointer to a vector of 3 doubles containing the
*                       coordinates of the vector x_car expressed in
*                       {r,alpha,delta}.
*                       x_radec should be pre-allocated.
*   (int) car2radec     Error code: 1 if the position or velocity of the
*                                     orbiting body is null,
*                                   2 if x_car is null
*                                   0 otherwise.
*
* NON-STANDARD LIBRARIES:
*   mathUtils
*
* REFERENCES:
*   (none)
*
* AUTHOR:
*   Camilla Colombo, 23/10/2007, MATLAB, car2radec.m
*
* PORTING:
*   Nicolas Croisard, 20/09/2008
* 
* CHANGELOG
* -------------------------------------------------------------------------
*/

int rth2tnh(const double x_rth[3], const double a, const double e, const double f, const double mu, double x_tnh[3]);
/* rth2tnh - Vector reference frame transformation from radial-transversal-h
*	reference frame to tangent-normal-h reference frame.
*
* DESCRIPTION:
*	This functions does not work for parabola.
*
*	rth reference frame: {r,t,h}
*		r-axis: direction of the orbit radius
*		h-axis: direction of angular momentum
*		t-axis: in the orbit plane, completes the reference frame (inward)
*	thn reference frame: {t,n,h}
*		t-axis: tangent to the motion
*		h-axis: direction of angular momentum
*		n-axis: inward normal to t, in the orbit plane
*
* PROTOTYPE:
*   int rth2tnh(const double x_rth[3], const double a, const double e,
*		const double f, const double mu, double x_tnh[3])
*
* INPUT:
*   (double) x_rth[3]   Pointer to a vector of 3 doubles containing the
*                       coordinates of the vector to be transformed,
*                       expressed in {r,t,h}.
*   (double) a          Double constaining the semi-major axis.
*   (double) e          Double constaining the eccentricity.
*   (double) f          Double constaining the true anomaly from the
*                       pericenter in rad.
*   (double) mu         Double constaining the gravitational constant of
*                       the central body.
*
* OUTPUT:
*   (double) x_tnh[3]   Pointer to a vector of 3 doubles containing the
*                       coordinates of the vector x_rth expressed in
*                       {t,n,h}.
*                       x_tnh should be pre-allocated.
*   (int) rth2tnh       Error code: 1 if e is within 1 +/- 1e-10 (parabola),
*                                   0 otherwise.
*
* NON-STANDARD LIBRARIES:
*   mathUtils
*
* REFERENCES:
*   (none)
*
* AUTHOR:
*   Camilla Colombo, 10/03/2006, MATLAB, rth2tnh.m
*
* PORTING:
*   Nicolas Croisard, 20/09/2008
* 
* CHANGELOG:
* -------------------------------------------------------------------------
*/

int tnh2rth(const double x_tnh[3], const double a, const double e, const double f, const double mu, double x_rth[3]);
/* tnh2rth - Vector reference frame transformation from radial-transversal-h
*	reference frame to tangent-normal-h reference frame.
*
* DESCRIPTION:
*	rth reference frame: {r,t,h}
*		r-axis: direction of the orbit radius
*		h-axis: direction of angular momentum
*		t-axis: in the orbit plane, completes the reference frame (inward)
*	thn reference frame: {t,n,h}
*		t-axis: tangent to the motion
*		h-axis: direction of angular momentum
*		n-axis: inward normal to t, in the orbit plane
*
* PROTOTYPE:
*   int tnh2rth(const double x_tnh[3], const double a, const double e,
*		const double f, const double mu, double x_rth[3])
*
* INPUT:
*   (double) x_tnh[3]   Pointer to a vector of 3 doubles containing the
*                       coordinates of the vector to be transformed,
*                       expressed in {t,n,h}.
*   (double) a          Double constaining the semi-major axis.
*   (double) e          Double constaining the eccentricity.
*   (double) f          Double constaining the true anomaly from the
*                       pericenter in rad.
*   (double) mu         Double constaining the gravitational constant of
*                       the central body.
*
* OUTPUT:
*   (double) x_rth[3]   Pointer to a vector of 3 doubles containing the
*                       coordinates of the vector x_rth expressed in
*                       {r,t,h}.
*                       x_rth should be pre-allocated.
*   (int) tnh2rth       Error code: 1 if a is zero,
*                                   0 otherwise.
*
* NON-STANDARD LIBRARIES:
*   mathUtils
*
* REFERENCES:
*   (none)
*
* AUTHOR:
*   Camilla Colombo, 23/02/2006, MATLAB, tnh2rth.m
*
* PORTING:
*   Nicolas Croisard, 20/09/2008
* 
* CHANGELOG:
* -------------------------------------------------------------------------
*/








/******************************************************************************
 *                    CALENDAR CONVERSION FUNCTIONS                           *
 ******************************************************************************/

double hms2fracday(const double hrs, const double mn, const double sec);
/* hms2fracday - Convert hours, minutes, and seconds into a fraction of day.
* 
* PROTOTYPE:
*   double hms2fracday(const double hrs, const double mn,
*                      const double sec)
*
* INPUT:
*   (double) hrs            Number of hours
*   (double) mn             Number of minutes
*   (double) sec            Number of seconds
* 
* OUTPUT:
*   (double) hms2fracday    A double representing the time expressed in
*                           day
* 
* REFERENCES:
*   (none)
*
* AUTHOR:
*   Nicolas Croisard, 16/02/2008, MATLAB, hms2fracday.m
*
* PORTING:
*   Nicolas Croisard, 22/09/2008
* 
* CHANGELOG:
* -------------------------------------------------------------------------
*/


int fracday2hms(const double fracday, double hms[3]);
/* FRACDAY2HMS - Convert a fraction of day into hours, minutes, and seconds.
* 
* PROTOTYPE:
*   int fracday2hms(const double fracday, double hms[3])
*
* INPUT:
*   (double) fracday        A double representing the time expressed in
*                           day
* 
* OUTPUT:
*   (double) hms[3]         Array of 3 doubles.
*                               hms[0] = Number of hours
*                               hms[1] = Number of minutes
*                               hms[2] = Number of seconds
*   (int) fracday2hms       Etrror code(0)
*
* REFERENCES:
*   (none)
*
* AUTHOR:
*   Nicolas Croisard, 16/02/2008, MATLAB, fracday2hms.m
*
* PORTING:
*   Nicolas Croisard, 23/09/2008
* 
* CHANGELOG:
* -------------------------------------------------------------------------
*/

double jd2mjd(const double jd);
/* jd2mjd - Modified Julian day number from Julian day number.
* 
* DESCRIPTION:
*	This function returns the modified Julian day number corresponding to
*	the given Julian day number.
* 
* PROTOTYPE:
*   double jd2mjd(const double jd)
*
* INPUT:
*   (double) jd         Double containing the date in Julian Day.
*                       The JD (Julian day) count is from 0 at 12:00
*                       noon, 1 January -4712 (4713 BC), Julian
*                       proleptic calendar. The corresponding date in
*                       Gregorian calendar is 12:00 noon, 24 November
*                       -4713.
* 
* OUTPUT:
*   (double) jd2mjd     Double containing the date in Modified Julian
*                       Day.
*                       The MJD count is from 00:00 midnight at the
*                       beginning of Wednesday November 17, 1858.
* 
* REFERENCES:
*   (none)
* 
* AUTHOR:
*   Nicolas Croisard, 16/02/2008, MATLAB, jd2mjd.m
*
* PORTING:
*   Nicolas Croisard, 23/09/2008
* 
* CHANGELOG:
* -------------------------------------------------------------------------
*/

double mjd2jd(const double mjd);
/* mjd2jd - Julian day number from modified Julian day number.
*
* DESCRIPTION:
*	This function returns the Julian day number corresponding to the given
*	modified Julian day number.
* 
* PROTOTYPE:
*   double mjd2jd(const double mjd)
*
* INPUT:
*   (double) mjd        Double containing the date in Modified Julian
*                       Day.
*                       The MJD count is from 00:00 midnight at the
*                       beginning of Wednesday November 17, 1858.
* 
* OUTPUT:
*   (double) mjd2jd     Double containing the date in Julian Day.
*                       The JD (Julian day) count is from 0 at 12:00
*                       noon, 1 January -4712 (4713 BC), Julian
*                       proleptic calendar. The corresponding date in
*                       Gregorian calendar is 12:00 noon, 24 November
*                       -4713.
* 
* REFERENCES:
*   (none)
* 
* AUTHOR:
*   Nicolas Croisard, 16/02/2008, MATLAB, mjd2jd.m
*
* PORTING:
*   Nicolas Croisard, 23/09/2008
* 
* CHANGELOG:
*   20/02/2008, Revised by Matteo Ceriotti
*   23/09/2008, Nicolas Croisard: Conversion of the M-file into C
* 
* ------------------------ - SpaceART Toolbox - ------------------------ */

double jd2mjd2000(const double jd);
/* jd2mjd2000 - Modified Julian day 2000 number from Julian day number.
* 
* DESCRIPTION:
*	This function returns the modified Julian day 2000 number corresponding to
*	the given Julian day number.
* 
* PROTOTYPE:
*   double jd2mjd2000(const double jd)
*
* INPUT:
*   (double) jd         Double containing the date in Julian Day.
*                       The JD (Julian day) count is from 0 at 12:00
*                       noon, 1 January -4712 (4713 BC), Julian
*                       proleptic calendar. The corresponding date in
*                       Gregorian calendar is 12:00 noon, 24 November
*                       -4713.
* 
* OUTPUT:
*   (double) jd2mjd2000 Double containing the date in Modified Julian
*                       Day 2000.
*                       MJD2000 is defined as the number of days since
*                       01-01-2000, 12:00 noon.
* 
* REFERENCES:
*   (none)
* 
* AUTHOR:
*   Nicolas Croisard, 16/02/2008, MATLAB, jd2mjd2000.m
*
* PORTING:
*   Nicolas Croisard, 24/09/2008
* 
* CHANGELOG:
* -------------------------------------------------------------------------
*/

double mjd20002jd(const double mjd2000);
/* mjd20002jd - Julian day number from modified Julian day 2000 number.
* 
* DESCRIPTION:
*	This function returns the Julian day number corresponding to the given
*	modified Julian day 2000 number.
* 
* PROTOTYPE:
*   double mjd20002jd(const double mjd2000)
* 
* INPUT:
*   (double) mjd2000    Double containing the date in Modified Julian
*                       Day 2000.
*                       MJD2000 is defined as the number of days since
*                       01-01-2000, 12:00 noon.
* 
* OUTPUT:
*   (double) mjd20002jd Double containing the date in Julian Day.
*                       The JD (Julian day) count is from 0 at 12:00
*                       noon, 1 January -4712 (4713 BC), Julian
*                       proleptic calendar. The corresponding date in
*                       Gregorian calendar is 12:00 noon, 24 November
*                       -4713.
* 
* REFERENCES:
*   (none)
* 
* AUTHOR:
*   Nicolas Croisard, 16/02/2008, MATLAB, mjd20002jd.m
*
* PORTING:
*   Nicolas Croisard, 24/09/2008
* 
* CHANGELOG:
* -------------------------------------------------------------------------
*/

double mjd20002mjd(const double mjd2000);
/* mjd20002mjd - Modified Julian day number from modified Julian day 2000
*	number.
*
* DESCRIPTION:
*	This function returns the modified Julian day number corresponding to
*	the given modified Julian day 2000 number.
* 
* PROTOTYPE:
*   double mjd20002mjd(const double mjd2000)
* 
* INPUT:
*   (double) mjd2000        Double containing the date in Modified
*                           Julian Day 2000.
*                           MJD2000 is defined as the number of days
*                           since 01-01-2000, 12:00 noon.
* 
* OUTPUT:
*   (double) mjd20002mjd    Double containing the date in Modified
*                           Julian Day.
*                           The MJD count is from 00:00 midnight at the
*                           beginning of Wednesday November 17, 1858.
* 
* REFERENCES:
*   (none)
* 
* AUTHOR:
*   Nicolas Croisard, 16/02/2008, MATLAB, mjd20002mjd.m
*
* PORTING:
*   Nicolas Croisard, 24/09/2008
* 
* CHANGELOG:
* -------------------------------------------------------------------------
*/

double mjd2mjd2000(const double mjd);
/* mjd20002mjd - Modified Julian day 2000 number from modified Julian day
*             number.
* 
* DESCRIPTION:
*	This function returns the modified Julian day 2000 number corresponding
*	to the given modified Julian day number.
* 
* PROTOTYPE:
*   double mjd2mjd2000(const double mjd)
* 
* INPUT:
*   (double) mjd            Double containing the date in Modified
*                           Julian Day.
*                           The MJD count is from 00:00 midnight at the
*                           beginning of Wednesday November 17, 1858.
* 
* OUTPUT:
*   (double) mjd2mjd2000    Double containing the date in Modified
*                           Julian Day 2000.
*                           MJD2000 is defined as the number of days
*                           since 01-01-2000, 12:00 noon.
* 
* REFERENCES:
*   (none)
* 
* AUTHOR:
*   Nicolas Croisard, 16/02/2008, MATLAB, mjd2mjd2000.m
*
* PORTING:
*   Nicolas Croisard, 24/09/2008
* 
* CHANGELOG:
* -------------------------------------------------------------------------
*/

int date2jd(const double date[6], double *jd);
/* date2jd - Julian day number from Gregorian date.
* 
* DESCRIPTION:
*	This function computes the Julian day number of the given date
*	(Gregorian calendar) plus a fractional part depending on the time of day
* 
*	Note:
*		1/ The function is valid for the whole range of dates since 12:00
*			noon 24 November -4713, Gregorian calendar. (This bound is set
*			in order to have symmetry with the inverse function JD2DATE)
*		2/ The inputs must be feasible (i.e. the date must exist!). If an
*			unfeasible date is inputed, wrong results are given because no
*			check is done on that.
*		3/ For dates before 1582, the resulting date components are valid
*			only in the Gregorian proleptic calendar. This is based on the
*			Gregorian calendar but extended to cover dates before its
*			introduction.
* 
* PROTOTYPE:
*   int date2jd(const double date[6], double *jd)
*
* INPUT:
* 	(double) date[6]    Pointer to a vector of 6 doubles containing the
*                       date in the Gregorian calendar.
*                       date = [year, month, day, hour, minute, second]
* 
* OUTPUT:
*   (double) jd         Double containing the date in Julian Day.
*                       The JD (Julian day) count is from 0 at 12:00
*                       noon, 1 January -4712 (4713 BC), Julian
*                       proleptic calendar. The corresponding date in
*                       Gregorian calendar is 12:00 noon, 24 November
*                       -4713.
*   (int) date2jd       Error code: 1 if the date is invalid (jd<0),
*                                   0 otherwise.
* 
* REFERENCES:
*   Formula from http://scienceworld.wolfram.com/astronomy/JulianDate.html
*   (last visited 15/02/2008)
*   Compared to http://pdc.ro.nu/mjd.cgi for a few dates, the same results
*   were found
* 
* AUTHOR:
*   Nicolas Croisard, 16/02/2008, MATLAB, date2jd.m
*
* PORTING:
*   Nicolas Croisard, 22/09/2008
* 
* CHANGELOG:
*   03/03/2008, Revised by Camilla Colombo
*   22/09/2008, Nicolas Croisard: Conversion of the M-file into C
* 
* ------------------------ - SpaceART Toolbox - ------------------------ */

int jd2date(const double jd, double date[6]);
/* jd2date - Gregorian date from Julian day number.
*
* DESCRIPTION:
*	This function computes the given date (Gregorian calendar) from the
*	Julian day number.
* 
*	Note:
*		1/ jd must be a non-negative real. This means that the function is
*			valid for the whole range of dates since 12:00 noon 24 November
*			-4713, Gregorian calendar.
*		2/ For dates before 1582, the resulting date components are valid
*			only in the Gregorian proleptic calendar. This is based on the
*			Gregorian calendar but extended to cover dates before its
*			introduction.
* 
* PROTOTYPE:
*   int jd2date(const double jd, double date[6])
*
* INPUT:
*   (double) jd         Double containing the date in Julian Day.
*                       The JD (Julian day) count is from 0 at 12:00
*                       noon, 1 January -4712 (4713 BC), Julian
*                       proleptic calendar. The corresponding date in
*                       Gregorian calendar is 12:00 noon, 24 November
*                       -4713.
* 
* OUTPUT:
* 	(double) date[6]    Pointer to a vector of 6 doubles containing the
*                       date in the Gregorian calendar.
*                       date = [year, month, day, hour, minute, second]
*                       date must be pre-allocated.
*   (int) jd2date       Error code: 1 if the date is invalid (jd<0),
*                                   0 otherwise.
* 
* REFERENCES:
*   Formula from http://scienceworld.wolfram.com/astronomy/JulianDate.html
*   (last visited 15/02/2008)
*   Compared to http://pdc.ro.nu/mjd.cgi for a few dates, the same results
*   were found
* 
* AUTHOR:
*   Nicolas Croisard, 16/02/2008, MATLAB, jd2date.m
*
* PORTING:
*   Nicolas Croisard, 23/09/2008
* 
* CHANGELOG:
*   03/03/2008, Revised by Camilla Colombo
*   23/09/2008, Nicolas Croisard: Conversion of the M-file into C
* 
* ------------------------ - SpaceART Toolbox - ------------------------ */

int date2mjd(const double date[6], double *mjd);
/* DATE2MJD - Modified Julian day number from Gregorian date.
*
* DESCRIPTION:
*	This function computes the modified Julian day number of the given date
*	(Gregorian calendar) plus a fractional part depending on the time of day
* 
*	Note:
*		1/ The function is valid for the whole range of dates since 12:00
*			noon 24 November -4713, Gregorian calendar.
*		2/ The inputs must be feasible (i.e. the date must exist!). If an
*			unfeasible date is inputed, wrong results are given because no
*			check is done on that.
*		3/ For dates before 1582, the resulting date components are valid
*			only in the Gregorian proleptic calendar. This is based on the
*			Gregorian calendar but extended to cover dates before its
*			introduction.
* 
* PROTOTYPE:
*   int date2mjd(const double date[6], double *mjd)
* 
* INPUT:
* 	(double) date[6]    Pointer to a vector of 6 doubles containing the
*                       date in the Gregorian calendar.
*                       date = [year, month, day, hour, minute, second]
* 
* OUTPUT:
*   (double) mjd        Double containing the date in Modified Julian
*                       Day.
*                       The MJD count is from 00:00 midnight at the
*                       beginning of Wednesday November 17, 1858.
*   (int) date2mjd      Error code:
*							1 if the date is invalid (corresponding jd<0),
*                           0 otherwise.
*
* REFERENCES:
*   (none)
* 
* AUTHOR:
*   Nicolas Croisard, 16/02/2008, MATLAB, date2mjd.m
*
* PORTING:
*   Nicolas Croisard, 23/09/2008
* 
* CHANGELOG:
*   03/03/2008, Revised by Camilla Colombo
*   23/09/2008, Nicolas Croisard: Conversion of the M-file into C
* 
* ------------------------ - SpaceART Toolbox - ------------------------ */

int mjd2date(const double mjd, double date[6]);
/* mjd2date - Gregorian date from modified Julian day number.
*
* DESCRIPTION:
*	This function computes the given date (Gregorian calendar) from the
*	modified Julian day number.
* 
*	Note:
*		1/ mjd must be a valid (ie the corresponding jd>=0). This means that
*			the function is valid for the whole range of dates since 12:00
*			noon 24 November -4713, Gregorian calendar.
*		2/ For dates before 1582, the resulting date components are valid
*			only in the Gregorian proleptic calendar. This is based on the
*			Gregorian calendar but extended to cover dates before its
*			introduction.
* 
* PROTOTYPE:
*   int mjd2date(const double mjd, double date[6])
*
* INPUT:
*   (double) mjd        Double containing the date in Modified Julian
*                       Day.
*                       The MJD count is from 00:00 midnight at the
*                       beginning of Wednesday November 17, 1858.
* 
* OUTPUT:
* 	(double) date[6]    Pointer to a vector of 6 doubles containing the
*                       date in the Gregorian calendar.
*                       date = [year, month, day, hour, minute, second]
*                       date must be preallocated.
*   (int) mjd2date      Error code:
*							1 if the date is invalid (corresponding jd<0),
*                           0 otherwise.
* 
* REFERENCES:
*   (none)
* 
* AUTHOR:
*   Nicolas Croisard, 16/02/2008, MATLAB, mjd2date.m
*
* PORTING:
*   Nicolas Croisard, 23/09/2008
* 
* CHANGELOG:
* -------------------------------------------------------------------------
*/

int date2mjd2000(const double date[6], double *mjd2000);
/* date2mjd2000 - Modified Julian day 2000 number from Gregorian date.
* 
* DESCRIPTION:
*	This function computes the modified Julian day 2000 number of the given
*	date (Gregorian calendar) plus a fractional part depending on the time
*	of day.
* 
*	Note:
*		1/ The function is valid for the whole range of dates since 12:00
*			noon 24 November -4713, Gregorian calendar.
*		2/ The inputs must be feasible (i.e. the date must exist!). If an
*			unfeasible date is inputed, wrong results are given because no
*			check is done on that.
*		3/ For dates before 1582, the resulting date components are valid
*			only in the Gregorian proleptic calendar. This is based on the
*			Gregorian calendar but extended to cover dates before its
*			introduction.
* 
* PROTOTYPE:
*   int date2mjd2000(const double date[6], double *mjd2000)
* 
* INPUT:
* 	(double) date[6]    Pointer to a vector of 6 doubles containing the
*                       date in the Gregorian calendar.
*                       date = [year, month, day, hour, minute, second]
* 
* OUTPUT:
*   (double) mjd2000    Double containing the date in Modified Julian
*                       Day 2000.
*                       MJD2000 is defined as the number of days since
*                       01-01-2000, 12:00 noon.
*   (int) date2mjd2000  Error code:
*							1 if the date is invalid (corresponding jd<0),
*                           0 otherwise.
* 
* REFERENCES:
*   (none)
* 
* AUTHOR:
*   Nicolas Croisard, 16/02/2008, MATLAB, date2mjd2000.m
*
* PORTING:
*   Nicolas Croisard, 24/09/2008
* 
* CHANGELOG:
* -------------------------------------------------------------------------
*/

int mjd20002date(const double mjd2000, double date[6]);
/* mjd20002date - Gregorian date from modified Julian day 2000 number.
* 
* DESCRIPTION:
*	This function computes the given date (Gregorian calendar) from the
*	modified Julian day 2000 number.
* 
*	Note:
*		1/ mjd2000 must be a valid (ie the corresponding jd>=0). This means
*			that the function is valid for the whole range of dates since
*			12:00 noon 24 November -4713, Gregorian calendar.
*		2/ For dates before 1582, the resulting date components are valid
*			only in the Gregorian proleptic calendar. This is based on the
*			Gregorian calendar but extended to cover dates before its
*			introduction.
* 
* PROTOTYPE:
*   int mjd20002date(const double mjd2000, double date[6])
*
* INPUT:
*   (double) mjd2000    Double containing the date in Modified Julian
*                       Day 2000.
*                       MJD2000 is defined as the number of days since
*                       01-01-2000, 12:00 noon.
* 
* OUTPUT:
* 	(double *) date[6]  Pointer to a vector of 6 doubles containing the
*                       date in the Gregorian calendar.
*                       date = [year, month, day, hour, minute, second]
*                       date must be preallocated.
*   int mjd20002date    Error code:
*							1 if the date is invalid (corresponding jd<0),
*                           0 otherwise.
* 
* REFERENCES:
*   (none)
* 
* AUTHOR:
*   Nicolas Croisard, 16/02/2008, MATLAB, mjd20002date.m
*
* PORTING:
*   Nicolas Croisard, 24/09/2008
* 
* CHANGELOG:
* -------------------------------------------------------------------------
*/

#endif /* CONVERSION_H */
