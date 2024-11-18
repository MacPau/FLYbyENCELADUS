/******************************************************************************
*                              ephemerides.h                           *
*                           Ephemerides functions                      *
*                                                                      *
*                              Space Toolbox                           *
******************************************************************************/ 

#ifndef EPHEMERIDES_H
#define EPHEMERIDES_H

/* Include libraries
 * ----------------- */
#include "mathUtils.h"         /* In-house C mathematical utilities */
#include "astroConstants.h"    /* In-house astronautical constants */
#include "conversion.h"        /* In-house orbital element conversions */
#include "ephNEO.h"			   /* In-house database of NEOs */

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

/* Constants 
 * --------- */


/* Prototype
 * --------- */
int uplanet(const double mjd2000, const int ibody, double kep[]);
/* uplanet - Analytical ephemerides for planets - P. Dysli, 1977.
*
* PROTOTYPE:
*   int uplanet(const double mjd2000, const int ibody, double kep[6])
*
* DESCRIPTION:
*	Planetary orbital elements restituted in a Sun-Centred ecliptic system.
* 
* INPUT:
*   (double) mjd2000	Date in MJD2000
*						(Modified Julian Day since 01/01/2000, 12:00
*						MJD2000 = MJD-51544.5)
*   (int) ibody         Number of the celestial body (< 11)
*                           1:   Mercury
*                           2:   Venus
*                           3:   Earth
*                           4:   Mars
*                           5:   Jupiter
*                           6:   Saturn
*                           7:   Uranus
*                           8:   Neptune
*                           9:   Pluto
*                           10:  Sun
*
* OUTPUT:
*   (double) kep[6]     Pointer to a vector of 6 doubles containing the
*                       mean keplerian elements at date.
*                       kep should be pre-allocated.
*                           kep[0] = semimajor axis in km
*                           kep[1] = eccentricity
*                           kep[2] = inclination in rad
*                           kep[3] = right ascension of the ascending
*                                    node in rad
*                           kep[4] = argument of perigee in rad
*                           kep[5] = true anomaly in rad
*   (int) uplanet       Error code: 0 if 1 <= ibody <= 10
*                                   1 if 11 (MOON),
*                                   2 if ibody >= 12.
*
* NON-STANDARD LIBRARIES:
*   mathUtils
*
* REFERENCES:
*   (none)
*
* AUTHOR:
*   P. Dysli, 1977, MATLAB, uplanet.m
*
* PORTING:
*   Nicolas Croisard, 16/09/2008, from MATLAB, uplanet.m
*
* CHANGELOG:
*
* -------------------------------------------------------------------------
*/

int ephMoon_pos(const double mjd2000, double xp[]);
/* ephMoon_pos - Cartesian position of the Moon.
*
* PROTOTYPE:
*   int ephMoon_pos(const double mjd2000, double xp[])
*
* DESCRIPTION:
*	This function gives the position of the Moon at a given epoch in the
*	Geocentric Equatorial Reference Frame (IAU-76/FK5 J2000, mean equator,
*	mean equinox frame). This frame {x,y,z} is characterised by:
*		x-axis: on the equatorial plane, along the direction of the gamma
*			point
*		z-axis: direction of the north pole
*		y-axis: on the equatorial plane, completes the reference frame
*
* INPUT:
*   (double) mjd2000    Date in MJD2000
*                       (Modified Julian Day since 01/01/2000, 12:00
*                        MJD2000 = MJD-51544.5)
*
* OUTPUT:
*   (double) xp[3]      Pointer to a vector of 3 doubles containing the
*                       position vector of the Moon in cartesian
*                       coordinates, expressed in the Geocentric
*                       Equatorial Reference Frame [km].
*                       xp should be pre-allocated.
*   (int) ephMoon_pos   Error code (0).
*
* NON-STANDARD LIBRARIES:
*   mathUtils
*
* REFERENCES:
*   D. A. Vallado, "Fundamentals of Astrodynamics and Applications"
*   (3rd edition), p.290 (algorithm 31).
*
* AUTHOR:
*   Daniel Novak, 04/03/2008
*
* PORTING:
*   Nicolas Croisard, 16/09/2008
*
* CHANGELOG:
*
* -------------------------------------------------------------------------
*/

int ephMoon(const double mjd2000, double xp[], double vp[]);
/* ephMoon - Cartesian position and velocity of the Moon.
*
* PROTOTYPE:
*   int ephMoon(const double mjd2000, double xp[], double vp[])
*
* DESCRIPTION:
*	This function gives the position and velocity of the Moon at a given
*	epoch in the Geocentric Equatorial Reference Frame (IAU-76/FK5 J2000,
*	mean equator, mean equinox frame). This frame {x,y,z} is characterised
*	by:
*		x-axis: on the equatorial plane, along the direction of the gamma
*			point
*		z-axis: direction of the north pole
*		y-axis: on the equatorial plane, completes the reference frame
*
* INPUT:
*   (double) mjd2000    Date in MJD2000
*                       (Modified Julian Day since 01/01/2000, 12:00
*                        MJD2000 = MJD-51544.5)
*
* OUTPUT:
*   (double) xp[3]      Pointer to a vector of 3 doubles containing the
*                       position vector of the Moon in cartesian
*                       coordinates, expressed in the Geocentric
*                       Equatorial Reference Frame [km].
*                       xp should be pre-allocated.
*   (double) vp[3]      Pointer to a vector of 3 doubles containing the
*                       velocity vector of the Moon in cartesian
*                       coordinates, expressed in the Geocentric
*                       Equatorial Reference Frame [km/s].
*                       vp should be pre-allocated.
*   (int) ephMoon		Error code (0).
*
* NON-STANDARD LIBRARIES:
*   mathUtils
*
* REFERENCES:
*   D. A. Vallado, "Fundamentals of Astrodynamics and Applications"
*   (3rd edition), p.290 (algorithm 31).
*
* AUTHOR:
*   Daniel Novak, 04/03/2008
*
* PORTING:
*   Nicolas Croisard, 16/09/2008
*
* CHANGELOG:
*
* -------------------------------------------------------------------------
*/

int ephSS_car(const int ibody, const double mjd2000, double r[], double v[]);
/* ephSS_car - Cartesian ephemerides of the solar system.
*
* PROTOTYPE:
*   int EphSS_car(const int ibody, const double mjd2000,
*             double r[3], double v[3])
*
* DESCRIPTION:
*	This function returns the ephemeris of the solar system bodies. It uses
*	UPLANET for planets, ephMoon for the Moon, and ephNEO for
*	asteroid ephemerides. It returns the cartesian position and velocity of
*	the body, centered in the Sun for all the bodies but the Moon (for which
*	a cartesian Earth-centered reference frame is chosen).
*
* INPUT:
*   (int) ibody         Number of the celestial body
*                           1 to 9: Planets
*                           10: Sun
*                           11: Moon
*                           >=12: NEOs
*   (double) mjd2000    Date in MJD2000
*                       (Modified Julian Day since 01/01/2000, 12:00
*                        MJD2000 = MJD-51544.5)
*
* OUTPUT:
*   (double) r[3]       Pointer to a vector of 3 doubles containing the
*                       position vector of the body in cartesian
*                       coordinates, expressed in the Heliocentric
*                       Ecliptic Reference Frame (Geocentric Equatorial
*                       for the Moon). [km]
*                       r should be pre-allocated.
*   (double) v[3]       Pointer to a vector of 3 doubles containing the
*                       velocity vector of the body in cartesian
*                       coordinates, expressed in the Heliocentric
*                       Ecliptic Reference Frame (Geocentric Equatorial
*                       for the Moon). [km/s]
*                       v should be pre-allocated.
*   (int) ephSS_car     Error code: 0 if no error
*                                   1 if an error occurs in uplanet
*                                   2 if an error occurs in ephNEO
*
* NON-STANDARD LIBRARIES:
*   mathUtils, astroConstants, conversion, ephNEO
*
* REFERENCES:
*   (none)
*
* AUTHOR:
*   Matteo Ceriotti, 10/01/2007, MATLAB, EphSS.m
*
* PORTING:
*   Nicolas Croisard, 16/09/2008
*
* CHANGELOG:
*
* -------------------------------------------------------------------------
*/


int ephSS_kep(const int ibody, const double mjd2000, double kep[]);
/* ephSS_kep - Ephemerides of the solar system in Keplerian parameters.
*
* PROTOTYPE:
*   int ephSS_kep(const int ibody, const double mjd2000, double kep[6])
*
* DESCRIPTION:
*	This function returns the ephemeris of the solar system bodies. It uses
*	uplanet for planets, ephMoon for the Moon, and ephNEO for
*	asteroid ephemerides. It returns the keplerian elementsof the body,
*	centered in the Sun for all the bodies but the Moon (for which a
*	cartesian Earth-centered reference frame is chosen).
*
* INPUT:
*   (int) ibody         Number of the celestial body
*                           1 to 9: Planets
*                           10: Sun
*                           11: Moon
*                           >=12: NEOs
*   (double) mjd2000    Date in MJD2000
*                       (Modified Julian Day since 01/01/2000, 12:00
*                        MJD2000 = MJD-51544.5)
*
* OUTPUT:
*   (double *) kep      Pointer to a vector of 6 doubles containing the
*                       mean keplerian elements at date.
*                       kep should be pre-allocated.
*                           kep[0] = semimajor axis in km
*                           kep[1] = eccentricity
*                           kep[2] = inclination in rad
*                           kep[3] = right ascension of the ascending
*                                    node in rad
*                           kep[4] = argument of perigee in rad
*                           kep[5] = true anomaly in rad
*   (int) EphSS_kep     Error code: 0 if no error
*                                   1 if an error occurs in UPLANET
*                                   2 if an error occurs in NEOEPHEMERIS
*
* NON-STANDARD LIBRARIES:
*   mathUtils, astroConstants, conversions, ephNEO
*
* REFERENCES:
*   (none)
*
* AUTHOR:
*   Nicolas Croisard, 03/05/2008, MATLAB, EphSS_kep.m
*
* PORTING:
*   Nicolas Croisard, 16/09/2008
*
* CHANGELOG:
*
* -------------------------------------------------------------------------
*/

#endif /* EPHEMERIDES_H */
