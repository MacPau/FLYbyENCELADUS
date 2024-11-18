/******************************************************************************
*                             anomConv.h                                      *
*                 Functions for converting anomalies                          *
*                                                                             *
*                                                       Matteo Ceriotti, 2009 *
*******************************************************************************/
 
#ifndef ANOMCONV_H
#define ANOMCONV_H

#include <math.h>

#define MEANAN2ECCAN_MAXITER 5

int meanAn2eccAn(const double m, const double e, double *psi);
/* meanAn2eccAn - Mean anomaly to eccentric anomaly.
* 
* PROTOTYPE:
*   int meanAn2eccAn(const double m, const double e, double *psi)
* 
* DESCRIPTION:
*   Works for ellipses only.
*
* INPUT:
*   (double) m          Mean anomaly, radians.
*   (double) e          Eccentricity.
*     
* OUTPUT:
*   (double *) psi      Eccentric anomaly, radians.
*   (int) meanAn2eccAn  Error code: 0: no error. 1: not ellipse.
*
* ORIGINAL VERSION:
*   uplanet.m
*   
* AUTHOR:
*   Matteo Ceriotti, 08/12/2009
*
* CHANGELOG:
*
* -------------------------------------------------------------------------
*/

int eccAn2trueAn(const double psi, const double e, double *th);
/* eccAn2trueAn - Eccentric anomaly to true anomaly.
* 
* PROTOTYPE:
*   int eccAn2trueAn(const double psi, const double e, double *th)
* 
* DESCRIPTION:
*   Works for ellipses only.
*
* INPUT:
*   (double) psi        Eccentric anomaly, radians.
*   (double) e          Eccentricity.
*     
* OUTPUT:
*   (double *) th       True anomaly, radians.
*   (int) eccAn2trueAn  Error code: 0: no error. 1: not ellipse.
*
* AUTHOR:
*   Matteo Ceriotti, 08/12/2009
*
* CHANGELOG:
*
* -------------------------------------------------------------------------
*/

#endif /* ANOMCONV_H */
