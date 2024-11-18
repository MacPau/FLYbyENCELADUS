/******************************************************************************
*                             anomConv.c                                      *
*                 Functions for converting anomalies                          *
*                                                                             *
*                                                       Matteo Ceriotti, 2009 *
*******************************************************************************/

#include "anomConv.h"

int meanAn2eccAn(const double m, const double e, double *psi){
    double g, gd;
    int i;
    
    if(e >= 1)          /* Works only for ellipses */
        return 1;
    
    *psi = m;             /* psi is the eccentric anomaly, uses m as a first guess */
    
    for (i=0;i<MEANAN2ECCAN_MAXITER;i++)
    {
        g       = m - (*psi - e*sin(*psi));
        gd      = -1. + e*cos(*psi);
        *psi    = *psi - g/gd;               /* Computes the eccentric anomaly kep */
    }
    return 0;
}

int eccAn2trueAn(const double psi, const double e, double *th){
    if(e >= 1)          /* Works only for ellipses */
        return 1;
    
    (*th) = 2. * atan(sqrt((1. + e) / (1. - e)) * tan(psi/2.)); /* In radians */
    return 0;
}
