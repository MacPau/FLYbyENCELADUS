/******************************************************************************
*                                  dynamics.c                                 *
*                             Functions for dynamics                          *
*                                                                             *
*                                Space Toolbox                                *
******************************************************************************/

#include "dynamics.h"

int dyn_2BP(const double t, const double s[6], const double mu, double out[6])
{
    /* Copy in out the velocity */
    dcopy(s+3,out,3);
    
    /* Compute the acceleration */
    dscalarmult(s, -mu/DCUBE(dnorm(s,3)), 3, out+3);
    
    return 0;
}

int dyn_2BP_u(const double t, const double s[6], const double u[3], const double mu, double out[6])
{
    /* Copy in out the velocity */
    dcopy(s+3,out,3);
    
    /* Compute the acceleration */
    dscalarmult(s, -mu/DCUBE(dnorm(s,3)), 3, out+3);
    out[3] += u[0];
    out[4] += u[1];
    out[5] += u[2];
    
    return 0;
}
