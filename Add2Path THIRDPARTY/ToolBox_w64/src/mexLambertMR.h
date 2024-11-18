/******************************************************************************
 *                  Mex Gateway routine for lambertMR                         *
 *                                                                            *
 *                                                 Matteo Ceriotti, 2008      *
 ******************************************************************************/

#ifndef MEXLAMBERTMR_H
#define MEXLAMBERTMR_H

/* Include */
#include "mex.h"
#include "transfer_highThrust.h"

/* Constants */

/* Prototypes */
void cpdouble(const double *in, double *out, int n);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* MEXLAMBERTMR_H */
