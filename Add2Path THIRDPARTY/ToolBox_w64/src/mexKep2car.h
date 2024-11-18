/******************************************************************************
 *                      Mex Gateway routine for kep2car                       *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

/* To compile the mex file, use the following command:
mex mexKep2cart.c conversions.c math_utils.c -output kep2cart -DMEXCOMPILE
*/

#ifndef MEXKEP2CAR_H
#define MEXKEP2CAR_H

/* Include */
#include "mex.h"
#include "conversion.h"

/* Prototype */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* MEXKEP2CAR_H */
