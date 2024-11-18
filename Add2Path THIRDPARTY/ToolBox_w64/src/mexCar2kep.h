/******************************************************************************
 *                      Mex Gateway routine for car2kep                      *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

/* To compile the mex file, use the following command:
mex mexCart2kep.c conversions.c math_utils.c -output cart2kep -DMEXCOMPILE
*/

#ifndef MEXCAR2KEP_H
#define MEXCAR2KEP_H

/* Include */
#include "mex.h"
#include "conversion.h"

/* Prototype */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* MEXCAR2KEP_H */
