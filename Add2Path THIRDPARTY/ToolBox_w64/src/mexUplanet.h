/******************************************************************************
 *                      Mex Gateway routine for uplanet                       *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

/* To compile the mex file, use the following command:
mex mexUplanet.c ephemerides.c conversions.c astro_constants.c math_utils.c -output uplanet -DMEXCOMPILE
*/

#ifndef MEXUPLANET_H
#define MEXUPLANET_H

/* Include */
#include "mex.h"
#include "ephemerides.h"

/* Prototype */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* MEXUPLANET_H */
