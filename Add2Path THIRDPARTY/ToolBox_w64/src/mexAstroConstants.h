/******************************************************************************
 *                  Mex Gateway routine for astroConstants                    *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

/* To compile the mex file, use the following command:
mex mexAstro_constants.c astro_constants.c math_utils.c -output astro_constants -DMEXCOMPILE
*/

#ifndef MEXASTRO_CONSTANTS_H
#define MEXASTRO_CONSTANTS_H

/* Include */
#include "mex.h"
#include "astroConstants.h"
#include "mathUtils.h"

/* Prototype */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* MEXASTRO_CONSTANTS_H */
