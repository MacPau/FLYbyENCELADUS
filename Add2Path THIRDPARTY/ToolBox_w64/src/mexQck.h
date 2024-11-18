/******************************************************************************
 *                        Mex Gateway routine for qck                         *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

/* To compile the mex file, use the following command:
mex mexQck.c math_utils.c -output qck -DMEXCOMPILE
*/

#ifndef MEXQCK_H
#define MEXQCK_H

/* Include */
#include "mex.h"
#include "mathUtils.h"

/* Prototype */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* MEXQCK_H */
