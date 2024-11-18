/******************************************************************************
 *                     Mex Gateway routine for mjd20002mjd                    *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

/* To compile the mex file, use the following command:
mex mexMjd20002mjd.c conversions.c math_utils.c -output mjd20002mjd -DMEXCOMPILE
*/

#ifndef MEXMJD20002MJD_H
#define MEXMJD20002MJD_H

/* Include */
#include "mex.h"
#include "conversion.h"

/* Prototype */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* MEXMJD20002MJD_H */
