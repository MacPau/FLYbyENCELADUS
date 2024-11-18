/******************************************************************************
 *                     Mex Gateway routine for mjd2mjd2000                    *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

/* To compile the mex file, use the following command:
mex mexMjd2mjd2000.c conversions.c math_utils.c -output mjd2mjd2000 -DMEXCOMPILE
*/

#ifndef MEXMJD2MJD2000_H
#define MEXMJD2MJD2000_H

/* Include */
#include "mex.h"
#include "conversion.h"

/* Prototype */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* MEXMJD2MJD2000_H */
