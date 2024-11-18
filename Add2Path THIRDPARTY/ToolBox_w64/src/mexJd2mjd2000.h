/******************************************************************************
 *                      Mex Gateway routine for jd2mjd2000                    *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

/* To compile the mex file, use the following command:
mex mexJd2mjd2000.c conversions.c math_utils.c -output jd2mjd2000 -DMEXCOMPILE
*/

#ifndef MEXJD2MJD2000_H
#define MEXJD2MJD2000_H

/* Include */
#include "mex.h"
#include "conversion.h"

/* Prototype */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* MEXJD2MJD2000_H */
