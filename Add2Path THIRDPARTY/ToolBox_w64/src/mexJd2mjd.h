/******************************************************************************
 *                        Mex Gateway routine for jd2mjd                      *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

/* To compile the mex file, use the following command:
mex mexJd2mjd.c conversions.c math_utils.c -output jd2mjd -DMEXCOMPILE
*/

#ifndef MEXJD2MJD_H
#define MEXJD2MJD_H

/* Include */
#include "mex.h"
#include "conversion.h"

/* Prototype */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* MEXJD2MJD_H */
