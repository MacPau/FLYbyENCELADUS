/******************************************************************************
 *                     Mex Gateway routine for date2mjd2000                   *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

/* To compile the mex file, use the following command:
mex mexDate2mjd2000.c conversions.c math_utils.c -output date2mjd2000 -DMEXCOMPILE
*/

#ifndef MEXDATE2MJD2000_H
#define MEXDATE2MJD2000_H

/* Include */
#include "mex.h"
#include "conversion.h"

/* Prototype */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* MEXDATE2MJD2000_H */
