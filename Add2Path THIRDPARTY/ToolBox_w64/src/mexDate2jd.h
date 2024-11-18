/******************************************************************************
 *                       Mex Gateway routine for date2jd                      *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

/* To compile the mex file, use the following command:
mex mexDate2jd.c conversions.c math_utils.c -output date2jd -DMEXCOMPILE
*/

#ifndef MEXDATE2JD_H
#define MEXDATE2JD_H

/* Include */
#include "mex.h"
#include "conversion.h"

/* Prototype */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* MEXDATE2JD_H */
