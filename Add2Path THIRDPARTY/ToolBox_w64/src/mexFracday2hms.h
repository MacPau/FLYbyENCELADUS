/******************************************************************************
 *                     Mex Gateway routine for fracday2hms                    *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

/* To compile the mex file, use the following command:
mex mexFracday2hms.c conversions.c math_utils.c -output fracday2hms -DMEXCOMPILE
*/

#ifndef MEXFRACDAY2HMS_H
#define MEXFRACDAY2HMS_H

/* Include */
#include "mex.h"
#include "conversion.h"

/* Prototype */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* MEXFRACDAY2HMS_H */
