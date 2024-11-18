/******************************************************************************
 *                     Mex Gateway routine for hms2fracday                    *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

/* To compile the mex file, use the following command:
mex mexHms2fracday.c conversions.c math_utils.c -output hms2fracday -DMEXCOMPILE
*/

#ifndef MEXHMS2FRACDAY_H
#define MEXHMS2FRACDAY_H

/* Include */
#include "mex.h"
#include "conversion.h"

/* Prototype */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* MEXHMS2FRACDAY_H */
