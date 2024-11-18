/******************************************************************************
 *                  Mex Gateway routine for kepPro                            *
 *                                                                            *
 *                                                 Matteo Ceriotti, 2008      *
 ******************************************************************************/

/* To compile the mex function, use the following command from MATLAB:
   mex keppro3.c mexKeppro3.c -output keppro3 -DMEXCOMPILE
*/
   
#ifndef MEXKEPPRO_H
#define MEXKEPPRO_H

/* Include */
#include "mex.h"
#include "transfer_highThrust.h"

/* Define */

/* Constants */

/* Prototypes */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* MEXKEPPRO_H */
