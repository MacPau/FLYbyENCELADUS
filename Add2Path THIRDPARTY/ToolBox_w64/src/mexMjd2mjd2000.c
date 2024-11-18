/******************************************************************************
 *                     Mex Gateway routine for mjd2mjd2000                    *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

#include "mexMjd2mjd2000.h"

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {

    /* Declaration */
    double mjd;         /* Date in MJD (input) */
    double mjd2000;     /* Date in MJD2000 (output) */
    
    
    /* Check the number of inputs and outputs */
    if(nrhs != 1)
        mexErrMsgTxt("MJD2MJD2000 requires 1 and only 1 input.");
    if(nlhs > 1)
        mexErrMsgTxt("Too many outputs. MJD2MJD2000 returns only 1 output.");
    
    
    /* Get values of input */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])
                             || mxGetM(prhs[0])!=1
                             || mxGetN(prhs[0])!=1
                             || mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("The only input of MJD2MJD2000 must be a real scalar!");
    mjd = *mxGetPr(prhs[0]);
    
    
    /* Compute the output */
    mjd2000 = mjd2mjd2000(mjd);
    
    
    /* Return the output */
    plhs[0] = mxCreateDoubleScalar(mjd2000);
    

    /* Exit the function */
    return;
}
