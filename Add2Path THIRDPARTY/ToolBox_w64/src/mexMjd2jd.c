/******************************************************************************
 *                        Mex Gateway routine for mjd2jd                      *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

#include "mexMjd2jd.h"

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {

    /* Declaration */
    double mjd;         /* Date in MJD (input) */
    double jd;          /* Date in JD (output) */
    
    
    /* Check the number of inputs and outputs */
    if(nrhs != 1)
        mexErrMsgTxt("MJD2JD requires 1 and only 1 input.");
    if(nlhs > 1)
        mexErrMsgTxt("Too many outputs. MJD2JD returns only 1 output.");
    
    
    /* Get values of input */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])
                             || mxGetM(prhs[0])!=1
                             || mxGetN(prhs[0])!=1
                             || mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("The only input of MJD2JD must be a real scalar!");
    mjd = *mxGetPr(prhs[0]);
    
    
    /* Compute the output */
    jd = mjd2jd(mjd);
    
    
    /* Return the output */
    plhs[0] = mxCreateDoubleScalar(jd);
    

    /* Exit the function */
    return;
}
