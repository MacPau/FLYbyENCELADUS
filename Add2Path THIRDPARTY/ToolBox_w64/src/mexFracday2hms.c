/******************************************************************************
 *                     Mex Gateway routine for fracday2hms                    *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

#include "mexFracday2hms.h"

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {

    /* Declaration */
    double in;          /* Fraction of day (input) */
    double hms[3];      /* Array containing the number of hours, minutes and seconds (output) */
    
    
    /* Check the number of inputs and outputs */
    if(nrhs != 1)
        mexErrMsgTxt("FRACDAY2HMS requires 1 and only 1 input.");
    if(nlhs > 3)
        mexErrMsgTxt("Too many outputs. FRACDAY2HMS can return up to 3 outputs.");
    
    
    /* Get values of input */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])
                             || mxGetM(prhs[0])!=1
                             || mxGetN(prhs[0])!=1
                             || mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("The only input of FRACDAY2HMS must be a real scalar!");
    in = *mxGetPr(prhs[0]);
    
    
    /* Compute the output */
    fracday2hms(in,hms);
    
    
    /* Return the outputs */
    plhs[0] = mxCreateDoubleScalar(hms[0]);
    if (nlhs>1)
    {
        plhs[1] = mxCreateDoubleScalar(hms[1]);
        
        if (nlhs>2)
            plhs[2] = mxCreateDoubleScalar(hms[2]);
    }
    

    /* Exit the function */
    return;
}
