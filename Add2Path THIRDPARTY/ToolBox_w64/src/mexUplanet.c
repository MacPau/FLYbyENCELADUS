/******************************************************************************
 *                      Mex Gateway routine for uplanet                       *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

#include "mexUplanet.h"

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {

    /* Declaration */
    double mjd2000;
    int ibody;
    double *kep;
    double *ksun;
    int error;
    
    /* Check the inputs and outputs */
    if(nrhs != 2)
        mexErrMsgTxt("UPLANET requires 2 and only 2 inputs.");
    if(nlhs > 2)
        mexErrMsgTxt("Too many outputs. Uplanet can return 1 or 2 outputs.");

    /* Get values of inputs */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])
                             || mxGetN(prhs[0])*mxGetM(prhs[0]) != 1
                             || mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("UPLANET: The first input (mjd2000) must be a scalar.");
    mjd2000 = mxGetScalar(prhs[0]);
    
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])
                             || mxGetN(prhs[1])*mxGetM(prhs[1]) != 1
                             || mxGetNumberOfDimensions(prhs[1]) != 2)
        mexErrMsgTxt("UPLANET: The second input (ibody) must be an integer scalar.");
    ibody = roundd2i(mxGetScalar(prhs[1]));
    
    
    /*Allocate memory and assign output pointer */
    plhs[0] = mxCreateDoubleMatrix(1,6,mxREAL);
    
    /* Pointer to first output */
    kep = mxGetPr(plhs[0]);

    /* Call the uplanet function */
    error = uplanet(mjd2000, ibody, kep);
    if (error == 1)
        mexErrMsgTxt("UPLANET: The ephemeris of the MOON shall be computed via MOON_EPH (by calling EPHSS).");
    if (error == 2)
        mexErrMsgTxt("UPLANET: The second input (id) is out of bound!");

    /* If 2 outputs are requested, return also the Sun planetary constant */
    if(nlhs==2)
    {
        plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
        ksun    = mxGetPr(plhs[1]);
        ksun[0] = MU_SUN;
    }
    
    /* Exit the function */
    return;
}
