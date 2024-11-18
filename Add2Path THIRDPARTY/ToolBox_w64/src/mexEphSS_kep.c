/******************************************************************************
 *                     Mex Gateway routine for ephSS_kep                      *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

#include "mexEphSS_kep.h"

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {

    /* Declaration */
    double mjd2000;         /* Date in MJD2000 (input) */
    int ibody;              /* Identifier of the body (input) */
    double *kep;            /* array to hold the Keplerian elements (output) */
    int error;
    
    
    /* Check the number of inputs and outputs */
    if(nrhs != 2)
        mexErrMsgTxt("EPHSS_KEP requires 2 and only 2 inputs.");
    if(nlhs > 1)
        mexErrMsgTxt("Too many outputs. EPHSS_KEP returns only 1 output.");
    
    
    /* Get values of inputs */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])
                             || mxGetN(prhs[0])*mxGetM(prhs[0]) != 1
                             || mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("EPHSS_KEP: The first input (ibody) must be an integer scalar.");
    ibody = roundd2i(mxGetScalar(prhs[0]));
    
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])
                             || mxGetN(prhs[1])*mxGetM(prhs[1]) != 1
                             || mxGetNumberOfDimensions(prhs[1]) != 2)
        mexErrMsgTxt("EPHSS_KEP: The second input (mjd2000) must be a scalar.");
    mjd2000 = mxGetScalar(prhs[1]);
    
    
    /*Allocate memory and assign output pointer */
    plhs[0] = mxCreateDoubleMatrix(1,6,mxREAL);
    
    /* Pointer to first output */
    kep = mxGetPr(plhs[0]);
    
    
    /* Call the uplanet function */
    error = ephSS_kep(ibody, mjd2000, kep);

	if (error){
		if (error == 1)
			mexErrMsgTxt("EPHSS_KEP: Error in calling uplanet!");
		else if (error == 2)
			mexErrMsgTxt("EPHSS_KEP: Error in calling ephNEO!");
		else
			mexErrMsgTxt("EPHSS_KEP: Unknown error!");
	}

    /* Exit the function */
    return;
}
