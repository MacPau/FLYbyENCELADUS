/******************************************************************************
 *                       Mex Gateway routine for ephSS_car                    *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

#include "mexEphSS_car.h"

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {

    /* Declaration */
    int ibody;              /* Identifier of the body (input) */
    double mjd2000;         /* Date in MJD2000 (input) */
    double *pXP, *pVP;      /* pointer to Cartesian position and velocity (outputs) */
    double xp[3], vp[3];	/* Cartesian position and velocity */
    int error;              /* error flag returned */
    
    /* Check the inputs and outputs */
    if(nrhs != 2)
        mexErrMsgTxt("EPHSS_CAR requires 2 and only 2 inputs.");
    if(nlhs > 3)
        mexErrMsgTxt("Too many outputs. EPHSS_CAR can return 1 or 2 outputs. output 3 is error code.");

    /* Get values of inputs */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])
                             || mxGetN(prhs[0])*mxGetM(prhs[0]) != 1
                             || mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("EPHSS_CAR: The first input (ibody) must be an integer scalar.");
    ibody = roundd2i(mxGetScalar(prhs[0]));
    
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])
                             || mxGetN(prhs[1])*mxGetM(prhs[1]) != 1
                             || mxGetNumberOfDimensions(prhs[1]) != 2)
        mexErrMsgTxt("EPHSS_CAR: The second input (mjd2000) must be a scalar.");
    mjd2000  = mxGetScalar(prhs[1]);
    

    /* Call the uplanet function */
    error = ephSS_car(ibody, mjd2000, xp, vp);
	if (error){
		if (error == 1)
			mexErrMsgTxt("EPHSS_CAR: Error in calling uplanet!");
		else if (error == 2)
			mexErrMsgTxt("EPHSS_CAR: Error in calling ephNEO!");
		else
			mexErrMsgTxt("EPHSS_CAR: Unknown error!");
	}

    /*Allocate memory and assign output pointer */
    plhs[0] = mxCreateDoubleMatrix(1,3,mxREAL);
    
    /* Pointer to first output */
    pXP    = mxGetPr(plhs[0]);
    pXP[0] = xp[0];
    pXP[1] = xp[1];
    pXP[2] = xp[2];

    /* If 2 outputs are requested, return also the velocity vector */
    if(nlhs>1)
    {
        plhs[1] = mxCreateDoubleMatrix(1,3,mxREAL);
        pVP     = mxGetPr(plhs[1]);
        pVP[0]  = vp[0];
        pVP[1]  = vp[1];
        pVP[2]  = vp[2];
        
        if(nlhs>2)      /* error requested */
            plhs[2] = mxCreateDoubleScalar((double) error);

    }
    
    /* Exit the function */
    return;
}
