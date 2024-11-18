/******************************************************************************
 *                  Mex Gateway routine for euler_axis_angle                  *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

#include "mexEulerAxisAngle.h"

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {

    /* Declaration */
    double *v;      /* array to hold the coordinates of the vector to be rotated
                       (input) */
    double *n;      /* array to hold the coordinates of the axis of rotation
                       (input) */
    double theta;   /* Angle of rotation in radians (input) */
    double *v1;     /* Arrays to hold coordinates of the rotated vector
                       (output) */
    int error;      /* Error code */
    mwSize nr, nc;      /* number of rows and columns of inputs */
    
    
    /* Check the number of inputs and outputs */
    if(nrhs != 3)
        mexErrMsgTxt("EULERAXISANGLE requires 3 and only 3 inputs!");
    if(nlhs > 2)
        mexErrMsgTxt("Too many outputs. EULERAXISANGLE can return up to 2 outputs!");
    
    
    /* Get values of inputs */
    nr = mxGetM(prhs[1]);
    nc = mxGetN(prhs[1]);
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])
                             || (!(nr==3 && nc==1) && !(nr==1 && nc==3))
                             || mxGetNumberOfDimensions(prhs[1]) != 2)
        mexErrMsgTxt("The second input of EULERAXISANGLE must be a 3-elements vector of reals!");
    n = mxGetPr(prhs[1]);
    
    nr = mxGetM(prhs[0]);   /* Kept to return the same shape vector as output */
    nc = mxGetN(prhs[0]);   /* Kept to return the same shape vector as output */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])
                             || (!(nr==3 && nc==1) && !(nr==1 && nc==3))
                             || mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("The first input of EULERAXISANGLE must be a 3-elements vector of reals!");
    v = mxGetPr(prhs[0]);
    
    if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])
                             || mxGetM(prhs[2])*mxGetN(prhs[2]) != 1
                             || mxGetNumberOfDimensions(prhs[2]) != 2)
        mexErrMsgTxt("The third input of EULERAXISANGLE must be a real scalar!");
    theta = (double) mxGetScalar(prhs[2]);
    
    
    /*Allocate memory and assign output pointer */
    plhs[0] = mxCreateDoubleMatrix(nr,nc,mxREAL);

    /* Pointer to first output */
    v1 = mxGetPr(plhs[0]);
    
    
    /* Compute the outputs */
    error = eulerAxisAngle(v, n, theta, v1);

    /* Assigne the extra outputs if necessary */
    if(nlhs>1)      /* Error code requested */
        plhs[1] = mxCreateDoubleScalar((double) error);

    /* Exit the function */
    return;
}
