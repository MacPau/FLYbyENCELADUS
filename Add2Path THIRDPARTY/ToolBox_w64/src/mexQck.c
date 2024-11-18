/******************************************************************************
 *                        Mex Gateway routine for qck                         *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

#include "mexQck.h"

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {

    /* Declarations */
    double *pangle; /* pointer to the angle(s) to be reduced (input) */
    double *out;    /* pointer to the reduced angles (output) */
    mwSize n,m;     /* Dimension of the input array */
    int i;          /* Counter */
    
    /* Check the number of inputs and outputs */
    if(nrhs != 1)
        mexErrMsgTxt("QCK requires 1 and only 1 inputs.");
    if(nlhs > 1)
        mexErrMsgTxt("Too many outputs. QCK returns only 1 output.");
    
    
    /* Get values of inputs */
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])
                             || mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("The only input of QCK must be a real scalar, vector or matrix!");
    pangle = mxGetPr(prhs[0]);
    
    /* Allocate memory and assign output pointer */
    plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);
    out = mxGetPr(plhs[0]);
                     
    /* Compute the output */
    for(i=0;i<m*n;i++)
        out[i] = QCK(pangle[i]);       

    /* Exit the function */
    return;
}
