/******************************************************************************
 *                  Mex Gateway routine for astroConstants                    *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

#include "mexAstroConstants.h"


void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {

    /* Declarations */
    int *in;            /* Pointer to input */
    double *out;        /* Pointer to output */
    int error;          /* error flag returned */
    mwSize n_query;     /* Number of elements of input vector */
    int i;              /* Counter */

    /* Check the number of inputs and outputs */
    if(nrhs != 1)
        mexErrMsgTxt("ASTROCONSTANT requires 1 and only 1 input.");
    if(nlhs > 2)
        mexErrMsgTxt("ASTROCONSTANT returns up to 2 outputs.");

    /* Number of queries */
    if (mxGetN(prhs[0])==1)
        n_query = mxGetM(prhs[0]);
    else {
        if (mxGetM(prhs[0])==1)
            n_query = mxGetN(prhs[0]);
        else
            mexErrMsgTxt("ASTROCONSTANTS accepts a scalar or a vector as input.");
    }
        
    /* Get pointer to inputs */
    if(!(in = mxMalloc(n_query*sizeof(int))))
        mexErrMsgTxt("Out of Memory.");
        
    for(i=0;i<n_query;i++)
    {
        in[i] = roundd2i(*(mxGetPr(prhs[0]) + i));
/*        in[i] = (*(mxGetPr(prhs[0]) + i) > 0) ?
                     (int)(*(mxGetPr(prhs[0]) + i) + 0.5)
                   : (int)(*(mxGetPr(prhs[0]) + i) - 0.5); /* Rounding to int */

    }

    /* Allocate memory fot the output */
    if (mxGetN(prhs[0])==1)
        plhs[0] = mxCreateDoubleMatrix(n_query,1,mxREAL);
    else
        plhs[0] = mxCreateDoubleMatrix(1,n_query,mxREAL);
    
    /* Pointer to output */
    out = mxGetPr(plhs[0]);

    /* Call the astro_constants function */
    error = astroConstants(n_query, in, out);
    
    if(nlhs>1)      /* error requested */
        plhs[1] = mxCreateDoubleScalar((double) error);

    
    /* Exit the function */
    mxFree(in);
    return;
}
