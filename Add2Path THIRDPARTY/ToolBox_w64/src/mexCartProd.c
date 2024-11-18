/******************************************************************************
 *                      Mex Gateway routine for cartprod                      *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

#include "mexCartProd.h"

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {

    /* Declaration */
    double *in;     /* array to hold the number of indexes for each variable
                           (input) */
    int *n;         /* array to hold the number of indexes for each variable */
    int *index;     /* Array to hold the cartesian product matrix */
    double *out;    /* arry to be returned to matlab (output) */
    int npar;       /* Munber of parameters */
    int nindex;     /* Number of index combinations */
    int i,j;        /* Counter */
    
    
    
    /* Check the number of inputs and outputs */
    if(nrhs != 1)
        mexErrMsgTxt("CARTPROD requires 1 and only 1 inputs.");
    if(nlhs > 1)
        mexErrMsgTxt("Too many outputs. CARTPROD returns only 1 output.");
    
    
    /* Get values of inputs */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])
                             || (mxGetM(prhs[0])!=1 && mxGetN(prhs[0])!=1)
                             || mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("The only input of CARTPROD must be a n-elements vector of integers!");
    npar = mxGetM(prhs[0])*mxGetN(prhs[0]);
    n = mxMalloc(npar*sizeof(int));
    in = mxGetPr(prhs[0]);
    
    /* Round the input to integers */
    for(i=0;i<npar;i++)
        n[i] = roundd2i(in[i]);
    
    
    /* Computes nindex, which is the product of the number of the bpas for each
     * parameter */
    nindex=1; /* Counter for rows in index matrix */
    for(i=0;i<npar;i++)
        nindex = nindex * n[i];
    
    /* Allocates memory for index[nindex][npar] */
    index = (int *)mxMalloc(nindex*npar*sizeof(int));
    

    /* Compute the outputs */
    if (cartProd_index_prealloc(npar, n, nindex, index))
        mexErrMsgTxt("CARTPROD: Out of memory!");
    
    
    /*Allocate memory and assign output pointer */
    plhs[0] = mxCreateDoubleMatrix(nindex,npar,mxREAL);
    out = mxGetPr(plhs[0]);
    
    /* Copy the values of index to put, and converting into doubles */
    for(i=0;i<nindex;i++)
        for(j=0;j<npar;j++)
            out[j*nindex+i] = (double) index[i*npar+j]+1;
    
    /* Free the allocated memory */
    mxFree(index);

    /* Exit the function */
    return;
}
