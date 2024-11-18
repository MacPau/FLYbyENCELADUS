/******************************************************************************
 *                       Mex Gateway routine for date2jd                      *
 *                                                                            *
 *                                                     Nicolas Croisard, 2008 *
 ******************************************************************************/

#include "mexDate2jd.h"

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {

    /* Declaration */
    double *date;       /* 6-element array of double containing the date in
                           [Y M D h m s] to be converted (input) */
    double jd;          /* Date expressed in JD (output) */
    int error;          /* error code */
    mwSize nr, nc;      /* number of rows and columns of input */

    
    /* Check the number of inputs and outputs */
    if(nrhs != 1)
        mexErrMsgTxt("DATE2JD requires 1 and only 1 inputs.");
    if(nlhs > 1)
        mexErrMsgTxt("Too many outputs. DATE2JD returns only 1 output.");
    
    
    /* Get values of inputs */
    nr = mxGetM(prhs[0]);
    nc = mxGetN(prhs[0]);
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])
                             || (!(nr==6 && nc==1) && !(nr==1 && nc==6))
                             || mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("The input of DATE2JD must be a 6-elements vector of reals!");
    date = mxGetPr(prhs[0]);
    
    
    /* Compute the output */
    error = date2jd(date,&jd);
    
    
    /* Check for error */
    if (error>0)
        mexErrMsgTxt("DATE2JD: The function is valid for dates after since 12:00 noon 24 November -4713, Gregorian calendar");
    
    
    /* Return the output */
    plhs[0] = mxCreateDoubleScalar(jd);

    

    /* Exit the function */
    return;
}
