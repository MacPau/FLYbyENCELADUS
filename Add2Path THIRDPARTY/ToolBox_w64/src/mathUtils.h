/*******************************************************************************
*                                mathUtils.h                                   *
*                          C Mathematical utilities                            *
*                                                                              *
*                    Matteo Ceriotti, Nicolas Croisard, 2008                   *
*                                                                              *
*                              Space Toolbox                                   *
*******************************************************************************/

/* List of functions available:
*
*      1/ Constants
*         ---------
*   PI:  3.1415926535897931
*   PI2: 2*PI
*   RAD2DEG: 180.0/PI
*   DEG2RAD: PI/180.0
*
*
*      2/ Functions for scalar
*         --------------------
*   DSQR, FSQR: Square of double or float
*   DCUBE, FCUBE: Cube of a double or float
*   DMAX, FMAX, LMAX, IMAX: Maximum of two doubles, floats, longs or intergers
*   DMIN, FMIN, LMIN, IMIN: Minimum of two doubles, floats, longs or intergers
*   roundd2i, roundf2i: Round a double or float to the closest integer
*   fsign, dsign: Float, double signum function
*   DASINH, FASINH: Inverse of the hyperbolic sinus for a double, float
*   DACOSH, FACOSH: Inverse of the hyperbolic cossinus for a double, float
*   DATANH, FATANH: Inverse of the hyperbolic tangent for a double, float
*
*
*      3/ Functions for vector
*         --------------------
*   in_vector: Find the first occurence of an int in an int vector.
*   f2d: Copies a vector of floats into a vector of doubles.
*   d2f: Copies a vector of doubles into a vector of floats.
*   fscalarmult, dscalarmult: Multiplies all elements of a float, double vector
*           by a float, double scalar
*   fscalardiv, dscalardiv: Divides all elements of a float, double vector by a
*           float, double scalar.
*   max_ivector, max_fvector, max_dvector: Finds the maximum value among
*           elements of a integer, float or double vector
*   min_ivector, min_fvector, min_dvector: Finds the minimum value among
*           elements of a integer, float or double vector
*   fdot, ddot: Float, double vector dot product
*   fcross, dcross: Float, double vector cross product
*   fnorm, dnorm: Float, double vector norm
*
*
*      4/ Functions for matrix
*         --------------------
*   fmatmult, dmatmult: Float, double matrix multiply
*   fmattransp, dmattransp: Float, double matrix transposition
*
*
*      5/ Functions for vertor or matrix
*         ------------------------------
*   fcopy, dcopy: Copy matrices or vectors, float or double version
*   fsum, dsum: Sum of float, double elements
*   fmatsum, dmatsum: Elementwise float, double vector/matrix summation
*   fmatdiff, dmatdiff: Float, double elementwise vector/matrix subtraction
*   print_dmatrix, print_imatrix, print_fmatrix, print_cmatrix: Print a matrix
*           of double, integer, float, null-terminated strings on the screen
*
*
*      6/ Solvers
*         -------
*   fnewton, dnewton: - Newton solver with numerical derivative, float or
*           double version
*
*
*      7/ Miscellaneous
*         -------------
*   freeall: Free the allocated memory pointed by each pointer in argument list
*   cart_prod_index, cart_prod_index_prealloc: Generates the index matrix for
*                                              the cartesian product
*   
*
*/

#ifndef MATH_UTILS_H
#define	MATH_UTILS_H

/* Standard include */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <malloc.h>
#include <stdlib.h>

/* Personal include */

/* Mex compiler switch
* ------------------- */
#ifdef MEXCOMPILE
#include "matrix.h"
#include "mex.h"
#define printf mexPrintf
#define malloc mxMalloc
#define calloc mxCalloc
#define free mxFree
#endif
#ifndef MEXCOMPILE
#include <stdio.h>
#endif

/* Constants */
#define PI 3.1415926535897931
#define PI2 6.2831853071795862
#define RAD2DEG 180.0/PI
#define DEG2RAD PI/180.0

/* Inline definitions
   From: Press et al., "Numerical recipes in C, the art of scientific computing
   (2nd edition)" Appendix B (non-protected material).
*/

/* Square of a float */
static float fsqrarg;
#define FSQR(a) ((fsqrarg=(a)) == 0.0 ? 0.0 : (fsqrarg*fsqrarg))

/* Square of a double */
static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : (dsqrarg*dsqrarg))

/* Maximum of two doubles */
static double dmaxargl,dmaxarg2;
#define DMAX(a,b) (dmaxargl=(a),dmaxarg2=(b),(dmaxargl) > (dmaxarg2) ? (dmaxargl) : (dmaxarg2))

/* Minimum of two doubles */
static double dminargl,dminarg2;
#define DMIN(a,b) (dminargl=(a),dminarg2=(b),(dminargl) < (dminarg2) ?(dminargl) : (dminarg2))

/* Maximum of two floats */
static float maxargl,maxarg2;
#define FMAX(a,b) (maxargl=(a),maxarg2=(b),(maxargl) > (maxarg2) ? (maxargl) : (maxarg2))

/* Minimum of two floats */
static float minargl,minarg2;
#define FMIN(a,b) (minargl=(a),minarg2=(b),(minargl) < (minarg2) ? (minargl) : (minarg2))

/* Maximum of two long integers */
static long lmaxargl,Imaxarg2;
#define LMAX(a,b) (lmaxargl=(a),lmaxarg2=(b),(lmaxargl) > (Imaxarg2) ? (lmaxargl) : (Imaxarg2))

/* Minimum of two long integers */
static long lminargl,Iminarg2;
#define LMIN(a,b) (lminargl=(a),lminarg2=(b),(lminargl) < (Iminarg2) ? (lminargl) : (Iminarg2))

/* Maximum of two integers */
static int imaxargl,imaxarg2;
#define IMAX(a,b) (imaxargl=(a),imaxarg2=(b),(imaxargl) > (imaxarg2) ? (imaxargl) : (imaxarg2)) 

/* Minimum of two integers */
static int iminargl,iminarg2; 
#define IMIN(a,b) (iminargl=(a),iminarg2=(b),(iminargl) < (iminarg2) ? (iminargl) : (iminarg2))


/* Personal Inline definitions
*/

/* Cube of a float */
static float fcubearg;
#define FCUBE(a) ((fcubearg=(a)) == 0.0 ? 0.0 : (fcubearg*fcubearg*fcubearg))

/* Cube of a double */
static double dcubearg;
#define DCUBE(a) ((dcubearg=(a)) == 0.0 ? 0.0 : (dcubearg*dcubearg*dcubearg))

/* Inverse of the hyperbolic sinus for a double */
static double dasinharg;
#define DASINH(a) (dasinharg=(a), log(dasinharg+sqrt(DSQR(dasinharg)+1.0)))

/* Inverse of the hyperbolic sinus for a float */
static float fasinharg;
#define FASINH(a) (fasinharg=(a), log(fasinharg+sqrt(FSQR(fasinharg)+1.0)))

/* Inverse of the hyperbolic cosinus for a double */
static double dacosharg;
#define DACOSH(a) (dacosharg=(a), dacosharg<1.0 ? -HUGE_VAL : log(dacosharg+sqrt(DSQR(dacosharg)-1.0)))

/* Inverse of the hyperbolic cosinus for a float */
static float facosharg;
#define FACOSH(a) (facosharg=(a), facosharg<1.0 ? -HUGE_VAL : log(facosharg+sqrt(FSQR(facosharg)-1.0)))

/* Inverse of the hyperbolic tangent for a double */
static double datanharg;
#define DATANH(a) (datanharg=(a), datanharg>=-1.0 ? -HUGE_VAL : (datanharg<=1.0 ? HUGE_VAL : 0.5*log((1.0+datanharg)/(1.0+datanharg))))

/* Inverse of the hyperbolic tangent for a float */
static float fatanharg;
#define FATANH(a) (fatanharg=(a), fatanharg>=-1.0 ? -HUGE_VAL : (fatanharg<=1.0 ? HUGE_VAL : 0.5*log((1.0+fatanharg)/(1.0+fatanharg))))

/* Reduce an angle to the interval [0 2*PI[ */
static double qckarg1;
#define QCK(a) (qckarg1=(a), qckarg1<0 ? qckarg1 - PI2*(ceil(qckarg1/PI2)-1) : qckarg1 - PI2*floor(qckarg1/PI2))

/* Function PROTOTYPE:s */

int roundd2i(const double n);
/* roundd2i - Round a double to the closest integer
*
* DESCRIPTION:
*	No overflow check: the double must be in the integer range.
* 
* PROTOTYPE:
*   int roundd2i(const double n)
* 
* INPUT:
*   (double) n      Value to round.
*
* OUTPUT:
*   (int) roundd2i  Rounded value.
* 
* AUTHOR:
*   Matteo Ceriotti, 04/09/2008
* 
* CHANGELOG:
*/

int roundf2i(const float n);
/* roundf2i - Round a float to the closest integer
*
* DESCRIPTION:
*	No overflow check: the float must be in the integer range.
* 
* PROTOTYPE:
*   int roundf2i(const float n)
* 
* INPUT:
*   (float) n       Value to round.
*
* OUTPUT:
*   (int) roundf2i  Rounded value.
* 
* AUTHOR:
*   Matteo Ceriotti, 04/09/2008
* 
* CHANGELOG:
*/

int print_dmatrix(const double*mat, const int nrows, const int ncols);
/* print_dmatrix - Print a matrix of doubles on the screen.
* 
* PROTOTYPE:
*   int print_dmatrix(const double*mat, const int nrows, const int ncols)
* 
* INPUT:
*   (double*) mat           Pointer to the matrix.
*   (int) nrows             Number of rows.
*   (int) ncols             Number of columns.
*   
* OUTPUT:
*   (int) print_dmatrix     Error code (0).
* 
* AUTHOR:
*   Matteo Ceriotti, 04/09/2008
* 
* CHANGELOG:
*/

int print_imatrix(const int*mat, const int nrows, const int ncols);
/* print_imatrix - Print a matrix of integers on the screen.
* 
* PROTOTYPE:
*   int print_imatrix(const int*mat, const int nrows, const int ncols)
* 
* INPUT:
*   (int*) mat              Pointer to the matrix.
*   (int) nrows             Number of rows.
*   (int) ncols             Number of columns.
*   
* OUTPUT:
*   (int) print_imatrix     Error code (0).
* 
* AUTHOR:
*   Matteo Ceriotti, 04/09/2008
* 
* CHANGELOG:
*/

int print_fmatrix(const float*mat, const int nrows, const int ncols);
/* print_fmatrix - Print a matrix of floats on the screen.
* 
* PROTOTYPE:
*   int print_fmatrix(const float*mat, const int nrows, const int ncols)
* 
* INPUT:
*   (float*) mat            Pointer to the matrix.
*   (int) nrows             Number of rows.
*   (int) ncols             Number of columns.
*   
* OUTPUT:
*   (int) print_fmatrix     Error code (0).
* 
* AUTHOR:
*   Matteo Ceriotti, 04/09/2008
* 
* CHANGELOG:
*/

int print_cmatrix(const char*mat, const int nrows, const int ncols);
/* print_cmatrix - Print a matrix of null-terminated strings on the screen.
* 
* PROTOTYPE:
*   int print_cmatrix(const char*mat, const int nrows, const int ncols)
* 
* INPUT:
*   (char*) mat             Pointer to the matrix.
*   (int) nrows             Number of rows (number of strings).
*   (int) ncols             Number of columns (max length of each string).
*   
* OUTPUT:
*   (int) print_cmatrix     Error code (0).
* 
* AUTHOR:
*   Matteo Ceriotti, 04/09/2008
* 
* CHANGELOG:
*/

int in_vector(int a, int*v, int n);
/* in_vector - Finds the first occurence of int a in int vector v.
* 
* PROTOTYPE:
*   int in_vector(int a, int*v, int n)
* 
* INPUT:
*   (int) a            Number to look for.
*   (int*) v           Vector to be searched.
*   (int) n            Length of vector v.
*   
* OUTPUT:
*   (int) in_vector    Index of the first occurrence of a in v.
* 
* AUTHOR:
*   Matteo Ceriotti, 04/09/2008
* 
* CHANGELOG:
*/

float max_fvector(const float*v, const int n);
/* max_fvector - Finds the maximum value among elements of float vector v.
* 
* PROTOTYPE:
*   float max_fvector(const float*v, const int n)
* 
* INPUT:
*   (float*) v      Pointer to the vector of floats.
*   (int) n         Length of vector v.
*   
* OUTPUT:
*   (float) max_fvector     Maximum.
* 
* AUTHOR:
*   Matteo Ceriotti, 04/09/2008
* 
* CHANGELOG:
*/

float min_fvector(const float*v, int n);
/* min_fvector - Finds the minimum value among elements of float vector v.
*  
* PROTOTYPE:
*   float min_fvector(const float*v, const int n)
* 
* INPUT:
*   (float*) v      Pointer to the vector of floats.
*   (int) n         Length of vector v.
*   
* OUTPUT:
*   (float) min_fvector     Minimum.
* 
* AUTHOR:
*   Matteo Ceriotti, 04/09/2008
* 
* CHANGELOG:
*/

double max_dvector(const double*v, const int n);
/* max_dvector - Finds the maximum value among elements of double vector v.
* 
* PROTOTYPE:
*   double max_fvector(const double*v, const int n)
* 
* INPUT:
*   (double*) v     Pointer to the vector of doubles.
*   (int) n         Length of vector v.
*   
* OUTPUT:
*   (double) max_dvector     Maximum.
* 
* AUTHOR:
*   Matteo Ceriotti, 04/09/2008
* 
* CHANGELOG:
*/

double min_dvector(const double*v, int n);
/* min_dvector - Finds the minimum value among elements of double vector v.
* 
* PROTOTYPE:
*   double min_dvector(const double*v, const int n)
* 
* INPUT:
*   (double*) v     Pointer to the vector of doubles.
*   (int) n         Length of vector v.
*   
* OUTPUT:
*   (double) min_dvector     Minimum.
* 
* AUTHOR:
*   Matteo Ceriotti, 04/09/2008
* 
* CHANGELOG:
*/

int max_ivector(const int*v, const int n);
/* max_ivector - Finds the maximum value among elements of integer vector v.
* 
* PROTOTYPE:
*   int max_fvector(const int*v, const int n)
* 
* INPUT:
*   (int*) v        Pointer to the vector of integers.
*   (int) n         Length of vector v.
*   
* OUTPUT:
*   (itn) max_ivector     Maximum.
* 
* AUTHOR:
*   Matteo Ceriotti, 04/09/2008
* 
* CHANGELOG:
*/

int min_ivector(const int*v, int n);
/* min_ivector - Finds the minimum value among elements of integer vector v.
* 
* PROTOTYPE:
*   int min_dvector(const int*v, const int n)
* 
* INPUT:
*   (int*) v        Pointer to the vector of integers.
*   (int) n         Length of vector v.
*   
* OUTPUT:
*   (int) min_ivector     Minimum.
* 
* AUTHOR:
*   Matteo Ceriotti, 04/09/2008
* 
* CHANGELOG:
*/

void freeall(int num, ...);
/* freeall - Free the allocated memory pointed by each pointer in argument
*   list.
* 
* PROTOTYPE:
*   void freeall(int num, ...)
* 
* INPUT:
*   (int) num       Number of variable arguments.
*   ...             Pointers to allocated memory.
* 
* OUTPUT:
*   none.
* 
* AUTHOR:
*   Matteo Ceriotti, 04/09/2008
* 
* CHANGELOG:
*/

int f2d(const float*a, double*b, int size);
/* f2d - Copies a vector of float into a vector of doubles.
* 
* PROTOTYPE:
*   void f2d(const float*a, double*b, int size)
* 
* INPUT:
*   (float*) v		Pointer to the vector of floats (source).
*   (double*) b     Pointer to the vector of doubles (destination).
*   (int) size      Length of vectors v and b. Must be preallocated.
*   
* OUTPUT:
*   (int) f2d       Error code (0).
* 
* AUTHOR:
*   Matteo Ceriotti, 04/09/2008
* 
* CHANGELOG:
*/

int d2f(const double*a, float*b, int size);
/* f2d - Copies a vector of float into a vector of doubles.
* 
* PROTOTYPE:
*   int d2f(const double*a, float*b, int size)
* 
* INPUT:
*   (double*) v     Pointer to the vector of doubles (source).
*   (float*) b      Pointer to the vector of floats (destination).
*   (int) size      Length of vectors v and b. Must be preallocated.
*   
* OUTPUT:
*   (int) d2f       Error code (0).
* 
* AUTHOR:
*   Matteo Ceriotti, 04/09/2008
* 
* CHANGELOG:
*/

/*************************************************/

void fscalarmult(const float*a, const float b, const int size, float*c);
/* fscalarmult - Multiplies all elements of a float vector by a float scalar.
* 
* PROTOTYPE:
*   void fscalarmult(const float*a, const float b, const int size,
*        float*c)
*
* INPUT:
*   float*a     vector to be multiplied.
*   float b     scalar to be multiplied.
*   int size    size of vectors a and c.
*
* OUTPUT:
*   float*c    resulting vector (pre-allocated memory).
* 
* AUTHOR:
*   Matteo Ceriotti, 2006
*/

int fscalardiv(const float*a, const float b, const int size, float*c);
/* fscalardiv - Divides all elements of a float vector by a float scalar.
* 
* PROTOTYPE:
*   int fscalardiv(const float*a, const float b, const int size, float*c)
*
* INPUT:
*   float*a     vector to be divided.
*   float b     scalar which divides.
*   int size    size of vectors a and c.
*
* OUTPUT:
*   float*c          resulting vector (pre-allocated memory).
*   int scalardiv    0 if worked correctly, else 1.
*
* AUTHOR:
*   Matteo Ceriotti, 2006
*/

void fmatmult(const float*a, const int ax, const int ay, const int by, const float*b, float*c);
/* fmatmult - Float matrix multiply
*
* DESCRIPTION:
*	Multiplies a vector or matrix a (size ax*ay) by a vector or matrix b (size
*	ay*by) giving the vector or matrix c (size ax*by).
*
* PROTOTYPE:
*   void fmatmult(const float*a, const int ax, const int ay, const int by,
*        const float*b, float*c)
* INPUT:
*   float*a    first factor.
*   int ax     number of rows of a.
*   int ay     number of columns of a, and rows of b.
*   int by     number of columns of b.
*   float*b    second factor.
*
* OUTPUT:
*   float*c          resulting vector or matrix (size ax* by, pre-allocated).
*   int scalardiv    0 if worked correctly, else 1.
*
* AUTHOR:
*   Matteo Ceriotti, 2006
*/

void fmatsum(const float*a, const float*b, const int size, float*c);
/* fmatsum - Elementwise float vector/matrix summation
*
* DESCRIPTION:
*	Sums two matrices or vectors elementwise: c = a + b.
*
* PROTOTYPE:
*   void fmatsum(const float*a, const float*b, const int size, float*c)
*
* INPUT:
*   float*a,*b    matrices to be summed.
*   int size      size of matrices a and b (number of rows* number of cols).
*
* OUTPUT:
*   float*c       resulting matrix (pre-allocated memory).
*
* AUTHOR:
*   Matteo Ceriotti, 2006
*/

void fmatdiff(const float*a, const float*b, const int size, float*c);
/* fmatdiff - Float elementwise vector/matrix subtraction
*
* DESCRIPTION:
*	Subtracts two matrices or vectors elementwise: c = a - b.
*
* PROTOTYPE:
*   void fmatdiff(const float*a, const float*b, const int size, float*c)
*
* INPUT:
*   float*a,*b    matrices to be subtracted.
*   int size      size of matrices a and b (number of rows* number of cols).
*
* OUTPUT:
*   float*c       resulting matrix (pre-allocated memory).
*
* AUTHOR:
*   Matteo Ceriotti, 2006
*/

void fmattransp(const float*a, const int x, const int y, float*b);
/* fmattransp - Float matrix transposition.
*
* DESCRIPTION:
*	Returns the transposed matrix: b = a'.
*
* PROTOTYPE:
*   void fmattransp(const float*a, const int x, const int y, float*b)
*
* INPUT:
*   float*a     matrix to be transposed (size x*y).
*   int x, y    dimensions of matrix a.
*
* OUTPUT:
*   float*b     transposed matrix (size y*x), memory pre-allocated.
*
* AUTHOR:
*   Matteo Ceriotti, 2006
*/

void fcross(const float*a, const float*b, float*c);
/* fcross - Float vector cross product.
*
* DESCRIPTION:
*	Returns the cross product of the vectors a and b.
*	That is, c = a x b. a and b must be 3 element vectors.
* 
* PROTOTYPE:
*   void fcross(const float*a, const float*b, float*c)
*
* INPUT:
*   float*a,*b    input vectors (each of size 3).
*
* OUTPUT:
*   float*c       resulting vector (c = a X b).
*
* AUTHOR:
*   Matteo Ceriotti, 2006
*/

float fdot(const float*a, const float*b, const int size);
/* fdot - Float vector dot product.
*
* DESCRIPTION:
*	Returns the scalar product of the vectors a and b.
*	a and b must be vectors both with size size.
*
* PROTOTYPE:
*   float fdot(const float*a, const float*b, const int size)
*
* INPUT:
*   float*a,*b    vectors.
*   int size      size of both a and b.
*
* OUTPUT:
*   float dot    dot product.
*
* AUTHOR:
*   Matteo Ceriotti, 2006
*/

int fsign(const float n);
/* fsign - Float signum function
*
* DESCRIPTION:
*	Returns 1 if the input is greater than zero, 0 if it equals zero and -1 if
*	it is less than zero.
*
* PROTOTYPE:
*   int fsign(const float n)
*
* INPUT:
*   float n    Number.
*
* OUTPUT:
*   int sign   Sign.
*
* AUTHOR:
*   Matteo Ceriotti, 2006
*/

float fnorm(const float*a, const int size);
/* fnorm - Float vector norm
*
* DESCRIPTION:
*	Returns sum(a.^2)^(1/2).
*
* PROTOTYPE:
*   float fnorm(const float*a, const int size)
*
* INPUT:
*   float*a     Vector.
*   int size    size of a.
*
* OUTPUT:
*   float norm	Norm.
*
* AUTHOR:
*   Matteo Ceriotti, 2006
*
* CHANGELOG:
*	18/02/2009, Matteo Ceriotti: Added explicit cast to float to avoid
*		warnings.
*/

int fnewton(const float x0, const int nmax, const float toll,
    const float deltax, int (*fun)(const float x, void**pointers, float*out),
    void**pointers, float*xout);
/* fnewton - Newton solver with numerical derivative, float version
*
* DESCRIPTION:
*	Finds the roots x for the scalar equation f(x) = 0, using a Newton method
*	and evaluating the derivative of fun numerically, with an incremental ratio.
*
* PROTOTYPE:
*   int fnewton(const float x0, const int nmax, const float toll,
*       const float deltax, int (*fun)(const float x, void**pointers,
*       float*out), void**pointers, float*xout)
*
* INPUT:
*   float x0        starting point.
*   int nmax        maximum number of iterations.
*   float toll      tolerance required.
*   float deltax	increment used to compute the first derivative
*                   numerically.
*   int (*fun)      handle to function fun.
*
* OUTPUT:
*   float*xout		root of fun.
*   int newton		0 if worked succesfully;
*					1 if derivative becomes 0;
*					2 if not converged in nmax iterations;
*					3 if error in evaluating fun.
*
* AUTHOR:
*   Matteo Ceriotti, 2006
*
* CHANGELOG:
*	18/02/2009, Matteo Ceriotti: Added explicit cast to float to avoid
*		warnings.
*/

void fcopy(const float*a, float*b, const int size);
/* fcopy - Copy matrices or vectors, float version
*
* DESCRIPTION:
*	Copy the values in array a in array b. That is the matlab assignment b = a.
*
* PROTOTYPE:
*   void fcopy(const float*a, float*b, const int size)
*
* INPUT:
*   float*a		vector or matrix to be copied.
*   int size    size of a and b.
*
* OUTPUT:
*   float*b		vector or matrix that will be equal to a[] (preallocated).
*
* AUTHOR:
*   Matteo Ceriotti, 2006
*/

float fsum(const float*a, const int size);
/* fsum - Sum of float elements
*
* DESCRIPTION:
*	Sums all the elements of a vector or matrix a.
* 
* PROTOTYPE:
*   float fsum(const float*a, const int size)
*
* INPUT:
*   float*a		vector or matrix whose elements will be summed.
*   int size    size of a.
*
* OUTPUT:
*   float*sum	value of the sum.
*
* AUTHOR:
*   Matteo Ceriotti, 2006
*/

/****************************************************/

void dscalarmult(const double*a, const double b, const int size, double*c);
/* dscalarmult - Multiplies all elements of a double vector by a double scalar.
* 
* PROTOTYPE:
*   void dscalarmult(const double*a, const double b, const int size,
*        double*c)
*
* INPUT:
*   double*a    vector to be multiplied.
*   double b    scalar to be multiplied.
*   int size    size of vectors a and c.
*
* OUTPUT:
*   double*c    resulting vector (pre-allocated memory).
* 
* AUTHOR:
*   Matteo Ceriotti, 2006
*/


int dscalardiv(const double*a, const double b, const int size, double*c);
/* dscalardiv - Divides all elements of a double vector by a double scalar.
* 
* PROTOTYPE:
*   int dscalardiv(const double*a, const double b, const int size, double*c)
*
* INPUT:
*   double*a    vector to be divided.
*   double b    scalar which divides.
*   int size    size of vectors a and c.
*
* OUTPUT:
*   double*c		resulting vector (pre-allocated memory).
*   int scalardiv   0 if worked correctly, else 1.
*
* AUTHOR:
*   Matteo Ceriotti, 2006
*/


void dmatmult(const double*a, const int ax, const int ay, const int by, const double*b, double*c);
/* dmatmult - Double matrix multiply
*
* DESCRIPTION:
*	Multiplies a vector or matrix a (size ax*ay) by a vector or matrix b (size
*	ay*by) giving the vector or matrix c (size ax*by).
*
* PROTOTYPE:
*   void dmatmult(const double*a, const int ax, const int ay, const int by,
*        const double*b, double*c)
*
* INPUT:
*   double*a	first factor.
*   int ax		number of rows of a.
*   int ay		number of columns of a, and rows of b.
*   int by		number of columns of b.
*   double*b    second factor.
*
* OUTPUT:
*   double*c		resulting vector or matrix (size ax* by, pre-allocated).
*   int scalardiv   0 if worked correctly, else 1.
*
* AUTHOR:
*   Matteo Ceriotti, 2006
*/


void dmatsum(const double*a, const double*b, const int size, double*c);
/* dmatsum - Elementwise double vector/matrix summation
*
* DESCRIPTION:
*	Sums two matrices or vectors elementwise: c = a + b.
*
* PROTOTYPE:
*   void dmatsum(const double*a, const double*b, const int size, double*c)
*
* INPUT:
*   double*a,*b		matrices to be summed.
*   int size		size of matrices a and b (number of rows* number of cols).
*
* OUTPUT:
*   double*c		resulting matrix (pre-allocated memory).
*
* AUTHOR:
*   Matteo Ceriotti, 2006
*/


void dmatdiff(const double*a, const double*b, const int size, double*c);
/* dmatdiff - Double elementwise vector/matrix subtraction
*
* DESCRIPTION:
*	Subtracts two matrices or vectors elementwise: c = a - b.
*
* PROTOTYPE:
*   void dmatdiff(const double*a, const double*b, const int size, double*c)
*
* INPUT:
*   double*a,*b		matrices to be subtracted.
*   int size		size of matrices a and b (number of rows* number of cols).
*
* OUTPUT:
*   double*c		resulting matrix (pre-allocated memory).
*
* AUTHOR:
*   Matteo Ceriotti, 2006
*/


void dmattransp(const double*a, const int x, const int y, double*b);
/* dmattransp - Double matrix transposition.
*
* DESCRIPTION:
*	Returns the transposed matrix: b = a'.
*
* PROTOTYPE:
*   void dmattransp(const double*a, const int x, const int y, double*b)
*
* INPUT:
*   double*a    matrix to be transposed (size x*y).
*   int x, y    dimensions of matrix a.
*
* OUTPUT:
*   double*b    transposed matrix (size y*x), memory pre-allocated.
*
* AUTHOR:
*   Matteo Ceriotti, 2006
*/


void dcross(const double*a, const double*b, double*c);
/* dcross - Double vector cross product.
*
* DESCRIPTION:
*	Returns the cross product of the vectors a and b.
*	That is, c = a x b. a and b must be 3 element vectors.
* 
* PROTOTYPE:
*   void dcross(const double*a, const double*b, double*c)
*
* INPUT:
*   double*a,*b		input vectors (each of size 3).
*
* OUTPUT:
*   double*c		resulting vector (c = a X b).
*
* AUTHOR:
*   Matteo Ceriotti, 2006
*/


double ddot(const double*a, const double*b, const int size);
/* ddot - Double vector dot product.
*
* DESCRIPTION:
*	Returns the scalar product of the vectors a and b.
*	a and b must be vectors both with size size.
*
* PROTOTYPE:
*   double ddot(const double*a, const double*b, const int size)
*
* INPUT:
*   double*a,*b		Vectors.
*   int size		size of both a and b.
*
* OUTPUT:
*   double dot		Dot product.
*
* AUTHOR:
*   Matteo Ceriotti, 2006
*/


int dsign(const double n);
/* dsign - Double signum function
*
* DESCRIPTION:
*	Returns 1 if the input is greater than zero, 0 if it equals zero and -1 if
*	it is less than zero.
*
* PROTOTYPE:
*   int dsign(const double n)
*
* INPUT:
*   double n    Number.
*
* OUTPUT:
*   int sign    Sing.
*
* AUTHOR:
*   Matteo Ceriotti, 2006
*/


double dnorm(const double*a, const int size);
/* dnorm - Double vector norm
*
* DESCRIPTION:
*	Returns sum(a.^2)^(1/2).
*
* PROTOTYPE:
*   double dnorm(const double*a, const int size)
*
* INPUT:
*   double*a		Vector.
*   int size		size of a.
*
* OUTPUT:
*   double norm		Norm.
*
* AUTHOR:
*   Matteo Ceriotti, 2006
*/


int dnewton(const double x0, const int nmax, const double toll,
    const double deltax, int (*fun)(const double x, void**pointers, double*out),
    void**pointers, double*xout);
/* dnewton - Newton solver with numerical derivative, double version
*
* DESCRIPTION:
*	Finds the roots x for the scalar equation f(x) = 0, using a Newton method
*	and evaluating the derivative of fun numerically, with an incremental ratio.
*
* PROTOTYPE:
*   int dnewton(const double x0, const int nmax, const double toll,
*       const double deltax, int (*fun)(const double x, void**pointers,
*       double*out), void**pointers, double*xout)
*
* INPUT:
*   double x0		starting point.
*   int nmax		maximum number of iterations.
*   double toll		tolerance required.
*   double deltax   increment used to compute the first derivative
*					numerically.
*   int (*fun)		handle to function fun.
*
* OUTPUT:
*   double*xout		root of fun.
*   int newton		0 if worked succesfully;
*					1 if derivative becomes 0;
*					2 if not converged in nmax iterations;
*					3 if error in evaluating fun.
*
* AUTHOR:
*   Matteo Ceriotti, 2006
*/


void dcopy(const double*a, double*b, const int size);
/* dcopy - Copy matrices or vectors, double version
*
* DESCRIPTION:
*	Copy the values in array a in array b. That is the matlab assignment b = a.
*
* PROTOTYPE:
*   void dcopy(const double*a, double*b, const int size)
*
* INPUT:
*   double*a    vector or matrix to be copied.
*   int size    size of a and b.
*
* OUTPUT:
*   double*b    vector or matrix that will be equal to a[] (preallocated).
*
* AUTHOR:
*   Matteo Ceriotti, 2006
*/


double dsum(const double*a, const int size);
/* dsum - Sum of double elements
*
* DESCRIPTION:
*	Sums all the elements of a vector or matrix a.
* 
* PROTOTYPE:
*   double dsum(const double*a, const int size)
*
* INPUT:
*   double*a		vector or matrix whose elements will be summed.
*   int size		size of a.
*
* OUTPUT:
*   double*sum		value of the sum.
*
* AUTHOR:
*   Matteo Ceriotti, 2006
*/

int cartProd_index(int npar, int*n, int**pp_index, int*nindex);
/* cartProd_index - Generates the index matrix for the cartesian product
*
* DESCRIPTION:
*	This function dynamically allocates memory needed for the output matrix
*	index. This memory should be freed somewhere outside.
*
* PROTOTYPE:
*   int cartProd_index(int npar, int*n, int**pp_index, int*nindex)
*
* INPUT:
*   (int) npar              Number of parameters (length of vector n), and
*                           number of columns of OUTPUT:*pp_index matrix.
*   (int*) n                Pointer to a vector containing the number of
*                           elements given for each parameter.
*
* OUTPUT:
*   (int**) pp_index        Pointer to pointer to the index matrix for the
*                           cartesian product; the matrix is dynamically
*                           allocated inside the function. The function
*                           modifies the second pointer to point to the
*                           correct memory area. The size of this matrix is
*                           (*nindex)*npar.
*                           For example, if:
*                              npar = 3; n[3] = {2, 2, 3};
*                           The output is:
*                                          0 0 0
*                                          0 0 1
*                                          0 0 2
*                                          0 1 0
*                                          0 1 1
*                             *pp_index = 0 1 2      *nindex = 12
*                                          1 0 0
*                                          1 0 1
*                                          1 0 2
*                                          1 1 0
*                                          1 1 1
*                                          1 1 2
*   (int*) nindex           Pointer to nindex, which is the number of rows
*                           of the matrix*pp_index.
*   (int) cartProd_index    Error code: 1 if out of memory,
*                                       0 otherwise
*
* AUTHOR:
*   Matteo Ceriotti, 2005
*/

int cartProd_index_prealloc(const int npar, const int*n, const int nindex, int index[]);
/* cartProd_index_prealloc - Generates the index matrix for the cartesian
*                          product
*
* DESCRIPTION:
*	This function assume a preallocated output (on the contrary of
*	cart_prod_index)
*
* PROTOTYPE:
*   int cartProd_index_prealloc(int npar, int*n, int**pp_index, int*nindex)
*
* INPUT:
*   (int) npar      Number of parameters (length of vector n), and number of
*                   columns of OUTPUT:*pp_index matrix.
*   (int*) n        Pointer to a vector containing the number of elements
*                   given for each parameter.
*   (int) nindex    Integer  the number of rows of the matrix*pp_index.
*
* OUTPUT:
*   (int*) index    Pointer to the index matrix for the cartesian product.
*                   index must be preallocated. It size is nindex*npar.
*                       For example, if: npar = 3; n[3] = {2, 2, 3};
*                       The output is:
*                                          0 0 0
*                                          0 0 1
*                                          0 0 2
*                                          0 1 0
*                                          0 1 1
*                             *pp_index = 0 1 2      *nindex = 12
*                                          1 0 0
*                                          1 0 1
*                                          1 0 2
*                                          1 1 0
*                                          1 1 1
*                                          1 1 2
*
*   (int) cartProd_index_prealloc		Error code (0)
*
* AUTHOR:
*   Matteo Ceriotti, 2005
*/

#endif /* MATH_UTILS_H*/
