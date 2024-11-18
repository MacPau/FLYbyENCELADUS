/*******************************************************************************
*                                mathUtils.h                                   *
*                          C Mathematical utilities                            *
*                                                                              *
*                    Matteo Ceriotti, Nicolas Croisard, 2008                   *
*                                                                              *
*                              Space Toolbox                                   *
*******************************************************************************/

#include "mathUtils.h"

int roundd2i(const double n){
    return (n > 0) ? (int)(n + 0.5) : (int)(n - 0.5); /* Rounding to int */
}

int roundf2i(const float n){
    return (n > 0) ? (int)(n + 0.5) : (int)(n - 0.5); /* Rounding to int */
}

int print_dmatrix(const double *mat, const int nrows, const int ncols){
	int i, j;
    for(i=0; i<nrows; i++){
        for(j=0; j<ncols; j++){
            printf(" %.15e", *(mat+ncols*i+j));
        }
        printf("\n");
    }
    return 0;
}

int print_imatrix(const int *mat, const int nrows, const int ncols){
	int i, j;
    for(i=0; i<nrows; i++){
        for(j=0; j<ncols; j++){
            printf(" %4.0f",(float)*(mat+ncols*i+j));
        }
        printf("\n");
    }
    return 0;
}

int print_fmatrix(const float *mat, const int nrows, const int ncols){
	int i,j;
    for(i=0;i<nrows;i++){
        for(j=0;j<ncols;j++){
            printf(" %.8f", *(mat+ncols*i+j));
        }
		printf("\n");
    }
    return 0;
}

int print_cmatrix(const char *mat, const int nrows, const int ncols){
	int i;
	for(i=0;i<nrows;i++){
		printf("%s ",mat+ncols*i);
	}
	printf("\n");
	return 0;
}

int in_vector(int a, int *v, int n){
	int i;
	for(i=0;i<n;i++){
		if(v[i]==a){
			return i;
		}
	}
	return -1;
}

float max_fvector(const float *v, const int n){
	int i;
	float max=v[0];
	for(i=1;i<n;i++){
		if(v[i]>max) max=v[i];
	}
	return max;
}

float min_fvector(const float *v, int n){
	int i;
	float min=v[0];
	for(i=1;i<n;i++){
		if(v[i]<min) min=v[i];
	}
	return min;
}

double max_dvector(const double *v, const int n){
	int i;
	double max=v[0];
	for(i=1;i<n;i++){
		if(v[i]>max) max=v[i];
	}
	return max;
}

double min_dvector(const double *v, int n){
	int i;
	double min=v[0];
	for(i=1;i<n;i++){
		if(v[i]<min) min=v[i];
	}
	return min;
}

int max_ivector(const int *v, const int n){
	int i;
	int max=v[0];
	for(i=1;i<n;i++){
		if(v[i]>max) max=v[i];
	}
	return max;
}

int min_ivector(const int *v, int n){
	int i;
	int min=v[0];
	for(i=1;i<n;i++){
		if(v[i]<min) min=v[i];
	}
	return min;
}

void freeall(int num, ...){
	va_list argptr;
	
	va_start(argptr,num);
	
	for(;num;num--){
		free( va_arg(argptr, void *) );
	}
	
	va_end(argptr);
}

int f2d(const float *a, double *b, int size){
    for(;size;size--){
        b[size-1]=(double)a[size-1];
    }
    return 0;
}

int d2f(const double *a, float *b, int size){
    for(;size;size--){
        b[size-1]=(float)a[size-1];
    }
    return 0;
}

/**************************/

void fscalarmult(const float *a, const float b, const int size, float *c){
    int i;
    for(i=0;i<size;i++){
        c[i]=a[i]*b;
    }
}

int fscalardiv(const float *a, const float b, const int size, float *c){
    int i;
    if(b==0.){
        return 1;
    }
    for(i=0;i<size;i++){
        c[i]=a[i]/b;
    }
    return 0;
}

void fmatmult(const float *a, const int ax, const int ay, const int by, const float *b, float *c){
    int i,j,k;
    for(i=0;i<ax;i++){
        for(j=0;j<by;j++){
            c[by*i+j]=0;
            for(k=0;k<ay;k++){
                c[by*i+j]+=a[ay*i+k]*b[by*k+j];
            }
        }
    }
}

void fmatsum(const float *a, const float *b, const int size, float *c){
    int i;
    for(i=0;i<size;i++){
        c[i]=a[i]+b[i];
    }
}

void fmatdiff(const float *a, const float *b, const int size, float *c){
    int i;
    for(i=0;i<size;i++){
        c[i]=a[i]-b[i];
    }
}

void fmattransp(const float *a, const int x, const int y, float *b){
    int i,j;
    for(i=0;i<x;i++){
        for(j=0;j<y;j++){
            b[j*x+i]=a[i*y+j];
        }
    }
}

void fcross(const float *a, const float *b, float *c){
    c[0]=a[1]*b[2]-a[2]*b[1];
    c[1]=a[2]*b[0]-a[0]*b[2];
    c[2]=a[0]*b[1]-a[1]*b[0];
    return;
}

float fdot(const float *a, const float *b, const int size){
    int i;
    float dot=0;
    for(i=0;i<size;i++){
        dot+=a[i]*b[i];
    }
    return dot;
}

int fsign(const float n){
    if(n>0) return 1;
    else if(n==0) return 0;
    else return -1;
}

float fnorm(const float *a, const int size){
    int i;
    float out=0;
    for(i=0;i<size;i++){
        out+=(float)DSQR(a[i]);
    }
    out=(float)sqrt(out);
    return out;
}

int fnewton(const float x0, const int nmax, const float toll,
    const float deltax, int (*fun)(const float x, void **pointers, float *out),
    void **pointers, float *xout)
{
    float err,fx,dfx,fx1,fx2,xn;
    float x,x1,x2;
    int nit;
    x=x0;
    if(fun(x,pointers,&fx)) return 3; /* Evaluates fx=fun(x) */
    err=toll+1;
    for(nit=1;nit<=nmax;nit++){
        if(err<=toll){
            *xout=x;
            return 0;
        }
        x1=x-deltax/2;
        x2=x+deltax/2;
        
        if(fun(x1,pointers,&fx1)) return 3; /* Evaluates fx1=fun(x1) */
        if(fun(x2,pointers,&fx2)) return 3; /* Evaluates fx2=fun(x2) */
        dfx=(fx2-fx1)/deltax;
        if(dfx==0){
            return 1;
        }
        else{
            xn=x-fx/dfx;
            err=(float)fabs(xn-x);
            x=xn;
            if(fun(x,pointers,&fx)) return 3; /* Evaluates fx=fun(x) */
        }
    }
    /* Too many itarations */
    return 2;
}

void fcopy(const float *a, float *b, const int size){
    int i;
    for(i=0;i<size;i++){
        b[i]=a[i];
    }
}

float fsum(const float *a, const int size){
    int i;
    float s=0;
    for(i=0;i<size;i++){
        s+=a[i];
    }
    return s;
}

/**************************/

void dscalarmult(const double *a, const double b, const int size, double *c){
    int i;
    for(i=0;i<size;i++){
        c[i]=a[i]*b;
    }
}

int dscalardiv(const double *a, const double b, const int size, double *c){
    int i;
    if(b==0.){
        return 1;
    }
    for(i=0;i<size;i++){
        c[i]=a[i]/b;
    }
    return 0;
}

void dmatmult(const double *a, const int ax, const int ay, const int by, const double *b, double *c){
    int i,j,k;
    for(i=0;i<ax;i++){
        for(j=0;j<by;j++){
            c[by*i+j]=0;
            for(k=0;k<ay;k++){
                c[by*i+j]+=a[ay*i+k]*b[by*k+j];
            }
        }
    }
}

void dmatsum(const double *a, const double *b, const int size, double *c){
    int i;
    for(i=0;i<size;i++){
        c[i]=a[i]+b[i];
    }
}

void dmatdiff(const double *a, const double *b, const int size, double *c){
    int i;
    for(i=0;i<size;i++){
        c[i]=a[i]-b[i];
    }
}

void dmattransp(const double *a, const int x, const int y, double *b){
    int i,j;
    for(i=0;i<x;i++){
        for(j=0;j<y;j++){
            b[j*x+i]=a[i*y+j];
        }
    }
}

void dcross(const double *a, const double *b, double *c){
    c[0]=a[1]*b[2]-a[2]*b[1];
    c[1]=a[2]*b[0]-a[0]*b[2];
    c[2]=a[0]*b[1]-a[1]*b[0];
    return;
}

double ddot(const double *a, const double *b, const int size){
    int i;
    double dot=0;
    for(i=0;i<size;i++){
        dot+=a[i]*b[i];
    }
    return dot;
}

int dsign(const double n){
    if(n>0) return 1;
    else if(n==0) return 0;
    else return -1;
}

double dnorm(const double *a, const int size){
    int i;
    double out=0;
    for(i=0;i<size;i++){
        out+=pow(a[i],2);
    }
    out=sqrt(out);
    return out;
}

int dnewton(const double x0, const int nmax, const double toll,
    const double deltax, int (*fun)(const double x, void **pointers, double *out),
    void **pointers, double *xout)
{
    double err,fx,dfx,fx1,fx2,xn;
    double x,x1,x2;
    int nit;
    x=x0;
    if(fun(x,pointers,&fx)) return 3; /* Evaluates fx=fun(x) */
    err=toll+1;
    for(nit=1;nit<=nmax;nit++){
        if(err<=toll){
            *xout=x;
            return 0;
        }
        x1=x-deltax/2;
        x2=x+deltax/2;
        
        if(fun(x1,pointers,&fx1)) return 3; /* Evaluates fx1=fun(x1) */
        if(fun(x2,pointers,&fx2)) return 3; /* Evaluates fx2=fun(x2) */
        dfx=(fx2-fx1)/deltax;
        if(dfx==0){
            return 1;
        }
        else{
            xn=x-fx/dfx;
            err=fabs(xn-x);
            x=xn;
            if(fun(x,pointers,&fx)) return 3; /* Evaluates fx=fun(x) */
        }
    }
    /* Too many itarations */
    return 2;
}

void dcopy(const double *a, double *b, const int size){
    int i;
    for(i=0;i<size;i++){
        b[i]=a[i];
    }
}

double dsum(const double *a, const int size){
    int i;
    double s=0;
    for(i=0;i<size;i++){
        s+=a[i];
    }
    return s;
}

int cartProd_index(int npar, int *n, int **pp_index, int *nindex)
{
    /* Declarations */
    int i,k;            /* Counters */
    int j;              /* Counter (can be negative) */
    int indexcount=-1;  /*=0 in matlab */ /* Counter to build index matrix */
    int *index;
    
    index=*pp_index; /* Used to simplify pointer algebra */
    
    /* Computes nindex, which is the product of the number of the bpas for each
     * parameter */
    *nindex=1; /* Counter for rows in index matrix */
    for(i=0;i<npar;i++)
        *nindex = *nindex * n[i];
    
    /* Allocates memory for index[nindex][npar], and initializes it to 0. */
    if(!(index=(int *)calloc((*nindex)*npar,sizeof(int))))
        return 1;       /* out of memory */
    
    /* Creates index matrix */
    for(i=0;i<*nindex;i++)  /* Counter for the rows of the index matrix, which are nindex */
    {
        /* The first row of index must not be increased */
        if(i)
        {
            index[npar*i+npar-1] = indexcount+1; /* Increases the last element
                                                    of the i-th row of the index
                                                    matrix (except for the first
                                                    row) */
        }   /* if(i) */
        
        for(j=npar-1;j>=0;j--)  /* Cycles on columns of index matrix, which are
                                   npar, from the last to the first */
        {
            if(index[i*npar+j] >= n[j]) /* If index if j-th parameter is major
                                           than the corresponding number of bpas
                                           n[j]... */
            {
                /* Sets to 0 all the indexes in the column j, for all the rows
                   from i (the current) to the last one */
                for(k=i;k<*nindex;k++)
                    index[k*npar+j]=0;
                /* Increases all the indexes in the column j-1 (the previous of
                   the current), for all the rows from i (the current) to the
                   last one */
                for(k=i;k<*nindex;k++)
                    index[k*npar+j-1] = index[k*npar+j-1]+1;
                    
            } /* if(index[i*npar+j] >= n[j]) */
        } /* Next j */
        
        indexcount = index[i*npar+npar-1];
        
    } /* Next i */
    
    /* Modifies the pointer to point to the allocated matrix mc */
    *pp_index=index;
    
    
    return 0;
}

int cartProd_index_prealloc(const int npar, const int *n, const int nindex, int index[])
{
    /* Declarations */
    int i,k;            /* Counters */
    int j;              /* Counter (can be negative) */
    int indexcount=-1;  /*=0 in matlab */ /* Counter to build index matrix */
    
    /* Initialise the index array elements to 0 */
    for(i=0;i<nindex*npar;i++)
        index[i] = 0;
        
    /* Creates index matrix */
    for(i=0;i<nindex;i++)  /* Counter for the rows of the index matrix, which are nindex */
    {
        /* The first row of index must not be increased */
        if(i)
        {
            index[npar*i+npar-1] = indexcount+1; /* Increases the last element
                                                    of the i-th row of the index
                                                    matrix (except for the first
                                                    row) */
        }   /* if(i) */
        
        for(j=npar-1;j>=0;j--)  /* Cycles on columns of index matrix, which are
                                   npar, from the last to the first */
        {
            if(index[i*npar+j] >= n[j]) /* If index if j-th parameter is major
                                           than the corresponding number of bpas
                                           n[j]... */
            {
                /* Sets to 0 all the indexes in the column j, for all the rows
                   from i (the current) to the last one */
                for(k=i;k<nindex;k++)
                    index[k*npar+j]=0;
                /* Increases all the indexes in the column j-1 (the previous of
                   the current), for all the rows from i (the current) to the
                   last one */
                for(k=i;k<nindex;k++)
                    index[k*npar+j-1] = index[k*npar+j-1]+1;
                    
            } /* if(index[i*npar+j] >= n[j]) */
        } /* Next j */
        
        indexcount = index[i*npar+npar-1];
        
    } /* Next i */
    
    
    return 0;
}
