/******************************************************************************
*                              keplerianMotion.c                              *
*                       Functions for Keplerian motion                        *
*                                                                             *
*                                Space Toolbox                                *
******************************************************************************/

#include "keplerianMotion.h"

int kepPro(double out[6], const double in[6], const double dt, const double mu, const int npar, double tol, int maxiter)
{
	/* Declarations */
    double r0[3], v0[3];
	double nr0, inv_mu, alpha, r0v0, chi0, h2, p, s, w, a, chin, chin1, sqrt_mu, numer;
	double psi, sqrt_psi, c2, c3, r0v0_sqrt_mu, r, nr, chin2, gd, fd, f, g;
	double h[3], r_vect[3], v_vect[3];
	int i;
	
	/* Code */
	/* Optional inputs check */
	if(npar < 2){
        maxiter = KEPPRO_DEFAULTMAXITER;
        if(npar < 1){
            tol = KEPPRO_DEFAULTTOL;
        }
    }
       
	if(fabs(dt) <= KEPPRO_TINYNUMBER){
		for(i = 0; i < 6; i++){
			out[i] = in[i];
        }
		return 0;
    }
	
	r0[0] = in[0]; r0[1] = in[1]; r0[2] = in[2];
	v0[0] = in[3]; v0[1] = in[4]; v0[2] = in[5];
	
	nr0 = sqrt(r0[0]*r0[0] + r0[1]*r0[1] + r0[2]*r0[2]);
	inv_mu = 1./mu;
	sqrt_mu = sqrt(mu);
	
	alpha = -(v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2]) * inv_mu + 2./nr0;
	r0v0 = r0[0]*v0[0] + r0[1]*v0[1] + r0[2]*v0[2];
	
	if(alpha > tol){
		chi0 = sqrt_mu*dt*alpha;
    
#ifdef DEBUG
		if(fabs(alpha - 1.0) < tol)
			printf("Warning: kepPro3:alpha == 1.\n");
#endif
    } else if (fabs(alpha) < tol){
		
		h[0] = r0[1]*v0[2] - r0[2]*v0[1];
		h[1] = r0[2]*v0[0] - r0[0]*v0[2];
		h[2] = r0[0]*v0[1] - r0[1]*v0[0];
		
		h2 = h[0]*h[0] + h[1]*h[1] + h[2]*h[2];
		p = h2 * inv_mu;
		s = 0.5 * atan(1 / (3.*sqrt(mu/pow(p,3))*dt));
		w = atan(pow(tan(s), 1./3.));
		chi0 = sqrt(p) * 2 / tan(2*w);
	} else if (alpha < -tol){
		a = 1./alpha;
		if (dt >= 0.)
			chi0 =  sqrt(-a) * log(-2.0*mu*alpha*dt / (r0v0 + sqrt(-mu*a) * (1.0 - nr0*alpha)));
		else
			chi0 = -sqrt(-a) * log(-2.0*mu*alpha*dt / (r0v0 - sqrt(-mu*a) * (1.0 - nr0*alpha)));
	}
	
	i = 0;
	chin = chi0;
	chin1 = chi0 + tol + 1;
	
	numer = sqrt_mu * (tol+1);
	
	while ((fabs(numer / sqrt_mu) > tol) & (i <= maxiter))
	{
		i = i + 1;
		chin = chin1;
		psi = chin*chin*alpha;
	 	
		if (psi > tol){
			sqrt_psi = sqrt(psi);
			c2 = (1.0 - cos(sqrt_psi)) / psi;
			c3 = (sqrt_psi-sin(sqrt_psi)) / (psi*sqrt_psi);
		} else if (psi < -tol){
			sqrt_psi = sqrt(-psi);
			c2 = (1.0 - cosh(sqrt_psi)) / psi;
			c3 = (sinh(sqrt_psi) - sqrt_psi) / (-psi*sqrt_psi);
		} else{
			c2 = 0.5;
			c3 = 1./6.;
		}
		if(c2 > KEPPRO_HUGENUMBER || c3 > KEPPRO_HUGENUMBER){
			return 3;
		}
		
		r0v0_sqrt_mu = r0v0 / sqrt_mu;
		r = chin*chin * c2 + r0v0_sqrt_mu * chin * (1.-psi*c3) + nr0*(1.-psi*c2);
		nr = fabs(r);
		numer = sqrt_mu*dt - chin*chin*chin*c3 - r0v0_sqrt_mu*chin*chin*c2 - nr0*chin*(1.-psi*c3);
		chin1 = chin+numer/nr;
	}
	
	if (i > maxiter){
		return 2;
	}
	
	chin2 = chin*chin;
	f = 1. - chin2/nr0*c2;
	g = dt - (chin2*chin) / sqrt_mu * c3;
	r_vect[0] = f*r0[0] + g*v0[0];
	r_vect[1] = f*r0[1] + g*v0[1];
	r_vect[2] = f*r0[2] + g*v0[2];
	nr = sqrt(r_vect[0]*r_vect[0]+r_vect[1]*r_vect[1]+r_vect[2]*r_vect[2]);
	gd = 1. - chin2/nr*c2;
	fd = sqrt_mu/nr/nr0 * chin * (psi*c3-1.);
	v_vect[0] = fd*r0[0] + gd*v0[0];
	v_vect[1] = fd*r0[1] + gd*v0[1];
	v_vect[2] = fd*r0[2] + gd*v0[2];
	
	out[0]=r_vect[0]; out[1]=r_vect[1]; out[2]=r_vect[2];
	out[3]=v_vect[0]; out[4]=v_vect[1]; out[5]=v_vect[2];
	
	if (fabs(f*gd - fd*g - 1.) > 0.1)
		return 1;
	return 0;
}

int dvInsertion(const int ibody, const double t, const double vf[3], const double rp, const double ecc, double *dv)
{
    /* Declarations */
    double mu_Planets[10] = MU_PLANETS;     /* array of the planetary constants
                                               of the planets in [km^3/s^2] */
    double mu_P;    /* Planetary constants of the target Planet [km^3/s^2] */
    double xp2[3],vp2[3];   /* Cartesian position and velocity of the target
                               planet [km] and [km/s] */
    double vinf[3], nvinf;  /* Velocity (and its norm) relative to the planet
                               (at infinite on the hyperbola) */
    double nvp1;   /* Norm of the velocity at the pericentre of the hyperbola */
    double nvp2;   /* Norm of velocity at the pericentre of the target orbit */
    
    /* Check the target planet index */
    if(ibody<1 || ibody>9)
        return 1;
        
    /* Planetary constants of the target Planet [km^3/s^2] */
    mu_P = mu_Planets[ibody-1];
    
    /* Cartesian position and velocity of the target planet [km] and [km/s] */
    if(ephSS_car(ibody,t,xp2,vp2)>0)
        return 2;
    
    /* Velocity (and its norm) relative to the planet
       (at infinite on the hyperbola) */
    dmatdiff(vf,vp2,3,vinf);
    nvinf = dnorm(vinf,3);
    
    /* Velocity at the pericentre of the hyperbola */
    nvp1 = sqrt(DSQR(nvinf)+2*mu_P/rp);

    /* Velocity at the pericentre of the target orbit */
    nvp2 = sqrt(mu_P/rp*(1+ecc));
    
    /* Delta-v at the pericentre */
    (*dv) = fabs(nvp1-nvp2);
    
    return 0;
}

int kepEq_t(double f, const double a, const double ecc, const double mu,
        double f0, const double t0, double *t)
{
    /* Declarations */
    long int k, k0;     /* 2PI moduli of resp. f and f0 */
    double E, E0;       /* Eccentric anomalies at resp. t and t0 */
    
    /* Eccentricity check */
    if(ecc>=1)
        return 1;

    /* Transform f to have f in ]-pi, pi] */
    k = (long int) (f/PI2);
    f = f - k*PI2;
    if(f>PI){
        f -= PI2;
        k++;
    }
    else if(f<=-PI){
        f += PI2;
        k--;
    }
        
    
    /* Transform f0 to have f0 in ]-pi, pi] */
    k0 = (long int) (f0/PI2);
    f0 = f0 - k0*PI2;
    if(f0>PI){
        f0 -= PI2;
        k0++;
    }
    else if(f0<=-PI){
        f0 += PI2;
        k0--;
    }
    
    /* Eccentric anomaly at t */
    E = 2.*atan(sqrt((1-ecc)/(1+ecc))*tan(f/2.));
    /* Eccentric anomaly at t */
    E0 = 2.*atan(sqrt((1-ecc)/(1+ecc))*tan(f0/2.));
    
    /* Time corresponding to true anomaly f */
    (*t) = sqrt(DCUBE(a)/mu)*(E-E0-ecc*(sin(E)-sin(E0)) + PI2*(k-k0)) + t0;
    
    return 0;
}

int kepEq_f(const double t, const double a, const double ecc, const double mu,
            const double f0, const double t0, int const imax, double const tol, 
            double *f)
{
    /* Declarations */
    double n;       /* Mean motion [rad/T] */
    double M, M0;      /* Mean anomalies at resp. t and t0 */
    double E, E0;       /* Eccentric anomalies at resp. t and t0 */
    long int nrev;      /* Number of revolution */
    double Mk, Mkold, Mkoldold, Eold;
    int niter;
    int error;      /* Error code */
    
    /* Eccentricity check */
    if(ecc>=1)
        return 2;
    
    /* Mean motion */
    n = sqrt(mu/DCUBE(a));  /* [rad/T] */
    
    /* Calculation of M0
     * ----------------- */
    E0 = acos((ecc+cos(f0))/(1+ecc*cos(f0)));
    if(sin(f0)<0)
        E0 = -E0;
        
    M0 = E0-ecc*sin(E0);
    M = n*(t-t0)+M0;
    nrev = (long int) (n*(t-t0)/PI2);

    /* Newton loop from Battin p. 217 */
    E = M + ecc*sin(M)/(1-sin(M+ecc)+sin(M)); /* Initial condition from Battin
                                               p. 194, Eq. (5.4) */
    /* Initialise the variables used in the loop */
    Mkold = E-ecc*sin(E)+1;     /* To have it for sure different ... */
    Mk = E-ecc*sin(E)+1;        /* ... than Mk at iteration 1 */
    niter = 0;
    Eold = E + 2*tol;           /* To force at least one iteration */
    while((fabs(Eold-E)>tol) && (niter<imax))
    {
        niter += 1;
        Eold = E;
        Mkoldold = Mkold;
        Mkold = Mk;
        Mk = E-ecc*sin(E);
        if(Mk==Mkoldold)    /* Newton loop is in a limit cycle */
            Mk = (Mk + Mkold)/2;   /* Set Mk in the middle of the limit cycle */
        E = E + (M-Mk)/(1-ecc*cos(E));
    }   /* while((fabs(Eold-E)>tol) && (niter<imax)) */
    
    /* Error code depending on the Newton loop convergence */
    if(fabs(Eold-E)>tol)
        error = 1;
    else
        error = 0;
    
    /* True anomaly */
    (*f) = 2*atan(sqrt((1+ecc)/(1-ecc)) * tan(E/2));   /* [rad] [-pi pi] */
    
    if(t>=t0)   /* Forward in time */
    {
        /* f must follow f0 */
        if((*f)<f0)
            (*f) += PI2;
    }
    else /* Backward in time */
    {
        /* f must preceed f0 */
        if((*f)>f0)
            (*f) -= PI2;
    } /* if(t>=t0) */
    
    /* Multiple revolutions */
    (*f) += nrev*PI2;
    
    return error;
}

int swingby(const double v1[3], const double rp, const double mu, const double gamma, const double n_r[3], double v2[3])
{
    /* Declaration */
    double mu_rp;           /* Planetary constant of the planet divided by the
                               radius of pericentre of the hyperbola */
    double theta_inf;
    double delta;           /* Deflection angle */
    double n_pi[3];         /* Rotated n_r around v1 */
    
    /* Compute mu_rp */
    mu_rp = mu/rp;
    
    /* Compute theta_inf (cf. Kaplan "Modern spacecraft dynamics and control",
     * pag. 93
     */
    theta_inf = acos(-mu_rp/(DSQR(dnorm(v1,3))+mu_rp));
    
    /* Compute the deflection angle */
    delta = 2.0*theta_inf - PI;
    
    /* Compute n_pi */
    if (eulerAxisAngle(n_r,v1,gamma,n_pi)>0)
        return 1;
    
    /* Compute v2 */
    if (eulerAxisAngle(v1,n_pi,delta,v2)>0)
        return 2;
    
    
    return 0;
}
