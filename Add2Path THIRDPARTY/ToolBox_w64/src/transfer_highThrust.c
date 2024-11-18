/*******************************************************************************
*                          transfer_highThrust.h                               *
*                   Functions for high-thrust transfers                        *
*                                                                              *
*                              Space Toolbox                                   *
*******************************************************************************/

#include "transfer_highThrust.h"


int lambertMR(const double ri[3], const double rf[3], double tof, double mu, int orbittype, int nrev, int ncase, int optionslmr,
	double *a, double *p, double *e, int *error, double vi[3], double vf[3], double *tpar, double *theta){
	
	double rim, rim2, rfm, rfm2, cth, cr[3], sth, c, s, beta, pmin, tmin, lambda, w, r1, s1, l, m, l1;
	double y = 1., x0, x = -1e8, eta, delta, u, sigma, c1, gamma, denom, h1, h2, qr, xpll, lp2xp1;
	double b, u0, ku, r11, s11, t11, constant, h3, Erss, Erss_old, h, dh;
	int b1, checkfeas = 0, n1 = 0, n = 0, m1, i8, i, checknconvrss, checknconvoss, n3;
	int nnew;
	
#ifdef DEBUG
	printf("ri=%f %f %f\nrf=%f %f %f\ntof=%f\nmu=%f\norb=%d\nnrev=%d\nncase=%d\n",ri[0],ri[1],ri[2],rf[0],rf[1],rf[2],tof,mu,orbittype,nrev,ncase);
#endif

	/* Code */
	vi[0] = 0; vi[1] = 0; vi[2] = 0;
	vf[0] = 0; vf[1] = 0; vf[2] = 0;
	*a = 0; *p = 0; *e = 0; *tpar = 0;
	
	rim2 = ri[0]*ri[0] + ri[1]*ri[1] + ri[2]*ri[2];
	rim = sqrt(rim2);
	rfm2 = rf[0]*rf[0] + rf[1]*rf[1] + rf[2]*rf[2];
	rfm = sqrt(rfm2);
	cth = (ri[0]*rf[0] + ri[1]*rf[1] + ri[2]*rf[2]) / (rim*rfm);
	cr[0] = ri[1]*rf[2] - ri[2]*rf[1];
	cr[1] = ri[2]*rf[0] - ri[0]*rf[2];
	cr[2] = ri[0]*rf[1] - ri[1]*rf[0];
	sth = sqrt(cr[0]*cr[0] + cr[1]*cr[1] + cr[2]*cr[2]) / (rim*rfm);
	
	switch(orbittype){
		case 0: /* direct tranfer */
			if(cr[2] < 0) sth = -sth;
			break;
		case 1: /* retrograde transfer */
			if(cr[2] > 0) sth = -sth;
			break;
		default:
			*error = 1;
			printf("%d is not an allowed orbitType\n", orbittype);
			return 1;
	}
	
	*theta = QCK(atan2(sth, cth));
	
	if(((*theta) == PI2) | ((*theta) == 0)){
		*error = 2;
		return 1;
	}
	
	if(sth >= 0)
		b1 = 1;
	else
		b1 = -1;
	
	c = sqrt(rim2 + rfm2 - 2.*rim*rfm*cth);
	s = (rim + rfm +c)/2.;
	beta = 2.*asin(sqrt((s-c)/s));
	pmin = PI2 * sqrt(s*s*s/(8.*mu));
	tmin = pmin * (PI - beta + sin(beta)) / PI2;
	lambda = (double)b1 * sqrt((s-c)/s);
	
	if( (4.*tof*lambda == 0)  || (fabs((s-c)/s) < LAMBERTMR_TOL) ){
		*error = -1;
		return 1;
	}

	if((*theta)*180./PI <= 5.){
		w = atan(pow(rfm/rim, 0.25)) - PI/4;
		r1 = DSQR(sin((*theta)/4.));
		s1 = DSQR(tan(2.*w));
		l = (r1+s1) / (r1 + s1 + cos((*theta)/2));
	} else
		l = DSQR((1 - lambda) / (1 + lambda));
	
	m = 8. * mu * tof*tof / (s*s*s * pow(1+lambda, 6));
	*tpar = (sqrt(2./mu) / 3) * (pow(s, 1.5) - (double)b1 * pow(s-c, 1.5));
	l1 = (1-l) / 2.;
	
	if(nrev == 0){
		
		*error = 0;
		
		if(tof-(*tpar) <= 1e-3)
			x0 = 0;
		else
			x0 = l;
		
		while((fabs(x0-x) >= fabs(x)*LAMBERTMR_TOL+LAMBERTMR_TOL) & (n <= LAMBERTMR_NITERMAX)){
			n++;
			x = x0;
			eta = x / DSQR(sqrt(1+x) + 1);
			checkfeas = 1;
			
			delta = 1.;
			u = 1.;
			sigma = 1.;
			m1 = 0;
			
			while((fabs(u) > LAMBERTMR_TOL) & (m1 <= LAMBERTMR_NITERMAX)){
				m1 = m1 + 1;
				gamma = DSQR(m1 + 3.) / (4. * DSQR(m1 + 3.) - 1.);
				delta = 1. / (1 + gamma*eta*delta);
				u = u * (delta-1.);
				sigma = sigma + u;
			}
			
			c1 = 8. * (sqrt(1+x) + 1) / (3 + 1. / (5 + eta + (9.*eta/7) * sigma));
			
			if(n == 1){
				denom = (1 + 2.*x + l) * (3.*c1 + x*c1 + 4.*x);
				h1 = DSQR(l+x) * (c1 + 1. + 3.*x) / denom;
				h2 = m * (c1 + x - l) / denom;
			} else{
				qr = sqrt(l1*l1 + m/(y*y));
				xpll = qr - l1;
				lp2xp1 = 2.*qr;
				denom = lp2xp1 * (3.*c1 + x*c1 + 4.*x);
				h1 = ((xpll*xpll) * (c1 + 1 + 3.*x)) / denom;
				h2 = m * (c1 + x - l) / denom;
			}
			
			b = 27. * h2 / (4. * DCUBE(1.+h1));
			u = -b / (2. * (sqrt(b+1) + 1));
			
			delta = 1;
			u0 = 1;
			sigma = 1;
			n1 = 0;
			
			while((n1 < LAMBERTMR_NITERMAX) & (fabs(u0) >= LAMBERTMR_TOL)){
				if(n1 == 0){
					gamma = 4./27;
					delta = 1./(1 - gamma*u*delta);
					u0 = u0 * (delta-1);
					sigma = sigma + u0;
				} else {
					for(i8 = 1; i8 < 3; i8++){
						if(i8 == 1)
							gamma = 2.*(3.*(double)n1+1)*(6.*(double)n1-1)/(9.*(4.*(double)n1 - 1)*(4.*(double)n1+1));
						else
							gamma = 2.*(3.*(double)n1+2)*(6.*(double)n1+1)/(9.*(4.*(double)n1 + 1)*(4.*(double)n1+3));
						delta = 1. / (1-gamma*u*delta);
						u0 = u0 * (delta-1);
						sigma = sigma + u0;
					}
				}
				n1++;
			} /* while(n1 < LAMBERTMR_NITERMAX & fabs(u0) >= tol) */
			ku = DSQR(sigma/3.);
			y = ((1.+h1)/3.)*(2.+sqrt(b+1.)/(1.-2.*u*ku));
			x0 = sqrt(DSQR((1.-l)/2.) + m/(y*y))-(1.+l)/2.;
		} /* while(fabs(x0-x) >= LAMBERTMR_TOL & n <= LAMBERTMR_NITERMAX) */
	} else if((nrev > 0) & (4.*tof*lambda != 0)){ /* if(nrev == 0) */
		
		checknconvrss = 1;
		checknconvoss = 1;
		n3 = 1;
		
		while(n3 < 3){
			if((ncase == 0) | (checknconvrss == 0)){
				y = 1;
				n = 0;
				n1 = 0;
				*error = 0;
				if(checknconvoss == 0){
					x0 = 2.*x0;
					checknconvoss = 1;
				} else if(checknconvrss == 0)
					;
				else
					x0 = l;
				x = -1e8;
				
				while((fabs(x0-x) >= fabs(x)*LAMBERTMR_TOL+LAMBERTMR_TOL) & (n <= LAMBERTMR_NITERMAX)){
					n++;
#ifdef DEBUG
					printf("n=%d\n",n);
#endif
					x = x0;
					eta = x / DSQR(sqrt(1+x) + 1);
					checkfeas = 1;
					delta = 1;
					u = 1;
					sigma = 1;
					m1 = 0;
					
					while((fabs(u) > LAMBERTMR_TOL) & (m1 <= LAMBERTMR_NITERMAX)){
						m1++;
						gamma = DSQR(m1+3) / (4.*DSQR(m1+3) - 1);
						delta = 1. / (1+gamma*eta*delta);
						u = u*(delta - 1);
						sigma = sigma + u;
					}
					
					c1 = 8.*(sqrt(1+x)+1) / (3. + 1 / (5.+eta+(9.*eta/7.)*sigma));

					if(n == 1){
						denom = (1 + 2.*x + l) * (3.*c1 + x*c1 + 4.*x);
						h1 = DSQR(l+x) * (c1 + 1. + 3.*x) / denom;
						h2 = m * (c1 + x - l) / denom;
					} else{
						qr = sqrt(l1*l1 + m/(y*y));
						xpll = qr - l1;
						lp2xp1 = 2.*qr;
						denom = lp2xp1 * (3.*c1 + x*c1 + 4.*x);
						h1 = ((xpll*xpll) * (c1 + 1 + 3.*x)) / denom;
						h2 = m * (c1 + x - l) / denom;
					}
					
					h3 = m*(double)nrev*PI/(4*x*sqrt(x));
					h2 = h3 + h2;
					
					b = 27. * h2 / (4. * DCUBE(1.+h1));
					u = -b / (2. * (sqrt(b+1) + 1));
		 			
					delta = 1;
					u0 = 1;
					sigma = 1;
					n1 = 0;
					while((n1 < LAMBERTMR_NITERMAX) & (fabs(u0) >= LAMBERTMR_TOL)){
						if(n1 == 0){
							gamma = 4./27.;
							delta = 1./(1-gamma*u*delta);
							u0 = u0*(delta-1);
							sigma = sigma + u0;
						} else{
							for(i8 = 1; i8 < 3; i8++){
								if(i8 == 1)
									gamma = 2.*(3.*(double)n1+1)*(6.*(double)n1-1)/(9.*(4.*(double)n1 - 1)*(4.*(double)n1+1));
								else
									gamma = 2.*(3.*(double)n1+2)*(6.*(double)n1+1)/(9.*(4.*(double)n1 + 1)*(4.*(double)n1+3));
								delta = 1. / (1-gamma*u*delta);
								u0 = u0 * (delta-1);
								sigma = sigma + u0;
							}
						}
						n1++;
					}
					ku = DSQR(sigma/3.);
					y = ((1.+h1)/3.)*(2.+sqrt(b+1.)/(1.-2.*u*ku));
					
					if(y > sqrt(m/l)){
						if(optionslmr == 2)
							printf("lambertMR:SuccessiveSubstitutionDiverged, Original Successive Substitution is diverging\n-> Reverse Successive Substitution used to find the proper X0.\n");
						checknconvoss = 0;
						break;
					}
					
					x0 = sqrt( DSQR((1-l)/2.)+ m/(y*y)) - (1+l)/2.;
#ifdef DEBUG
					printf("n:%d  x0: %.16f\n",n,x0);
#endif
				}
				if(n >= LAMBERTMR_NITERMAX){ /* Checks if previous loop ended due to maximum number of iterations */
					if(optionslmr == 2)
						printf("lambertMR:SuccessiveSubstitutionExceedMaxIter, Original Successive Substitution exceeded max number of iterations\n-> Reverse Successive Substitution used to find the proper X0.\n");
					checknconvoss = 0;
				}
			}
			
			if(((ncase == 1) | (checknconvoss == 0)) & !((checknconvrss == 0) & (checknconvoss == 0))){
				n = 0;
				n1 = 0;
				*error = 0;
				if(checknconvrss == 0){
					x0 = x0/2;
					checknconvrss = 1;
				} else if(checknconvoss == 0)
					;
				else
					x0 = l;
				
				x = -1e8;
				
				while((fabs(x0-x) >= fabs(x)*LAMBERTMR_TOL+LAMBERTMR_TOL) & (n <= LAMBERTMR_NITERMAX)){
					n++;
					x = x0;
					checkfeas = 1;
					y = sqrt(m/((l+x)*(1+x)));
					
					if(y < 1){
						if(optionslmr == 2)
							printf("lambertMR:SuccessiveSubstitutionDiverged, Reverse Successive Substitution is diverging\n-> Original Successive Substitution used to find the proper X0.\n");
						checknconvrss = 0;
						break;
					}
					
					/* The following Newton-Raphson method should always
					converge, given the previous first guess choice,
					according to the paper. Therefore, the condition on
					number of iterations should not be neccesary. It could be
					necessary for the case tof < tofmin. */

					/* To assure the Newton-Raphson method to be convergent */
					Erss = 2.*atan(sqrt(x));
					while(h_E(Erss,y,m,nrev) < 0)
						Erss = Erss/2.;
					
					nnew = 1;
					Erss_old = -1.e8;
					while((fabs(Erss-Erss_old) >= fabs(Erss)*LAMBERTMR_TOL+LAMBERTMR_TOL) & (nnew <= LAMBERTMR_NITERMAX)){
						nnew++;
						h = h_E2(Erss,y,m,nrev,&dh);
						Erss_old = Erss;
						Erss = Erss - h/dh;
#ifdef DEBUG
						printf("Nnew= %d Erss= %.16f h(Erss)= %.16f\n",nnew,Erss,h);
#endif
					}
					if(nnew >= LAMBERTMR_NITERMAX)
						if(optionslmr != 0)
							printf("lambertMR:NewtonRaphsonIterExceeded\nNewton-Raphson exceeded max iterations.\n");
					x0 = DSQR(tan(Erss/2));
				}
			}
			
			if((checknconvoss == 1) & (checknconvrss == 1))
				break;
			
			if((checknconvrss == 0) & (checknconvoss == 0)){
				if(optionslmr != 0)
					printf("lambertMR:SuccessiveSubstitutionDiverged, Both Original Successive Substitution and Reverse Successive Substitution diverge because Nrev > NrevMAX.\nWork in progress to calculate NrevMAX.\n");
				*error = 3;
				return 1;
			}
			
			n3++;
		}
		
		if(n3 == 3){
			if(optionslmr != 0)
				printf("lambertMR:SuccessiveSubstitutionDiverged, Either Original Successive Substitution or Reverse Successive Substitution is always diverging\nbecause Nrev > NrevMAX or because large-a solution = small-a solution (limit case).\nWork in progress to calculate NrevMAX.\n");
			*error = 3;
			return 1;
			
		}
	}
	
	if(checkfeas == 0){
		*error = 1;
		return 1;
	}
	if((n1 >= LAMBERTMR_NITERMAX) | (n >= LAMBERTMR_NITERMAX)){
		*error = 4;
		if(optionslmr != 0){
			printf("Lambert algorithm has not converged\n");
		}
		return 1;
	}
	
	constant = m * s * DSQR(1+lambda);
	*a = constant / (8*x0*y*y);
	
	r11 = DSQR(1 + lambda) / (4.*tof*lambda);
	s11 = y*(1+x0);
	t11 = (m*s*DSQR(1+lambda)) / s11;
	
	for(i=0; i<3; i++){
		vi[i] = -r11*(s11*(ri[i]-rf[i]) - t11*ri[i]/rim);
		vf[i] = -r11*(s11*(ri[i]-rf[i]) + t11*rf[i]/rfm);
	}
	
	*p = (2.*rim*rfm*y*y*DSQR(1+x0)*DSQR(sin((*theta)/2.))) / constant;
	*e = sqrt(1 - (*p)/(*a));
	
	return 0;
}

double h_E(const double E, const double y, const double m, const int nrev){
    return (nrev*PI + E - sin(E)) / DCUBE(tan(E/2.)) - 4./m * (DCUBE(y) - DSQR(y));
}

double h_E2(const double E, const double y, const double m, const int nrev, double *dh){
	double tanE2;
	tanE2 = tan(E/2.); /* For speed */
	*dh = (1-cos(E))/DCUBE(tanE2) - 3./2.*(nrev*PI+E-sin(E))*DSQR(1./cos(E/2.)) / (DSQR(tanE2)*DSQR(tanE2));
	return (nrev*PI + E - sin(E)) / DCUBE(tanE2) - 4./m * (DCUBE(y) - DSQR(y));
}
