/*******************************************************************************
*                               conversion.c                                   *
*                Functions for time and frame conversion                       *
*                                                                              *
*                              Space Toolbox                                   *
*******************************************************************************/

#include "conversion.h"

int kep2car(const double kep[], const double mu, const double p_input, double out[])
{
    /* Declaration */
    double a   = kep[0];    /* In case of parabola, a corresponds to rp */
    double e   = kep[1];
    double i   = kep[2];
    double Om  = kep[3];
    double om  = kep[4];
    double tho = kep[5];
    
    double rotmat[3][3];
    
    double p, r, wom_dot, r_dot;
    double xp, yp, vxp, vyp;
    
    int orbitType = 1;      /* Type of orbit: 0 for circular orbit
                                              1 for elliptic orbit
                                              2 for parabola
                                              3 for hyperbola */

    /* Rotation matrix */
    rotmat[0][0] = cos(om)*cos(Om)-sin(om)*cos(i)*sin(Om);
    rotmat[1][0] = cos(om)*sin(Om)+sin(om)*cos(i)*cos(Om);
    rotmat[2][0] = sin(om)*sin(i);

    rotmat[0][1] = -sin(om)*cos(Om)-cos(om)*cos(i)*sin(Om);
    rotmat[1][1] = -sin(om)*sin(Om)+cos(om)*cos(i)*cos(Om);
    rotmat[2][1] = cos(om)*sin(i);

    rotmat[0][2] = sin(i)*sin(Om);
    rotmat[1][2] = -sin(i)*cos(Om);
    rotmat[2][2] = cos(i);

    /* In plane Parameters */
    if ((e<(1+KEP2CAR_ELIMIT_PAR)) & (e>(1-KEP2CAR_ELIMIT_PAR))){ /* Parabola */
        /* In this case, we assume p is given and meaningful */
        p = p_input;
        orbitType = 2;
    }
    else
    {
        p = a*(1.0-DSQR(e));
        
        if (e<KEP2CAR_ELIMIT_CIR)     /* Circular orbit */
            orbitType = 0;
        else if (e>=(1+KEP2CAR_ELIMIT_PAR))     /* hyperbola */
            orbitType = 3;
    }
    
    r       = p/(1.0+e*cos(tho));
    xp      = r*cos(tho);
    yp      = r*sin(tho);
    wom_dot = sqrt(mu*p)/DSQR(r);
    r_dot   = sqrt(mu/p)*e*sin(tho);
    vxp     = r_dot*cos(tho)-r*sin(tho)*wom_dot;
    vyp     = r_dot*sin(tho)+r*cos(tho)*wom_dot;

    /* 3D cartesian vector */
    out[0] = rotmat[0][0]*xp + rotmat[0][1]*yp;
    out[1] = rotmat[1][0]*xp + rotmat[1][1]*yp;
    out[2] = rotmat[2][0]*xp + rotmat[2][1]*yp;

    out[3] = rotmat[0][0]*vxp + rotmat[0][1]*vyp;
    out[4] = rotmat[1][0]*vxp + rotmat[1][1]*vyp;
    out[5] = rotmat[2][0]*vxp + rotmat[2][1]*vyp;

    return orbitType;
}

int car2kep(const double in[], const double mu, double kep[], double *eccAnom, double *meanAnom, double *dt, double *p)
{
   /* Declaration */
    double r[3], v[3];      /* Cartesian position and velocity */
    double nr;              /* Norm of r */
    double h[3];            /* cross(r,v): Angular momentum vector */
    double nh;              /* Norm of h */
    double plocal;          /* local value of p */
    double ev[3];           /* Eccentricity vector */
    double e;               /* Eccentricity */
    double ne;              /* Fictitious eccentricity */
    double a;               /* Semi-major axis */
    double inc;             /* inclination */
    double n[3];            /* mean motion */
    double nn;              /* norm of mean motion */
    double Om;              /* Argument of the ascending node */
    double om;              /* Argument of the pericentre */
    double th;              /* True anomaly */
    
    double eccAnom_local;   /* local value of Eccentric anomaly, hyperbolic anomaly or parabolic anomaly */
    double sinE, cosE;      /* Sinus and cosinus of the eccentric (or hyperbolic) anomaly */
    double meanAnom_local;  /* local value of the mean anomaly */
    double hyperAnom;               /* Hyperbolic anomaly */
    double sinhH;           /* Sinus of the hyperbolic anomaly */
    double parabAnom;               /* Parabolic anomaly */
    
    double temp3[3];        /* Temporary 3-component vector */
    
    int orbitType = 1;      /* Type of orbit: 0 for circular orbit
                                              1 for elliptic orbit
                                              2 for parabola
                                              3 for hyperbola */
    
    dcopy(in,r,3);
    dcopy(in+3,v,3);
    
    /* Norm of r */
    nr = dnorm(r,3);
    
	/* cross(r,v): Angular momentum vector */
	dcross(r,v,h);
    
    /* Norm of h */
    nh = dnorm(h,3);
    
    /* Inclination */
    inc = acos(h[2]/nh);    
    
    /* Line of nodes vector */
    if (inc != 0.0 && inc != acos(-1)) /* Orbit is out of xy plane */
    {
        /* n=cross([0 0 1],h); n=n/norm(n); */
        nn = dnorm(h,2);   /* 2 is correct even if h has 3 elements */
        n[0] = -h[1] / nn;
        n[1] =  h[0] / nn;
        n[2] = 0.0;
    }
    else /* Orbit is in xy plane: n is not defined */
    {
        /* Arbitrary choice */
        n[0] = 1.0;
        n[1] = 0.0;
        n[2] = 0.0;
        
#ifdef MEXCOMPILE
        mexWarnMsgIdAndTxt("spaceToolboxC:car2kep:planarOrbit","Planar orbit. Arbitrary choice of Omega = 0.");
#endif
    }
    
    /* Argument of the ascending node */
    Om = acos(n[0]);
    if (n[1]<0)
        Om = fmod(PI2 - Om, PI2);
    /* ---> CL: 22/10/2009, Matteo Ceriotti, Camilla Colombo: added mod(Om,2*pi)
                 to avoid numerical approximation error and hence Om=2*pi. */
    
    /* Parameter */
    plocal = DSQR(nh) / mu;
    
    /* Eccentricity vector */
    ev[0] = (v[1]*h[2]-v[2]*h[1]) / mu - r[0]/nr;
    ev[1] = (v[2]*h[0]-v[0]*h[2]) / mu - r[1]/nr;
    ev[2] = (v[0]*h[1]-v[1]*h[0]) / mu - r[2]/nr;
    
    /* Eccentricity */
    e = dnorm(ev,3);
    
    /* Argument of pericentre */
    if (e < CAR2KEP_ELIMIT_CIR){        /* Circular orbit */
        ev[0] = n[0]; ev[1] = n[1]; ev[2] = n[2];
        ne = 1;
#ifdef MEXCOMPILE
        mexWarnMsgIdAndTxt("spaceToolboxC:car2kep:circularOrbit","Circular orbit. Arbitrary choice of omega = 0.");
#endif
    } else {                        /* Non circular orbit */
        ne = e;
    }
    
    om = acos( DMIN(DMAX(ddot(n,ev,3)/ne, -1.), 1.) );
    /* In the circular case om = 0
    * ---> CL: 22/10/2009, Matteo Ceriotti, Camilla Colombo: numerical roundoff
    *           error in the calculation of theta to avoid argument of acos to
    *           be >1 or <(-1) */
    
    dcross(n,ev,temp3);
    if (ddot(h,temp3,3) < 0){
        om = fmod(PI2-om,PI2);
    }
    /*
    * ---> CL: 22/10/2009, Matteo Ceriotti coded modifications by
    *           Jeannette Heiligers, Camilla Colombo: this condition works
    *           for prograde/retrograde, planar/inclined case.
    * ---> CL: 22/10/2009, Matteo Ceriotti, Camilla Colombo: added mod(om,2*pi)
    *           to avoid numerical approximation error and hence om=2*pi. */
    
    /* Semi-major axis */
    if ((e<(1+CAR2KEP_ELIMIT_PAR)) & (e>(1-CAR2KEP_ELIMIT_PAR))){        /* Parabola */
        a = HUGE_VAL;
#ifdef DEBUG
        printf("a = %e\n", a);
        printf("CAR2KEP Warning: Parabola. Semi-major axis is Inf.\n");
#endif
    }
    else
        a = plocal / (1-DSQR(e));
    
    /* True anomaly */
    th = acos(DMIN(DMAX(ddot(ev,r,3)/ne/nr,-1.0),1.0));
    
    /*
    * ---> CL: 06/11/2007, Daniel Novak: numerical roundoff error in the
    *           calculation of theta to avoid argument of acos to be >1 or
    *           <(-1)
    * ---> CL: 16/10/2009, Matteo Ceriotti coded modifications by
    *           Jeannette Heiligers, Camilla Colombo:
    *           Note: do not substitute ne with e in the previous
    *           formula, because we want to have substitute [ne = norm(ev) = 1]
    *           in the circular case and we want e to be the real value of the
    *           eccentricity */
    
    dcross(ev,r,temp3); 
    if(ddot(h,temp3,3)<0)
        /* the condition dot(r,cross(h,ev)) < 0 works in the same way */
        th = fmod(PI2 - th, PI2);
    /*
    * ---> CL: 22/10/2009, Matteo Ceriotti coded modifications by
    *           Jeannette Heiligers, Camilla Colombo: this condition
    *           works for prograde/retrograde, planar/inclined case
    * ---> CL: 22/10/2009, Matteo Ceriotti, Camilla Colombo: added mod(th,2*pi)
    *           to avoid numerical approximation error and hence th=2*pi. */
    
    /* Assign the elements of the output kep */
    kep[0] = a;
    kep[1] = e;
    kep[2] = inc;
    kep[3] = Om;
    kep[4] = om;
    kep[5] = th;
    
    
    /* The following formulas have been found in Vallado chapter 2.
     * ------------------------------------------------------------ */

    /* eccAnom_local  = eccentric or hyperbolic anomaly
     * meanAnom_local  = mean anomaly
     * dt = time from the pericentre passage
     */
    if (e>=CAR2KEP_ELIMIT_CIR && e<=(1-CAR2KEP_ELIMIT_PAR)) /* Ellipse */
    {
        sinE          = sin(th) * sqrt(1-DSQR(e)) / (1+ e*cos(th));
        cosE          = (e+cos(th)) / (1+e*cos(th));
        eccAnom_local = atan2(sinE,cosE);
        if(eccAnom_local<0)
            eccAnom_local += PI2;
        
        meanAnom_local = eccAnom_local - e*sinE;
        
        nn = sqrt(mu/DCUBE(a));
    }
    
    else if (e<CAR2KEP_ELIMIT_CIR) /* Circumference */
    {
        eccAnom_local  = th;
        meanAnom_local = eccAnom_local;
        
        orbitType = 0;
    }
    
    else if(e>=(1+CAR2KEP_ELIMIT_PAR)) /* Hyperbola */
    {
        sinhH          = sin(th) * sqrt(DSQR(e)-1.0) / (1.0 + e*cos(th));
        /* coshH=(e+cos(th))/(1+e*cos(th)); : Not needed */
        hyperAnom      = DASINH(sinhH);
        eccAnom_local  = hyperAnom;
        meanAnom_local = e*sinhH - hyperAnom;
        nn             = sqrt(-mu/DCUBE(a));
        
        orbitType = 3;
    }
    else /* Parabola */
    {
        parabAnom      = tan(th/2); /* Parabolic anomaly */
        eccAnom_local  = parabAnom;
        meanAnom_local = parabAnom + DCUBE(parabAnom)/3.0;
        nn             = 2 * sqrt(mu/DCUBE(plocal));
        
        orbitType = 2;
    }
    
    
    /* Modify the value of the scalar outputs */
    (*eccAnom)  = eccAnom_local;
    (*meanAnom) = meanAnom_local;
    (*p)        = plocal;
    (*dt)       = meanAnom_local / nn;
    
    
    return orbitType;
}

int car2bpl(const double x_car[], const double U_car[], const double vp_car[], double x_bpl[])
{
    double nn[3];       /* normalised vector of U_cart */
    double ee[3];       
    double cc[3];
    double mat[3*3];    /* Transformation matrix */
    double temp[3];    /* temporal vector */
    

    /* nn */
    if(dscalardiv(U_car,dnorm(U_car,3),3,nn))
        return 1;
    
    /* ee */
    dcross(vp_car,nn,temp);
    if(dscalardiv(temp,dnorm(temp,3),3,ee))
        return 1;
        
    /* cc */
    dcross(ee,nn,cc);

    /* Transformation matrix */
    dcopy(ee,mat,3);
    dcopy(nn,mat+3,3);
    dcopy(cc,mat+6,3);
    
    /* Apply transformation */
    dmatmult(mat, 3, 3, 1, x_car, x_bpl);
    
    return 0;
}

int eulerAxisAngle(const double v[],const double n[], const double theta, double v1[])
{
    double nn[3];       /* normalised vector of n */
    double rot[3*3];    /* rotation matrix */
    
    
    /* Normalise n */
    if(dscalardiv(n,dnorm(n,3),3,nn))
        return 1;
    
    /* Define the rotation matrix by column */
    rot[0] = cos(theta)+(1-cos(theta))*DSQR(nn[0]);
    rot[1] = (1-cos(theta))*nn[0]*nn[1]+sin(theta)*nn[2];
    rot[2] = (1-cos(theta))*nn[0]*nn[2]-sin(theta)*nn[1];
    
    rot[3] = (1-cos(theta))*nn[0]*nn[1]-sin(theta)*nn[2];
    rot[4] = cos(theta)+(1-cos(theta))*DSQR(nn[1]);
    rot[5] = (1-cos(theta))*nn[1]*nn[2]+sin(theta)*nn[0];
    
    rot[6] = (1-cos(theta))*nn[0]*nn[2]+sin(theta)*nn[1];
    rot[7] = (1-cos(theta))*nn[1]*nn[2]-sin(theta)*nn[0];
    rot[8] = cos(theta)+(1-cos(theta))*DSQR(nn[2]);
    
    /* Rotation of v */
    dmatmult(v,1,3,3,rot,v1);
    
    return 0;
}

int rth2car(const double x_rth[3], const double s_car[6], double x_car[3])
{
    double r[3],v[3];           /* Position and velocity vector */
    double h[3];
    double rn[3],hn[3],tn[3];   /* normalised vectors */
    double mat[3*3];            /* transformation matrix */
    
    /* Copy the position and velocity vector */
    dcopy(s_car,r,3);
    dcopy(s_car+3,v,3);
    
    /* Compute rn */
    if(dscalardiv(r,dnorm(r,3),3,rn))
        return 1;
    
    /* Compute h and hn */
    dcross(r,v,h);
    if(dscalardiv(h,dnorm(h,3),3,hn))
        return 1;
    
    /* Compute tn */
    dcross(hn,rn,tn);
    
    /* Transformation matrix */
    dcopy(rn,mat,3);
    dcopy(tn,mat+3,3);
    dcopy(hn,mat+6,3);
    
    /* Apply transformation */
    dmatmult(x_rth, 1, 3, 3, mat, x_car);
    
    return 0;
}

int car2rth(const double x_car[3], const double s_car[6], double x_rth[3])
{
    double r[3],v[3];           /* Position and velocity vector */
    double h[3];
    double rn[3],hn[3],tn[3];   /* normalised vectors */
    double mat[3*3];            /* transformation matrix */
    
    /* Copy the position and velocity vector */
    dcopy(s_car,r,3);
    dcopy(s_car+3,v,3);
    
    /* Compute rn */
    if(dscalardiv(r,dnorm(r,3),3,rn))
        return 1;
    
    /* Compute h and hn */
    dcross(r,v,h);
    if(dscalardiv(h,dnorm(h,3),3,hn))
        return 1;
    
    /* Compute tn */
    dcross(hn,rn,tn);
    
    /* Transformation matrix */
    dcopy(rn,mat,3);
    dcopy(tn,mat+3,3);
    dcopy(hn,mat+6,3);
    
    /* Apply transformation */
    dmatmult(mat, 3, 3, 1, x_car, x_rth);
    
    return 0;
}

int tnh2car(const double x_tnh[3], const double s_car[6], double x_car[3])
{
    double r[3],v[3];           /* Position and velocity vector */
    double h[3];
    double tn[3],nn[3],hn[3];   /* normalised vectors */
    double mat[3*3];            /* transformation matrix */
    
    /* Copy the position and velocity vector */
    dcopy(s_car,r,3);
    dcopy(s_car+3,v,3);
    
    /* Compute tn */
    if(dscalardiv(v,dnorm(v,3),3,tn))
        return 1;
    
    /* Compute h and hn */
    dcross(r,v,h);
    if(dscalardiv(h,dnorm(h,3),3,hn))
        return 1;
    
    /* Compute nn */
    dcross(hn,tn,nn);
    
    /* Transformation matrix */
    dcopy(tn,mat,3);
    dcopy(nn,mat+3,3);
    dcopy(hn,mat+6,3);
    
    /* Apply transformation */
    dmatmult(x_tnh, 1, 3, 3, mat, x_car);
    
    return 0;
}

int car2tnh(const double x_car[3], const double s_car[6], double x_tnh[3])
{
    double r[3],v[3];           /* Position and velocity vector */
    double h[3];
    double tn[3],nn[3],hn[3];   /* normalised vectors */
    double mat[3*3];            /* transformation matrix */
    
    /* Copy the position and velocity vector */
    dcopy(s_car,r,3);
    dcopy(s_car+3,v,3);
    
    /* Compute tn */
    if(dscalardiv(v,dnorm(v,3),3,tn))
        return 1;
    
    /* Compute h and hn */
    dcross(r,v,h);
    if(dscalardiv(h,dnorm(h,3),3,hn))
        return 1;
    
    /* Compute nn */
    dcross(hn,tn,nn);
    
    /* Transformation matrix */
    dcopy(tn,mat,3);
    dcopy(nn,mat+3,3);
    dcopy(hn,mat+6,3);
    
    /* Apply transformation */
    dmatmult(mat, 3, 3, 1, x_car, x_tnh);
    
    return 0;
}

int radec2tnh(const double x_radec[3], double x_tnh[3])
{
    x_tnh[0] = x_radec[0] * cos(x_radec[2])*cos(x_radec[1]);
    x_tnh[1] = x_radec[0] * cos(x_radec[2])*sin(x_radec[1]);
    x_tnh[2] = x_radec[0] * sin(x_radec[2]);
    
    return 0;
}

int tnh2radec(const double x_tnh[3], double x_radec[3])
{
    x_radec[0] = dnorm(x_tnh,3);
    
    if (x_radec[0] > 0)
    {
        x_radec[1] = atan2(x_tnh[1],x_tnh[0]);
        x_radec[2] = asin(x_tnh[2]/x_radec[0]);
    }
    else
    {
        x_radec[1] = 0;
        x_radec[2] = 0;
    }
    return 0;
}

int radec2car(const double x_radec[3], const double s_car[6], double x_car[3])
{
    double x_tnh[3];
    int error;
    
    radec2tnh(x_radec, x_tnh);
    error = tnh2car(x_tnh, s_car, x_car);
    
    return error;
}

int car2radec(const double x_car[3], const double s_car[6], double x_radec[3])
{
    double x_tnh[3];
    int error;
    
    if ((error = car2tnh(x_car, s_car, x_tnh)==0))
    {
        error = tnh2radec(x_tnh, x_radec);
        if (error>0)
            error++;
    }
    
    return error;
}

int rth2tnh(const double x_rth[3], const double a, const double e, const double f, const double mu, double x_tnh[3])
{
    double p, n, h, r, v;
    double sinb, cosb;
    double rot[3*3];
    
    if ((e>(1-RTH2TNH_TOL_ECC)) & (e<(1+RTH2TNH_TOL_ECC)))
        return 1;

    /* Parameter */
    p = a*(1-DSQR(e));
    n = sqrt(mu/DCUBE(a));
    h = n*DSQR(a);
    h = h*sqrt(1-DSQR(e));
    r = p/(1+e*cos(f));
    v = sqrt(2*mu/r - mu/a);
    
    /* Sinus and cosinus */
    sinb = h*e/(p*v)*sin(f);
    cosb = h/(p*v)*(1+e*cos(f));
    
    /* Rotation matrix  specified by column*/
    rot[0] = sinb;
    rot[1] = -cosb;
    rot[2] = 0;
    
    rot[3] = cosb;
    rot[4] = sinb;
    rot[5] = 0;
    
    rot[6] = 0;
    rot[7] = 0;
    rot[8] = 1;
    
    /* Apply the rotation matrix */
    dmatmult(x_rth,1,3,3,rot,x_tnh);
    
    return 0;
}

int tnh2rth(const double x_tnh[3], const double a, const double e, const double f, const double mu, double x_rth[3])
{
    double p, n, h, r, v;
    double sinb, cosb;
    double rot[3*3];
    
    if ((e>(1-TNH2RTH_TOL_ECC)) & (e<(1+TNH2RTH_TOL_ECC)))
        return 1;

    /* Parameter */
    p = a*(1-DSQR(e));
    n = sqrt(mu/DCUBE(a));
    h = n*DSQR(a);
    h = h*sqrt(1-DSQR(e));
    r = p/(1+e*cos(f));
    v = sqrt(2*mu/r - mu/a);
    
    /* Sinus and cosinus */
    sinb = h*e/(p*v)*sin(f);
    cosb = h/(p*v)*(1+e*cos(f));
    
    /* Rotation matrix  specified by column*/
    rot[0] = sinb;
    rot[1] = -cosb;
    rot[2] = 0;
    
    rot[3] = cosb;
    rot[4] = sinb;
    rot[5] = 0;
    
    rot[6] = 0;
    rot[7] = 0;
    rot[8] = 1;
    
    /* Apply the rotation matrix */
    dmatmult(rot,3,3,1,x_tnh,x_rth);
    
    return 0;
}


double hms2fracday(const double hrs, const double mn, const double sec)
{
    return (hrs + (mn + sec/60)/60) / 24;
}

int fracday2hms(const double fracday, double hms[3])
{
    double temp;
    
    temp = fracday*24;
    
    /* Hours */
    if(temp>=0)
        hms[0] = floor(temp);
    else
        hms[0] = floor(temp)+1;
        
    /* Minutes */
    temp -= hms[0];
    if(temp>=0)
        hms[1] = floor(temp*60);
    else
        hms[1] = floor(temp*60)+1;
    
    /* Seconds */
    hms[2] = (temp-hms[1]/60)*3600;
    
    return 0;
}

double jd2mjd(const double jd)
{
    return jd - MJD_ZERO_IN_JD;
}

double mjd2jd(const double mjd)
{
    return mjd + MJD_ZERO_IN_JD;
}

double jd2mjd2000(const double jd)
{
    /*                      mjd2000   
     *    |-----------------------------------------|
     *             mjd
     *    |-------------------|                       */
    return jd - MJD_ZERO_IN_JD - MJD2000_ZERO_IN_MJD;
}

double mjd20002jd(const double mjd2000)
{
    /*                           jd   
     *    |----------------------------------------------|
     *                  mjd
     *    |-----------------------------|                  */
    return mjd2000 + MJD2000_ZERO_IN_MJD + MJD_ZERO_IN_JD;
}

double mjd20002mjd(const double mjd2000)
{
    return mjd2000 + MJD2000_ZERO_IN_MJD;
}

double mjd2mjd2000(const double mjd)
{
    return mjd - MJD2000_ZERO_IN_MJD;
}

int date2jd(const double date[6], double *jd)
{
    double Y,M,D,hrs,mn,sec;
    double jd_local;
    
    /* Manage the input */
    Y   = date[0];
    M   = date[1];
    D   = date[2];
    hrs = date[3];
    mn  = date[4];
    sec = date[5];
    
    
    /* Formula converting Gregorian date into JD */
    jd_local = 367*Y - floor(7*(Y+floor((M+9)/12))/4)
                     - floor(3*floor((Y+(M-9)/7)/100+1)/4)
                     + floor(275*M/9)
                     + D + 1721028.5 + hms2fracday(hrs,mn,sec);
    
    /* Check if the input date was valid */
    if (jd_local<0)
        return 1;
    else
    {
        (*jd) = jd_local;
        return 0;
    }
}

int jd2date(const double jd, double date[6])
{
    
    double j, g, dg, c, dc, b, db, a, da, y, m, d;
    
    /* Check if the date is valid */
    if(jd<0)
        return 1;
    
    /* Calculations */
    j = floor(jd+0.5) + 32044;
    g = floor(j/146097);
    dg = fmod(j,146097);
    c = floor((floor(dg/36524)+1) * 3/4);
    dc = dg - c*36524;
    b = floor(dc/1461);
    db = fmod(dc,1461);
    a = floor((floor(db/365)+1) * 3/4);
    da = db - a*365;
    y = g*400 + c*100 + b*4 + a;
    m = floor((da*5 + 308)/153) - 2;
    d = da - floor((m+4)*153/5) + 122;
    
    /* Year, Month and Day */
    date[0] = y-4800 + floor((m+2)/12);
    date[1] = fmod((m+2),12) + 1;
    date[2] = floor(d+1);
    
    /* Hour, Minute and Second */
    fracday2hms(fmod(jd+0.5,floor(jd+0.5)),(date+3));
    
    return 0;   
}

int date2mjd(const double date[6], double *mjd)
{
    double jd;
    int error;
    
    error = date2jd(date, &jd);
    
    if(error==0)
        (*mjd) = jd2mjd(jd);
    
    return error;
}

int mjd2date(const double mjd, double date[6])
{
    return jd2date(mjd2jd(mjd),date);
}

int date2mjd2000(const double date[6], double *mjd2000)
{
    double jd;
    int error;
    
    /* Converte date to JD */
    error = date2jd(date,&jd);
    if(error == 0)
        (*mjd2000) = jd - MJD_ZERO_IN_JD - MJD2000_ZERO_IN_MJD;
    
    return error;
}

int mjd20002date(const double mjd2000, double date[6])
{
    return jd2date(mjd20002jd(mjd2000),date);
}
