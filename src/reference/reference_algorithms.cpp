#include <Eigen/Core>

#include "GNSS_utils/satellite.h"
#include "GNSS_utils/gtime.h"
#include "GNSS_utils/datetime.h"
#include "GNSS_utils/test_common.h"
#include "GNSS_utils/wgs84.h"

using namespace Eigen;

#define SQR(x) (x*x)

#define dbg(x) std::cout << #x": " << x << std::endl;

void eph2pos(const GTime& t, const eph_t* eph, Vector3d& pos, double* dts)
{
    // From RTKLIB eph2pos() in ephemeris.c
    static const int MAX_ITER_KEPLER = 30;
    static const double RTOL_KEPLER = 1E-13;
    static const double MU_GPS = 3.9860050E14;
    static const double OMGE = 7.2921151467E-5;
    static const double CLIGHT = 299792458.0;

    double tk,M,E,Ek,sinE,cosE,u,r,i,O,sin2u,cos2u,x,y,sinO,cosO,cosi,mu,omge;
    double xg,yg,zg,sino,coso;
    int n,sys,prn;

//    trace(4,"eph2pos : time=%s sat=%2d\n",time_str(time,3),eph->sat);

//    if (eph->A<=0.0) {
//        rs[0]=rs[1]=rs[2]=*dts=*var=0.0;
//        return;
//    }
    tk= (t -eph->toe).toSec();

//    switch ((sys=satsys(eph->sat,&prn))) {
//    case SYS_GAL: mu=MU_GAL; omge=OMGE_GAL; break;
//    case SYS_CMP: mu=MU_CMP; omge=OMGE_CMP; break;
//    default:      mu=MU_GPS; omge=OMGE;     break;
//    }
    mu = MU_GPS;
    omge = OMGE;
    M=eph->M0+(sqrt(mu/(eph->A*eph->A*eph->A))+eph->deln)*tk;

    for (n=0,E=M,Ek=0.0;fabs(E-Ek)>RTOL_KEPLER&&n<MAX_ITER_KEPLER;n++) {
        Ek=E; E-=(E-eph->e*sin(E)-M)/(1.0-eph->e*cos(E));
    }
    if (n>=MAX_ITER_KEPLER) {
//        trace(2,"eph2pos: kepler iteration overflow sat=%2d\n",eph->sat);
        return;
    }
    sinE=sin(E); cosE=cos(E);


//    trace(4,"kepler: sat=%2d e=%8.5f n=%2d del=%10.3e\n",eph->sat,eph->e,n,E-Ek);

    u=atan2(sqrt(1.0-eph->e*eph->e)*sinE,cosE-eph->e)+eph->omg;

    r=eph->A*(1.0-eph->e*cosE);
    i=eph->i0+eph->idot*tk;
    sin2u=sin(2.0*u); cos2u=cos(2.0*u);
    u+=eph->cus*sin2u+eph->cuc*cos2u;
    r+=eph->crs*sin2u+eph->crc*cos2u;
    i+=eph->cis*sin2u+eph->cic*cos2u;
    x=r*cos(u); y=r*sin(u); cosi=cos(i);

    /* beidou geo satellite */
//    if (sys==SYS_CMP&&(eph->flag==2||(eph->flag==0&&prn<=5))) {
//        O=eph->OMG0+eph->OMGd*tk-omge*eph->toes;
//        sinO=sin(O); cosO=cos(O);
//        xg=x*cosO-y*cosi*sinO;
//        yg=x*sinO+y*cosi*cosO;
//        zg=y*sin(i);
//        sino=sin(omge*tk); coso=cos(omge*tk);
//        rs[0]= xg*coso+yg*sino*COS_5+zg*sino*SIN_5;
//        rs[1]=-xg*sino+yg*coso*COS_5+zg*coso*SIN_5;
//        rs[2]=-yg*SIN_5+zg*COS_5;
//    }
//    else {
        O=eph->OMG0+(eph->OMGd-omge)*tk-omge*eph->toes;
        sinO=sin(O); cosO=cos(O);
        pos[0]=x*cosO-y*cosi*sinO;
        pos[1]=x*sinO+y*cosi*cosO;
        pos[2]=y*sin(i);
//    }

    tk= (t - eph->toc).toSec();
    *dts=eph->f0+eph->f1*tk+eph->f2*tk*tk;

    /* relativity correction */
    *dts-=2.0*sqrt(mu*eph->A)*eph->e*sinE/SQR(CLIGHT);
}
/* inner product ---------------------------------------------------------------
* inner product of vectors
* args   : double *a,*b     I   vector a,b (n x 1)
*          int    n         I   size of vector a,b
* return : a'*b
*-----------------------------------------------------------------------------*/
extern double dot(const double *a, const double *b, int n)
{
    double c=0.0;

    while (--n>=0) c+=a[n]*b[n];
    return c;
}
/* transform ecef to geodetic postion ------------------------------------------
* transform ecef position to geodetic position
* args   : double *r        I   ecef position {x,y,z} (m)
*          double *pos      O   geodetic position {lat,lon,h} (rad,m)
* return : none
* notes  : WGS84, ellipsoidal height
*-----------------------------------------------------------------------------*/
static void ecef2pos(const double *r, double *pos)
{
    static int RE_WGS84 = 6378137.0;           /* earth semimajor axis (WGS84) (m) */
    static int FE_WGS84 = (1.0/298.257223563); /* earth flattening (WGS84) */
    double e2=FE_WGS84*(2.0-FE_WGS84),r2=dot(r,r,2),z,zk,v=RE_WGS84,sinp;

    for (z=r[2],zk=0.0;fabs(z-zk)>=1E-4;) {
        zk=z;
        sinp=z/sqrt(r2+z*z);
        v=RE_WGS84/sqrt(1.0-e2*sinp*sinp);
        z=r[2]+v*e2*sinp;
    }
    pos[0]=r2>1E-12?atan(z/sqrt(r2)):(r[2]>0.0?M_PI/2.0:-M_PI/2.0);
    pos[1]=r2>1E-12?atan2(r[1],r[0]):0.0;
    pos[2]=sqrt(r2+z*z)-v;
}

/* ecef to local coordinate transfromation matrix ------------------------------
* compute ecef to local coordinate transfromation matrix
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *E        O   ecef to local coord transformation matrix (3x3)
* return : none
* notes  : matirix stored by column-major order (fortran convention)
*-----------------------------------------------------------------------------*/
static void xyz2enu(const double *pos, double *E)
{
    double sinp=sin(pos[0]),cosp=cos(pos[0]),sinl=sin(pos[1]),cosl=cos(pos[1]);

    E[0]=-sinl;      E[3]=cosl;       E[6]=0.0;
    E[1]=-sinp*cosl; E[4]=-sinp*sinl; E[7]=cosp;
    E[2]=cosp*cosl;  E[5]=cosp*sinl;  E[8]=sinp;
}

/* multiply matrix -----------------------------------------------------------*/
static void matmul(const char *tr, int n, int k, int m, double alpha,
                   const double *A, const double *B, double beta, double *C)
{
    double d;
    int i,j,x,f=tr[0]=='N'?(tr[1]=='N'?1:2):(tr[1]=='N'?3:4);

    for (i=0;i<n;i++) for (j=0;j<k;j++) {
        d=0.0;
        switch (f) {
            case 1: for (x=0;x<m;x++) d+=A[i+x*n]*B[x+j*m]; break;
            case 2: for (x=0;x<m;x++) d+=A[i+x*n]*B[j+x*k]; break;
            case 3: for (x=0;x<m;x++) d+=A[x+i*m]*B[x+j*m]; break;
            case 4: for (x=0;x<m;x++) d+=A[x+i*m]*B[j+x*k]; break;
        }
        if (beta==0.0) C[i+j*n]=alpha*d; else C[i+j*n]=alpha*d+beta*C[i+j*n];
    }
}

// Doppler Calculation
// Adatped from resdop in RTKLIB pntpos.c
double doppler(GTime& t, Satellite &sat, const Eigen::Vector3d& rec_pos,
               const Eigen::Vector3d& rec_vel)
{
//    double rate, pos[3], E[9], a[3], e[3], vs[3];
////    int i, j
////            , nv = 0;

////    trace(3, "resdop  : n=%d\n", n);

    //// THIS IS  Super complicated (and error-prone) way to calculate the los....
//    ecef2pos(rec_pos.data(), pos);
//    xyz2enu(pos, E);

////    for (i = 0; i<n&&i<MAXOBS; i++) {


////        if (obs[i].D[0] == 0.0 || lam == 0.0 || !vsat[i] || norm(rs + 3 + i * 6, 3) <= 0.0)
////        {
////            continue;
////        }
//        Vector2d azel = sat.azimuthElevation(t, rec_pos);
//        /* line-of-sight vector in ecef */
//        double cosel = cos(azel[1]);
//        a[0] = sin(azel[0])*cosel;
//        a[1] = cos(azel[0])*cosel;
//        a[2] = sin(azel[1]);
//        matmul("TN", 3, 1, 3, 1.0, E, a, 0.0, e);
        sat.update(t);
        Vector3d e = (sat.pos - rec_pos).normalized();

        /* satellite velocity relative to receiver in ecef */
        double vs[3];
        for (int j = 0; j<3; j++)
            vs[j] = sat.vel[j] - rec_vel[j];

        std::cout << "e: " << e[0] << ", " << e[1] << ", " << e[2] << std::endl;
        /* range rate with earth rotation correction */
        double rate = dot(vs, e.data(), 3);
//                + Satellite::OMEGA_EARTH / Satellite::C_LIGHT
//                  * (sat.vel[1] * rec_pos[0] + sat.pos[1] * rec_vel[0] - sat.vel[0] * rec_pos[1] - sat.pos[0] * rec_vel[1]);

        /* doppler residual */
        return rate - Satellite::C_LIGHT * sat.clk[1];

        /* design matrix */
//        for (j = 0; j<4; j++) H[j + nv * 4] = j<3 ? -e[j] : 1.0;

//        nv++;
//    }
//    return nv;
}

/* ionosphere model ------------------------------------------------------------
* From RTKLIB rtkcmn.c
* compute ionospheric delay by broadcast ionosphere model (klobuchar model)
* args   : gtime_t t        I   time (gpst)
*          double *ion      I   iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} (JJ: ALWAYS USE DEFAULT)
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
* return : ionospheric delay (L1) (m)
*-----------------------------------------------------------------------------*/
double ionmodel(const GTime &t,const double *pos, const double *azel)
{

    static const double CLIGHT = 299792458.0;
    const double ion_default[]={ /* 2004/1/1 */
        0.1118E-07,-0.7451E-08,-0.5961E-07, 0.1192E-06,
        0.1167E+06,-0.2294E+06,-0.1311E+06, 0.1049E+07
    };
    double tt,f,psi,phi,lam,amp,per,x;
    int week;

//    if (pos[2]<-1E3||azel[1]<=0) return 0.0;
//    if (norm(ion,8)<=0.0)
    const double* ion=ion_default;

    /* earth centered angle (semi-circle) */
    psi=0.0137/(azel[1]/M_PI+0.11)-0.022;

    /* subionospheric latitude/longitude (semi-circle) */
    phi=pos[0]/M_PI+psi*cos(azel[0]);
    if      (phi> 0.416) phi= 0.416;
    else if (phi<-0.416) phi=-0.416;
    lam=pos[1]/M_PI+psi*sin(azel[0])/cos(phi*M_PI);
    /* geomagnetic latitude (semi-circle) */
    phi+=0.064*cos((lam-1.617)*M_PI);

    /* local time (s) */
    tt=43200.0*lam+t.tow_sec;
    tt-=floor(tt/86400.0)*86400.0; /* 0<=tt<86400 */

    /* slant factor */
    f=1.0+16.0*pow(0.53-azel[1]/M_PI,3.0);

    /* ionospheric delay */
    amp=ion[0]+phi*(ion[1]+phi*(ion[2]+phi*ion[3]));
    per=ion[4]+phi*(ion[5]+phi*(ion[6]+phi*ion[7]));
    amp=amp<    0.0?    0.0:amp;
    per=per<72000.0?72000.0:per;
    x=2.0*M_PI*(tt-50400.0)/per;

    return CLIGHT*f*(fabs(x)<1.57 ? 5E-9+amp*(1.0+x*x*(-0.5+x*x/24.0)) : 5E-9);
}


/*! \brief Compute range between a satellite and the receiver
 * FROM gps-sdr-sim gpssim.c https://github.com/osqzss/gps-sdr-sim
 *  \param[out] rho The computed range
 *  \param[in] eph Ephemeris data of the satellite
 *  \param[in] g GPS time at time of receiving the signal
 *  \param[in] xyz position of the receiver
 *
 */
static constexpr double SECONDS_IN_WEEK = 604800.0;
static constexpr double SECONDS_IN_HALF_WEEK = 302400.0;
static constexpr double SECONDS_IN_DAY = 86400.0;
static constexpr double SECONDS_IN_HOUR = 3600.0;
static constexpr double SECONDS_IN_MINUTE = 60.0;
static constexpr double LEAP_SECONDS = 18;
static const double SPEED_OF_LIGHT = 299792458.0;
static const double OMEGA_EARTH = 7.2921151467e-5;
static constexpr double PI = M_PI;


double ionosphericDelay(const ionoutc_t *ionoutc, GTime g, double *llh, double *azel)
{
    double iono_delay = 0.0;
    double E,phi_u,lam_u,F;

    if (ionoutc->enable==false)
        return (0.0); // No ionospheric delay

    E = azel[1]/PI;
    phi_u = llh[0]/PI;
    lam_u = llh[1]/PI;

    // Obliquity factor
    F = 1.0 + 16.0*pow((0.53 - E),3.0);

    if (ionoutc->vflg==false)
        iono_delay = F*5.0e-9*SPEED_OF_LIGHT;
    else
    {
        double t,psi,phi_i,lam_i,phi_m,phi_m2,phi_m3;
        double AMP,PER,X,X2,X4;

        // Earth's central angle between the user position and the earth projection of
        // ionospheric intersection point (semi-circles)
        psi = 0.0137/(E + 0.11) - 0.022;

        // Geodetic latitude of the earth projection of the ionospheric intersection point
        // (semi-circles)
        phi_i = phi_u + psi*cos(azel[0]);
        if(phi_i>0.416)
            phi_i = 0.416;
        else if(phi_i<-0.416)
            phi_i = -0.416;

        // Geodetic longitude of the earth projection of the ionospheric intersection point
        // (semi-circles)
        lam_i = lam_u + psi*sin(azel[0])/cos(phi_i*PI);

        // Geomagnetic latitude of the earth projection of the ionospheric intersection
        // point (mean ionospheric height assumed 350 km) (semi-circles)
        phi_m = phi_i + 0.064*cos((lam_i - 1.617)*PI);
        phi_m2 = phi_m*phi_m;
        phi_m3 = phi_m2*phi_m;

        AMP = ionoutc->alpha0 + ionoutc->alpha1*phi_m
            + ionoutc->alpha2*phi_m2 + ionoutc->alpha3*phi_m3;
        if (AMP<0.0)
            AMP = 0.0;

        PER = ionoutc->beta0 + ionoutc->beta1*phi_m
            + ionoutc->beta2*phi_m2 + ionoutc->beta3*phi_m3;
        if (PER<72000.0)
            PER = 72000.0;

        // Local time (sec)
        t = SECONDS_IN_DAY/2.0*lam_i + g.tow_sec;
        while(t>=SECONDS_IN_DAY)
            t -= SECONDS_IN_DAY;
        while(t<0)
            t += SECONDS_IN_DAY;

        // Phase (radians)
        X = 2.0*PI*(t - 50400.0)/PER;

        if(fabs(X)<1.57)
        {
            X2 = X*X;
            X4 = X2*X2;
            iono_delay = F*(5.0e-9 + AMP*(1.0 - X2/2.0 + X4/24.0))*SPEED_OF_LIGHT;
        }
        else
            iono_delay = F*5.0e-9*SPEED_OF_LIGHT;
    }

    return (iono_delay);
}

//#define DBG(x) printf(#x": %6.6f\n", x); std::cout << std::flush;

void computeRange(range_t *rho, Satellite& eph, ionoutc_t *ionoutc, GTime g, Vector3d& xyz, Vector3d& dxyz)
{

    Vector3d pos,vel,los;
    Vector2d clk;
    double tau;
    double range,rate;
    double xrot,yrot;

    Vector3d lla,ned;
    double tmat[3][3];

    // SV position at time of the pseudorange observation.
//    satpos(eph, g, pos, vel, clk);
    eph.computePositionVelocityClock(g, pos, vel, clk);


    // Receiver to satellite vector and light-time.
    los = pos - xyz;
    tau = los.norm()/SPEED_OF_LIGHT;

    // Extrapolate the satellite position backwards to the transmission time.
    pos[0] -= vel[0]*tau;
    pos[1] -= vel[1]*tau;
    pos[2] -= vel[2]*tau;

    // Earth rotation correction. The change in velocity can be neglected.
    xrot = pos[0] + pos[1]*OMEGA_EARTH*tau;
    yrot = pos[1] - pos[0]*OMEGA_EARTH*tau;
    pos[0] = xrot;
    pos[1] = yrot;

    // New observer to satellite vector and satellite range.
    los =  pos - xyz;
    range = los.norm();
    rho->d = range;
//    DBG(rho->d);

    // Pseudorange.
    rho->range = range - SPEED_OF_LIGHT*clk[0];

    // Relative velocity of SV and receiver.
    rate = (vel - dxyz).dot(los)/range;

    // Pseudorange rate.
    rho->rate = rate - SPEED_OF_LIGHT*clk[1];

    // Time of application.
    rho->g = g;

    // Azimuth and elevation angles.
    Vector2d az_el;
    eph.los2azimuthElevation(xyz, los, az_el);
    rho->azel[0] = az_el[0];
    rho->azel[1] = az_el[1];
    lla = WGS84::ecef2lla(xyz);

    // Add ionospheric delay
    rho->iono_delay = ionosphericDelay(ionoutc, g, lla.data(), rho->azel);
    rho->range += rho->iono_delay;


    return;
}
