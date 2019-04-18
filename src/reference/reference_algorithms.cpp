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

void computeRange(range_t *rho, Satellite& eph, ionoutc_t *ionoutc, GTime g, Vector3d& xyz)
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
    rate = vel.dot(los)/range;

    // Pseudorange rate.
    rho->rate = rate - SPEED_OF_LIGHT*clk[1];

    // Time of application.
    rho->g = g;

    // Azimuth and elevation angles.
    Vector2d az_el;
    eph.los2azimuthElevation(xyz, los, az_el);
    rho->azel[0] = az_el[0];
    rho->azel[1] = az_el[1];
    lla = WSG84::ecef2lla(xyz);

    // Add ionospheric delay
    rho->iono_delay = ionosphericDelay(ionoutc, g, lla.data(), rho->azel);
    rho->range += rho->iono_delay;


    return;
}
