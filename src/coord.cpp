#include "coord.h"
#include "math_core.h"
#include <cmath>
#include <cassert>
#include <stdexcept>

namespace coord{

ProlSph::ProlSph(double Delta) :
    Delta2(Delta*Delta)
{
    if(Delta<=0)
        throw std::invalid_argument("Invalid parameters for Prolate Spheroidal coordinate system");
}

//--------  angular momentum functions --------//

template<> double Ltotal(const PosVelCar& p) {
    return sqrt(pow_2(p.y*p.vz-p.z*p.vy) + pow_2(p.z*p.vx-p.x*p.vz) + pow_2(p.x*p.vy-p.y*p.vx));
}
template<> double Ltotal(const PosMomCar& p) {
    return sqrt(pow_2(p.y*p.pz-p.z*p.py) + pow_2(p.z*p.px-p.x*p.pz) + pow_2(p.x*p.py-p.y*p.px));
}
template<> double Ltotal(const PosVelCyl& p) {
    return sqrt((pow_2(p.R) + pow_2(p.z)) * pow_2(p.vphi) + pow_2(p.R*p.vz-p.z*p.vR));
}
template<> double Ltotal(const PosMomCyl& p) {
	return p.R>0? sqrt((1 + pow_2(p.z)/pow_2(p.R)) * pow_2(p.pphi) + pow_2(p.R*p.pz-p.z*p.pR))
			: fabs(p.R*p.pz-p.z*p.pR);
}
template<> double Ltotal(const PosVelSph& p) {
    return sqrt(pow_2(p.vtheta) + pow_2(p.vphi)) * p.r;
}
template<> double Ltotal(const PosMomSph& p) {
    return sqrt(pow_2(p.ptheta) + pow_2(p.pphi));
}
template<> double Ltotal(const PosVelAxi& p) {
    return Ltotal(toPosVelCyl(p));
}

template<> double Lz(const PosVelCar& p) { return p.x * p.vy - p.y * p.vx; }
template<> double Lz(const PosVelCyl& p) { return p.R * p.vphi; }
template<> double Lz(const PosVelSph& p) {
    double sintheta, costheta;
    math::sincos(p.theta, sintheta, costheta);  // this gives exactly sintheta=0 for theta=M_PI
    return p.r * sintheta * p.vphi;
}
template<> double Lz(const PosVelAxi& p) {
    double chi = p.cs.Delta2>=0 ? p.rho : sqrt(pow_2(p.rho) - p.cs.Delta2);
    double sinnu = 1 / sqrt(1 + pow_2(p.cotnu));
    return p.vphi * chi * sinnu;
}

// multiply two numbers, replacing {anything including INFINITY} * 0 with 0;
// the same result may be achieved by nan2num(x*y), but with two comparisons instead of one
inline double mul(double x, double y) { return y==0 ? 0 : x*y; }

//--------  position conversion functions ---------//

template<> PosCar toPos(const PosCyl& p, const Car) {
    double sinphi, cosphi;
    math::sincos(p.phi, sinphi, cosphi);
    return PosCar(mul(p.R, cosphi), mul(p.R, sinphi), p.z);
}
template<> PosCar toPos(const PosSph& p, const Car) {
    double sintheta, costheta, sinphi, cosphi;
    math::sincos(p.theta, sintheta, costheta);
    math::sincos(p.phi, sinphi, cosphi);
    return PosCar(mul(p.r, sintheta*cosphi), mul(p.r, sintheta*sinphi), mul(p.r, costheta));
}
template<> PosCyl toPos(const PosCar& p, const Cyl) {
    return PosCyl(sqrt(pow_2(p.x) + pow_2(p.y)), p.z, math::atan2(p.y, p.x));
}
template<> PosCyl toPos(const PosSph& p, const Cyl) {
    double sintheta, costheta;
    math::sincos(p.theta, sintheta, costheta);
    return PosCyl(mul(p.r, sintheta), mul(p.r, costheta), p.phi);
}
template<> PosSph toPos(const PosCar& p, const Sph) {
    return PosSph(sqrt(pow_2(p.x)+pow_2(p.y)+pow_2(p.z)),
        math::atan2(sqrt(pow_2(p.x) + pow_2(p.y)), p.z), math::atan2(p.y, p.x));
}
template<> PosSph toPos(const PosCyl& p, const Sph) {
    return PosSph(sqrt(pow_2(p.R) + pow_2(p.z)), math::atan2(p.R, p.z), p.phi);
}
template<> PosCyl toPos(const PosProlSph& p, const Cyl) {
    if(fabs(p.nu)>p.coordsys.Delta2 || p.lambda<p.coordsys.Delta2)
        throw std::invalid_argument("Incorrect ProlSph coordinates");
    const double R = sqrt( (p.lambda-p.coordsys.Delta2) * (1 - fabs(p.nu) / p.coordsys.Delta2) );
    const double z = sqrt( p.lambda * fabs(p.nu) / p.coordsys.Delta2) * (p.nu>=0 ? 1 : -1);
    return PosCyl(R, z, p.phi);
}
// declare an instantiation which will be defined later
template<> PosProlSph toPosDeriv(const PosCyl& from,
    PosDerivT<Cyl, ProlSph>* derivs, PosDeriv2T<Cyl, ProlSph>* derivs2, const ProlSph cs);
template<> PosProlSph toPos(const PosCyl& from, const ProlSph cs) {
    return toPosDeriv<Cyl,ProlSph>(from, NULL, NULL, cs);
}


//-------- position conversion with derivatives --------//

template<>
PosCyl toPosDeriv(const PosCar& p, PosDerivT<Car, Cyl>* deriv, PosDeriv2T<Car, Cyl>* deriv2, const Cyl)
{
    const double R2=pow_2(p.x)+pow_2(p.y), R=sqrt(R2);
    if(R==0) {
        // degenerate case, but provide something meaningful nevertheless,
        // assuming that these numbers will be multiplied by 0 anyway
        if(deriv!=NULL)
            deriv->dRdx=deriv->dRdy=deriv->dphidx=deriv->dphidy=1.;
        if(deriv2!=NULL)
            deriv2->d2Rdx2=deriv2->d2Rdy2=deriv2->d2Rdxdy=
            deriv2->d2phidx2=deriv2->d2phidy2=deriv2->d2phidxdy=1.;
        return PosCyl(0, p.z, 0);
    }
    const double cosphi=p.x/R, sinphi=p.y/R;
    if(deriv!=NULL) {
        deriv->dRdx=cosphi;
        deriv->dRdy=sinphi;
        deriv->dphidx=-sinphi/R;
        deriv->dphidy=cosphi/R;
    }
    if(deriv2!=NULL) {
        deriv2->d2Rdx2 =pow_2(sinphi)/R;
        deriv2->d2Rdy2 =pow_2(cosphi)/R;
        deriv2->d2Rdxdy=-sinphi*cosphi/R;
        deriv2->d2phidx2 =2*sinphi*cosphi/R2;
        deriv2->d2phidy2 =-deriv2->d2phidx2;
        deriv2->d2phidxdy=(pow_2(sinphi)-pow_2(cosphi))/R2;
    }
    return PosCyl(R, p.z, math::atan2(p.y, p.x));
}

template<>
PosSph toPosDeriv(const PosCar& p, PosDerivT<Car, Sph>* deriv, PosDeriv2T<Car, Sph>* deriv2, const Sph)
{
    const double x2=pow_2(p.x), y2=pow_2(p.y), z2=pow_2(p.z);
    const double R2=x2+y2, R=sqrt(R2);
    const double r2=R2+z2, r=sqrt(r2), invr=1/r;
    if(deriv!=NULL) {
        deriv->drdx=p.x*invr;
        deriv->drdy=p.y*invr;
        deriv->drdz=p.z*invr;
        const double temp=p.z/(R*r2);
        deriv->dthetadx=p.x*temp;
        deriv->dthetady=p.y*temp;
        deriv->dthetadz=-R/r2;
        deriv->dphidx=-p.y/R2;
        deriv->dphidy=p.x/R2;
    }
    if(deriv2!=NULL) {
        const double invr3=invr/r2;
        deriv2->d2rdx2=(r2-x2)*invr3;
        deriv2->d2rdy2=(r2-y2)*invr3;
        deriv2->d2rdz2=R2*invr3;
        deriv2->d2rdxdy=-p.x*p.y*invr3;
        deriv2->d2rdxdz=-p.x*p.z*invr3;
        deriv2->d2rdydz=-p.y*p.z*invr3;
        const double invr4=1/(r2*r2);
        const double temp=p.z*invr4/(R*R2);
        deriv2->d2thetadx2=(r2*y2-2*R2*x2)*temp;
        deriv2->d2thetady2=(r2*x2-2*R2*y2)*temp;
        deriv2->d2thetadz2=2*R*p.z*invr4;
        deriv2->d2thetadxdy=-p.x*p.y*(r2+2*R2)*temp;
        const double temp2=(R2-z2)*invr4/R;
        deriv2->d2thetadxdz=p.x*temp2;
        deriv2->d2thetadydz=p.y*temp2;
        deriv2->d2phidx2=2*p.x*p.y/pow_2(R2);
        deriv2->d2phidy2=-deriv2->d2phidx2;
        deriv2->d2phidxdy=(y2-x2)/pow_2(R2);
    }
    return PosSph(r, math::atan2(R, p.z), math::atan2(p.y, p.x));
}

template<>
PosCar toPosDeriv(const PosCyl& p, PosDerivT<Cyl, Car>* deriv, PosDeriv2T<Cyl, Car>* deriv2, const Car)
{
    double sinphi, cosphi;
    math::sincos(p.phi, sinphi, cosphi);
    const double x=mul(p.R, cosphi), y=mul(p.R, sinphi);
    if(deriv!=NULL) {
        deriv->dxdR=cosphi;
        deriv->dydR=sinphi;
        deriv->dxdphi=-y;
        deriv->dydphi= x;
    }
    if(deriv2!=NULL) {
        deriv2->d2xdRdphi=-sinphi;
        deriv2->d2ydRdphi=cosphi;
        deriv2->d2xdphi2=-x;
        deriv2->d2ydphi2=-y;
    }
    return PosCar(x, y, p.z);
}

template<>
PosSph toPosDeriv(const PosCyl& p, PosDerivT<Cyl, Sph>* deriv, PosDeriv2T<Cyl, Sph>* deriv2, const Sph)
{
    const double r = sqrt(pow_2(p.R) + pow_2(p.z));
    const double rinv= 1./r;
    const double costheta=p.z*rinv, sintheta=p.R*rinv;
    if(deriv!=NULL) {
        deriv->drdR=sintheta;
        deriv->drdz=costheta;
        deriv->dthetadR=costheta*rinv;
        deriv->dthetadz=-sintheta*rinv;
    }
    if(deriv2!=NULL) {
        deriv2->d2rdR2=pow_2(costheta)*rinv;
        deriv2->d2rdz2=pow_2(sintheta)*rinv;
        deriv2->d2rdRdz=-costheta*sintheta*rinv;
        deriv2->d2thetadR2=-2*costheta*sintheta*pow_2(rinv);
        deriv2->d2thetadz2=-deriv2->d2thetadR2;
        deriv2->d2thetadRdz=(pow_2(sintheta)-pow_2(costheta))*pow_2(rinv);
    }
    return PosSph(r, math::atan2(p.R, p.z), p.phi);
}

template<>
PosCar toPosDeriv(const PosSph& p, PosDerivT<Sph, Car>* deriv, PosDeriv2T<Sph, Car>* deriv2, const Car)
{
    double sintheta, costheta, sinphi, cosphi;
    math::sincos(p.theta, sintheta, costheta);
    math::sincos(p.phi, sinphi, cosphi);
    const double R=mul(p.r, sintheta), x=mul(R, cosphi), y=mul(R, sinphi), z=mul(p.r, costheta);
    if(deriv!=NULL) {
        deriv->dxdr=sintheta*cosphi;
        deriv->dydr=sintheta*sinphi;
        deriv->dzdr=costheta;
        deriv->dxdtheta=z*cosphi;
        deriv->dydtheta=z*sinphi;
        deriv->dzdtheta=-R;
        deriv->dxdphi=-y;
        deriv->dydphi= x;
    }
    if(deriv2!=NULL) {
        deriv2->d2xdrdtheta=costheta*cosphi;
        deriv2->d2ydrdtheta=costheta*sinphi;
        deriv2->d2zdrdtheta=-sintheta;
        deriv2->d2xdrdphi=-sintheta*sinphi;
        deriv2->d2ydrdphi= sintheta*cosphi;
        deriv2->d2xdtheta2=-x;
        deriv2->d2ydtheta2=-y;
        deriv2->d2zdtheta2=-z;
        deriv2->d2xdthetadphi=-z*sinphi;
        deriv2->d2ydthetadphi= z*cosphi;
        deriv2->d2xdphi2=-x;
        deriv2->d2ydphi2=-y;
    }
    return PosCar(x, y, z);
}

template<>
PosCyl toPosDeriv(const PosSph& p, PosDerivT<Sph, Cyl>* deriv, PosDeriv2T<Sph, Cyl>* deriv2, const Cyl)
{
    double sintheta, costheta;
    math::sincos(p.theta, sintheta, costheta);
    const double R=mul(p.r, sintheta), z=mul(p.r, costheta);
    if(deriv!=NULL) {
        deriv->dRdr=sintheta;
        deriv->dRdtheta=z;
        deriv->dzdr=costheta;
        deriv->dzdtheta=-R;
    }
    if(deriv2!=NULL) {
        deriv2->d2Rdrdtheta=costheta;
        deriv2->d2Rdtheta2=-p.r*sintheta;
        deriv2->d2zdrdtheta=-sintheta;
        deriv2->d2zdtheta2=-p.r*costheta;
    }
    return PosCyl(R, z, p.phi);
}

template<>
PosCyl toPosDeriv(const PosProlSph& p, PosDerivT<ProlSph, Cyl>* deriv, PosDeriv2T<ProlSph, Cyl>* deriv2, const Cyl)
{
    const double absnu = fabs(p.nu);
    const double sign = p.nu>=0 ? 1 : -1;
    const double lminusd = p.lambda-p.coordsys.Delta2;
    const double nminusd = absnu-p.coordsys.Delta2;  // note: |nu|<=Delta^2
    if(nminusd>0 || lminusd<0)
        throw std::invalid_argument("Incorrect ProlSph coordinates");
    const double R = sqrt( lminusd * (1 - absnu / p.coordsys.Delta2) );
    const double z = sqrt( p.lambda * absnu /  p.coordsys.Delta2 ) * (p.nu>=0 ? 1 : -1);
    if(deriv!=NULL) {
        deriv->dRdlambda = 0.5*R/lminusd;
        deriv->dRdnu     = 0.5*R/nminusd * sign;
        deriv->dzdlambda = 0.5*z/p.lambda;
        deriv->dzdnu     = 0.5*z/p.nu;
    }
    if(deriv2!=NULL) {
        deriv2->d2Rdlambda2   = -0.25*R / pow_2(lminusd);
        deriv2->d2Rdnu2       = -0.25*R / pow_2(nminusd);
        deriv2->d2Rdlambdadnu = -0.25*R / (lminusd * nminusd * sign);
        deriv2->d2zdlambda2   = -0.25*z / pow_2(p.lambda);
        deriv2->d2zdnu2       = -0.25*z / pow_2(p.nu);
        deriv2->d2zdlambdadnu = -0.25*z / (p.lambda * p.nu);
    }
    return PosCyl(R, z, p.phi);
}

template<>
PosProlSph toPosDeriv(const PosCyl& from, PosDerivT<Cyl, ProlSph>* deriv, PosDeriv2T<Cyl, ProlSph>* deriv2, const ProlSph cs)
{
    // lambda and nu are roots "t" of equation  R^2/(t-Delta^2) + z^2/t = 1
    double R2     = pow_2(from.R), z2 = pow_2(from.z);
    double signz  = from.z>=0 ? 1 : -1;   // nu will have the same sign as z
    double sum    = R2+z2+cs.Delta2;
    double dif    = R2+z2-cs.Delta2;
    double sqD    = sqrt(pow_2(dif) + 4*R2*cs.Delta2);   // determinant is always non-negative
    if(z2==0) sqD = sum;
    if(R2==0) sqD = fabs(dif);
    double lmd, dmn;  // lambda-Delta^2, Delta^2-|nu| - separately from lambda and nu, to avoid roundoffs
    if(dif >= 0) {
        lmd       = 0.5 * (sqD + dif);
        dmn       = R2>0 ? cs.Delta2 * R2 / lmd : 0;
    } else {
        dmn       = 0.5 * (sqD - dif);
        lmd       = cs.Delta2 * R2 / dmn;
    }
    double lambda = cs.Delta2 + lmd;
    double absnu  = 2 * cs.Delta2 / (sum + sqD) * z2;
    if(absnu*2 > cs.Delta2)             // compare |nu| and Delta^2-|nu|
        absnu     = cs.Delta2 - dmn;    // avoid roundoff errors when Delta^2-|nu| is small
    else
        dmn       = cs.Delta2 - absnu;  // same in the opposite case, when |nu| is small
    if(deriv!=NULL || deriv2!=NULL) {
        if(sqD==0)
            throw std::runtime_error("Error in coordinate conversion Cyl=>ProlSph: "
                "the special case lambda = nu = Delta^2 is not implemented");
        if(deriv!=NULL) {  // accurate expressions valid for arbitrary large/small values (no cancellations)
            deriv->dlambdadR = from.R * 2*lambda / sqD;
            deriv->dlambdadz = from.z * 2*lmd    / sqD;
            deriv->dnudR     = from.R * 2*-absnu / sqD * signz;
            deriv->dnudz     = from.z * 2*dmn    / sqD * signz;
        }
        if(deriv2!=NULL) {  // here no attempts were made to avoid cancellation errors
            double common = 8 * cs.Delta2 * R2 * z2 / pow_3(sqD);
            deriv2->d2lambdadR2 = 1 + sum/sqD - common;
            deriv2->d2lambdadz2 = 1 + dif/sqD + common;
            deriv2->d2nudR2     =(1 - sum/sqD + common) * signz;
            deriv2->d2nudz2     =(1 - dif/sqD - common) * signz;
            deriv2->d2lambdadRdz= 2 * from.R * from.z * (1 - sum * dif / pow_2(sqD)) / sqD;
            deriv2->d2nudRdz    = -deriv2->d2lambdadRdz * signz;
        }
    }
    return PosProlSph(lambda, absnu*signz, from.phi, cs);
}

template<>
PosCyl toPosDeriv(const PosAxi& p,
    PosDerivT<Axi, Cyl> *deriv, PosDeriv2T<Axi, Cyl> *deriv2, const Cyl)
{
    double
        eta      = sqrt(pow_2(p.rho) + fabs(p.cs.Delta2)),
        chi      = p.cs.Delta2 >= 0  ?  p.rho  :    eta,
        psi      = p.cs.Delta2 >= 0  ?    eta  :  p.rho,
        sinnu    = 1 / sqrt(1 + pow_2(p.cotnu)),
        cosnu    = sinnu!=0 ? p.cotnu * sinnu : (p.cotnu>0 ? 1 : -1),
        R        = mul(chi, sinnu),
        z        = mul(psi, cosnu),
        dchidrho = p.cs.Delta2 >= 0  ?  1  :  p.rho / chi,
        dpsidrho = p.cs.Delta2 <= 0  ?  1  :  p.rho / psi;
    if(deriv) {
        deriv->dRdrho = dchidrho * sinnu;
        deriv->dRdnu  = chi * cosnu;
        deriv->dzdrho = dpsidrho * cosnu;
        deriv->dzdnu  =-psi * sinnu;
    }
    if(deriv2) {
        deriv2->d2Rdrho2   = p.cs.Delta2 >= 0  ?  0  : -p.cs.Delta2 / pow_3(chi) * sinnu;
        deriv2->d2Rdrhodnu = dchidrho * cosnu;
        deriv2->d2Rdnu2    = -R;
        deriv2->d2zdrho2   = p.cs.Delta2 <= 0  ?  0  :  p.cs.Delta2 / pow_3(psi) * cosnu;
        deriv2->d2zdrhodnu = -dpsidrho * sinnu;
        deriv2->d2zdnu2    = -z;
    }
    return coord::PosCyl(R, z, p.phi);
}

// common fragment shared between toPosDeriv<Cyl, Axi> and toPosVel<Cyl, Axi>
inline void getPosAxi(const PosCyl& p, const Axi cs,
    double& chi, double& psi, double& cosnu, double& sinnu)
{
    double r2 = pow_2(p.R) + pow_2(p.z);
    double sum = 0.5 * (r2 + cs.Delta2);
    double dif = 0.5 * (r2 - cs.Delta2);
    // these branches select the more accurate way of computing quantities without cancellation,
    // but formally any choice is mathematically correct
    double det = cs.Delta2 >= 0 ?
        sqrt(pow_2(dif) + mul(pow_2(p.R), cs.Delta2)):
        sqrt(pow_2(sum) - mul(pow_2(p.z), cs.Delta2));
    // 2*det = (chi * cosnu)^2 + (psi * sinnu)^2 = chi^2 + D^2 sinnu^2 = psi^2 - D^2 cosnu^2
    if(det == INFINITY) {
        chi = psi = INFINITY;
        if(p.R == INFINITY) { cosnu = 0; sinnu = 1; }
        else { cosnu = p.z > 0 ? 1 : -1; sinnu = 0; }
        return;
    }
    if(sum >= 0) {
        psi   = sqrt(det + sum);
        cosnu = p.z!=0 ? p.z / psi : 0;
    } else {  // implies Delta^2 < 0, r < |Delta|
        cosnu = sqrt((det - sum) / -cs.Delta2) * (p.z>=0 ? 1 : -1);
        psi   = p.z / cosnu;
    }
    if(dif >= 0) {
        chi   = sqrt(det + dif);
        sinnu = p.R!=0 ? p.R / chi : 0;
    } else {  // implies Delta^2 > 0, r < Delta
        sinnu = sqrt((det - dif) / cs.Delta2);
        chi   = p.R / sinnu;
    }
}

template<>
PosAxi toPosDeriv(const PosCyl& p,
    PosDerivT<Cyl, Axi> *deriv, PosDeriv2T<Cyl, Axi> *deriv2, const Axi cs)
{
    double rho, eta, chi, psi, cosnu, sinnu;
    getPosAxi(p, cs, chi, psi, cosnu, sinnu);
    // this branch, by contrast, critically distinguishes between prolate and oblate cases;
    // rho = min(chi, psi) and eta = max(chi, psi)
    if(cs.Delta2 >= 0) {
        rho = chi;
        eta = psi;
    } else {
        eta = chi;
        rho = psi;
    }
    double M = 1 / (pow_2(chi * cosnu) + pow_2(psi * sinnu));  // = 1 / (2*det)
    if(deriv) {
        deriv->drhodR = M * psi * sinnu * eta;
        deriv->drhodz = M * chi * cosnu * eta;
        deriv->dnudR  = M * chi * cosnu;
        deriv->dnudz  =-M * psi * sinnu;
    }
    if(deriv2) {
        double common = pow_2(chi * chi * cosnu) + pow_2(psi * psi * sinnu) - 3 * pow_2(chi * psi);
        deriv2->d2rhodR2 = pow_2(M) * psi * (
            chi * eta * pow_2(cosnu) * (1 - 4*M * cs.Delta2 * pow_2(sinnu)) +
            (cs.Delta2 >= 0 ? 0 : cs.Delta2 * pow_2(sinnu) ) );
        deriv2->d2rhodz2 = pow_2(M) * chi * (
            psi * eta * pow_2(sinnu) * (1 + 4*M * cs.Delta2 * pow_2(cosnu)) -
            (cs.Delta2 <= 0 ? 0 : cs.Delta2 * pow_2(cosnu) ) );
        deriv2->d2rhodRdz= pow_2(M) * cosnu * sinnu * eta * (pow_2(rho) + M * common);
        deriv2->d2nudR2  = pow_3(M) * cosnu * sinnu * common;
        deriv2->d2nudz2  = -deriv2->d2nudR2;
        deriv2->d2nudRdz = pow_3(M) * psi * chi *
            (pow_2(psi * sinnu) * (3 - 2*pow_2(sinnu)) - pow_2(chi * cosnu) * (3 - 2*pow_2(cosnu)));
    }
    return PosAxi(rho, cosnu / sinnu, p.phi, cs);
}

// shortcuts for coordinate conversions without derivatives
template<> PosCyl toPos(const PosAxi& from, const Cyl) {
    return toPosDeriv<Axi, Cyl>(from, NULL);
}

template<> PosAxi toPos(const PosCyl& from, const Axi cs) {
    return toPosDeriv<Cyl, Axi>(from, NULL, NULL, cs);
}

template<>
PosCar toPosDeriv(const PosEll& p,
    PosDerivT<Ell, Car> *deriv, PosDeriv2T<Ell, Car> *deriv2, const Car)
{
    double lambda=fabs(p.lambda),mu=fabs(p.mu),nu=fabs(p.nu);
    int sgnl=math::sign(p.lambda),sgnmu=math::sign(p.mu),sgnnu=math::sign(p.nu);
    double x2=(lambda+p.ell.alpha)*(mu+p.ell.alpha)*(nu+p.ell.alpha)/((p.ell.alpha-p.ell.beta)*(p.ell.alpha-p.ell.gamma));
    double y2=(lambda+p.ell.beta)*(mu+p.ell.beta)*(nu+p.ell.beta)/((p.ell.beta-p.ell.alpha)*(p.ell.beta-p.ell.gamma));
    double z2=(lambda+p.ell.gamma)*(mu+p.ell.gamma)*(nu+p.ell.gamma)/((p.ell.gamma-p.ell.beta)*(p.ell.gamma-p.ell.alpha));
    double x=(x2>0)?sgnl*sqrt(x2):0;
    double y=(y2>0)?sgnmu*sqrt(y2):0;
    double z=(z2>0)?sgnnu*sqrt(z2):0;

    if(deriv) {
        deriv->dxdlambda=x/(2*(lambda+p.ell.alpha));
        deriv->dxdnu=x/(2*(nu+p.ell.alpha));
        deriv->dxdmu = x/(2*(mu+p.ell.alpha));
        deriv->dydlambda=y/(2*(lambda+p.ell.beta));
        deriv->dydnu=y/(2*(nu+p.ell.beta));
        deriv->dydmu = y/(2*(mu+p.ell.beta));
        deriv->dzdlambda=z/(2*(lambda+p.ell.gamma));
        deriv->dzdnu=z/(2*(nu+p.ell.gamma));
        deriv->dzdmu = z/(2*(mu+p.ell.gamma));
    }
    if(deriv2) {
       deriv2->d2xdlambda2=x/(4*pow_2(lambda+p.ell.alpha));
       deriv2->d2xdlambdadmu=x/(4*(lambda+p.ell.alpha)*(mu+p.ell.alpha));
       deriv2->d2xdlambdadnu=x/(4*(lambda+p.ell.alpha)*(nu+p.ell.alpha));
       deriv2->d2xdmu2=x/(4*pow_2(mu+p.ell.alpha));
       deriv2->d2xdmudnu=x/(4*(mu+p.ell.alpha)*(nu+p.ell.alpha));
       deriv2->d2xdnu2=x/(4*pow_2(nu+p.ell.alpha));

       deriv2->d2ydlambda2=y/(4*pow_2(lambda+p.ell.beta));
       deriv2->d2ydlambdadmu=y/(4*(lambda+p.ell.beta)*(mu+p.ell.beta));
       deriv2->d2ydlambdadnu=y/(4*(lambda+p.ell.beta)*(nu+p.ell.beta));
       deriv2->d2ydmu2=y/(4*pow_2(mu+p.ell.beta));
       deriv2->d2ydmudnu=y/(4*(mu+p.ell.beta)*(nu+p.ell.beta));
       deriv2->d2ydnu2=y/(4*pow_2(nu+p.ell.beta));

       deriv2->d2zdlambda2=z/(4*pow_2(lambda+p.ell.gamma));
       deriv2->d2zdlambdadmu=z/(4*(lambda+p.ell.gamma)*(mu+p.ell.gamma));
       deriv2->d2zdlambdadnu=z/(4*(lambda+p.ell.gamma)*(nu+p.ell.gamma));
       deriv2->d2zdmu2=z/(4*pow_2(mu+p.ell.gamma));
       deriv2->d2zdmudnu=z/(4*(mu+p.ell.gamma)*(nu+p.ell.gamma));
       deriv2->d2zdnu2=z/(4*pow_2(nu+p.ell.gamma));
    }
    return coord::PosCar(x,y,z);
}

inline void getPosEll(const coord::PosCar &p,const Ell ell, double& lambda, double& mu, double& nu, double& b, double& c, double &d){
    const double alpha=ell.alpha,beta=ell.beta,gamma=ell.gamma;
    b=-(pow_2(p.x)+pow_2(p.y)+pow_2(p.z)-(alpha+beta+gamma));
    c=alpha*gamma+beta*gamma+alpha*beta-(beta+gamma)*pow_2(p.x)-(alpha+gamma)*pow_2(p.y)-(alpha+beta)*pow_2(p.z);
    d= alpha*beta*gamma-alpha*beta*pow_2(p.z)-alpha*gamma*pow_2(p.y)-beta*gamma*pow_2(p.x);
    double q=c/3-b*b/9,r=(b*c-3*d)/6-pow_3(b)/27;
    double phi1=acos(r/pow(-q,1.5))/3,phi2=phi1-2*M_PI/3;
    nu=2*sqrt(-q)*cos(phi1+2*M_PI/3)-b/3;
    mu=2*sqrt(-q)*cos(phi2)-b/3;
    lambda=(-b-(nu+mu));
}
template<>
PosEll toPosDeriv(const PosCar& p,
    PosDerivT<Car, Ell> *deriv, PosDeriv2T<Car, Ell> *deriv2, const Ell ell)
{
    double lambda,mu,nu,a2,a1,a0;
    getPosEll(p,ell,lambda,mu,nu,a2,a1,a0);
    //0=t^3+a2*t^2+a1*t+a0
    //a2=-(x^2+y^2+z^2-alpha-beta-gamma)
    //a1=alpha*gamma+beta*gamma+alpha*beta-(beta+gamma)*x^2-(alpha+gamma)*y^2-(alpha+beta)*z^2
    //a0=alpha*beta*gamma-alpha*beta*z^2-alpha*gamma*y^2-beta*gamma*x^2
    //(3*t^2+2*a2*t+a1)*dtdx=-da2dx*t^2-da1dx*t-da0dx
    double Al=(deriv||deriv2)?1/(3*pow_2(lambda)+2*a2*lambda+a1):0;
    double An=(deriv||deriv2)?1/(3*pow_2(nu)+2*a2*nu+a1):0;
    double Am=(deriv||deriv2)?1/(3*pow_2(mu)+2*a2*mu+a1):0;
    PosDerivT<Car,Ell> D;
    if(deriv||deriv2) {
        D.dlambdadx=2*Al*p.x*(pow_2(lambda)+(ell.beta+ell.gamma)*lambda+ell.beta*ell.gamma);
        D.dlambdady=2*Al*p.y*(pow_2(lambda)+(ell.alpha+ell.gamma)*lambda+ell.alpha*ell.gamma);
        D.dlambdadz=2*Al*p.z*(pow_2(lambda)+(ell.alpha+ell.beta)*lambda+ell.beta*ell.alpha);
        
        D.dmudx=2*Am*p.x*(pow_2(mu)+(ell.beta+ell.gamma)*mu+ell.beta*ell.gamma);
        D.dmudy=2*Am*p.y*(pow_2(mu)+(ell.alpha+ell.gamma)*mu+ell.alpha*ell.gamma);
        D.dmudz=2*Am*p.z*(pow_2(mu)+(ell.alpha+ell.beta)*mu+ell.beta*ell.alpha);

        D.dnudx=2*An*p.x*(pow_2(nu)+(ell.beta+ell.gamma)*nu+ell.beta*ell.gamma);
        D.dnudy=2*An*p.y*(pow_2(nu)+(ell.alpha+ell.gamma)*nu+ell.alpha*ell.gamma);
        D.dnudz=2*An*p.z*(pow_2(nu)+(ell.alpha+ell.beta)*nu+ell.beta*ell.alpha);
        if(deriv)*deriv=D;
    }
    if(deriv2) {
        double dAldx=-((6*lambda+2*a2)*D.dlambdadx-4*lambda*p.x-2*(ell.beta+ell.gamma)*p.x)*pow_2(Al);
        double dAldy=-((6*lambda+2*a2)*D.dlambdady-4*lambda*p.y-2*(ell.alpha+ell.gamma)*p.y)*pow_2(Al);
        double dAldz=-((6*lambda+2*a2)*D.dlambdadz-4*lambda*p.z-2*(ell.alpha+ell.beta)*p.z)*pow_2(Al);
        deriv2->d2lambdadx2=2*(Al+dAldx*p.x)*(pow_2(lambda)+(ell.beta+ell.gamma)*lambda+ell.beta*ell.gamma)
        +2*Al*p.x*(2*lambda+ell.beta+ell.gamma)*D.dlambdadx;
        deriv2->d2lambdadxdy=2*dAldy*p.x*(pow_2(lambda)+(ell.beta+ell.gamma)*lambda+ell.beta*ell.gamma)
        +2*Al*p.x*(2*lambda+ell.beta+ell.gamma)*D.dlambdady;
        deriv2->d2lambdadxdz=2*dAldz*p.x*(pow_2(lambda)+(ell.beta+ell.gamma)*lambda+ell.beta*ell.gamma)
        +2*Al*p.x*(2*lambda+ell.beta+ell.gamma)*D.dlambdadz;
        deriv2->d2lambdady2=2*(Al+dAldy*p.y)*(pow_2(lambda)+(ell.alpha+ell.gamma)*lambda+ell.alpha*ell.gamma)
        +2*Al*p.y*(2*lambda+ell.alpha+ell.gamma)*D.dlambdady;
        deriv2->d2lambdadydz=2*(dAldz*p.y)*(pow_2(lambda)+(ell.alpha+ell.gamma)*lambda+ell.alpha*ell.gamma)
        +2*Al*p.y*(2*lambda+ell.alpha+ell.gamma)*D.dlambdadz;
        deriv2->d2lambdadz2=2*(Al+dAldz*p.z)*(pow_2(lambda)+(ell.alpha+ell.beta)*lambda+ell.alpha*ell.beta)
        +2*Al*p.z*(2*lambda+ell.alpha+ell.beta)*D.dlambdadz;

        double dAmdx=-((6*mu+2*a2)*D.dmudx-4*mu*p.x-2*(ell.beta+ell.gamma)*p.x)*pow_2(Am);
        double dAmdy=-((6*mu+2*a2)*D.dmudy-4*mu*p.y-2*(ell.alpha+ell.gamma)*p.y)*pow_2(Am);
        double dAmdz=-((6*mu+2*a2)*D.dmudz-4*mu*p.z-2*(ell.alpha+ell.beta)*p.z)*pow_2(Am);
        deriv2->d2mudx2=2*(Am+dAmdx*p.x)*(pow_2(mu)+(ell.beta+ell.gamma)*mu+ell.beta*ell.gamma)
        +2*Am*p.x*(2*mu+ell.beta+ell.gamma)*D.dmudx;
        deriv2->d2mudxdy=2*dAmdy*p.x*(pow_2(mu)+(ell.beta+ell.gamma)*mu+ell.beta*ell.gamma)
        +2*Am*p.x*(2*mu+ell.beta+ell.gamma)*D.dmudy;
        deriv2->d2mudxdz=2*dAmdz*p.x*(pow_2(mu)+(ell.beta+ell.gamma)*mu+ell.beta*ell.gamma)
        +2*Am*p.x*(2*mu+ell.beta+ell.gamma)*D.dmudz;
        deriv2->d2mudy2=2*(Am+dAmdy*p.y)*(pow_2(mu)+(ell.alpha+ell.gamma)*mu+ell.alpha*ell.gamma)
        +2*Am*p.y*(2*mu+ell.alpha+ell.gamma)*D.dmudy;
        deriv2->d2mudydz=2*(dAmdz*p.y)*(pow_2(mu)+(ell.alpha+ell.gamma)*mu+ell.alpha*ell.gamma)
        +2*Am*p.y*(2*mu+ell.alpha+ell.gamma)*D.dmudz;
        deriv2->d2mudz2=2*(Am+dAmdz*p.z)*(pow_2(mu)+(ell.alpha+ell.beta)*mu+ell.alpha*ell.beta)
        +2*Am*p.z*(2*mu+ell.alpha+ell.beta)*D.dmudz;

        double dAndx=-((6*nu+2*a2)*D.dnudx-4*nu*p.x-2*(ell.beta+ell.gamma)*p.x)*pow_2(An);
        double dAndy=-((6*nu+2*a2)*D.dnudy-4*nu*p.y-2*(ell.alpha+ell.gamma)*p.y)*pow_2(An);
        double dAndz=-((6*nu+2*a2)*D.dnudz-4*nu*p.z-2*(ell.alpha+ell.beta)*p.z)*pow_2(An);
        deriv2->d2nudx2=2*(An+dAndx*p.x)*(pow_2(nu)+(ell.beta+ell.gamma)*nu+ell.beta*ell.gamma)
        +2*An*p.x*(2*nu+ell.beta+ell.gamma)*D.dnudx;
        deriv2->d2nudxdy=2*dAndy*p.x*(pow_2(nu)+(ell.beta+ell.gamma)*nu+ell.beta*ell.gamma)
        +2*An*p.x*(2*nu+ell.beta+ell.gamma)*D.dnudy;
        deriv2->d2nudxdz=2*dAndz*p.x*(pow_2(nu)+(ell.beta+ell.gamma)*nu+ell.beta*ell.gamma)
        +2*An*p.x*(2*nu+ell.beta+ell.gamma)*D.dnudz;
        deriv2->d2nudy2=2*(An+dAndy*p.y)*(pow_2(nu)+(ell.alpha+ell.gamma)*nu+ell.alpha*ell.gamma)
        +2*An*p.y*(2*nu+ell.alpha+ell.gamma)*D.dnudy;
        deriv2->d2nudydz=2*(dAndz*p.y)*(pow_2(nu)+(ell.alpha+ell.gamma)*nu+ell.alpha*ell.gamma)
        +2*An*p.y*(2*nu+ell.alpha+ell.gamma)*D.dnudz;
        deriv2->d2nudz2=2*(An+dAndz*p.z)*(pow_2(nu)+(ell.alpha+ell.beta)*nu+ell.alpha*ell.beta)
        +2*An*p.z*(2*nu+ell.alpha+ell.beta)*D.dnudz;
    }
    return PosEll(lambda,mu,nu,ell);
}

template<> PosCar toPos(const PosEll& from, const Car) {
    return toPosDeriv<Ell, Car>(from, NULL);
}

template<> PosEll toPos(const PosCar& from, const Ell ell) {
    return toPosDeriv<Car, Ell>(from, NULL, NULL, ell);
}

template<>
PosCar toPosDeriv(const PosEls& p,
    PosDerivT<Els, Car> *deriv, PosDeriv2T<Els, Car> *deriv2, const Car)
{
    double snchi=1/sqrt(1+pow_2(p.cotchi));
    double cschi=p.cotchi*snchi;
    double csphi=cos(p.phi),snphi=sin(p.phi);
    double E1=(p.els.Deltaz2-p.els.Deltay2)/p.els.Deltaz2;
    double E2=p.els.Deltay2/p.els.Deltaz2;
    double A1=sqrt(1-E1*pow_2(cschi));
    double A2=sqrt(1-E2*pow_2(snphi));
    double x=p.rho*csphi*A1;
    double y=sqrt(pow_2(p.rho)+p.els.Deltay2)*snphi*snchi;
    double z=sqrt(pow_2(p.rho)+p.els.Deltaz2)*cschi*A2;
    if(deriv) {
        deriv->dxdrho=csphi*A1;
        deriv->dxdchi=p.rho*csphi/A1*E1*cschi*snchi;
        deriv->dxdphi=-p.rho*snphi*A1;
        deriv->dydrho=p.rho/sqrt(pow_2(p.rho)+p.els.Deltay2)*snphi*snchi;
        deriv->dydchi=sqrt(pow_2(p.rho)+p.els.Deltay2)*snphi*cschi;
        deriv->dydphi=sqrt(pow_2(p.rho)+p.els.Deltay2)*csphi*snchi;
        deriv->dzdrho=p.rho/sqrt(pow_2(p.rho)+p.els.Deltaz2)*cschi*A2;
        deriv->dzdchi=sqrt(pow_2(p.rho)+p.els.Deltaz2)*(-snchi)*A2;
        deriv->dzdphi=sqrt(pow_2(p.rho)+p.els.Deltaz2)*cschi*(-E2*snphi*csphi)/A2;
    }
    if(deriv2) {
       deriv2->d2xdrho2=0;
       deriv2->d2xdrhodchi=csphi/A1*(E1*cschi*snchi);
       deriv2->d2xdrhodphi=-A1*snphi;
       deriv2->d2xdchi2=p.rho*csphi*E1/pow_3(A1)*(pow_2(cschi)-pow_2(snchi)-E1*pow_2(cschi*cschi));
       deriv2->d2xdchidphi=-E1*p.rho*cschi*snchi*snphi/A1;
       deriv2->d2xdchi2=-x;

       deriv2->d2ydrho2=p.els.Deltay2*snchi*snphi/pow_3(sqrt(pow_2(p.rho)+p.els.Deltay2));
       deriv2->d2ydrhodchi=p.rho/sqrt(pow_2(p.rho)+p.els.Deltay2)*snphi*cschi;
       deriv2->d2ydrhodphi=p.rho/sqrt(pow_2(p.rho)+p.els.Deltay2)*csphi*snchi;
       deriv2->d2ydchi2=-y;
       deriv2->d2ydchidphi=sqrt(pow_2(p.rho)+p.els.Deltay2)*csphi*cschi;
       deriv2->d2ydphi2=-y;

       deriv2->d2zdrho2=p.els.Deltaz2/pow_3(sqrt(pow_2(p.rho)+p.els.Deltaz2))*cschi*A2;
       deriv2->d2zdrhodchi=-p.rho/sqrt(pow_2(p.rho)+p.els.Deltaz2)*A2*snchi;
       deriv2->d2zdrhodphi=p.rho/sqrt(pow_2(p.rho)+p.els.Deltaz2)*cschi*(-E2*cschi*snchi)/A2;
       deriv2->d2zdchi2=-z;
       deriv2->d2zdchidphi=E2*sqrt(pow_2(p.rho)+p.els.Deltaz2)*csphi*snphi*snchi/A2;
       deriv2->d2zdphi2=-E2*sqrt(pow_2(p.rho)+p.els.Deltaz2)*cschi*(pow_2(csphi)-pow_2(snphi)+
            E2*pow_2(pow_2(snphi)))/pow_3(A2);
    }
    return coord::PosCar(x,y,z);
}
inline void getPosEls(const coord::PosCar &p,const Els els, double& rho, double& cschi, double &snchi, double& phi, double& a2, double& a1, double &a0){
    a2 = -(pow_2(p.x)+pow_2(p.y)+pow_2(p.z)-els.Deltay2-els.Deltaz2);
    a1 = els.Deltay2*els.Deltaz2-(els.Deltay2+els.Deltaz2)*pow_2(p.x)-els.Deltaz2*pow_2(p.y)-els.Deltay2*pow_2(p.z);
    a0 = -els.Deltay2*els.Deltaz2*pow_2(p.x);
    double q=a1/3-a2*a2/9,r=(a1*a2-3*a0)/6-pow_3(a2)/27;
    double phi1=acos(r/pow(-q,1.5))/3;
    double lambda=2*sqrt(-q)*cos(phi1)-a2/3;
    if(std::isinf(-q))lambda=INFINITY;
    rho=(lambda)>0?sqrt(lambda):1e-6;
    double E2=(els.Deltay2>0)?els.Deltay2/els.Deltaz2:0;
    double b2=1-E2;
    double b1=-(1+pow_2(p.z)/(pow_2(rho)+els.Deltaz2)*(1-E2));
    if(p.x!=0)b1+=E2*pow_2(p.x/rho);
    double b0=pow_2(p.z)/(pow_2(rho)+els.Deltaz2);
    double Dn=b1*b1-4*b0*b2;
    double cschi2=0;
    if(b2==0)cschi2=-b0/b1;
    else if(Dn>0)cschi2=(-b1-sqrt(Dn))/(2*b2);
    else cschi2=-b1/(2*b2);
    if(cschi2<0)cschi2=0;
    cschi=math::sign(p.z)*sqrt(cschi2);
    if(fabs(cschi2)>1)printf("(%f,%f,%f) %f\n",p.x,p.y,p.z,cschi2);
    if(fabs(cschi)>1)cschi=math::sign(cschi)*1.;
    snchi=sqrt(1-cschi*cschi);
    phi=0;
    if(cschi2>0.5&&rho>0){
        double X1=fabs(p.x/(rho*sqrt(1-(1-E2)*cschi2)));
        if(X1>1){
            printf("nox\n");
            X1=1;
        }
        phi=acos(X1);
    }
    else if(els.Deltay2>0||rho>0){
        double Y1=fabs(p.y/(sqrt(rho*rho+els.Deltay2)*snchi));
        if(fabs(Y1)>1)Y1=1;
        if(Y1>1){
            printf("noy\n");
            Y1=1;
        }
        phi=asin(Y1);
    }
    
    if(p.x<0)phi=M_PI-math::sign(p.y)*phi;
    else if(p.y<0)phi=2*M_PI-phi;
    //printf("chi:%f\n",snchi);
    //printf("phi:%f\n",phi);
}
template<>
PosEls toPosDeriv(const PosCar& p,
    PosDerivT<Car, Els> *deriv, PosDeriv2T<Car, Els> *deriv2, const Els els)
{
    double rho,snchi,cschi,phi,a2,a1,a0;
    getPosEls(p,els,rho,cschi,snchi,phi,a2,a1,a0);
    double cotchi=cschi/snchi;
    double cschi2=pow_2(cschi);
    double csphi=cos(phi),csphi2=pow_2(csphi),snphi=sin(phi);
    double Eyz=(els.Deltay2==0)?0:els.Deltay2/els.Deltaz2;
    double a0n=els.Deltaz2+pow_2(p.x)+pow_2(p.y)+2*pow_2(p.z)-Eyz*(pow_2(p.x)+pow_2(p.z)+els.Deltaz2);
    double a1n=-2*(1-Eyz)*(pow_2(p.x)+pow_2(p.y)+pow_2(p.z)+2*els.Deltaz2-els.Deltay2);
    double a2n=3*(els.Deltaz2-els.Deltay2)*(1-Eyz);
    double a0m=pow_2(p.x)+pow_2(p.y)+Eyz*(-els.Deltaz2+pow_2(p.x)+pow_2(p.z));
    double a1m=2*Eyz*(els.Deltay2+els.Deltaz2-pow_2(p.x)-pow_2(p.y)-pow_2(p.z));
    double lambda=rho*rho;
    //0=t^3+a2*t^2+a1*t+a0
    //a2=-(x^2+y^2+z^2-Deltay^2-Deltaz^2)
    //a1=Deltay^2*Deltaz^2-(Deltay^2+Deltaz^2)*x^2-Deltaz2*y^2-Deltay^2*z^2
    //a0=-Deltay^2*Deltaz^2*x^2
    //rho=sqrt(lambda),
    //(3*t^2+2*a2*t+a1)*dtdx=-da2dx*t^2-da1dx*t-da0dx
    double K12=1-(1-Eyz)*cschi2;
    double K1=(K12>0)?sqrt(K12):0;
    double ratK=(K1!=0)?snchi/K1:1;
    if(K1>1)K1=1;
    double Ar=(deriv||deriv2)?1/((3*pow_2(lambda)+2*a2*lambda+a1)):0;
    double Am=(deriv||deriv2)?1/(-3*Eyz*els.Deltay2*pow_2(csphi2)+a1m*csphi2+a0m):0;
    double An=(deriv||deriv2)?1/(a2n*pow_2(cschi2)+a1n*cschi2+a0n):0;
    double Rny=sqrt(pow_2(rho)+els.Deltay2);
    double Rnz=sqrt(pow_2(rho)+els.Deltaz2);
    double ratRny=(els.Deltay2>0)?rho/Rny:0;
    double ratRnz=(els.Deltaz2>0)?rho/Rnz:0;
    double N2=(1-(1-Eyz)*cschi2);
    double N=(N2>0)?sqrt(N2):0;
    double snphiy=sqrt(1-Eyz*csphi2);
    if(N>1)N=1;
    PosDerivT<Car,Els> D;
    if(deriv){
        D.drhodx=Ar*K1*csphi*(pow_2(lambda)+(els.Deltay2+els.Deltaz2)*lambda+els.Deltay2*els.Deltaz2);
        D.drhody=Ar*rho*p.y*(lambda+els.Deltaz2);
        D.drhodz=Ar*rho*p.z*(lambda+els.Deltay2);
        D.dchidx=An*p.x*(1-Eyz)*snchi*cschi;
        D.dchidy=An*Rny*snphi*cschi*(1-(1-Eyz)*cschi2);
        D.dchidz=-An*Rnz*snphiy*snchi*(1-cschi2+Eyz*cschi2);
        D.dphidx=-Am*rho*N*snphi*(1-Eyz*csphi2);
        D.dphidy=Am*csphi*(1-Eyz*csphi2)*Rny*snchi;
        D.dphidz=Am*Eyz*snphi*csphi*p.z;
        if(deriv)*deriv=D;
    }
    if(deriv2){
        deriv2->d2rhodx2=-D.drhodx*Ar*((6*lambda+2*a2)*2*rho*D.drhodx-4*p.x*lambda-2*(els.Deltay2+els.Deltaz2)*p.x)
                +Ar*((1-Eyz)*cschi*ratK*D.dchidx*csphi-snphi*D.dphidx*K1)*(pow_2(lambda)+(els.Deltay2+els.Deltaz2)*lambda+els.Deltay2*els.Deltaz2)
                +Ar*K1*csphi*(2*lambda+els.Deltay2+els.Deltaz2)*D.drhodx*2*rho;
        deriv2->d2rhodxdy=-Ar*((6*lambda+2*a2)*2*rho*D.drhody-4*p.y*lambda-2*els.Deltaz2*p.y)*D.drhodx
                +Ar*((1-Eyz)*cschi*ratK*D.dchidy*csphi-snphi*D.dphidy*K1)*(pow_2(lambda)+(els.Deltay2+els.Deltaz2)*lambda+els.Deltay2*els.Deltaz2)
                +Ar*K1*csphi*(2*lambda+els.Deltay2+els.Deltaz2)*D.drhody*2*rho;
        deriv2->d2rhodxdz=-Ar*((6*lambda+2*a2)*2*rho*D.drhodz-4*p.z*lambda-2*els.Deltay2*p.z)*D.drhodx
                +Ar*((1-Eyz)*cschi*ratK*D.dchidz*csphi-snphi*D.dphidz*K1)*(pow_2(lambda)+(els.Deltay2+els.Deltaz2)*lambda+els.Deltay2*els.Deltaz2)
                +Ar*K1*csphi*(2*lambda+els.Deltay2+els.Deltaz2)*D.drhodz*2*rho;
        deriv2->d2rhody2=-Ar*D.drhody*((6*lambda+2*a2)*2*rho*D.drhody-4*p.y*lambda-2*p.y*els.Deltaz2)
                +Ar*(D.drhody*p.y*(lambda+els.Deltaz2)+rho*(lambda+els.Deltaz2)+rho*p.y*2*rho*D.drhody);
        deriv2->d2rhodydz=-Ar*D.drhody*((6*lambda+2*a2)*2*rho*D.drhodz-4*p.z*lambda-2*p.z*els.Deltay2)
                +Ar*(D.drhodz*p.y*(lambda+els.Deltaz2)+rho*p.y*2*rho*D.drhodz);
        deriv2->d2rhodz2=-Ar*D.drhodz*((6*lambda+2*a2)*2*rho*D.drhodz-4*p.z*lambda-2*p.z*els.Deltay2)
                +Ar*(D.drhodz*p.z*(lambda+els.Deltay2)+rho*(lambda+els.Deltay2)+rho*p.z*2*rho*D.drhodz);

        deriv2->d2chidx2=-An*D.dchidx*((2*a2n*cschi2+a1n)*(-2*cschi*snchi)*D.dchidx-4*(1-Eyz)*p.x*cschi2+2*p.x*(1-Eyz))
                +An*(1-Eyz)*(snchi*cschi+p.x*(2*cschi2-1)*D.dchidx);
        deriv2->d2chidxdy=-An*D.dchidx*((2*a2n*cschi2+a1n)*(-2*cschi*snchi)*D.dchidy-4*(1-Eyz)*p.y*cschi2+2*p.y)
                +An*(1-Eyz)*p.x*(2*cschi2-1)*D.dchidy;
        deriv2->d2chidxdz=-An*D.dchidx*((2*a2n*cschi2+a1n)*(-2*cschi*snchi)*D.dchidz-4*(1-Eyz)*p.z*cschi2+2*p.z*(2-Eyz))
                +An*(1-Eyz)*p.x*(2*cschi2-1)*D.dchidz;
        deriv2->d2chidy2=-An*D.dchidy*((2*a2n*cschi2+a1n)*(-2*cschi*snchi)*D.dchidy-4*(1-Eyz)*p.y*cschi2+2*p.y)
                +An*((ratRny*D.drhody*snphi+Rny*csphi*D.dphidy)*cschi*(1-(1-Eyz)*cschi2)
                +Rny*snphi*(1-3*(1-Eyz)*cschi2)*(-snchi)*D.dchidy);
        deriv2->d2chidydz=-An*D.dchidy*((2*a2n*cschi2+a1n)*(-2*cschi*snchi)*D.dchidz-4*(1-Eyz)*p.z*cschi2+2*p.z*(2-Eyz))
                +An*((ratRny*D.drhodz*snphi+Rny*csphi*D.dphidz)*cschi*(1-(1-Eyz)*cschi2)
                +Rny*snphi*(1-3*(1-Eyz)*cschi2)*(-snchi)*D.dchidz);
        deriv2->d2chidz2=-An*D.dchidz*((2*a2n*cschi2+a1n)*(-2*cschi*snchi)*D.dchidz-4*(1-Eyz)*p.z*cschi2+2*p.z*(2-Eyz))
                -An*((ratRnz*D.drhodz*snphiy+Rnz*Eyz*csphi*snphi/snphiy*D.dphidz)*snchi*(1-cschi2+Eyz*cschi2)
                +Rnz*snphiy*(3*pow_2(snchi)*cschi+Eyz*(cschi2*cschi-2*cschi*pow_2(snchi)))*D.dchidz);


        double dAmdy=-Am*((-6*Eyz*els.Deltay2*csphi2+a1m)*(-2*csphi*snphi)*D.dphidy-4*Eyz*p.y*csphi2+2*p.y);
        double dAmdz=-Am*((-6*Eyz*els.Deltay2*csphi2+a1m)*(-2*csphi*snphi)*D.dphidz-4*Eyz*p.z*csphi2+2*Eyz*p.z);
        deriv2->d2phidx2=-Am*D.dphidx*((-6*Eyz*els.Deltay2*csphi2+a1m)*(-2*csphi*snphi)*D.dphidx-4*Eyz*p.x*csphi2+2*p.x*(1+Eyz))
                -Am*((D.drhodx*N+rho/N*(1-Eyz)*cschi*snchi*D.dchidx)*snphi*(1-Eyz*csphi2)
                +rho*N*(csphi-Eyz*(csphi*csphi2-2*csphi*pow_2(snphi)))*D.dphidx);
        deriv2->d2phidxdy=D.dphidx*dAmdy
                -Am*((D.drhody*N+rho/N*(1-Eyz)*cschi*snchi*D.dchidy)*snphi*(1-Eyz*csphi2)
                +rho*N*(csphi-Eyz*(csphi*csphi2-2*csphi*pow_2(snphi)))*D.dphidy);
        deriv2->d2phidxdz=D.dphidx*dAmdz-Am*((D.drhodz*N+rho/N*(1-Eyz)*cschi*snchi*D.dchidz)*snphi*(1-Eyz*csphi2)
                +rho*N*(csphi-Eyz*(csphi*csphi2-2*csphi*pow_2(snphi)))*D.dphidz);
        deriv2->d2phidy2=D.dphidy*dAmdy
                +Am*((D.drhody*ratRny*snchi+Rny*cschi*D.dchidy)*(1-Eyz*csphi2)*csphi+Rny*snchi*(-1+3*Eyz*csphi2)*snphi*D.dphidy);
        deriv2->d2phidydz=D.dphidy*dAmdz
                +Am*((D.drhodz*ratRny*snchi+Rny*cschi*D.dchidz)*(1-Eyz*csphi2)*csphi+Rny*snchi*(-1+3*Eyz*csphi2)*snphi*D.dphidz);
        deriv2->d2phidz2=D.dphidz*dAmdz+Am*Eyz*(csphi*snphi+p.z*(2*csphi2-1)*D.dphidz);

    }
    return coord::PosEls(rho,cotchi,phi,els);
}

template<> PosCar toPos(const PosEls& from, const Car) {
    return toPosDeriv<Els, Car>(from, NULL);
}

template<> PosEls toPos(const PosCar& from, const Els els) {
    return toPosDeriv<Car, Els>(from, NULL, NULL, els);
}

//--------  position+velocity conversion functions  ---------//

template<> PosVelCar toPosVel(const PosVelCyl& p, const Car) {
    double sinphi, cosphi;
    math::sincos(p.phi, sinphi, cosphi);
    const double vx = p.vR * cosphi - p.vphi * sinphi;
    const double vy = p.vR * sinphi + p.vphi * cosphi;
    return PosVelCar(p.R * cosphi, p.R * sinphi, p.z, vx, vy, p.vz);
}

template<> PosVelCar toPosVel(const PosVelSph& p, const Car) {
    double sintheta, costheta, sinphi, cosphi;
    math::sincos(p.theta, sintheta, costheta);
    math::sincos(p.phi, sinphi, cosphi);
    const double R = p.r * sintheta, vR = p.vr * sintheta + p.vtheta * costheta;
    const double vx = vR * cosphi - p.vphi * sinphi;
    const double vy = vR * sinphi + p.vphi * cosphi;
    const double vz = p.vr * costheta - p.vtheta * sintheta;
    return PosVelCar(R * cosphi, R * sinphi, p.r * costheta, vx, vy, vz);
}

template<> PosVelCyl toPosVel(const PosVelCar& p, const Cyl) {
    const double R=sqrt(pow_2(p.x) + pow_2(p.y));
    if(R==0)  // determine phi from vy/vx rather than y/x
        return PosVelCyl(R, p.z, math::atan2(p.vy, p.vx), sqrt(pow_2(p.vx) + pow_2(p.vy)), p.vz, 0);
    const double cosphi = p.x / R, sinphi = p.y / R;
    const double vR   = p.vx * cosphi + p.vy * sinphi;
    const double vphi =-p.vx * sinphi + p.vy * cosphi;
    return PosVelCyl(R, p.z, math::atan2(p.y, p.x), vR, p.vz, vphi);
}

template<> PosVelCyl toPosVel(const PosVelSph& p, const Cyl) {
    double sintheta, costheta;
    math::sincos(p.theta, sintheta, costheta);
    const double R  = p.r  * sintheta, z = p.r * costheta;
    const double vR = p.vr * sintheta + p.vtheta * costheta;
    const double vz = p.vr * costheta - p.vtheta * sintheta;
    return PosVelCyl(R, z, p.phi, vR, vz, p.vphi);
}

template<> PosVelSph toPosVel(const PosVelCar& p, const Sph) {
    const double R2 = pow_2(p.x) + pow_2(p.y), R = sqrt(R2), invR = 1/R;
    const double r2 = R2 + pow_2(p.z), r = sqrt(r2), invr = 1/r;
    if(R==0) {  // point along the z axis - determine phi from velocity rather than position
        const double vR = sqrt(pow_2(p.vx) + pow_2(p.vy));
        const double phi = math::atan2(p.vy, p.vx);
        if(p.z==0)  // point at origin - an even more special case
            return PosVelSph(0, math::atan2(vR, p.vz), phi, sqrt(pow_2(vR) + pow_2(p.vz)), 0, 0);
        return PosVelSph(r, p.z>=0 ? 0 : M_PI, phi,
            p.vz * (p.z>=0 ? 1 : -1), vR * (p.z>=0 ? 1 : -1), 0);
    }
    const double temp   = p.x * p.vx + p.y * p.vy;
    const double vr     = (temp + p.z * p.vz) * invr;
    const double vtheta = (temp * p.z * invR - p.vz * R) * invr;
    const double vphi   = (p.x * p.vy - p.y * p.vx) * invR;
    return PosVelSph(r, math::atan2(R, p.z), math::atan2(p.y, p.x), vr, vtheta, vphi);
}

template<> PosVelSph toPosVel(const PosVelCyl& p, const Sph) {
    const double r=sqrt(pow_2(p.R) + pow_2(p.z));
    if(r==0) {
        return PosVelSph(0, math::atan2(p.vR, p.vz), p.phi, sqrt(pow_2(p.vR) + pow_2(p.vz)), 0, 0);
    }
    const double invr = 1./r;
    const double costheta = p.z * invr, sintheta = p.R * invr;
    const double vr = p.vR * sintheta + p.vz * costheta;
    const double vtheta = p.vR * costheta - p.vz * sintheta;
    return PosVelSph(r, math::atan2(p.R, p.z), p.phi, vr, vtheta, p.vphi);
}

template<>
PosVelCyl toPosVel(const PosMomSph& p) {
	double sintheta, costheta;
	math::sincos(p.theta, sintheta, costheta);
	const double R=p.r*sintheta, z=p.r*costheta;
	const double vR=p.pr*sintheta+p.ptheta/p.r*costheta;
	const double vz=p.pr*costheta-p.ptheta/p.r*sintheta;
	return PosVelCyl(R, z, p.phi, vR, vz, p.pphi/(p.r*sintheta));
}

template<> PosMomCyl toPosMom(const PosMomSph& p) {
    double sintheta, costheta;
    math::sincos(p.theta, sintheta, costheta);
	const double R  = p.r  * sintheta, z = p.r * costheta;
    const double pR = p.pr * sintheta + p.ptheta * costheta/p.r;
    const double pz = p.pr * costheta - p.ptheta * sintheta/p.r;
    return PosMomCyl(R, z, p.phi, pR, pz, p.pphi);
}
template<>
PosMomCar toPosMom(const PosMomCyl& p) {
	double csp=cos(p.phi), snp=sin(p.phi), vphi=p.pphi/p.R;
	return PosMomCar(p.R*csp, p.R*snp, p.z,
			 p.pR*csp-vphi*snp,p.pR*snp+vphi*csp,p.pz);
}

template<> PosMomCar toPosMom(const PosMomSph& p) {
	PosMomCyl pC = toPosMom<Sph,Cyl>(p);
	return toPosMom<Cyl, Car>(pC);
}


template<> PosVelProlSph toPosVel(const PosVelCyl& from, const ProlSph cs) {
    PosDerivT<Cyl, ProlSph> derivs;
    const PosProlSph pprol = toPosDeriv<Cyl, ProlSph> (from, &derivs, NULL, cs);
    double lambdadot = derivs.dlambdadR*from.vR + derivs.dlambdadz*from.vz;
    double nudot     = derivs.dnudR    *from.vR + derivs.dnudz    *from.vz;
    double phidot    = from.vphi!=0 ? from.vphi/from.R : 0;
    return PosVelProlSph(pprol, lambdadot, nudot, phidot);
}

template<> PosVelAxi toPosVel(const PosVelCyl& pc, const Axi cs) {
    if(cs.Delta2 == 0) {  // shortcut for the spherical case
        if(pc.R == 0 && pc.z == 0) {
            // degenerate case - cannot determine nu  from the position alone, use velocity instead
            double vel   = sqrt(pow_2(pc.vR) + pow_2(pc.vz));
            double cotnu = vel>0 ? pc.vz / pc.vR : 0;
            return PosVelAxi(PosAxi(0, cotnu, pc.phi, cs), VelAxi(vel, 0, pc.vphi));
        }
        double rho  = sqrt(pow_2(pc.R) + pow_2(pc.z));
        double cosnu= pc.z / rho, sinnu = pc.R / rho;
        double vrho = sinnu * pc.vR + cosnu * pc.vz;
        double vnu  = cosnu * pc.vR - sinnu * pc.vz;
        return PosVelAxi(PosAxi(rho, cosnu / sinnu, pc.phi, cs), VelAxi(vrho, vnu, pc.vphi));
    }
    double chi, psi, cosnu, sinnu;
    getPosAxi(pc, cs, chi, psi, cosnu, sinnu);
    double
    rho  = cs.Delta2 >= 0  ? chi : psi,
    den  = 1 / sqrt(pow_2(chi * cosnu) + pow_2(psi * sinnu)),
    sinxi= den != INFINITY ? den * psi * sinnu : 1.0,
    cosxi= den != INFINITY ? den * chi * cosnu : 0.0,
    vrho = sinxi * pc.vR + cosxi * pc.vz,
    vnu  = cosxi * pc.vR - sinxi * pc.vz;
    return PosVelAxi(PosAxi(rho, cosnu / sinnu, pc.phi, cs), VelAxi(vrho, vnu, pc.vphi));
}

template<> PosVelCyl toPosVel(const PosVelAxi& ps, const Cyl) {
    double
    sinnu = 1 / sqrt(1 + pow_2(ps.cotnu)),
    cosnu = sinnu!=0 ? ps.cotnu * sinnu : (ps.cotnu>0 ? 1 : -1);
    if(ps.cs.Delta2 == 0) {  // shortcut for the spherical case
        if(ps.rho == 0)  // assume that vnu = 0
            return PosVelCyl(0, 0, ps.phi, ps.vrho * sinnu, ps.vrho * cosnu, ps.vphi);
        double R  = ps.rho  * sinnu, z = ps.rho * cosnu;
        double vR = ps.vrho * sinnu + ps.vnu * cosnu;
        double vz = ps.vrho * cosnu - ps.vnu * sinnu;
        return PosVelCyl(R, z, ps.phi, vR, vz, ps.vphi);
    }
    double
    eta  = sqrt(pow_2(ps.rho) + fabs(ps.cs.Delta2)),
    chi  = ps.cs.Delta2 >= 0  ?  ps.rho  :     eta,
    psi  = ps.cs.Delta2 >= 0  ?     eta  :  ps.rho,
    R    = chi * sinnu,
    z    = psi * cosnu,
    den  = 1 / sqrt(pow_2(chi * cosnu) + pow_2(psi * sinnu)),
    sinxi= den != INFINITY ? den * psi * sinnu : 1.0,
    cosxi= den != INFINITY ? den * chi * cosnu : 0.0,
    vR   = sinxi * ps.vrho + cosxi * ps.vnu,
    vz   = cosxi * ps.vrho - sinxi * ps.vnu;
    return PosVelCyl(R, z, ps.phi, vR, vz, ps.vphi);
}
template<> PosVelEll toPosVel(const PosVelCar& from, const Ell ell) {
    PosDerivT<Car, Ell> derivs;
    const PosEll pell = toPosDeriv<Car, Ell> (from, &derivs, NULL, ell);
    double lambdadot = derivs.dlambdadx*from.vx + derivs.dlambdady*from.vy+derivs.dlambdadz*from.vz;
    double mudot = derivs.dmudx*from.vx + derivs.dmudy*from.vy+derivs.dmudz*from.vz;
    double nudot = derivs.dnudx*from.vx + derivs.dnudy*from.vy+derivs.dnudz*from.vz;
    return PosVelEll(pell, VelEll(lambdadot, mudot, nudot));
}

template<> PosVelCar toPosVel(const PosVelEll& ps, const Car) {
    PosDerivT<Ell, Car> derivs;
    const PosCar xyz = toPosDeriv<Ell, Car> (ps,&derivs,NULL);
    double vx = derivs.dxdlambda*ps.vlambda + derivs.dxdmu*ps.vmu+derivs.dxdnu*ps.vnu;
    double vy = derivs.dydlambda*ps.vlambda + derivs.dydmu*ps.vmu+derivs.dydnu*ps.vnu;
    double vz = derivs.dzdlambda*ps.vlambda + derivs.dzdmu*ps.vmu+derivs.dzdnu*ps.vnu;
    return PosVelCar(xyz,VelCar(vx,vy,vz));
}

template<> PosVelEls toPosVel(const PosVelCar& from, const Els els) {
    PosDerivT<Car, Els> derivs;
    const PosEls pell = toPosDeriv<Car, Els> (from, &derivs, NULL, els);
    double rhodot = derivs.drhodx*from.vx + derivs.drhody*from.vy+derivs.drhodz*from.vz;
    double chidot = derivs.dchidx*from.vx + derivs.dchidy*from.vy+derivs.dchidz*from.vz;
    double phidot = derivs.dphidx*from.vx + derivs.dphidy*from.vy+derivs.dphidz*from.vz;
    //double snchi=1/sqrt(1+pow_2(pell.cotchi));
    return PosVelEls(pell, VelEls(rhodot, chidot, phidot));
}

template<> PosVelCar toPosVel(const PosVelEls& ps, const Car) {
    PosDerivT<Els, Car> derivs;
    const PosCar xyz = toPosDeriv<Els, Car> (ps,&derivs,NULL);
    double vx = derivs.dxdrho*ps.rhodot + derivs.dxdchi*ps.chidot + derivs.dxdphi*ps.phidot;
    double vy = derivs.dydrho*ps.rhodot + derivs.dydchi*ps.chidot+derivs.dydphi*ps.phidot;
    double vz = derivs.dzdrho*ps.rhodot + derivs.dzdchi*ps.chidot+derivs.dzdphi*ps.phidot;
    return PosVelCar(xyz,VelCar(vx,vy,vz));
}

void PosVelSph::momenta(double& pr, double& ptheta, double& pphi) const
{
    pr     = vr;
    ptheta = vtheta*r;
    pphi   = Lz(*this);
}

void PosVelAxi::momenta(double& prho, double& pnu, double& pphi) const
{
    double
    sinnu= 1 / sqrt(1 + pow_2(cotnu)),
    eta  = sqrt(pow_2(rho) + fabs(cs.Delta2)),
    chi  = cs.Delta2 >= 0  ?  rho  :  eta,
    mul  = sqrt(pow_2(chi) + cs.Delta2 * pow_2(sinnu));
    prho = (cs.Delta2 == 0 ? 1 : mul / eta) * vrho;
    pnu  = mul * vnu;
    pphi = chi * sinnu * vphi;
}

//-------- implementations of functions that convert gradients --------//

template<>
GradCar toGrad(const GradCyl& src, const PosDerivT<Car, Cyl>& deriv) {
    GradCar dest;
    dest.dx = src.dR*deriv.dRdx + src.dphi*deriv.dphidx;
    dest.dy = src.dR*deriv.dRdy + src.dphi*deriv.dphidy;
    dest.dz = src.dz;
    return dest;
}

template<>
GradCar toGrad(const GradSph& src, const PosDerivT<Car, Sph>& deriv) {
    GradCar dest;
    dest.dx = src.dr*deriv.drdx + src.dtheta*deriv.dthetadx + src.dphi*deriv.dphidx;
    dest.dy = src.dr*deriv.drdy + src.dtheta*deriv.dthetady + src.dphi*deriv.dphidy;
    dest.dz = src.dr*deriv.drdz + src.dtheta*deriv.dthetadz;
    return dest;
}

template<>
GradCyl toGrad(const GradCar& src, const PosDerivT<Cyl, Car>& deriv) {
    GradCyl dest;
    dest.dR = src.dx*deriv.dxdR + src.dy*deriv.dydR;
    dest.dz = src.dz;
    dest.dphi = src.dx*deriv.dxdphi + src.dy*deriv.dydphi;
    return dest;
}

template<>
GradCyl toGrad(const GradSph& src, const PosDerivT<Cyl, Sph>& deriv) {
    GradCyl dest;
    dest.dR = src.dr*deriv.drdR + src.dtheta*deriv.dthetadR;
    dest.dz = src.dr*deriv.drdz + src.dtheta*deriv.dthetadz;
    dest.dphi = src.dphi;
    return dest;
}

template<>
GradSph toGrad(const GradCar& src, const PosDerivT<Sph, Car>& deriv) {
    GradSph dest;
    dest.dr     = src.dx*deriv.dxdr     + src.dy*deriv.dydr     + src.dz*deriv.dzdr;
    dest.dtheta = src.dx*deriv.dxdtheta + src.dy*deriv.dydtheta + src.dz*deriv.dzdtheta;
    dest.dphi   = src.dx*deriv.dxdphi   + src.dy*deriv.dydphi;
    return dest;
}

template<>
GradSph toGrad(const GradCyl& src, const PosDerivT<Sph, Cyl>& deriv) {
    GradSph dest;
    dest.dr     = src.dR*deriv.dRdr     + src.dz*deriv.dzdr;
    dest.dtheta = src.dR*deriv.dRdtheta + src.dz*deriv.dzdtheta;
    dest.dphi   = src.dphi;
    return dest;
}

template<>
GradCyl toGrad(const GradProlSph& src, const PosDerivT<Cyl, ProlSph>& deriv) {
    GradCyl dest;
    dest.dR   = src.dlambda*deriv.dlambdadR + src.dnu*deriv.dnudR;
    dest.dz   = src.dlambda*deriv.dlambdadz + src.dnu*deriv.dnudz;
    dest.dphi = src.dphi;
    return dest;
}

template<>
GradProlSph toGrad(const GradCyl& src, const PosDerivT<ProlSph, Cyl>& deriv) {
    GradProlSph dest;
    dest.dlambda = src.dR*deriv.dRdlambda + src.dz*deriv.dzdlambda;
    dest.dnu     = src.dR*deriv.dRdnu     + src.dz*deriv.dzdnu;
    dest.dphi    = src.dphi;
    return dest;
}

template<>
GradCyl toGrad(const GradAxi& src, const PosDerivT<Cyl, Axi>& deriv) {
    GradCyl dest;
    dest.dR   = src.drho*deriv.drhodR + src.dnu*deriv.dnudR;
    dest.dz   = src.drho*deriv.drhodz + src.dnu*deriv.dnudz;
    dest.dphi = src.dphi;
    return dest;
}

template<>
GradAxi toGrad(const GradCyl& src, const PosDerivT<Axi, Cyl>& deriv) {
    GradAxi dest;
    dest.drho = src.dR*deriv.dRdrho + src.dz*deriv.dzdrho;
    dest.dnu  = src.dR*deriv.dRdnu  + src.dz*deriv.dzdnu;
    dest.dphi = src.dphi;
    return dest;
}
template<>
GradCar toGrad(const GradEll& src, const PosDerivT<Car, Ell>& deriv) {
    GradCar dest;
    dest.dx   = src.dlambda*deriv.dlambdadx+src.dmu*deriv.dmudx+src.dnu*deriv.dnudx;
    dest.dy   = src.dlambda*deriv.dlambdady+src.dmu*deriv.dmudy+src.dnu*deriv.dnudy;
    dest.dz   = src.dlambda*deriv.dlambdadz+src.dmu*deriv.dmudz+src.dnu*deriv.dnudz;
    return dest;
}

template<>
GradEll toGrad(const GradCar& src, const PosDerivT<Ell, Car>& deriv) {
    GradEll dest;
    dest.dlambda = src.dx*deriv.dxdlambda+src.dy*deriv.dydlambda + src.dz*deriv.dzdlambda;
    dest.dmu  = src.dx*deriv.dxdmu+src.dy*deriv.dydmu  + src.dz*deriv.dzdmu;
    dest.dnu  = src.dx*deriv.dxdnu+src.dy*deriv.dydnu  + src.dz*deriv.dzdnu;
    return dest;
}

template<>
GradCar toGrad(const GradEls& src, const PosDerivT<Car, Els>& deriv) {
    GradCar dest;
    dest.dx   = src.drho*deriv.drhodx+src.dchi*deriv.dchidx+src.dphi*deriv.dphidx;
    dest.dy   = src.drho*deriv.drhody+src.dchi*deriv.dchidy+src.dphi*deriv.dphidy;
    dest.dz   = src.drho*deriv.drhodz+src.dchi*deriv.dchidz+src.dphi*deriv.dphidz;
    return dest;
}

template<>
GradEls toGrad(const GradCar& src, const PosDerivT<Els, Car>& deriv) {
    GradEls dest;
    dest.drho = src.dx*deriv.dxdrho+src.dy*deriv.dydrho + src.dz*deriv.dzdrho;
    dest.dchi  = src.dx*deriv.dxdchi+src.dy*deriv.dydchi  + src.dz*deriv.dzdchi;
    dest.dphi = src.dx*deriv.dxdphi+src.dy*deriv.dydphi  + src.dz*deriv.dzdphi;
    return dest;
}

//-------- implementations of functions that convert hessians --------//

template<>
HessCar toHess(const GradCyl& srcGrad, const HessCyl& srcHess,
    const PosDerivT<Car, Cyl>& deriv, const PosDeriv2T<Car, Cyl>& deriv2) {
    HessCar dest;
    dest.dx2 =
        (srcHess.dR2   *deriv.dRdx + srcHess.dRdphi*deriv.dphidx) * deriv.dRdx +
        (srcHess.dRdphi*deriv.dRdx + srcHess.dphi2 *deriv.dphidx) * deriv.dphidx +
        srcGrad.dR*deriv2.d2Rdx2   + srcGrad.dphi*deriv2.d2phidx2;
    dest.dxdy =
        (srcHess.dR2   *deriv.dRdy + srcHess.dRdphi*deriv.dphidy) * deriv.dRdx +
        (srcHess.dRdphi*deriv.dRdy + srcHess.dphi2 *deriv.dphidy) * deriv.dphidx +
        srcGrad.dR*deriv2.d2Rdxdy  + srcGrad.dphi*deriv2.d2phidxdy;
    dest.dy2 =
        (srcHess.dR2   *deriv.dRdy + srcHess.dRdphi*deriv.dphidy) * deriv.dRdy +
        (srcHess.dRdphi*deriv.dRdy + srcHess.dphi2 *deriv.dphidy) * deriv.dphidy +
        srcGrad.dR*deriv2.d2Rdy2   + srcGrad.dphi*deriv2.d2phidy2;
    dest.dxdz = srcHess.dRdz*deriv.dRdx + srcHess.dzdphi*deriv.dphidx;
    dest.dydz = srcHess.dRdz*deriv.dRdy + srcHess.dzdphi*deriv.dphidy;
    dest.dz2  = srcHess.dz2;
    return dest;
}

template<>
HessCar toHess(const GradSph& srcGrad, const HessSph& srcHess,
    const PosDerivT<Car, Sph>& deriv, const PosDeriv2T<Car, Sph>& deriv2) {
    HessCar dest;
    dest.dx2 =
        (srcHess.dr2     *deriv.drdx + srcHess.drdtheta  *deriv.dthetadx + srcHess.drdphi    *deriv.dphidx) * deriv.drdx +
        (srcHess.drdtheta*deriv.drdx + srcHess.dtheta2   *deriv.dthetadx + srcHess.dthetadphi*deriv.dphidx) * deriv.dthetadx +
        (srcHess.drdphi  *deriv.drdx + srcHess.dthetadphi*deriv.dthetadx + srcHess.dphi2     *deriv.dphidx) * deriv.dphidx +
        srcGrad.dr*deriv2.d2rdx2     + srcGrad.dtheta*deriv2.d2thetadx2  + srcGrad.dphi*deriv2.d2phidx2;
    dest.dxdy =
        (srcHess.dr2     *deriv.drdy + srcHess.drdtheta  *deriv.dthetady + srcHess.drdphi    *deriv.dphidy) * deriv.drdx +
        (srcHess.drdtheta*deriv.drdy + srcHess.dtheta2   *deriv.dthetady + srcHess.dthetadphi*deriv.dphidy) * deriv.dthetadx +
        (srcHess.drdphi  *deriv.drdy + srcHess.dthetadphi*deriv.dthetady + srcHess.dphi2     *deriv.dphidy) * deriv.dphidx +
        srcGrad.dr*deriv2.d2rdxdy    + srcGrad.dtheta*deriv2.d2thetadxdy + srcGrad.dphi*deriv2.d2phidxdy;
    dest.dxdz =
        (srcHess.dr2     *deriv.drdz + srcHess.drdtheta  *deriv.dthetadz) * deriv.drdx +
        (srcHess.drdtheta*deriv.drdz + srcHess.dtheta2   *deriv.dthetadz) * deriv.dthetadx +
        (srcHess.drdphi  *deriv.drdz + srcHess.dthetadphi*deriv.dthetadz) * deriv.dphidx +
        srcGrad.dr*deriv2.d2rdxdz    + srcGrad.dtheta*deriv2.d2thetadxdz;
    dest.dy2 =
        (srcHess.dr2     *deriv.drdy + srcHess.drdtheta  *deriv.dthetady + srcHess.drdphi    *deriv.dphidy) * deriv.drdy +
        (srcHess.drdtheta*deriv.drdy + srcHess.dtheta2   *deriv.dthetady + srcHess.dthetadphi*deriv.dphidy) * deriv.dthetady +
        (srcHess.drdphi  *deriv.drdy + srcHess.dthetadphi*deriv.dthetady + srcHess.dphi2     *deriv.dphidy) * deriv.dphidy +
        srcGrad.dr*deriv2.d2rdy2     + srcGrad.dtheta*deriv2.d2thetady2  + srcGrad.dphi*deriv2.d2phidy2;
    dest.dydz =
        (srcHess.dr2     *deriv.drdz + srcHess.drdtheta  *deriv.dthetadz) * deriv.drdy +
        (srcHess.drdtheta*deriv.drdz + srcHess.dtheta2   *deriv.dthetadz) * deriv.dthetady +
        (srcHess.drdphi  *deriv.drdz + srcHess.dthetadphi*deriv.dthetadz) * deriv.dphidy +
        srcGrad.dr*deriv2.d2rdydz    + srcGrad.dtheta*deriv2.d2thetadydz;
    dest.dz2 =
        (srcHess.dr2     *deriv.drdz + srcHess.drdtheta  *deriv.dthetadz) * deriv.drdz +
        (srcHess.drdtheta*deriv.drdz + srcHess.dtheta2   *deriv.dthetadz) * deriv.dthetadz +
        srcGrad.dr*deriv2.d2rdz2     + srcGrad.dtheta*deriv2.d2thetadz2;
    return dest;
}

template<>
HessCyl toHess(const GradCar& srcGrad, const HessCar& srcHess,
    const PosDerivT<Cyl, Car>& deriv, const PosDeriv2T<Cyl, Car>& deriv2) {
    HessCyl dest;
    dest.dR2 =
        (srcHess.dx2 *deriv.dxdR + srcHess.dxdy*deriv.dydR) * deriv.dxdR +
        (srcHess.dxdy*deriv.dxdR + srcHess.dy2 *deriv.dydR) * deriv.dydR;
    dest.dRdz = srcHess.dxdz*deriv.dxdR + srcHess.dydz*deriv.dydR;
    dest.dRdphi =
        (srcHess.dx2 *deriv.dxdphi + srcHess.dxdy*deriv.dydphi) * deriv.dxdR +
        (srcHess.dxdy*deriv.dxdphi + srcHess.dy2 *deriv.dydphi) * deriv.dydR +
        srcGrad.dx*deriv2.d2xdRdphi + srcGrad.dy*deriv2.d2ydRdphi;
    dest.dz2 = srcHess.dz2;
    dest.dzdphi = (srcHess.dxdz*deriv.dxdphi + srcHess.dydz*deriv.dydphi);
    dest.dphi2 =
        (srcHess.dx2 *deriv.dxdphi + srcHess.dxdy*deriv.dydphi) * deriv.dxdphi +
        (srcHess.dxdy*deriv.dxdphi + srcHess.dy2 *deriv.dydphi) * deriv.dydphi +
        srcGrad.dx*deriv2.d2xdphi2 + srcGrad.dy*deriv2.d2ydphi2;
    return dest;
}

template<>
HessCyl toHess(const GradSph& srcGrad, const HessSph& srcHess,
    const PosDerivT<Cyl, Sph>& deriv, const PosDeriv2T<Cyl, Sph>& deriv2) {
    HessCyl dest;
    dest.dR2 =
        (srcHess.dr2     *deriv.drdR + srcHess.drdtheta*deriv.dthetadR) * deriv.drdR +
        (srcHess.drdtheta*deriv.drdR + srcHess.dtheta2 *deriv.dthetadR) * deriv.dthetadR +
        srcGrad.dr*deriv2.d2rdR2 + srcGrad.dtheta*deriv2.d2thetadR2;
    dest.dRdz =
        (srcHess.dr2     *deriv.drdz + srcHess.drdtheta*deriv.dthetadz) * deriv.drdR +
        (srcHess.drdtheta*deriv.drdz + srcHess.dtheta2 *deriv.dthetadz) * deriv.dthetadR +
        srcGrad.dr*deriv2.d2rdRdz + srcGrad.dtheta*deriv2.d2thetadRdz;
    dest.dz2 =
        (srcHess.dr2     *deriv.drdz + srcHess.drdtheta*deriv.dthetadz) * deriv.drdz +
        (srcHess.drdtheta*deriv.drdz + srcHess.dtheta2 *deriv.dthetadz) * deriv.dthetadz +
        srcGrad.dr*deriv2.d2rdz2 + srcGrad.dtheta*deriv2.d2thetadz2;
    dest.dRdphi = srcHess.drdphi*deriv.drdR + srcHess.dthetadphi*deriv.dthetadR;
    dest.dzdphi = srcHess.drdphi*deriv.drdz + srcHess.dthetadphi*deriv.dthetadz;
    dest.dphi2  = srcHess.dphi2;
    return dest;
}

template<>
HessSph toHess(const GradCar& srcGrad, const HessCar& srcHess,
    const PosDerivT<Sph, Car>& deriv, const PosDeriv2T<Sph, Car>& deriv2) {
    HessSph dest;
    dest.dr2 =
        (srcHess.dx2 *deriv.dxdr + srcHess.dxdy*deriv.dydr + srcHess.dxdz*deriv.dzdr) * deriv.dxdr +
        (srcHess.dxdy*deriv.dxdr + srcHess.dy2 *deriv.dydr + srcHess.dydz*deriv.dzdr) * deriv.dydr +
        (srcHess.dxdz*deriv.dxdr + srcHess.dydz*deriv.dydr + srcHess.dz2 *deriv.dzdr) * deriv.dzdr;
    dest.drdtheta =
        (srcHess.dx2 *deriv.dxdtheta + srcHess.dxdy*deriv.dydtheta + srcHess.dxdz*deriv.dzdtheta) * deriv.dxdr +
        (srcHess.dxdy*deriv.dxdtheta + srcHess.dy2 *deriv.dydtheta + srcHess.dydz*deriv.dzdtheta) * deriv.dydr +
        (srcHess.dxdz*deriv.dxdtheta + srcHess.dydz*deriv.dydtheta + srcHess.dz2 *deriv.dzdtheta) * deriv.dzdr +
        srcGrad.dx*deriv2.d2xdrdtheta + srcGrad.dy*deriv2.d2ydrdtheta + srcGrad.dz*deriv2.d2zdrdtheta;
    dest.drdphi =
        (srcHess.dx2 *deriv.dxdphi + srcHess.dxdy*deriv.dydphi)*deriv.dxdr +
        (srcHess.dxdy*deriv.dxdphi + srcHess.dy2 *deriv.dydphi)*deriv.dydr +
        (srcHess.dxdz*deriv.dxdphi + srcHess.dydz*deriv.dydphi)*deriv.dzdr +
        srcGrad.dx*deriv2.d2xdrdphi + srcGrad.dy*deriv2.d2ydrdphi;
    dest.dtheta2 =
        (srcHess.dx2 *deriv.dxdtheta + srcHess.dxdy*deriv.dydtheta + srcHess.dxdz*deriv.dzdtheta) * deriv.dxdtheta +
        (srcHess.dxdy*deriv.dxdtheta + srcHess.dy2 *deriv.dydtheta + srcHess.dydz*deriv.dzdtheta) * deriv.dydtheta +
        (srcHess.dxdz*deriv.dxdtheta + srcHess.dydz*deriv.dydtheta + srcHess.dz2 *deriv.dzdtheta) * deriv.dzdtheta +
        srcGrad.dx*deriv2.d2xdtheta2 + srcGrad.dy*deriv2.d2ydtheta2 + srcGrad.dz*deriv2.d2zdtheta2;
    dest.dthetadphi =
        (srcHess.dx2 *deriv.dxdphi + srcHess.dxdy*deriv.dydphi) * deriv.dxdtheta +
        (srcHess.dxdy*deriv.dxdphi + srcHess.dy2 *deriv.dydphi) * deriv.dydtheta +
        (srcHess.dxdz*deriv.dxdphi + srcHess.dydz*deriv.dydphi) * deriv.dzdtheta +
        srcGrad.dx*deriv2.d2xdthetadphi + srcGrad.dy*deriv2.d2ydthetadphi;
    dest.dphi2 =
        (srcHess.dx2 *deriv.dxdphi + srcHess.dxdy*deriv.dydphi) * deriv.dxdphi +
        (srcHess.dxdy*deriv.dxdphi + srcHess.dy2 *deriv.dydphi) * deriv.dydphi +
        srcGrad.dx*deriv2.d2xdphi2 + srcGrad.dy*deriv2.d2ydphi2;
    return dest;
}

template<>
HessSph toHess(const GradCyl& srcGrad, const HessCyl& srcHess,
    const PosDerivT<Sph, Cyl>& deriv, const PosDeriv2T<Sph, Cyl>& deriv2) {
    HessSph dest;
    dest.dr2 =
        (srcHess.dR2 *deriv.dRdr + srcHess.dRdz*deriv.dzdr) * deriv.dRdr +
        (srcHess.dRdz*deriv.dRdr + srcHess.dz2 *deriv.dzdr) * deriv.dzdr;
    dest.drdtheta =
        (srcHess.dR2 *deriv.dRdtheta + srcHess.dRdz*deriv.dzdtheta) * deriv.dRdr +
        (srcHess.dRdz*deriv.dRdtheta + srcHess.dz2 *deriv.dzdtheta) * deriv.dzdr +
        srcGrad.dR*deriv2.d2Rdrdtheta + srcGrad.dz*deriv2.d2zdrdtheta;
    dest.dtheta2 =
        (srcHess.dR2 *deriv.dRdtheta + srcHess.dRdz*deriv.dzdtheta) * deriv.dRdtheta +
        (srcHess.dRdz*deriv.dRdtheta + srcHess.dz2 *deriv.dzdtheta) * deriv.dzdtheta +
        srcGrad.dR*deriv2.d2Rdtheta2 + srcGrad.dz*deriv2.d2zdtheta2;
    dest.drdphi     = srcHess.dRdphi*deriv.dRdr     + srcHess.dzdphi*deriv.dzdr;
    dest.dthetadphi = srcHess.dRdphi*deriv.dRdtheta + srcHess.dzdphi*deriv.dzdtheta;
    dest.dphi2      = srcHess.dphi2;
    return dest;
}

//TODO// remove
template<>
HessCyl toHess(const GradProlSph& srcGrad, const HessProlSph& srcHess,
    const PosDerivT<Cyl, ProlSph>& deriv, const PosDeriv2T<Cyl, ProlSph>& deriv2) {
    HessCyl dest;
    dest.dR2 =
        (srcHess.dlambda2*deriv.dlambdadR + srcHess.dlambdadnu*deriv.dnudR)*deriv.dlambdadR +
        (srcHess.dlambdadnu*deriv.dlambdadR + srcHess.dnu2*deriv.dnudR)*deriv.dnudR +
        srcGrad.dlambda*deriv2.d2lambdadR2 + srcGrad.dnu*deriv2.d2nudR2;
    dest.dRdz =
        (srcHess.dlambda2*deriv.dlambdadz + srcHess.dlambdadnu*deriv.dnudz)*deriv.dlambdadR +
        (srcHess.dlambdadnu*deriv.dlambdadz + srcHess.dnu2*deriv.dnudz)*deriv.dnudR +
        srcGrad.dlambda*deriv2.d2lambdadRdz + srcGrad.dnu*deriv2.d2nudRdz;
    dest.dz2 =
        (srcHess.dlambda2*deriv.dlambdadz + srcHess.dlambdadnu*deriv.dnudz)*deriv.dlambdadz +
        (srcHess.dlambdadnu*deriv.dlambdadz + srcHess.dnu2*deriv.dnudz)*deriv.dnudz +
        srcGrad.dlambda*deriv2.d2lambdadz2 + srcGrad.dnu*deriv2.d2nudz2;
    dest.dRdphi = dest.dzdphi = dest.dphi2 = 0;  //TODO// assuming no dependence on phi
    return dest;
}

template<>
HessCyl toHess(const GradAxi& srcGrad, const HessAxi& srcHess,
    const PosDerivT<Cyl, Axi>& deriv, const PosDeriv2T<Cyl, Axi>& deriv2) {
    HessCyl dest;
    dest.dR2 =
        (srcHess.drho2  *deriv.drhodR + srcHess.drhodnu*deriv.dnudR) * deriv.drhodR +
        (srcHess.drhodnu*deriv.drhodR + srcHess.dnu2   *deriv.dnudR) * deriv.dnudR  +
        srcGrad.drho *deriv2.d2rhodR2 + srcGrad.dnu * deriv2.d2nudR2;
    dest.dRdz =
        (srcHess.drho2  *deriv.drhodz + srcHess.drhodnu*deriv.dnudz) * deriv.drhodR +
        (srcHess.drhodnu*deriv.drhodz + srcHess.dnu2   *deriv.dnudz) * deriv.dnudR  +
        srcGrad.drho *deriv2.d2rhodRdz+ srcGrad.dnu * deriv2.d2nudRdz;
    dest.dz2 =
        (srcHess.drho2  *deriv.drhodz + srcHess.drhodnu*deriv.dnudz) * deriv.drhodz +
        (srcHess.drhodnu*deriv.drhodz + srcHess.dnu2   *deriv.dnudz) * deriv.dnudz  +
        srcGrad.drho *deriv2.d2rhodz2 + srcGrad.dnu * deriv2.d2nudz2;
    dest.dRdphi = srcHess.drhodphi*deriv.drhodR + srcHess.dnudphi*deriv.dnudR;
    dest.dzdphi = srcHess.drhodphi*deriv.drhodz + srcHess.dnudphi*deriv.dnudz;
    dest.dphi2  = srcHess.dphi2;
    return dest;
}

template<>
HessAxi toHess(const GradCyl& srcGrad, const HessCyl& srcHess,
    const PosDerivT<Axi, Cyl>& deriv, const PosDeriv2T<Axi, Cyl>& deriv2) {
    HessAxi dest;
    dest.drho2 =
        (srcHess.dR2 *deriv.dRdrho + srcHess.dRdz*deriv.dzdrho) * deriv.dRdrho +
        (srcHess.dRdz*deriv.dRdrho + srcHess.dz2 *deriv.dzdrho) * deriv.dzdrho +
        srcGrad.dR*deriv2.d2Rdrho2 + srcGrad.dz*deriv2.d2zdrho2;
    dest.drhodnu =
        (srcHess.dR2 *deriv.dRdnu + srcHess.dRdz*deriv.dzdnu) * deriv.dRdrho +
        (srcHess.dRdz*deriv.dRdnu + srcHess.dz2 *deriv.dzdnu) * deriv.dzdrho +
        srcGrad.dR*deriv2.d2Rdrhodnu + srcGrad.dz*deriv2.d2zdrhodnu;
    dest.dnu2 =
        (srcHess.dR2 *deriv.dRdnu + srcHess.dRdz*deriv.dzdnu) * deriv.dRdnu +
        (srcHess.dRdz*deriv.dRdnu + srcHess.dz2 *deriv.dzdnu) * deriv.dzdnu +
        srcGrad.dR*deriv2.d2Rdnu2 + srcGrad.dz*deriv2.d2zdnu2;
    dest.drhodphi = srcHess.dRdphi*deriv.dRdrho + srcHess.dzdphi*deriv.dzdrho;
    dest.dnudphi  = srcHess.dRdphi*deriv.dRdnu  + srcHess.dzdphi*deriv.dzdnu;
    dest.dphi2    = srcHess.dphi2;
    return dest;
}

template<>
HessCar toHess(const GradEll& srcGrad, const HessEll& srcHess,
    const PosDerivT<Car, Ell>& deriv, const PosDeriv2T<Car, Ell>& deriv2) {
    HessCar dest;
    dest.dx2 =
        (srcHess.dlambda2*deriv.dlambdadx+srcHess.dlambdadmu*deriv.dmudx+srcHess.dlambdadnu*deriv.dnudx)*deriv.dlambdadx+
        (srcHess.dlambdadmu*deriv.dlambdadx+srcHess.dmu2*deriv.dmudx+srcHess.dmudnu*deriv.dnudx)*deriv.dmudx+
        (srcHess.dlambdadnu*deriv.dlambdadx+srcHess.dmudnu*deriv.dmudx+srcHess.dnu2*deriv.dnudx)*deriv.dnudx+
        srcGrad.dlambda*deriv2.d2lambdadx2+srcGrad.dmu*deriv2.d2mudx2+srcGrad.dnu*deriv2.d2nudx2;
    dest.dxdy =
        (srcHess.dlambda2*deriv.dlambdady+srcHess.dlambdadmu*deriv.dmudy+srcHess.dlambdadnu*deriv.dnudy)*deriv.dlambdadx+
        (srcHess.dlambdadnu*deriv.dlambdady+srcHess.dmudnu*deriv.dmudy+srcHess.dnu2*deriv.dnudy)*deriv.dnudx+
        (srcHess.dlambdadmu*deriv.dlambdady+srcHess.dmudnu*deriv.dnudy+srcHess.dmu2*deriv.dmudy)*deriv.dmudx+
        srcGrad.dlambda*deriv2.d2lambdadxdy+srcGrad.dmu*deriv2.d2mudxdy+srcGrad.dnu*deriv2.d2nudxdy;
    dest.dxdz =
        (srcHess.dlambda2*deriv.dlambdadz+srcHess.dlambdadmu*deriv.dmudz+srcHess.dlambdadnu*deriv.dnudz)*deriv.dlambdadx+
        (srcHess.dlambdadnu*deriv.dlambdadz+srcHess.dmudnu*deriv.dmudz+srcHess.dnu2*deriv.dnudz)*deriv.dnudx+
        (srcHess.dlambdadmu*deriv.dlambdadz+srcHess.dmudnu*deriv.dnudz+srcHess.dmu2*deriv.dmudz)*deriv.dmudx+
        srcGrad.dlambda*deriv2.d2lambdadxdz+srcGrad.dmu*deriv2.d2mudxdz+srcGrad.dnu*deriv2.d2nudxdz;
    dest.dy2 =
        (srcHess.dlambda2*deriv.dlambdady+srcHess.dlambdadmu*deriv.dmudy+srcHess.dlambdadnu*deriv.dnudy)*deriv.dlambdady+
        (srcHess.dlambdadmu*deriv.dlambdady+srcHess.dmu2*deriv.dmudy+srcHess.dmudnu*deriv.dnudy)*deriv.dmudy+
        (srcHess.dlambdadnu*deriv.dlambdady+srcHess.dmudnu*deriv.dmudy+srcHess.dnu2*deriv.dnudy)*deriv.dnudy+
        srcGrad.dlambda*deriv2.d2lambdady2+srcGrad.dmu*deriv2.d2mudy2+srcGrad.dnu*deriv2.d2nudy2;
    dest.dydz =
        (srcHess.dlambda2*deriv.dlambdady+srcHess.dlambdadmu*deriv.dmudy+srcHess.dlambdadnu*deriv.dnudy)*deriv.dlambdadz+
        (srcHess.dlambdadnu*deriv.dlambdady+srcHess.dmudnu*deriv.dmudy+srcHess.dnu2*deriv.dnudy)*deriv.dnudz+
        (srcHess.dlambdadmu*deriv.dlambdady+srcHess.dmudnu*deriv.dnudy+srcHess.dmu2*deriv.dmudy)*deriv.dmudz+
        srcGrad.dlambda*deriv2.d2lambdadydz+srcGrad.dmu*deriv2.d2mudydz+srcGrad.dnu*deriv2.d2nudydz;
    dest.dz2 =
        (srcHess.dlambda2*deriv.dlambdadz+srcHess.dlambdadmu*deriv.dmudz+srcHess.dlambdadnu*deriv.dnudz)*deriv.dlambdadz+
        (srcHess.dlambdadmu*deriv.dlambdadz+srcHess.dmu2*deriv.dmudz+srcHess.dmudnu*deriv.dnudz)*deriv.dmudz+
        (srcHess.dlambdadnu*deriv.dlambdadz+srcHess.dmudnu*deriv.dmudz+srcHess.dnu2*deriv.dnudz)*deriv.dnudz+
        srcGrad.dlambda*deriv2.d2lambdadz2+srcGrad.dmu*deriv2.d2mudz2+srcGrad.dnu*deriv2.d2nudz2;
    return dest;
}
template<>
HessEll toHess(const GradCar& srcGrad, const HessCar& srcHess,
    const PosDerivT<Ell, Car>& deriv, const PosDeriv2T<Ell, Car>& deriv2) {
    HessEll dest;
    dest.dlambda2 =
        (srcHess.dx2 *deriv.dxdlambda+srcHess.dxdy*deriv.dydlambda+srcHess.dxdz*deriv.dzdlambda)*deriv.dxdlambda+
        (srcHess.dxdy*deriv.dxdlambda+srcHess.dy2*deriv.dydlambda+srcHess.dydz*deriv.dzdlambda)*deriv.dydlambda+
        (srcHess.dxdz*deriv.dxdlambda+srcHess.dydz*deriv.dydlambda+srcHess.dz2*deriv.dzdlambda)*deriv.dzdlambda+
        srcGrad.dx*deriv2.d2xdlambda2+srcGrad.dy*deriv2.d2ydlambda2+srcGrad.dz*deriv2.d2zdlambda2;
    dest.dlambdadmu =
        (srcHess.dx2 *deriv.dxdmu + srcHess.dxdy*deriv.dydmu + srcHess.dxdz*deriv.dzdmu) * deriv.dxdlambda +
        (srcHess.dxdy*deriv.dxdmu + srcHess.dy2*deriv.dydmu + srcHess.dydz *deriv.dzdmu) * deriv.dydlambda +
        (srcHess.dxdz*deriv.dxdmu + srcHess.dydz*deriv.dydmu + srcHess.dz2 *deriv.dzdmu) * deriv.dzdlambda +
        srcGrad.dx*deriv2.d2xdlambdadmu + srcGrad.dy*deriv2.d2ydlambdadmu+ srcGrad.dz*deriv2.d2zdlambdadmu;
    dest.dlambdadnu =
        (srcHess.dx2 *deriv.dxdnu + srcHess.dxdy*deriv.dydnu + srcHess.dxdz*deriv.dzdnu) * deriv.dxdlambda +
        (srcHess.dxdy*deriv.dxdnu + srcHess.dy2*deriv.dydnu + srcHess.dydz *deriv.dzdnu) * deriv.dydlambda +
        (srcHess.dxdz*deriv.dxdnu + srcHess.dydz*deriv.dydnu + srcHess.dz2 *deriv.dzdnu) * deriv.dzdlambda +
        srcGrad.dx*deriv2.d2xdlambdadnu + srcGrad.dy*deriv2.d2ydlambdadnu+ srcGrad.dz*deriv2.d2zdlambdadnu;
    dest.dmu2 =
        (srcHess.dx2 *deriv.dxdmu + srcHess.dxdy*deriv.dydmu + srcHess.dxdz*deriv.dzdmu) * deriv.dxdmu +
        (srcHess.dxdy*deriv.dxdmu + srcHess.dy2*deriv.dydmu + srcHess.dydz *deriv.dzdmu) * deriv.dydmu +
        (srcHess.dxdz*deriv.dxdmu + srcHess.dydz*deriv.dydmu + srcHess.dz2 *deriv.dzdmu) * deriv.dzdmu +
        srcGrad.dx*deriv2.d2xdmu2 + srcGrad.dy*deriv2.d2ydmu2+ srcGrad.dz*deriv2.d2zdmu2;
    dest.dmudnu =
        (srcHess.dx2 *deriv.dxdmu + srcHess.dxdy*deriv.dydmu + srcHess.dxdz*deriv.dzdmu) * deriv.dxdnu +
        (srcHess.dxdy*deriv.dxdmu + srcHess.dy2*deriv.dydmu + srcHess.dydz *deriv.dzdmu) * deriv.dydnu +
        (srcHess.dxdz*deriv.dxdmu + srcHess.dydz*deriv.dydmu + srcHess.dz2 *deriv.dzdmu) * deriv.dzdnu +
        srcGrad.dx*deriv2.d2xdmudnu + srcGrad.dy*deriv2.d2ydmudnu+ srcGrad.dz*deriv2.d2zdmudnu;
    dest.dnu2 =
        (srcHess.dx2 *deriv.dxdnu + srcHess.dxdy*deriv.dydnu + srcHess.dxdz*deriv.dzdnu) * deriv.dxdnu +
        (srcHess.dxdy*deriv.dxdnu + srcHess.dy2*deriv.dydnu + srcHess.dydz *deriv.dzdnu) * deriv.dydnu +
        (srcHess.dxdz*deriv.dxdnu + srcHess.dydz*deriv.dydnu + srcHess.dz2 *deriv.dzdnu) * deriv.dzdnu +
        srcGrad.dx*deriv2.d2xdnu2 + srcGrad.dy*deriv2.d2ydnu2+ srcGrad.dz*deriv2.d2zdnu2;
    return dest;
}

template<>
HessCar toHess(const GradEls& srcGrad, const HessEls& srcHess,
    const PosDerivT<Car, Els>& deriv, const PosDeriv2T<Car, Els>& deriv2) {
    HessCar dest;
    dest.dx2 =
        (srcHess.drho2*deriv.drhodx+srcHess.drhodchi*deriv.dchidx+srcHess.drhodphi*deriv.dphidx)*deriv.drhodx+
        (srcHess.drhodchi*deriv.drhodx+srcHess.dchi2*deriv.dchidx+srcHess.dchidphi*deriv.dphidx)*deriv.dchidx+
        (srcHess.drhodphi*deriv.drhodx+srcHess.dchidphi*deriv.dchidx+srcHess.dphi2*deriv.dphidx)*deriv.dphidx+
        srcGrad.drho*deriv2.d2rhodx2+srcGrad.dchi*deriv2.d2chidx2+srcGrad.dphi*deriv2.d2phidx2;
    dest.dxdy =
        (srcHess.drho2*deriv.drhody+srcHess.drhodchi*deriv.dchidy+srcHess.drhodphi*deriv.dphidy)*deriv.drhodx+
        (srcHess.drhodchi*deriv.drhody+srcHess.dchi2*deriv.dchidy+srcHess.dchidphi*deriv.dphidy)*deriv.dchidx+
        (srcHess.drhodphi*deriv.drhody+srcHess.dchidphi*deriv.dchidy+srcHess.dphi2*deriv.dphidy)*deriv.dphidx+
        srcGrad.drho*deriv2.d2rhodxdy+srcGrad.dchi*deriv2.d2chidxdy+srcGrad.dphi*deriv2.d2phidxdy;
    dest.dxdz =
        (srcHess.drho2*deriv.drhodz+srcHess.drhodchi*deriv.dchidz+srcHess.drhodphi*deriv.dphidz)*deriv.drhodx+
        (srcHess.drhodchi*deriv.drhodz+srcHess.dchidphi*deriv.dphidz+srcHess.dchi2*deriv.dchidz)*deriv.dchidx+
        (srcHess.drhodphi*deriv.drhodz+srcHess.dchidphi*deriv.dchidz+srcHess.dphi2*deriv.dphidz)*deriv.dphidx+
        srcGrad.drho*deriv2.d2rhodxdz+srcGrad.dchi*deriv2.d2chidxdz+srcGrad.dphi*deriv2.d2phidxdz;
    dest.dy2 =
        (srcHess.drho2*deriv.drhody+srcHess.drhodchi*deriv.dchidy+srcHess.drhodphi*deriv.dphidy)*deriv.drhody+
        (srcHess.drhodchi*deriv.drhody+srcHess.dchi2*deriv.dchidy+srcHess.dchidphi*deriv.dphidy)*deriv.dchidy+
        (srcHess.drhodphi*deriv.drhody+srcHess.dchidphi*deriv.dchidy+srcHess.dphi2*deriv.dphidy)*deriv.dphidy+
        srcGrad.drho*deriv2.d2rhody2+srcGrad.dchi*deriv2.d2chidy2+srcGrad.dphi*deriv2.d2phidy2;
    dest.dydz =
        (srcHess.drho2*deriv.drhody+srcHess.drhodchi*deriv.dchidy+srcHess.drhodphi*deriv.dphidy)*deriv.drhodz+
        (srcHess.drhodchi*deriv.drhody+srcHess.dchi2*deriv.dchidy+srcHess.dchidphi*deriv.dphidy)*deriv.dchidz+
        (srcHess.drhodphi*deriv.drhody+srcHess.dchidphi*deriv.dchidy+srcHess.dphi2*deriv.dphidy)*deriv.dphidz+
        srcGrad.drho*deriv2.d2rhodydz+srcGrad.dchi*deriv2.d2chidydz+srcGrad.dphi*deriv2.d2phidydz;
    dest.dz2 =
        (srcHess.drho2*deriv.drhodz+srcHess.drhodchi*deriv.dchidz+srcHess.drhodphi*deriv.dphidz)*deriv.drhodz+
        (srcHess.drhodchi*deriv.drhodz+srcHess.dchi2*deriv.dchidz+srcHess.dchidphi*deriv.dphidz)*deriv.dchidz+
        (srcHess.drhodphi*deriv.drhodz+srcHess.dchidphi*deriv.dchidz+srcHess.dphi2*deriv.dphidz)*deriv.dphidz+
        srcGrad.drho*deriv2.d2rhodz2+srcGrad.dchi*deriv2.d2chidz2+srcGrad.dphi*deriv2.d2phidz2;
    return dest;
}
template<>
HessEls toHess(const GradCar& srcGrad, const HessCar& srcHess,
    const PosDerivT<Els, Car>& deriv, const PosDeriv2T<Els, Car>& deriv2) {
    HessEls dest;
    dest.drho2 =
        (srcHess.dx2 *deriv.dxdrho+srcHess.dxdy*deriv.dydrho+srcHess.dxdz*deriv.dzdrho)*deriv.dxdrho+
        (srcHess.dxdy*deriv.dxdrho+srcHess.dy2*deriv.dydrho+srcHess.dydz*deriv.dzdrho)*deriv.dydrho+
        (srcHess.dxdz*deriv.dxdrho+srcHess.dydz*deriv.dydrho+srcHess.dz2*deriv.dzdrho)*deriv.dzdrho+
        srcGrad.dx*deriv2.d2xdrho2+srcGrad.dy*deriv2.d2ydrho2+srcGrad.dz*deriv2.d2zdrho2;
    dest.drhodchi=
        (srcHess.dx2 *deriv.dxdchi + srcHess.dxdy*deriv.dydchi + srcHess.dxdz*deriv.dzdchi) * deriv.dxdrho +
        (srcHess.dxdy*deriv.dxdchi + srcHess.dy2*deriv.dydchi + srcHess.dydz *deriv.dzdchi) * deriv.dydrho +
        (srcHess.dxdz*deriv.dxdchi + srcHess.dydz*deriv.dydchi + srcHess.dz2 *deriv.dzdchi) * deriv.dzdrho +
        srcGrad.dx*deriv2.d2xdrhodchi + srcGrad.dy*deriv2.d2ydrhodchi + srcGrad.dz*deriv2.d2zdrhodchi;
    dest.drhodphi =
        (srcHess.dx2 *deriv.dxdphi+srcHess.dxdy*deriv.dydphi+srcHess.dxdz*deriv.dzdphi)*deriv.dxdrho+
        (srcHess.dxdy*deriv.dxdphi+srcHess.dy2*deriv.dydphi+srcHess.dydz*deriv.dzdphi)*deriv.dydrho+
        (srcHess.dxdz*deriv.dxdphi+srcHess.dydz*deriv.dydphi+srcHess.dz2*deriv.dzdphi)*deriv.dzdrho+
        srcGrad.dx*deriv2.d2xdrhodphi+srcGrad.dy*deriv2.d2ydrhodphi+srcGrad.dz*deriv2.d2zdrhodphi;
    dest.drhodchi=
        (srcHess.dx2 *deriv.dxdchi + srcHess.dxdy*deriv.dydchi + srcHess.dxdz*deriv.dzdchi) * deriv.dxdrho +
        (srcHess.dxdy*deriv.dxdchi + srcHess.dy2*deriv.dydchi + srcHess.dydz *deriv.dzdchi) * deriv.dydrho +
        (srcHess.dxdz*deriv.dxdchi + srcHess.dydz*deriv.dydchi + srcHess.dz2 *deriv.dzdchi) * deriv.dzdrho +
        srcGrad.dx*deriv2.d2xdrhodchi + srcGrad.dy*deriv2.d2ydrhodchi + srcGrad.dz*deriv2.d2zdrhodchi;
    dest.dchi2 =
        (srcHess.dx2 *deriv.dxdchi + srcHess.dxdy*deriv.dydchi + srcHess.dxdz*deriv.dzdchi) * deriv.dxdchi +
        (srcHess.dxdy*deriv.dxdchi + srcHess.dy2*deriv.dydchi + srcHess.dydz *deriv.dzdchi) * deriv.dydchi +
        (srcHess.dxdz*deriv.dxdchi + srcHess.dydz*deriv.dydchi + srcHess.dz2 *deriv.dzdchi) * deriv.dzdchi +
        srcGrad.dx*deriv2.d2xdchi2 + srcGrad.dy*deriv2.d2ydchi2+ srcGrad.dz*deriv2.d2zdchi2;
    dest.dchidphi =
        (srcHess.dx2 *deriv.dxdchi + srcHess.dxdy*deriv.dydchi + srcHess.dxdz*deriv.dzdchi) * deriv.dxdphi +
        (srcHess.dxdy*deriv.dxdchi+ srcHess.dy2*deriv.dydchi + srcHess.dydz *deriv.dzdchi) * deriv.dydphi +
        (srcHess.dxdz*deriv.dxdchi + srcHess.dydz*deriv.dydchi + srcHess.dz2 *deriv.dzdchi) * deriv.dzdphi +
        srcGrad.dx*deriv2.d2xdchidphi+ srcGrad.dy*deriv2.d2ydchidphi+ srcGrad.dz*deriv2.d2zdchidphi;
    dest.dphi2 =
        (srcHess.dx2 *deriv.dxdphi + srcHess.dxdy*deriv.dydphi + srcHess.dxdz*deriv.dzdphi) * deriv.dxdphi +
        (srcHess.dxdy*deriv.dxdphi + srcHess.dy2*deriv.dydphi + srcHess.dydz *deriv.dzdphi) * deriv.dydphi +
        (srcHess.dxdz*deriv.dxdphi + srcHess.dydz*deriv.dydphi + srcHess.dz2 *deriv.dzdphi) * deriv.dzdphi +
        srcGrad.dx*deriv2.d2xdphi2 + srcGrad.dy*deriv2.d2ydphi2+ srcGrad.dz*deriv2.d2zdphi2;
    return dest;
}

//------ conversion of derivatives of f(r) into gradients/hessians in different coord.sys. ------//
template<>
void evalAndConvertSph(const math::IFunction& F,
    const PosCar& pos, double* value, GradCar* deriv, HessCar* deriv2)
{
    assert(F.numDerivs()>=2);
    const double r=sqrt(pow_2(pos.x)+pow_2(pos.y)+pow_2(pos.z));
    if(deriv==NULL && deriv2==NULL) {
        F.evalDeriv(r, value, NULL, NULL);
        return;
    }
    double der, der2;
    F.evalDeriv(r, value, &der, deriv2!=NULL ? &der2 : NULL);
    double x_over_r=pos.x/r, y_over_r=pos.y/r, z_over_r=pos.z/r;
    if(r==0) {
        x_over_r=y_over_r=z_over_r=0;
    }
    if(deriv) {
        deriv->dx = x_over_r*der;
        deriv->dy = y_over_r*der;
        deriv->dz = z_over_r*der;
    }
    if(deriv2) {
        double der_over_r=der/r, dd=der2-der_over_r;
        if(r==0) {
            dd=0;
            if(der==0) der_over_r=der2;
        }
        deriv2->dx2 = pow_2(x_over_r)*dd + der_over_r;
        deriv2->dy2 = pow_2(y_over_r)*dd + der_over_r;
        deriv2->dz2 = pow_2(z_over_r)*dd + der_over_r;
        deriv2->dxdy= x_over_r*y_over_r*dd;
        deriv2->dydz= y_over_r*z_over_r*dd;
        deriv2->dxdz= x_over_r*z_over_r*dd;
    }
}

template<>
void evalAndConvertSph(const math::IFunction& F,
    const PosCyl& pos, double* value, GradCyl* deriv, HessCyl* deriv2)
{
    assert(F.numDerivs()>=2);
    const double r=sqrt(pow_2(pos.R)+pow_2(pos.z));
    if(deriv==NULL && deriv2==NULL) {
        F.evalDeriv(r, value, NULL, NULL);
        return;
    }
    double der, der2;
    F.evalDeriv(r, value, &der, deriv2!=NULL ? &der2 : NULL);
    double R_over_r=pos.R/r, z_over_r=pos.z/r;
    if(r==0) {
        R_over_r=z_over_r=0;
    }
    if(deriv) {
        deriv->dR = R_over_r*der;
        deriv->dz = z_over_r*der;
        deriv->dphi = 0;
    }
    if(deriv2) {
        double der_over_r=der/r, dd=der2-der_over_r;
        if(r==0) {
            dd=0;
            if(der==0) der_over_r=der2;
        }
        deriv2->dR2 = pow_2(R_over_r)*dd + der_over_r;
        deriv2->dz2 = pow_2(z_over_r)*dd + der_over_r;
        deriv2->dRdz= R_over_r*z_over_r*dd;
        deriv2->dRdphi=deriv2->dzdphi=deriv2->dphi2=0;
    }
}

template<>
void evalAndConvertSph(const math::IFunction& F,
    const PosSph& pos, double* value, GradSph* deriv, HessSph* deriv2)
{
    assert(F.numDerivs()>=2);
    double der, der2;
    F.evalDeriv(pos.r, value, deriv!=NULL ? &der : NULL, deriv2!=NULL ? &der2 : NULL);
    if(deriv) {
        deriv->dr = der;
        deriv->dtheta = deriv->dphi = 0;
    }
    if(deriv2) {
        deriv2->dr2 = der2;
        deriv2->dtheta2 = deriv2->dphi2 = deriv2->drdtheta = deriv2->drdphi = deriv2->dthetadphi = 0;
    }
}


//------ 3x3 matrix representing a [passive] rotation specified by Euler angles ------//

Orientation::Orientation(double alpha, double beta, double gamma)
{
    double sa, ca, sb, cb, sc, cc;
    math::sincos(alpha, sa, ca);
    math::sincos(beta,  sb, cb);
    math::sincos(gamma, sc, cc);
    mat[0] =  ca * cc - sa * cb * sc;
    mat[1] =  sa * cc + ca * cb * sc;
    mat[2] =  sb * sc;
    mat[3] = -ca * sc - sa * cb * cc;
    mat[4] = -sa * sc + ca * cb * cc;
    mat[5] =  sb * cc;
    mat[6] =  sa * sb;
    mat[7] = -ca * sb;
    mat[8] =  cb;
}

void Orientation::toEulerAngles(double& alpha, double& beta, double& gamma) const
{
    // beta is between 0 and pi; sin(beta) >= 0,
    // and while beta=acos(mat[8]) is mathematically correct, it has poor accuracy
    // when beta is close to 0 or pi; the alternative formula with atan2 works in all cases.
    beta  = math::atan2(sqrt(pow_2(mat[2]) + pow_2(mat[5])), mat[8]);
    if(mat[2] == 0 && mat[5] == 0 && mat[6] == 0 && mat[7] == 0) {
        // degenerate case: beta=0 or beta=pi;
        // in this case we can determine only the sum or the difference of the other two angles
        alpha = 0;
        gamma = math::atan2(mat[1], mat[0]) * (mat[8]>0 ? 1 : -1);
    } else {
        gamma = math::atan2(mat[2], mat[5]);
        alpha = math::atan2(mat[6],-mat[7]);
    }
}

}  // namespace coord
