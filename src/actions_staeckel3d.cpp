#include "actions_staeckel3d.h"
#include "math_core.h"
#include "math_fit.h"
#include "utils.h"
#include <stdexcept>
#include <algorithm>
#include <cassert>
#include <cmath>

// debugging output
#include <fstream>

namespace actions{

namespace {  // internal routines

/** Accuracy of integrals for computing actions and angles
    is determined by the number of points in fixed-order Gauss-Legendre scheme with an order
    that depends on the orbit eccentricity, varying from the value below up to MAX_GL_ORDER */
static const unsigned int INTEGR_ORDER = 10;

/** order of Gauss-Legendre quadrature for actions, frequencies and angles:
    use a higher order for more eccentric orbits, as indicated by the ratio
    of pericenter to apocenter radii (R1/R2) */
inline unsigned int integrOrder(double R1overR2) {
    if(R1overR2==0)
        return math::MAX_GL_ORDER;
    int log2;  // base-2 logarithm of R1/R2
    frexp(R1overR2, &log2);
    return std::min<int>(math::MAX_GL_ORDER, INTEGR_ORDER - log2);
}

/** relative tolerance in determining the range of variables (nu,lambda) to integrate over */
static const double ACCURACY_RANGE = 1e-6;

/** minimum range of variation of nu, lambda that is considered to be non-zero */
static const double MINIMUM_RANGE = 1e-10;


// ------ Data structures for both Axisymmetric Staeckel and Fudge action-angle finders ------

/** integration intervals for actions and angles
    (shared between Staeckel and Fudge action finders). */
struct TriaxialIntLimits {
    double rho_min, rho_max, Xn_min, Xn_max,Xm_min,Xm_max;
    int integrOrder;   ///< order of GL quadrature that depends on the eccentricity
};

/** Derivatives of actions by integrals of motion */
struct TriaxialActionDerivatives {
    double dJrdE, dJrdI2, dJrdI3, dJzdE, dJzdI2,dJzdI3, dJphidE,dJphidI2,dJphidI3;
};

/** Derivatives of integrals of motion over actions (do not depend on angles).
    Note that dE/dJ are the frequencies of oscillations in three directions, 
    so that we reuse the `Frequencies` struct members (Omega***), 
    and add other derivatives with a somewhat different naming convention. */
struct TriaxialIntDerivatives : public Frequencies {
    double dI3dJr, dI3dJz, dI3dJphi, dI2dJr, dI2dJz, dI2dJphi;
};

/** Derivatives of generating function over integrals of motion (depend on angles) */
struct TriaxialGenFuncDerivatives {
    double dSdE, dSdI2, dSdI3;
};

/** aggregate class that contains the point in prolate spheroidal coordinates,
    integrals of motion, and reference to the potential,
    shared between Axisymmetric Staeckel  and Axisymmetric Fudge action finders;
    only the coordinates and the two classical integrals are in common between them.
    It also implements the IFunction interface, providing the "auxiliary function" 
    that is used in finding the integration limits and computing actions and angles;
    this function F(tau) is related to the canonical momentum p(tau)  as 
    \f$  p(tau)^2 = F(tau) / (2*(tau+alpha)^2*(tau+gamma))  \f$,
    and the actual implementation of this function is specific to each descendant class.
*/
class TriaxialFunctionBase: public math::IFunction {
public:
    coord::PosVelEls point; ///< position/derivative in prolate spheroidal coordinates
    const double E;             ///< total energy
    const double I2; 
    const double I3;            ///< third integral
    TriaxialFunctionBase(const coord::PosVelEls& _point, double _E, double _I2, double _I3) :
        point(_point), E(_E), I2(_I2), I3(_I3) {};
    virtual unsigned int numDerivs() const { return 2; }
};

// ------ SPECIALIZED functions for Staeckel action finder -------

/** parameters of potential, integrals of motion, and prolate spheroidal coordinates 
    SPECIALIZED for the Axisymmetric Staeckel action finder */
class TriaxialFunctionStaeckel: public  TriaxialFunctionBase{
public:
    const math::IFunction& fncG;  ///< single-variable function of a Staeckel potential
    TriaxialFunctionStaeckel(const coord::PosVelEls& _point, double _E, double _I2, double _I3,
        const math::IFunction& _fncG) :
        TriaxialFunctionBase(_point, _E, _I2, _I3), fncG(_fncG) {};

    /** auxiliary function that enters the definition of canonical momentum for 
        for the Staeckel potential: it is the numerator of eq.50 in de Zeeuw(1985);
        the argument tau is replaced by tau+gamma >= 0. */
    virtual void evalDeriv(const double tauplusgamma, 
        double* value=0, double* deriv=0, double* deriv2=0) const;
};
/** compute integrals of motion in the Staeckel potential of an oblate perfect ellipsoid, 
    together with the coordinates in its prolate spheroidal coordinate system 
*/
TriaxialFunctionStaeckel findIntegralsOfMotionPerfectEllipsoid(
    const potential::PerfectEllipsoidTriaxial& potential, 
    const coord::PosVelCar& point)
{
    double E = totalEnergy(potential, point);
    const coord::Els& els=potential.els();
    const coord::PosVelEls pe = coord::toPosVel<coord::Car, coord::Els>(point, els);
    double Glambda,Gmu,Gnu;
    double lambda=pow_2(pe.rho),mu=-pe.els.Deltay2*pow_2(cos(pe.phi));
    double snchi=1/sqrt(1+pow_2(pe.cotchi));
    double cschi=(!std::isinf(pe.cotchi))?snchi*pe.cotchi:1.;
    double Eyz=(els.Deltay2>0)?els.Deltay2/els.Deltaz2:0;
    double nu=-(1-(1-Eyz)*pow_2(cschi))*els.Deltaz2;
    //dmudphi=2*Dy2*csphi*snphi=2*Sqrt(-mu*(Dy2+mu))
    //cschi=sqrt((nu/Dz2+1)/(1-Eyz)),snchi=sqrt((-nu/Dz2-Eyz)/(1-Eyz))
    //dnudchi=-(1-Eyz)*2*cschi*snchi*Dz2=-2*sqrt((Dz2+nu)*(-Dy2-nu))
    potential.evalDeriv(lambda, &Glambda);
    potential.evalDeriv(mu, &Gmu);
    potential.evalDeriv(nu, &Gnu);
    double P2n=(lambda-mu)*(lambda-nu)/
        ((lambda+pe.els.Deltay2)*(lambda+els.Deltaz2));
    double Q2n=(mu-nu)*(lambda-mu)/(mu+pe.els.Deltaz2);
    //double R2n=(els.Deltay2>0)?(1-mu/nu)*(lambda-nu):lambda-nu;
    double X=.5*P2n*pow_2(pe.rhodot)-Glambda*lambda
            *(lambda+pe.els.Deltaz2)/((lambda-mu)*(lambda-nu));
    double Y=.5*Q2n*pow_2(pe.phidot)
    -Gmu*mu*(mu+els.Deltaz2)/((mu-lambda)*(mu-nu));
    //double Z1=.5*R2n*pow_2(pe.chidot)-Gnu*nu*(nu+els.Deltaz2)/((nu-lambda)*(nu-mu));
    double Z=E-X-Y;
    double J=(mu+nu)*X+(nu+lambda)*Y+(lambda+mu)*Z;
    double K=mu*nu*X+nu*lambda*Y+lambda*mu*Z;
    double I2=-K/els.Deltaz2;
    double I3=(pow_2(els.Deltaz2)*E+els.Deltaz2*J+K)/(els.Deltaz2);
    return TriaxialFunctionStaeckel(pe, E, I2, I3, potential);
}

/** auxiliary function that enters the definition of canonical momentum for 
    for the Staeckel potential: it is the numerator of eq.50 in de Zeeuw(1985);
    except that in our convention `tau` >= 0 is equivalent to `tau+gamma` from that paper. */
void TriaxialFunctionStaeckel::evalDeriv(const double tau, 
    double* val, double* der, double* der2) const
{
    //assert(tau>=0);
    double G, dG, d2G;
    fncG.evalDeriv(tau, &G, der || der2 ? &dG : NULL, der2 ? &d2G : NULL);
    const double taupDz=tau+point.els.Deltaz2;
    if(val) 
        *val = ( (E + G) * taupDz - I3 ) * tau - I2* taupDz;
    if(der){
        *der = (E + G) * (taupDz + tau) + dG * taupDz * tau - I3 - I2;
        if(std::isnan(*der))printf("no:%f %f %f %f\n",tau,G,dG);
    }
    if(der2)
        *der2 = 2 * (E + G) + 2 * dG * (taupDz + tau) + d2G * taupDz * tau;
}

// -------- COMMON routines for Staeckel and Fudge action finders --------
/** parameters for the function that computes actions and angles
    by integrating an auxiliary function "fnc" as follows:
    the canonical momentum is   p^2(tau) = fnc(tau) / [2 (tau+alpha)*(tau+beta)*(tau+gamma) ];
    the integrand is given by   p^n * (tau+alpha)^a*(tau+beta)^b*(tau+gamma)^c  if p^2>0, otherwise 0.
*/
class TriaxialIntegrand: public math::IFunctionNoDeriv {
public:
    const TriaxialFunctionBase& fnc;      ///< parameters of aux.fnc. (Staeckel or Fudge)
    enum { nplus1, nminus1 } n;         ///< power of p: +1 or -1
    enum { azero, aminus1 } a; ///< power of (tau+alpha): 0, -1
    enum { bzero, bminus1 } b;          ///< power of (tau+beta): 0 or -1
    enum { czero, cminus1 } c;          ///< power of (tau+gamma): 0 or -1
    double nu_max, dfdnu_at_nu_max;     ///< upper limit for nu and the fnc derivative at this point
    explicit TriaxialIntegrand(const TriaxialFunctionBase& d, double _nu_max=0) : fnc(d),nu_max(_nu_max)
    {
        fnc.evalDeriv(_nu_max, NULL, &dfdnu_at_nu_max);
    }

    /** integrand for the expressions for actions and their derivatives 
        (e.g.Sanders 2012, eqs. A1, A4-A12).  It uses the auxiliary function to compute momentum,
        and multiplies it by some powers of (tau-delta) and tau.
    */
    virtual double value(const double X) const {
        const coord::Els& els = fnc.point.els;
        double tau=0;
        double ft,p2;
        double result=0;
        //chi integral X+1=-abs(cschi), dtaudX=(Dz2-Dy2)*2*(X+1)
        if(X<-1){
            tau=-els.Deltaz2+(els.Deltaz2-els.Deltay2)*pow_2(X+1);
            ft = fnc(tau), p2 = 2*ft/(-tau*(1-pow_2(X+1)));
            if(p2<0) return 0;
            if(n==nminus1){
                double K=2*(tau*(tau+els.Deltay2)*(tau+els.Deltaz2))/ft;
                result=(K>0)?2*(X+1)*(els.Deltaz2-els.Deltay2)*sqrt(K):0;
            } 
            else result=(p2>0)?sqrt(p2):0;
        }
        //phi integral X=-abs(csphi)
        else if(X<0){
            tau=-els.Deltay2*X*X;
            ft=fnc(tau),p2=2*ft/(-(1-pow_2(X))*(els.Deltaz2+tau));
            if(p2<0) return 0;
            if(n==nminus1){
                double K=2*(tau*(tau+els.Deltay2)*(tau+els.Deltaz2))/ft;
                result=(K>0)?-2*X*els.Deltay2*sqrt(K):0;
            } 
            else result=(p2>0)?sqrt(p2):0;
        }
        //rho integal,X=rho
        else{
            tau=X*X;
            ft=fnc(tau),p2=2*ft/((tau+els.Deltay2)*(els.Deltaz2+tau));
            if(p2<0)return 0;
            if(n==nminus1){
                double K=2*(tau*(tau+els.Deltay2)*(tau+els.Deltaz2))/ft;
                result=(K>0)?2*X*sqrt(K):0;
            } 
            else result=(p2>0)?sqrt(p2):0;
        }
        double taupDy2=tau+els.Deltay2;
        double taupDz2=tau+els.Deltaz2;
        if(a==aminus1)
            result /= tau;
        if(b==bminus1) {
            result /=taupDy2;
            if(a==aminus1&&X<-1){
                double dtdX=(els.Deltaz2-els.Deltay2)*2*(X+1);
                result+=sqrt(2 * (nu_max+els.Deltaz2) / dfdnu_at_nu_max / (tau-nu_max)) /tau*(dtdX);
            }
            //*/
                // subtract the singular component that will be integrated analytically and added later
              //  result += sqrt(2 * nu_max / dfdnu_at_nu_max / (tau-nu_max)) / tauminusdelta;
        }
        if(c==cminus1)
            result /= taupDz2;
        if(!isFinite(result))
            result=0;  // ad hoc fix to avoid problems at the boundaries of integration interval
        return result;
    }

    /** limiting case of the integration interval collapsing to a single point tau,
        i.e. f(tau)~=0, f'(tau)~=0, and f''(tau)<0 (downward-curving parabola).
        In this case, if f(tau) is in the numerator, the integral is assumed to be zero,
        while if the integrand contains f(tau)^(-1/2), then the limiting value of the integral
        is computed from the second derivative of f(tau) at the (single) point.
    */
    double limitingIntegralValue(const double tau) const {
        if(n==nplus1)
            return 0;
        const coord::Els& els = fnc.point.els;
        assert(tau+els.Deltaz2>=0);
        const double taupDeltay2 = tau+els.Deltay2;
        const double taupDeltaz2= tau+els.Deltaz2;
        double fncder2;
        fnc.evalDeriv(tau, NULL, NULL, &fncder2);  // ignore f(tau) and f'(tau), only take f''(tau)
        double result = 2*M_PI * sqrt(-taupDeltaz2/fncder2) * fabs(tau);
        if(a==aminus1)
            result /= tau;
        else if(b==bminus1)
            result /= taupDeltay2;
        if(c==cminus1)
            result /= taupDeltaz2;
        return result;
    }
};
class Rootfnc: public math::IFunction{
public:
    const TriaxialFunctionBase& fnc;
    explicit Rootfnc(const TriaxialFunctionBase& d) : fnc(d) {};
    virtual unsigned int numDerivs() const { return fnc.numDerivs(); }
    virtual void evalDeriv(const double X, 
        double* val=0, double* der=0, double* der2=0) const
    {
        const coord::Els& els = fnc.point.els;
        double tau=0;
        double fder=0,f2der=0;
        if(X<-1){
            tau=-els.Deltaz2+(els.Deltaz2-els.Deltay2)*pow_2(X+1);
            fnc.evalDeriv(tau, val, der? &fder : NULL,der2?&f2der:NULL);
            double dtdX=(els.Deltaz2-els.Deltay2)*2*(X+1),d2tdX2=2*(els.Deltaz2-els.Deltay2);
            double fder0=fder;
            if(der)*der=fder*dtdX;
            if(der2)*der2=f2der*pow_2(dtdX)+fder0*d2tdX2;
        }
        else if(X<0){
            tau=-els.Deltay2*pow_2(X);
            fnc.evalDeriv(tau, val, der? &fder : NULL,der2?&f2der:NULL);
            double dtdX=(els.Deltay2)*2*X,d2tdX2=2*(-els.Deltay2);
            double fder0=fder;
            if(der)*der=fder*dtdX;
            if(der2)*der2=f2der=f2der*pow_2(dtdX)+fder0*d2tdX2;
        }
        else{
            tau=X*X;
            fnc.evalDeriv(tau, val, der? &fder : NULL,der2?&f2der:NULL);
            double dtdX=2*X,d2tdX2=2;
            double fder0=fder;
            if(der)*der=fder*dtdX;
            if(der2)*der2=f2der=f2der*pow_2(dtdX)+fder0*d2tdX2;
        }
    }
};
/** A simple function that facilitates locating the root of auxiliary function 
    on a semi-infinite interval for lambda: instead of F(tau) we consider 
    F1(tau)=F(tau)/tau^2, which tends to a finite negative limit as tau tends to infinity. */
class TriaxialScaledForRootfinder: public math::IFunction {
public:
    const Rootfnc& fnc;
    explicit TriaxialScaledForRootfinder(const Rootfnc& d) : fnc(d) {};
    virtual unsigned int numDerivs() const { return fnc.numDerivs(); }
    virtual void evalDeriv(const double X, 
        double* val=0, double* der=0, double* der2=0) const
    {
        assert(X>0);
        if(der2)*der2=NAN;
        coord::Els els=fnc.fnc.point.els;
        double tau=X*X;
        double fval,fder;
        if(/*!isFinite(tau)*/tau>1e100) {
            if(val)
                *val = fnc.fnc.E; // the asymptotic value
            if(der)
                *der = NAN;    // we don't know it
            return;
        }
        fnc.evalDeriv(X,val||der?&fval:NULL,der?&fder:NULL);
        if(val)
            *val = fval / pow_2(tau+els.Deltaz2);
        if(der)
            *der = (fder - 2*fval*2*X/(tau+els.Deltaz2)) / pow_2(tau+els.Deltaz2);
    }
};

/** Compute the intervals of tau for which p^2(tau)>=0, 
    where  0    = nu_min     <= tau <= nu_max     <= delta    is the interval for the "nu" branch,
    and  delta <= lambda_min <= tau <= lambda_max < infinity  is the interval for "lambda".
*/
TriaxialIntLimits findIntegrationLimitsTriaxi(const TriaxialFunctionBase& fnc)
{
    TriaxialIntLimits lim;
    
    // figure out the value of function at and around some important points
    double f_zero = fnc(0);
    Rootfnc fnc2(fnc);
    double f_lambda, df_lambda, d2f_lambda;
    fnc2.evalDeriv(fnc.point.rho, &f_lambda, &df_lambda, &d2f_lambda);
    lim.Xn_max = -1;
    lim.Xn_min=lim.Xm_max = lim.rho_min = lim.rho_max = NAN;  // means not yet determined
    //double Xn_upper=1,Xm_upper = 1;     // upper bound on the interval to locate the root for nu_max
    double rho_lower = 0; // lower bound on the interval for lambda_min

    /*if(fnc.Lz==0) {
        // special case: f(delta) = -0.5 Lz^2 = 0, may have either tube or box orbit in the meridional plane
        double deltaminus = delta*(1-1e-15), deltaplus = delta*(1+1e-15);
        if(fnc(deltaminus)<0) {  // box orbit: f<0 at some interval left of delta
            if(f_zero>0)         // there must be a range of nu where the function is positive
                nu_upper = deltaminus;
            else 
                lim.nu_max = lim.nu_min;
        } else
            lim.nu_max = delta;
        if(fnc(deltaplus)<0)     // tube orbit: f must be negative on some interval right of delta
            lambda_lower = deltaplus;
        else
            lim.lambda_min = delta;
    }
    //*/
    double cschi=1/sqrt(1+pow_2(fnc.point.cotchi));
    double Xn=(!std::isinf(fnc.point.cotchi))?(-fabs(cschi*fnc.point.cotchi)-1):-2;
    double Xm=-fabs(cos(fnc.point.phi));
    if(!isFinite(lim.Xn_min)) 
    {   // find range for J_nu (i.e. J_z) if it has not been determined at the previous stage
        double Xn0=(Xn<-1)?Xn:-1-1e-5;
        lim.Xn_min = math::findRoot(fnc2, -2, Xn0, ACCURACY_RANGE);
        if(!isFinite(lim.Xn_min))
            // means that the value f(nu) was just very slightly negative, or that f(0)<=0
            // i.e. this is a clear upper boundary of the range of allowed nu
            lim.Xn_min = Xn;
    }
    if(!isFinite(lim.Xm_max)) 
    {   // find range for J_mu (i.e. J_phi) if it has not been determined at the previous stage
        if(f_zero<0)lim.Xm_max=0;
        else lim.Xm_max = math::findRoot(fnc2, Xm, 0, ACCURACY_RANGE);
        if(!isFinite(lim.Xm_max))
            // means that the value f(nu) was just very slightly negative, or that f(0)<=0
            // i.e. this is a clear upper boundary of the range of allowed nu
            lim.Xm_max = Xm;
        //*/
        if(fnc2(-1)>0)lim.Xm_min=math::findRoot(fnc2,-1,Xm,ACCURACY_RANGE);
        else lim.Xm_min=-1;
    }
    // find the range for J_lambda (i.e. J_r).
    // We assume that the point lambda is inside or at the edge of the interval where f(lambda)>=0,
    // so that we will search for roots on the intervals (delta, lambda) and (lambda, infinity).
    // However, due to roundoff errors, it may actually happen that f(lambda) is negative,
    // or even positive but very small, or simply zero. In this case at least one or both intervals
    // must be modified so as to robustly bracket the point where f passes through zero.

    // this will be the point that is guaranteed to lie "well inside" the interval of positive f,
    // or, in other words, that the intervals [delta, lambda_pos] and [lambda_pos, infinity)
    // both firmly bracket the roots (respectively, lambda_min and lambda_max).
    double rho_pos = fnc.point.rho;
    // linear extrapolation to estimate the location of the nearest root for lambda
    double dxToRoot   = -f_lambda / df_lambda;
    if(f_lambda<=0 || fabs(dxToRoot) < fnc.point.rho * MINIMUM_RANGE)
    {   // we are at the endpoint of the interval where f is positive,
        // so at least one of the endpoints may be assigned immediately
        if(df_lambda>=0) {
            lim.rho_min = fnc.point.rho;
        } 
        if(df_lambda<=0) {
            lim.rho_max = fnc.point.rho;
            if(fnc.point.rho == 0)  // can't be lower than that! means that the range is zero
                lim.rho_min = 0;
        }

        // now it may also happen that we are at or very near the shell orbit,
        // i.e. both lambda_min and lambda_max are very close (or equal) to lambda.
        // This happens when d^2 f / d lambda^2 < 0, i.e. the function is a downward parabola,
        // and the distance between its roots is very small, or even it does not cross zero at all
        // (of course, this could only happen due to roundoff errors).
        // the second derivative must be negative, and the determinant either small or negative.
        bool nearShell = d2f_lambda < 0 &&
            pow_2(df_lambda) - 2*f_lambda*d2f_lambda -
            pow_2(fnc.point.rho * MINIMUM_RANGE * d2f_lambda) < 0;

        // However, the complication is that when lambda = delta, the second derivative is not defined.
        // Therefore, we first shift the point (lambda_pos) by a small amount in the direction of
        // increasing f, then recompute the second derivative, and then again test the condition.
        double safeOffset = fmin(0.5*fabs(df_lambda / d2f_lambda), fnc.point.rho * MINIMUM_RANGE);
        rho_pos += fmax(dxToRoot, safeOffset) * (df_lambda>=0?1:-1);
        double f_lampos, df_lampos, d2f_lampos;
        fnc.evalDeriv(rho_pos, &f_lampos, &df_lampos, &d2f_lampos);
        nearShell |= d2f_lampos < 0 &&
            pow_2(df_lampos) - 2*f_lampos*d2f_lampos -
            pow_2(fnc.point.rho * MINIMUM_RANGE * d2f_lampos) < 0;

        // unfortunately, f(lambda) is subject to such a severe cancellation error
        // that the above procedure does not always correctly identify a near-shell orbit.
        // therefore, we declare this to be the case even when the distance-between-roots condition
        // is not met, but the function is still negative, which would fail the root-finder anyway.
        if(nearShell || f_lampos < 0) {
            lim.rho_min = lim.rho_max = fnc.point.rho;
        }
    }
    if(!isFinite(lim.rho_min)) {  // not yet determined 
        if(f_zero<0)lim.rho_min = math::findRoot(fnc2, rho_lower, rho_pos, ACCURACY_RANGE);
        else lim.rho_min=0;
    }
    if(!isFinite(lim.rho_max)) {
        lim.rho_max = math::findRoot(TriaxialScaledForRootfinder(fnc2),
            math::ScalingSemiInf(rho_pos) /* find root on [lambda_pos..+inf) */, ACCURACY_RANGE);
    }

    // sanity check
    if(utils::verbosityLevel >= utils::VL_WARNING &&
        (!isFinite(lim.rho_min+lim.rho_max+lim.Xn_max+lim.Xn_min)))
        utils::msg(utils::VL_WARNING, "findIntegrationLimitsAxisym", "failed at lambda="+
            utils::toString(fnc.point.rho)+", nu="+utils::toString(fnc.point.cotchi)+", E="+
            utils::toString(fnc.E)+", I2="+utils::toString(fnc.I2)+", I3="+utils::toString(fnc.I3));

    // ignore extremely small intervals

    // choose the order of Gauss-Legendre integration depending on the approximate eccentricity
    double RperiOverRapo = (lim.rho_min) / (lim.rho_max);
    lim.integrOrder = integrOrder(RperiOverRapo);
    return lim;
}

/** Compute the derivatives of actions (Jr, Jz, Jphi) over integrals of motion (E, Lz, I3),
    using the expressions A4-A9 in Sanders(2012).
    TODO: improve efficiency by saving the values of potential taken along the same
    integration path when computing the actions, and re-using them in this routine.
*/
TriaxialActionDerivatives computeActionDerivatives(
    const TriaxialFunctionBase& fnc, const TriaxialIntLimits& lim)
{
    TriaxialActionDerivatives der;
    coord::Els els=fnc.point.els;
    double nu_max=-els.Deltaz2+(els.Deltaz2-els.Deltay2)*(pow_2(1+lim.Xn_min));
    TriaxialIntegrand integrand(fnc, nu_max);
    math::ScaledIntegrand<math::ScalingCub>
        transf_l(math::ScalingCub(lim.rho_min, lim.rho_max), integrand),
        transf_n(math::ScalingCub(lim.Xn_min,     lim.Xn_max),     integrand),
        transf_m(math::ScalingCub(lim.Xm_min,     lim.Xm_max),     integrand);
    integrand.n = TriaxialIntegrand::nminus1;  // momentum goes into the denominator
    integrand.b =TriaxialIntegrand::bminus1;
    // derivatives w.r.t. E
    integrand.a = TriaxialIntegrand::azero;
    integrand.c = TriaxialIntegrand::czero;
    der.dJrdE = (lim.rho_min==lim.rho_max ? 
        integrand.limitingIntegralValue(lim.rho_min) :
        math::integrateGL(transf_l, 0, 1, lim.integrOrder) ) / (4*M_PI);
    der.dJzdE = -math::integrateGL(transf_n, 0, 1, lim.integrOrder) / (2*M_PI);
    der.dJphidE = math::integrateGL(transf_m, 0, 1, lim.integrOrder) / (2*M_PI);
    // derivatives w.r.t. I3
    integrand.a = TriaxialIntegrand::azero;
    integrand.b =TriaxialIntegrand::bminus1;
    integrand.c = TriaxialIntegrand::cminus1;
    der.dJrdI3 = - (lim.rho_min==lim.rho_max ? 
        integrand.limitingIntegralValue(lim.rho_min) :
        math::integrateGL(transf_l, 0, 1, lim.integrOrder) ) / (4*M_PI);
    der.dJzdI3 = math::integrateGL(transf_n, 0, 1, lim.integrOrder) / (2*M_PI);
    der.dJphidI3 = -math::integrateGL(transf_m, 0, 1, lim.integrOrder) / (2*M_PI);
    // derivatives w.r.t. Lz
    integrand.a = TriaxialIntegrand::aminus1;
    integrand.b = TriaxialIntegrand::bminus1;
    integrand.c = TriaxialIntegrand::czero;
    der.dJrdI2 = -(lim.rho_min==lim.rho_max ? 
        integrand.limitingIntegralValue(lim.rho_min) :
        math::integrateGL(transf_l, 0, 1, lim.integrOrder) ) / (4*M_PI);
    // the following integral is split into the analytically computed singular part
    // and the remaining regular part integrated numerically
    //double delta_minus_nu_max = fnc.point.coordsys.Delta2 - lim.nu_max;
    double singpart = 2 * sqrt(-2 * (nu_max+els.Deltaz2) / integrand.dfdnu_at_nu_max /(-nu_max) ) *
        (math::atan(sqrt( (nu_max+els.Deltaz2) / (-nu_max))));
    double Lz=math::sign(fnc.point.phidot)*sqrt(2*fnc.I2);
    der.dJphidI2 =(els.Deltay2/els.Deltaz2<1e-8)?1/Lz:-math::sign(fnc.point.phidot)*(math::integrateGL(transf_m, 0, 1, lim.integrOrder)) / (2*M_PI);
    der.dJzdI2 = (math::integrateGL(transf_n, 0, 1, lim.integrOrder)-singpart)/(2*M_PI);
    return der;
}

/** Compute the derivatives of integrals of motion (E, Lz, I3) over actions (Jr, Jz, Jphi), inverting
    the matrix of action derivatives by integrals.  These quantities are independent of angles,
    and in particular, the derivatives of energy w.r.t. the three actions are the frequencies. */
TriaxialIntDerivatives computeIntDerivatives(
    const TriaxialFunctionBase& fnc, const TriaxialIntLimits& lim)
{
    TriaxialActionDerivatives dJ = computeActionDerivatives(fnc, lim);
    TriaxialIntDerivatives der;
    // invert the matrix of derivatives
    math::Matrix<double> dJdI(3,3);
    dJdI(0,0)=dJ.dJrdE;
    dJdI(0,1)=dJ.dJrdI2;
    dJdI(0,2)=dJ.dJrdI3;
    dJdI(1,0)=dJ.dJzdE;
    dJdI(1,1)=dJ.dJzdI2;
    dJdI(1,2)=dJ.dJzdI3;
    dJdI(2,0)=dJ.dJphidE;
    dJdI(2,1)=dJ.dJphidI2;
    dJdI(2,2)=dJ.dJphidI3;
    math::LUDecomp LU(dJdI);
    std::vector<double> dIdJr=LU.solve({1,0,0});
    der.Omegar=dIdJr[0];
    der.dI2dJr=dIdJr[1];
    der.dI3dJr=dIdJr[2];
    std::vector<double> dIdJz=LU.solve({0,1,0});
    der.Omegaz=dIdJz[0];
    der.dI2dJz=dIdJz[1];
    der.dI3dJz=dIdJz[2];
    std::vector<double> dIdJphi=LU.solve({0,0,1});
    der.Omegaphi=dIdJphi[0];
    der.dI2dJphi=dIdJphi[1];
    der.dI3dJphi=dIdJphi[2];
    /*double det  = dJ.dJrdE * dJ.dJzdI3 - dJ.dJrdI3 * dJ.dJzdE;
    if(lim.nu_min==lim.nu_max || det==0) {
        // special case z==0: motion in z is irrelevant, but we could not compute dJzdI3 which is not zero
        der.Omegar   = 1 / dJ.dJrdE;
        der.Omegaphi =-dJ.dJrdLz / dJ.dJrdE;
        der.dI3dJr   = der.dI3dJz = der.dI3dJphi = der.Omegaz = 0;
    } else {  // everything as normal
        der.Omegar   = dJ.dJzdI3 / det;  // dE/dJr
        der.Omegaz   =-dJ.dJrdI3 / det;  // dE/dJz
        der.Omegaphi = (dJ.dJrdI3 * dJ.dJzdLz - dJ.dJrdLz * dJ.dJzdI3) / det;  // dE/dJphi
        der.dI3dJr   =-dJ.dJzdE / det;
        der.dI3dJz   = dJ.dJrdE / det;
        der.dI3dJphi =-(dJ.dJrdE * dJ.dJzdLz - dJ.dJrdLz * dJ.dJzdE) / det;
    }
    der.dLzdJr  = 0;
    der.dLzdJz  = 0;
    der.dLzdJphi= 1;
    //*/
    return der;
}

/** Compute the derivatives of generating function S over integrals of motion (E, Lz, I3),
    using the expressions A10-A12 in Sanders(2012).  These quantities do depend on angles.
    TODO: improve efficiency by computing all three integrals for each of two directions at once,
    saving on repetitive potential evaluations along the same paths.
*/
TriaxialGenFuncDerivatives computeGenFuncDerivatives(
    const TriaxialFunctionBase& fnc, const TriaxialIntLimits& lim)
{
    const double signldot = fnc.point.rhodot >= 0 ? +1 : -1;
    const double signndot = -fnc.point.cotchi*fnc.point.chidot >= 0 ? +1 : -1;
    const double signmdot = //-(fnc.point.phi<.5*M_PI||(fnc.point.phi>M_PI&&fnc.point.phi<1.5*M_PI))?1:-1;
    fnc.point.phidot >= 0 ? +1 : -1;
    //printf("sgnl:%f %f %f",signldot,signndot,signmdot);
    TriaxialGenFuncDerivatives der;
    coord::Els els=fnc.point.els;
    double nu_max=-els.Deltaz2+(els.Deltaz2-els.Deltay2)*(pow_2(1+lim.Xn_min));
    TriaxialIntegrand integrand(fnc,nu_max);
    math::ScaledIntegrand<math::ScalingCub>
        transf_l(math::ScalingCub(lim.rho_min, lim.rho_max), integrand),
        transf_m(math::ScalingCub(lim.Xm_min,     lim.Xm_max),     integrand),
        transf_n(math::ScalingCub(lim.Xn_min,     lim.Xn_max),     integrand);
    const double yl = math::scale(transf_l.scaling, fnc.point.rho);
    double Xn=(!std::isinf(fnc.point.cotchi))?(-1-abs(fnc.point.cotchi)/sqrt(1+pow_2(fnc.point.cotchi))):-2;
    const double yn = lim.Xn_min==lim.Xn_max ? 0 : math::scale(transf_n.scaling, Xn);
    const double ym = lim.Xm_min==lim.Xm_max ? 0 : math::scale(transf_m.scaling, -fabs(cos(fnc.point.phi)));
    integrand.n = TriaxialIntegrand::nminus1;  // momentum goes into the denominator
    integrand.b =TriaxialIntegrand::bminus1;
    // derivatives w.r.t. E
    integrand.a = TriaxialIntegrand::azero;
    integrand.c = TriaxialIntegrand::czero;
    der.dSdE =
        signldot * math::integrateGL(transf_l, 0, yl, lim.integrOrder) / 4
      + signmdot * math::integrateGL(transf_m, 0, ym, lim.integrOrder) / 4
      + signndot * -math::integrateGL(transf_n, yn, 1, lim.integrOrder) / 4;
    // derivatives w.r.t. I3
    integrand.a = TriaxialIntegrand::azero;
    integrand.b = TriaxialIntegrand::bminus1;
    integrand.c = TriaxialIntegrand::cminus1;
    der.dSdI3 = 
        signldot * -math::integrateGL(transf_l, 0, yl, lim.integrOrder) / 4
      + signmdot * -math::integrateGL(transf_m, 0, ym, lim.integrOrder) / 4
      + signndot * math::integrateGL(transf_n, yn, 1, lim.integrOrder) / 4;
    // derivatives w.r.t. Lz
    integrand.a = TriaxialIntegrand::aminus1;
    integrand.b = TriaxialIntegrand::bminus1;
    integrand.c = TriaxialIntegrand::czero;
    // the integral over nu is split into the analytically computed singular part
    // and the remaining regular part integrated numerically
    double nu=-els.Deltaz2+(els.Deltaz2-els.Deltay2)*(pow_2(1+Xn));
    double singpart = 2 * sqrt(-2 * (nu_max+els.Deltaz2) / integrand.dfdnu_at_nu_max /(-nu_max) ) *
        (math::atan(sqrt( (nu_max+els.Deltaz2) / (-nu_max))) -
         math::atan(sqrt((nu_max - nu) /(-nu_max))));
        //*/
    int sgnphi=(fnc.point.phi<.5*M_PI||(fnc.point.phi>M_PI&&fnc.point.phi<1.5*M_PI))?1:-1;
    double Lz=signmdot*sqrt(2*fnc.I2);
    double dSdI2nophi=signldot * - math::integrateGL(transf_l, 0, yl, lim.integrOrder) / 4
      + signndot * (math::integrateGL(transf_n, yn, 1, lim.integrOrder)-singpart)/4;
    double phicomp=(els.Deltay2==0)?0:-signmdot*sgnphi*math::integrateGL(transf_m, 0, ym, lim.integrOrder) / 4;
    //needto do properly
    if(els.Deltay2/els.Deltaz2<1e-8){
        der.dSdI2 =fnc.point.phi/Lz+dSdI2nophi;
    }
    else der.dSdI2=phicomp+dSdI2nophi;
    //printf("phi:%f\n",fnc.point.phi);
    return der;
}

/** Compute actions by integrating the momentum over the range of tau on which it is positive,
    separately for the "nu" and "lambda" branches (equation A1 in Sanders 2012). */
Actions computeActions(const TriaxialFunctionBase& fnc, const TriaxialIntLimits& lim)
{
    Actions acts;
    TriaxialIntegrand integrand(fnc);
    math::ScaledIntegrand<math::ScalingCub>
        transf_l(math::ScalingCub(lim.rho_min, lim.rho_max), integrand),
        transf_m(math::ScalingCub(lim.Xm_min,     lim.Xm_max),     integrand),
        transf_n(math::ScalingCub(lim.Xn_min,     lim.Xn_max),     integrand);
    integrand.n = TriaxialIntegrand::nplus1;  // momentum goes into the numerator
    integrand.a = TriaxialIntegrand::azero;
    integrand.b = TriaxialIntegrand::bzero;
    integrand.c = TriaxialIntegrand::czero;
    acts.Jr = math::integrateGL(transf_l, 0, 1, lim.integrOrder) / M_PI;
    // factor of 2 in Jz because we only integrate over half of the orbit (z>=0)
    acts.Jphi = math::integrateGL(transf_m, 0, 1, lim.integrOrder) / M_PI * 2;
    acts.Jz = math::integrateGL(transf_n, 0, 1, lim.integrOrder) / M_PI * 2;
    return acts;
}

/** Compute angles from the derivatives of integrals of motion and the generating function
    (equation A3 in Sanders 2012). */
Angles computeAngles(const TriaxialIntDerivatives& derI, const TriaxialGenFuncDerivatives& derS,
    bool addPiToThetaZ)
{
    Angles angs;
    angs.thetar   = derS.dSdE*derI.Omegar   + derS.dSdI3*derI.dI3dJr   + derS.dSdI2*derI.dI2dJr;
    angs.thetaz   = derS.dSdE*derI.Omegaz   + derS.dSdI3*derI.dI3dJz   + derS.dSdI2*derI.dI2dJz;
    angs.thetaphi = derS.dSdE*derI.Omegaphi + derS.dSdI3*derI.dI3dJphi + derS.dSdI2*derI.dI2dJphi;
    angs.thetar   = math::wrapAngle(angs.thetar);
    angs.thetaz   = math::wrapAngle(angs.thetaz + M_PI*addPiToThetaZ);
    angs.thetaphi = math::wrapAngle(angs.thetaphi);
    return angs;
}

}  // internal namespace

// -------- THE DRIVER ROUTINES --------

void evalTriaxialStaeckel(
    const potential::PerfectEllipsoidTriaxial& potential, const coord::PosVelCar& point,
    Actions* act, Angles* ang, Frequencies* freq)
{
    TriaxialFunctionStaeckel fnc=findIntegralsOfMotionPerfectEllipsoid(potential, point);
    const TriaxialIntLimits lim = findIntegrationLimitsTriaxi(fnc);
    /*int N=lim.integrOrder;
    std::vector<double> x(3*lim.integrOrder);
    coord::Els els=fnc.point.els;
    for(int i=0;i<N;i++){
        double Xn=math::unscale(math::ScalingCub(lim.Xn_min,lim.Xn_max),math::GLPOINTS[N][N-1-i]);
        x[i]=-els.Deltaz2+(els.Deltaz2-els.Deltay2)*pow_2(Xn+1);
        double Xm=math::unscale(math::ScalingCub(lim.Xm_min,lim.Xm_max),math::GLPOINTS[N][i]);
        x[i+N]=-els.Deltay2*pow_2(Xm);
        x[i+2*N]=math::unscale(math::ScalingCub(lim.rho_min,lim.rho_max),math::GLPOINTS[N][i]);
    }
        //*/
    //potential::PerfectEllipsoidTriaxialInterp potinterp(potential,x);
    //TriaxialFunctionStaeckel fnc(fnc.point,fnc.E,fnc.I2,fnc.I3,potential);
    //*/
    if(act)
        *act = computeActions(fnc, lim);
    if(ang || freq) {
        TriaxialIntDerivatives derI = computeIntDerivatives(fnc, lim);
        if(freq)  // store frequencies which are the first row of the derivatives matrix
            *freq = derI;
        if(ang) {
            TriaxialGenFuncDerivatives derS = computeGenFuncDerivatives(fnc, lim);
            bool addPiToThetaZ = fnc.point.chidot>0 && fnc.I2!=0;
            *ang = computeAngles(derI, derS, addPiToThetaZ);
        }
    }
    //evalTriaxial(fnc, act, ang, freq);
}
}