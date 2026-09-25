/** \file potential_perfect_ellipsoid_triaxial.h
    \brief potential in the triaxial Perfect Ellipsoid model 
    \author Tom Wright
    \date 2025-2026
*/
//See Elliptical galaxies with separable potentials by T. de Zueew (1985) for more details about this potential 
#pragma once
#include "potential_base.h"
#include "math_specfunc.h"
#include "math_spline.h"
namespace potential{
class PerfectEllipsoidTriaxial:public potential::BasePotential, 
public coord::IScalarFunction<coord::Els>, public math::IFunction{
public:
    PerfectEllipsoidTriaxial(double _mass,double _a,double _b,double _c):mass(_mass),a(_a),b(_b),c(_c),a2(_a*_a),b2(_b*_b),
    c2(_c*_c),elsc(_a*_a-_c*_c,_a*_a-_b*_b)
    {
        if(b>a||c>b||b<0)
            throw std::invalid_argument("Error in PerfectEllipsoidTriaxial: minor axis must be positive and c<=b<=a.");
        l=acos(c/a), snl=sin(l);
        snm=sqrt(1-pow_2(b/a))/snl;
        if(snm>1)snm=1;
        Flm=math::ellintF(l,snm);
        Elm=math::ellintE(l,snm);
        double snm2=snm*snm,csm2=1-snm*snm;
        K=mass/(M_PI*a*snl);
        //formula 
        Ga=2*K*Elm;//value of G(a2)
        if(a/b-1>1e-6)dGa=-2*K/(3*a2*pow_2(snl*snm))*((1-pow_2(snm))*Flm+(2*pow_2(snm)-1)*Elm-b/a*pow_2(snm)*snl*c/a);//G'(a2)
        else dGa=2*mass/(M_PI*c)*(-1./3.);
        dGc=-2*K/(3*a2*pow_2(snl)*(1-snm*snm))*((1-2*pow_2(snm))*Elm-(1-pow_2(snm))*Flm
                +b/c*snl*(pow_2(b/c)-1+pow_2(snm)));//G'(c2)
        double A4=((2+snm2)*Flm-2*(1+snm2)*Elm+snm2*b*c/a2*snl)/(3*pow_2(snm2));
        double A6=(pow_3(snl)*c*b/a2+4*(1+snm2)*A4-3*(Flm-Elm)/snm2)/(5*snm2);
        if(a/b-1>1e-6)d2Ga=4*K/pow_2(a2-c2)*(A4-snm2*A6);//G''(a2)
        else d2Ga=2*mass/(M_PI*c)*(2./5.);
        d2Gb=4*K/(3*pow_2((a2-c2)*snm2*csm2))*(2*(2*snm2-1)*Elm+csm2*(2-3*snm2)*Flm
                -snm2*snl*c/b*(a2/b2*csm2+4*snm2-2));//G''(b2)
        double J4=(2*(1+csm2)*Elm-csm2*Flm+snl*b/c*(csm2*a2/c2-2-2*csm2))/(3*pow_2(csm2));
        double J6=(pow_3(snl)*b/c*pow_2(a2/c2)-3*(b/c*snl-Elm)/(csm2)-4*(1+csm2)*J4)/(5*csm2);
        d2Gc=4*K/(pow_2(a2-c2))*(csm2*J6+J4);//G''(c2)
    }

    virtual coord::SymmetryType symmetry() const { 
        return coord::ST_TRIAXIAL; }
    virtual void evalDeriv(double tau, double* G=NULL, double* Gderiv=NULL, double* Gderiv2=NULL) const;
    const coord::Els& els() const { return elsc; }
    virtual std::string name() const { return myName(); }
    static std::string myName() { return "PerfectEllipsoidTriaxial"; }
    virtual double totalMass() const { return mass; }
    virtual void evalScalar(const coord::PosEls& pos,
        double* value=NULL, coord::GradEls* deriv=NULL, coord::HessEls* deriv2=NULL,
        double time=0) const;
private:
    const double mass,a,b,c,a2,b2,c2;
    const coord::Els elsc;
    double l,snl,snm,Flm,Elm;
    double K,Ga,dGa,dGc,d2Ga,d2Gb,d2Gc;
    virtual void evalCar(const coord::PosCar &pos,
        double* potential, coord::GradCar* deriv, coord::HessCar* deriv2, double time) const {
        // no direct conversion exists, use two-step
        coord::evalAndConvert<coord::Els, coord::Car>
            (*this, pos, potential, deriv, deriv2, time, elsc);
    }
    virtual void evalCyl(const coord::PosCyl &pos,
        double* potential, coord::GradCyl* deriv, coord::HessCyl* deriv2, double time) const {
        coord::PosCar xyz(pos.R,0,pos.z);
        coord::evalAndConvertTwoStep<coord::Els, coord::Car, coord::Cyl>
            (*this, pos, elsc, potential, deriv, deriv2, time);
    }
    virtual void evalSph(const coord::PosSph &pos,
        double* potential, coord::GradSph* deriv, coord::HessSph* deriv2, double time) const {
        coord::evalAndConvertTwoStep<coord::Els, coord::Car, coord::Sph>
            (*this, pos, elsc, potential, deriv, deriv2, time);  // use two-step conversion
    }

    /** the function that does the actual computation in prolate spheroidal coordinates 
        (implements the coord::IScalarFunction<ProlSph> interface) */
    virtual unsigned int numDerivs() const { return 2; }
};
class PerfectEllipsoidTriaxialInterp:public math::IFunction{
public:
    PerfectEllipsoidTriaxialInterp(const PerfectEllipsoidTriaxial &pot,const std::vector<double> xvals){
        const int N=xvals.size();
        std::vector<double> yvals(N);//, dyvals(N),d2yvals(N);
        for(int i=0;i<N;i++)pot.evalDeriv(xvals[i],&yvals[i]);
        //q=math::QuinticSpline(xvals,yvals,dyvals,d2yvals);
        interp=math::LinearInterpolator(xvals,yvals);
    }
    virtual void evalDeriv(double tau, double* G=NULL, double* Gderiv=NULL, double* Gderiv2=NULL) const{
        interp.evalDeriv(tau,G,Gderiv,Gderiv2);
    }
    virtual unsigned int numDerivs() const { return 2; }
private:
    math::LinearInterpolator interp;
    //math::QuinticSpline q;
};
}