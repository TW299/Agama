#include "potential_perfect_ellipsoid_triaxial.h"
#include "math_core.h"
#include <cassert>
namespace potential{
class func:public math::IFunctionNoDeriv{
public:
   const double a2,b2,c2,t;
    func(double _a2,double _b2,double _c2,double _t):a2(_a2),b2(_b2),c2(_c2),t(_t){}
    virtual double value(double s) const{
        double s2=s*s;
        double sqr1=sqrt(c2+(b2-c2)*s2),sqr2=sqrt(c2+(a2-c2)*s2);
        return sqr1/(sqr2*(c2+s2*(t-c2)));
    }
};
class derivfunc:public math::IFunctionNoDeriv{
public:
   const double a2,b2,c2,t;
    derivfunc(double _a2,double _b2,double _c2,double _t):a2(_a2),b2(_b2),c2(_c2),t(_t){}
    virtual double value(double s) const{
        double s2=s*s;
        double sqr1=sqrt(c2+(b2-c2)*s2),sqr2=sqrt(c2+(a2-c2)*s2);
        return -sqr1*s2/(sqr2*pow_2(c2+s2*(t-c2)));
    }
};

class d2func:public math::IFunctionNoDeriv{
public:
   const double a2,b2,c2,t;
    d2func(double _a2,double _b2,double _c2,double _t):a2(_a2),b2(_b2),c2(_c2),t(_t){}
    virtual double value(double s) const{
        double s2=s*s;
        double sqr1=sqrt(c2+(b2-c2)*s2),sqr2=sqrt(c2+(a2-c2)*s2);
        return 2*sqr1*s2*s2/(sqr2*pow_3(c2+s2*(t-c2)));
    }
};
void PerfectEllipsoidTriaxial::evalScalar(const coord::PosEls& pos,
    double* val, coord::GradEls* deriv, coord::HessEls* deriv2, double) const
{
    assert(pos.els.Deltay2 == elsc.Deltay2);
    assert(pos.els.Deltaz2 ==elsc.Deltaz2);
    if(pos.rho<0)
        throw std::invalid_argument("Error in OblatePerfectEllipsoid: "
            "incorrect values of spheroidal coordinates");
    if(!(pos.rho/a < 1e10)) {
        if(val)
            *val = 0;
        if(deriv)
            deriv->drho = deriv->dchi = deriv->dphi = 0;
        if(deriv2)
            deriv2->drho2 = deriv2->dchi2 = deriv2->drhodchi=deriv2->drhodphi=deriv2->dchidphi=deriv2->dphi2=0;
        return;
    }
    double snchi=1/sqrt(1+pow_2(pos.cotchi));
    double cschi=!(std::isinf(pos.cotchi))?pos.cotchi*snchi:1.;
    double lambda=pow_2(pos.rho);
    double mu=-elsc.Deltay2*pow_2(cos(pos.phi));
    double Eyz=(elsc.Deltay2>0)?elsc.Deltay2/elsc.Deltaz2:0;
    double nu=-(1-(1-Eyz)*pow_2(cschi))*pos.els.Deltaz2;
    double dldr=2*pos.rho,d2ldr2=2.;
    double dmudphi=2*elsc.Deltay2*cos(pos.phi)*sin(pos.phi);
    double d2mudphi2=2*elsc.Deltay2*(2*pow_2(cos(pos.phi))-1);
    double dnudchi=-2*(1-Eyz)*cschi*snchi*pos.els.Deltaz2;
    double d2nudchi2=-2*(1-Eyz)*(2*pow_2(cschi)-1)*pos.els.Deltaz2;
    double Glambda, dGdlambda, d2Gdlambda2, Gmu,dGdmu,d2Gdmu2, Gnu, dGdnu, d2Gdnu2;
    // values and derivatives of G(lambda) and G(|nu|)
    evalDeriv(lambda, &Glambda, &dGdlambda, &d2Gdlambda2);
    evalDeriv(mu,      &Gmu,     &dGdmu,     &d2Gdmu2);
    evalDeriv(nu,      &Gnu,     &dGdnu,     &d2Gdnu2);
    double Fl=lambda*(lambda+pos.els.Deltaz2)*Glambda;
    double Fm=mu*(mu+pos.els.Deltaz2)*Gmu;
    double Fn=nu*(nu+pos.els.Deltaz2)*Gnu;
    double dFldl=(2*lambda+pos.els.Deltaz2)*Glambda+lambda*(lambda+pos.els.Deltaz2)*dGdlambda;
    double dFmdm=(2*mu+pos.els.Deltaz2)*Gmu+mu*(mu+pos.els.Deltaz2)*dGdmu;
    double dFndn=(2*nu+pos.els.Deltaz2)*Gnu+nu*(nu+pos.els.Deltaz2)*dGdnu;
    double dPhidl=0,dPhidm=0,dPhidn=0;
    if(val!=NULL) 
        *val = (mu-nu>1e-8*a2)?-(Fl/((lambda-mu)*(lambda-nu))+Fm/((mu-lambda)*(mu-nu))
                +Fn/((nu-mu)*(nu-lambda))):-(Fl/((lambda-mu)*(lambda-nu))+dFndn/((nu-lambda))-Fn/pow_2(lambda-nu));
    if(deriv!=NULL||deriv2!=NULL) {
        if(mu-nu>1e-8*a2){
            dPhidl=-(dFldl/((lambda-mu)*(lambda-nu))-Fl/pow_2((lambda-mu)*(lambda-nu))
                    *(2*lambda-nu-mu)+Fm/(pow_2(mu-lambda)*(mu-nu))+Fn/((nu-mu)*pow_2(nu-lambda)));
            dPhidm=-(Fl/(pow_2(lambda-mu)*(lambda-nu))+dFmdm/((mu-lambda)*(mu-nu))
                -Fm/pow_2((mu-lambda)*(mu-nu))*(2*mu-nu-lambda)+Fn/(pow_2(nu-mu)*(nu-lambda)));
            dPhidn=-(Fl/((lambda-mu)*pow_2(lambda-nu))+Fm/((mu-lambda)*pow_2(mu-nu))
                +dFndn/((nu-lambda)*(nu-mu))-Fn/pow_2((nu-lambda)*(nu-mu))*(2*nu-lambda-mu));
        }
        else{
            dPhidl=-(dFldl/((lambda-mu)*(lambda-nu))-Fl/pow_2((lambda-mu)*(lambda-nu))
                    *(2*lambda-nu-mu)+dFndn/(pow_2(mu-lambda)));
            double d2Fndn2=2*Gnu+2*(2*nu+pos.els.Deltaz2)*dGdnu+nu*
                (nu+pos.els.Deltaz2)*d2Gdnu2;
            dPhidm=-(Fl/(pow_2(lambda-mu)*(lambda-nu))+.5*d2Fndn2/((nu-lambda))+Fn/(pow_2(lambda-mu)*(nu-lambda))-dFndn/((lambda-mu)*(lambda-nu)));
        }
        if(deriv!=NULL){
            deriv->drho=dldr*dPhidl;
            deriv->dchi=dnudchi*dPhidn;
            deriv->dphi=dPhidm*dmudphi;
        }
    }
    if(deriv2!=NULL) {
        double d2Fldl2=(2*Glambda+2*(2*lambda+pos.els.Deltaz2)*dGdlambda+lambda*
            (lambda+pos.els.Deltaz2)*d2Gdlambda2);
        double d2Fmdm2=2*Gmu+2*(2*mu+pos.els.Deltaz2)*dGdmu+mu*
            (mu+pos.els.Deltaz2)*d2Gdmu2;
        double d2Fndn2=2*Gnu+2*(2*nu+pos.els.Deltaz2)*dGdnu+nu*
            (nu+pos.els.Deltaz2)*d2Gdnu2;
        deriv2->drho2 = -pow_2(dldr)*(d2Fldl2/((lambda-mu)*(lambda-nu))-2*dFldl/pow_2((lambda-mu)
                    *(lambda-nu))*(2*lambda-nu-mu)-Fl/pow_2((lambda-mu)*(lambda-nu))*
                    (2-2*pow_2(2*lambda-nu-mu)/((lambda-mu)*(lambda-nu)))
                    +2*Fm/(pow_3(mu-lambda)*(mu-nu))+2*Fn/((nu-mu)*pow_3(nu-lambda)))+dPhidl*d2ldr2;
            deriv2->drhodphi=-dmudphi*dldr*(dFldl/(pow_2(lambda-mu)*(lambda-nu))-Fl/pow_2((lambda-mu)*(lambda-nu))
                    *(2*(2*lambda-nu-mu)/(lambda-mu)-1)+dFmdm/(pow_2(mu-lambda)*(mu-nu))
                    -Fm/(pow_3(mu-lambda)*pow_2(mu-nu))*(3*mu-2*nu-lambda)
                    +Fn/pow_2((nu-mu)*(nu-lambda)));
            deriv2->drhodchi=-dnudchi*dldr*(dFldl/(pow_2(lambda-nu)*(lambda-mu))-Fl/pow_2((lambda-mu)*(lambda-nu))
                    *(2*(2*lambda-nu-mu)/(lambda-nu)-1)+dFndn/(pow_2(nu-lambda)*(nu-mu))
                    -Fn/(pow_3(nu-lambda)*pow_2(nu-mu))*(3*nu-2*mu-lambda)
                    +Fm/pow_2((mu-nu)*(mu-lambda)));
            deriv2->dphi2 = -pow_2(dmudphi)*(d2Fmdm2/((mu-nu)*(mu-lambda))-2*dFmdm/pow_2((mu-nu)
                    *(mu-lambda))*(2*mu-lambda-nu)-Fm/pow_2((mu-nu)*(mu-lambda))*
                    (2-2*pow_2(2*mu-lambda-nu)/((mu-nu)*(mu-lambda)))
                    +2*Fn/((nu-lambda)*pow_3(nu-mu))+2*Fl/((lambda-nu)*pow_3(lambda-mu)))+dPhidm*d2mudphi2;
            deriv2->dchidphi = -dmudphi*dnudchi*(Fl/pow_2((lambda-mu)*(lambda-nu))+dFmdm/((mu-lambda)*pow_2(mu-nu))
                -Fm/pow_2((mu-lambda)*(mu-nu))*(2*(2*mu-nu-lambda)/(mu-nu)-1)
                +dFndn/(pow_2(nu-mu)*(nu-lambda))
                +Fn/(pow_3(nu-mu)*pow_2(nu-lambda))*(2*lambda-3*nu+mu));
            deriv2->dchi2 = -pow_2(dnudchi)*(d2Fndn2/((nu-mu)*(nu-lambda))-2*dFndn/pow_2((nu-mu)
                    *(nu-lambda))*(2*nu-lambda-mu)-Fn/pow_2((nu-mu)*(nu-lambda))*
                    (2-2*pow_2(2*nu-lambda-mu)/((nu-mu)*(nu-lambda)))
                    +2*Fm/((mu-lambda)*pow_3(mu-nu))+2*Fl/((lambda-mu)*pow_3(lambda-nu)))+dPhidn*d2nudchi2;
    }
}
void PerfectEllipsoidTriaxial::evalDeriv(double tau, double* G, double* Gderiv, double* Gderiv2) const{
    double P=math::ellintP(l,snm,-tau/(c2-a2));
    if(G){
        if(fabs(tau/a2)<1e-6){
            *G=Ga+dGa*(tau)+.5*d2Ga*pow_2(tau);
        }
        else *G=2*K/tau*((tau+a2-b2)*P+(b2-a2)*Flm);

    }
    if(Gderiv){
        if(fabs((tau)/a2)<1e-6)*Gderiv=dGa+d2Ga*(tau);
        else if(fabs((tau+a2-c2)/a2)<1e-6)*Gderiv=dGc+d2Gc*(tau+a2-c2);
        else{
            double dG=K/(pow_2(tau)*(tau+a2)*(tau+a2-c2))*((a2-c2)*(tau+a2)*tau*Elm+
                (tau+a2)*(b2*(c2-2*(tau+a2))+c2*(tau+a2)+a2*(b2-2*c2+tau+a2))*Flm+
                +(-a2*b2+a2*c2-b2*c2+2*b2*(tau+a2)-pow_2(tau+a2))*(tau+a2)*P+.5*a*b*sin(2*l)*pow_2(tau));
                *Gderiv=dG;
        }
        if(!std::isnan(mass)&&std::isnan(*Gderiv)){
            printf("tau:%f %f Gderiv:%f\n",tau+a2-c2,tau,*Gderiv);
            printf("V:%f %f %f %f\n",snm,sin(l),Elm,K);
            printf("a:%f,b:%f,c:%f\n",a,b,c);
            printf("dG;%f %f",dGa,d2Ga);
            exit(0);
        }
    }
    if(Gderiv2){
        if(fabs((tau)/a2)<1e-6)*Gderiv2=d2Ga;
        else if(fabs((tau+a2-b2)/a2)<1e-6)*Gderiv2=d2Gb;
        else if(fabs((tau+a2-c2)/a2)<1e-6)*Gderiv2=d2Gc;
        else{
            double d2G=-K/(4*pow_3(tau)*pow_2((tau+a2)-c2)*pow_2((tau+a2))*((tau+a2)-b2))*
                (2*(a2-c2)*(-tau)*pow_2(tau+a2)*(3*b2*(2*(tau+a2)-c2)-(5*(tau+a2)-2*c2)*(tau+a2)
                +a2*(2*(tau+a2)+c2-3*b2))*Elm+2*(b2-(tau+a2))*pow_2(tau+a2)*(3*pow_2(a2)*(b2-c2)+c2*(2*c2-5*(tau+a2))*(tau+a2)
                +a2*(-5*pow_2(c2)+2*b2*(c2-4*(tau+a2))+14*c2*(tau+a2)-3*pow_2((tau+a2)))+b2*(3*pow_2(c2)-8*c2*(tau+a2)+8*(tau+a2)*(tau+a2)))*Flm
                +2*pow_2((tau+a2))*P*(-3*pow_2(c2*b2)+pow_2(a2)*(-3*pow_2(b2)+2*b2*c2+pow_2(c2))
                -2*a2*b2*c2*(b2-c2)+4*(pow_2(a2)*(b2-c2)+b2*c2*(2*b2+c2)+a2*(2*b2*b2-b2*c2-c2*c2))*(tau+a2)
                -2*(4*b2*b2+5*b2*c2+5*a2*(b2-c2))*pow_2((tau+a2))+12*b2*pow_3((tau+a2))-3*pow_2((tau+a2)*(tau+a2)))
                +sin(2*l)*a*b*(5*pow_3((tau+a2))*pow_2((tau+a2))+pow_2((tau+a2))*(tau+a2)*(3*b2*(c2-2*(tau+a2))-2*c2*(tau+a2))
                +pow_3(a2)*((c2-4*(tau+a2))*(tau+a2)+b2*(5*(tau+a2)-2*c2))+
                pow_2(a2)*(tau+a2)*(b2*(7*c2-16*(tau+a2))+(tau+a2)*(13*(tau+a2)-4*c2))+a2*(tau+a2)*(tau+a2)*
                (5*c2*(tau+a2)-14*(tau+a2)*(tau+a2)+b2*(17*(tau+a2)-8*c2))));
            *Gderiv2=d2G;
        }
    }
}
}