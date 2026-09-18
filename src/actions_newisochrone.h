/* AngleAction coords of the isochrone sphere. Provides H, PE Omegar and OmegaL
 * in addition to p22J, pq2aa and aa2pq.
 * aa2pq optionally computes derivatives wrt (i) J, (ii) theta, (iii) M and b
*/
#pragma once
#include "math_core.h"
#include "math_specfunc.h"
#include "actions_base.h"
#include "coord.h"

namespace actions {

class Isochrone {
	public:
		double Js, b;
		Isochrone(double _Js=1, double _b=1) :  Js(_Js), b(_b) {
			if(Js<=0 || b<=0){
				printf("Illegal Isochrone parameter: %g %g\n",Js,b);
			}
		}
		void reset(double _Js,double _b){
			Js=_Js; b=_b;
			if(Js<=0 || b<=0){
				printf("Illegal Isochrone parameter: %g %g\n",Js,b);
			}
		}
		double H(coord::PosMomSph& rp) const;
		double H(Actions J) const{
			const double L=J.Jz+fabs(J.Jphi);
			return -.5*pow_2(pow_2(Js)/b)/pow_2(J.Jr+.5*(L+sqrt(L*L+4*Js*Js)));
		}
		double Omegar(const double Jr,const double L) const{
			return pow(Js,4)/(pow_2(b)*pow(Jr+.5*(L+sqrt(L+L+4*Js*Js)),3));
		}
		double OmegaL(const double Jr,const double L) const{
			return .5*(1+L/sqrt(L*L+4*Js*Js))*Omegar(Jr,L);
		}
		double PE(coord::PosMomSph& rp) const;
		coord::PosMomSph aa2pq(const ActionAngles& aa, Frequencies* freqs=NULL,
				       DerivAct<coord::Sph>* dJ=NULL, DerivAng<coord::Sph>* dA=NULL) const;
		coord::PosMomSph aa2pq(const Actions& J, const Angles& theta, Frequencies* freqs=NULL,
				       DerivAct<coord::Sph>* dJ=NULL, DerivAng<coord::Sph>* dA=NULL) const{
			return aa2pq(ActionAngles(J,theta), freqs, dJ, dA);
		}
		coord::PosMomSph aa2pq(const ActionAngles& aa,
				       coord::PosMomSph& drdM, coord::PosMomSph& drdb) const;
		coord::PosMomSph aa2pq(const Actions& J, const Angles& theta,
				       coord::PosMomSph& drdM, coord::PosMomSph& drdb) const{
			return  aa2pq(ActionAngles(J,theta), drdM, drdb);
		}
		ActionAngles pq2aa(const coord::PosMomSph& rtheta, Frequencies* freqs=NULL) const;
		Actions pq2J(const coord::PosMomSph& rtheta) const;
		Isochrone& operator *= (const double a) {
			Js *= a; b *= a;
			return *this;
		}
		Isochrone& operator += (const Isochrone& I){
			Js += I.Js; b+= I.b;
			return *this;
		}
		const Isochrone operator * (const double a) const {
			Isochrone I2(Js*a,b*a);
			return I2;
		}
		const Isochrone operator + (const Isochrone I) const{
			Isochrone I2(Js+I.Js,b+I.b);
			return I2;
		}
};

Isochrone interpIsochrone(const double, const Isochrone&, const Isochrone&);

Isochrone IsFromMass(const double mass,const double radius);

/** Compute any combination of actions, angles and frequencies
    in a spherical Isochrone potential specified by its total mass and scale radius.
    \param[in]  isochroneMass   is the total mass associated with the potential.
    \param[in]  isochroneRadius is the scale radius of the potential.
    \param[in]  point  is the position/velocity point.
    \param[out] act    if not NULL, will contain computed actions (Jr=NAN if E>=0).
    \param[out] ang    if not NULL, will contain corresponding angles (NAN if E>=0).
    \param[out] freq   if not NULL, will contain corresponding frequencies (NAN if E>=0).
*/
void evalIsochrone(
    const double isochroneMass, const double isochroneRadius,
    const coord::PosVelCyl& point,
    Actions* act=NULL,
    Angles* ang=NULL,
    Frequencies* freq=NULL);

/** Compute position/velocity from actions/angles in a spherical Isochrone potential.
    \param[in]  isochroneMass   is the total mass associated with the potential.
    \param[in]  isochroneRadius is the scale radius of the potential.
    \param[in]  actAng  is the action/angle point
    \param[out] freq    if not NULL, store the frequencies for these actions.
    \return     position and velocity point; NAN if Jr<0 or Jz<0.
*/
coord::PosVelCyl mapIsochrone(
    const double isochroneMass, const double isochroneRadius,
    const ActionAngles& actAng,
    Frequencies* freq=NULL);

/** Class for performing transformations between action/angle and coordinate/momentum for
    an isochrone potential (a trivial wrapper for the corresponding standalone functions) */
class ActionFinderIsochrone: public BaseActionFinder, public BaseActionMapper {
public:
    ActionFinderIsochrone(double _mass, double _radius): iso(IsFromMass(_mass,_radius)) {}

    virtual std::string name() const;

    virtual void eval(const coord::PosVelCyl& point,
        Actions* act=NULL, Angles* ang=NULL, Frequencies* freq=NULL) const{ 
		ActionAngles aa=iso.pq2aa(coord::toPosMomSph(point),freq);
		if(act)*act=aa;
		if(ang)*ang=Angles(math::wrapAngle(aa.thetar),aa.thetaz,aa.thetaphi);
	}

    virtual coord::PosVelCyl map(const ActionAngles& actAng, Frequencies* freq=NULL) const{
		coord::PosVelCyl Rv=coord::toPosVelCyl(iso.aa2pq(actAng,freq));
		return Rv;
	}

private:
    const Isochrone iso;  ///< parameters of the isochrone potential
};

// this class performs the conversion in both directions
typedef ActionFinderIsochrone ActionMapperIsochrone;

}//namespace