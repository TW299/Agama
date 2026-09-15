#pragma once
#include "actions_base.h"
#include "math_core.h"
#include "math_base.h"
#include "math_fit.h"
#include "orbit.h"
#include "potential_base.h"
#include "potential_utils.h"
namespace potential{
/// return scaledE as a function of E and invPhi0 = 1/Phi(0)
inline double scaleE(const double E, const double invPhi0, /*output*/ double* dEdscaledE=NULL)
{
	double expX = invPhi0 - 1/E;
	if(dEdscaledE)
		*dEdscaledE = E * E * expX;
	return log(expX);
}
		//return E and optionally dE/d(scaledE) as a function of scaledE and invPhi0
inline double unscaleE(const double scaledE, const double invPhi0, /*output*/ double* dEdscaledE=NULL)
{
	double expX = exp(scaledE);
	double E = 1 / (invPhi0 - expX);
	if(dEdscaledE)
		*dEdscaledE = E * E * expX;
	return E;
}
    /* The main job of PolarInterpolator is to hold the curve Jz(Jf) along
 * which the box/loop transition lies. In adition it holds the values
 * of Delta(E) that cause the I3 centrifugal barrier to vanish on this
 * curve and te associaed I3(E), where I3 is computed from the velocity
 * of the transition orbit at (Rsh,0) 
*/
class  PolarInterpolator{
	private:
		//math::LinearInterpolator interpI3, interpFD, interpUmin;
		//std::vector<double> coeffsJz;
		math::LinearInterpolator interpJzJ,interpJzE;
		math::ScalingSemiInf Sc;
		bool isSph;
	public:
		PolarInterpolator():isSph(true){}
		PolarInterpolator(const BasePotential&);//, const PtrShellInterpolator);
/*		PolarInterpolator(const std::vector<double>& gridEscaled, const std::vector<double>& gridI3,
				  const std::vector<double>& gridFD, const std::vector<double>& gridJfScaled,
				  const std::vector<double>& gridJz) :
		    interpI3(gridEscaled,gridI3), interpFD(gridEscaled,gridFD),
		    coeffsJz(math::fitPoly(15,gridJfScaled,gridJz)) {}*/
		/*void getFDI3critUmin(const double E,const double invPhi0,
			       double& Delta, double& I3, double& Umin) const{
			double scaledE = math::clip(scaleE(E, invPhi0),
				interpFD.xmin(), interpFD.xmax());
			Delta = interpFD.value(scaledE);
			I3 = interpI3.value(scaledE);
			Umin = interpUmin.value(scaledE);
		}			
		double getI3crit(const double E,const double invPhi0) const{
			const double scaledE = math::clip(scaleE(E, invPhi0),
				interpI3.xmin(), interpI3.xmax());
			return interpI3.value(scaledE);
		}
		double getFDcrit(const double E,const double invPhi0) const{
			const double scaledE = math::clip(scaleE(E, invPhi0),
				interpFD.xmin(), interpFD.xmax());
			return interpFD.value(scaledE);
		}
		double getUmin(const double E,const double invPhi0) const{
			const double scaledE = math::clip(scaleE(E, invPhi0),
				interpUmin.xmin(), interpUmin.xmax());
			return interpUmin.value(scaledE);
		}
        //*/
		double getJzcrit(const double E,const double invPhi0) const{
			if(isSph||interpJzE.numValues()==0)return 0;
			else{
				const double scaledE = math::clip(scaleE(E, invPhi0),
				interpJzE.xmin(), interpJzE.xmax());
				return interpJzE(scaledE);
			}
		}
		double getJzcrit(const double Jf) const{
			if(isSph||interpJzJ.numValues()==0)return 0;
			else{
				return interpJzJ(scale(Sc,Jf));
			}
			//return math::evalPoly(coeffsJz, scale(Sc,Jf));
		}
};

const extern PolarInterpolator PolarInterpolatorSph;
}