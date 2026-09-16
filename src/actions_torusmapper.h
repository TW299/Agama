#pragma once
#include "actions_newtorus.h"
#include "coord.h"
namespace actions{
	class ActionMapperTorus:public BaseActionMapper{
	public:
		ActionMapperTorus(const potential::BasePotential& pot,const double tol = 1e-9);
		virtual std::string name() const{return "Torus action mapper";}
		virtual coord::PosVelCyl map(const ActionAngles& actAng, Frequencies* freq=NULL) const;
	private:
		class Impl;
		Impl *impl;
	};
}