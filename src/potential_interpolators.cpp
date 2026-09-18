#include "potential_interpolators.h"
#include "actions_focal_distance_finder.h"
#include "utils.h"
#ifdef _OPENMP
#include <omp.h>
#endif
namespace potential{
/// number of sampling points for a shell orbit (equally spaced in time)
static const unsigned int NUM_STEPS_TRAJ = 64;
/// accuracy of root-finding for the radius of thin (shell) orbit
static const double ACCURACY_RSHELL = 1e-6;
/// accuracy of orbit integration for shell orbit
static const double ACCURACY_INTEGR = 1e-8;
/// accuracy parameter determining the spacing of the interpolation grid along the energy axis
static const double ACCURACY_INTERP2 = 1e-5;
/// upper limit on the number of timesteps in ODE solver (should be enough to track half of the orbit)
static const unsigned int MAX_NUM_STEPS_ODE = 2000;

namespace{
/** find the best-fit value of focal distance for a shell orbit.
    \param[in] traj  contains the trajectory of this orbit in R-z plane,
    \return  the parameter `delta` of a prolate spheroidal coordinate system which minimizes
    the variation of `lambda` coordinate for this orbit
    If the best-fit value is negative, it is replaced with zero.
*/
double fitFocalDistanceShellOrbit(const std::vector<coord::PosVelCyl>& traj)
{
	if(traj.size()==0)
		throw std::invalid_argument("Error in finding focal distance for a shell orbit: empty array");
	math::Matrix<double> coefs(traj.size(), 2);
	std::vector<double> rhs(traj.size());
	std::vector<double> result;  // regression parameters:  lambda/(lambda-delta), lambda
	for(unsigned int i=0; i<traj.size(); i++) {
		coefs(i, 0) = -pow_2(traj[i].R);
		coefs(i, 1) = 1.;
		rhs[i] = pow_2(traj[i].z);
	}
	// We want (R/b)^2+(z/a)^2=1 so result[1]+result[0](-R^2)=z^2 
	// so result[0]=(a/b)^2 and result[1]=a^2
	// -> FD^2=a^2-b^2=result[1]*(1-result[0])
	//If 1>result[0]>1 an ellipse has been fitted, while
	//a hyperbola is signalled by result[0]<0 
	math::linearMultiFit(coefs, rhs, NULL, result);
	//printf("R(%g %g %g)",result[0]-1,result[1],result[1] * (1. - result[0]));
	return sqrt(fmax(0,result[1] * (1 - result[0])));
//	return result[0]>1? sqrt( fmax( result[1] * (1 - 1/result[0]), 0) ) :
//			sqrt(fabs(result[1]));
}

    /* The next three classs are used to compute Jr & Jz at the box/loop
 * transition.
 * pzf.evalDeriv returns pz for motion along the z axis
 */
class pzf : public math::IFunction {
	private:
		const potential::BasePotential& pot;
		const double E,x0;
	public:
		pzf(const potential::BasePotential& _pot,double _E, double _x0=0) : pot(_pot), E(_E),x0(_x0) {};
		virtual void evalDeriv(const double z,
				       double* value = NULL, double* deriv = NULL,
				       double* deriv2 = NULL) const {
			coord::PosCar X(x0, 0, z);
			double Phi;
			coord::GradCar dx;
			pot.eval(X, &Phi, &dx);
				//std::cout << "Value:" << x << " " << pot.value(X) << "\n";
			double p2 = 2 * (E - Phi);
			if (p2 <= 0) {
				if (value)*value = 0;
			}
			else {
				double p = sqrt(p2);
				if (value)*value = p;
				if (deriv)*deriv = -dx.dz / p;
			}
		};
		virtual unsigned int numDerivs() const { return 1; }
};
/* Now some code used to compute (Jr,Jz) at the box/loop transition
 * when Jphi=0. We compute Jf=2*Jr+Jz by integrating along the
 * (unstable) short-axis orbit and then get Jr from the area of the SoS
 * when the instability of the orbit is triggered
 */
//pxf.evalDeriv returns px for motion along x axis
class pxf : public math::IFunction {
	private:
		const potential::BasePotential& pot;
		const double E;
	public:
		pxf(const potential::BasePotential& _pot, double _E) :pot(_pot), E(_E){};
		virtual void evalDeriv(const double x,
				       double* value = NULL, double* deriv = NULL,
				       double* deriv2 = NULL) const {
			coord::PosCar X(x, 0, 0);
			double Phi;
			coord::GradCar dx;
			pot.eval(X, &Phi, &dx);
				//std::cout << "Value:" << x << " " << pot.value(X) << "\n";
			double p2 = 2 * (E - Phi);
			if (p2 <= 0) {
				if (value) *value = 0;
			}
			else {
				double p = sqrt(p2);
				if (value) *value = p;
				if (deriv) *deriv = -dx.dx / p;
			}
		};
		virtual unsigned int numDerivs() const { return 1; }
};
//finds the z1 s.t. Iz=int pzdz between 0 and z1, pz is the z momentum at that energy.
class Ifind : public math::IFunction {
	private:
		const double I;
		pzf pzfunc;
	public:
		Ifind(double _Iz, const potential::BasePotential& _pot, double _E,double _x0=0)
				:I(_Iz), pzfunc(_pot,_E,_x0){
		};
		virtual void evalDeriv(const double z,
				       double* value = 0, double* deriv = 0, double* deriv2 = 0)  const {
			if (value)*value = math::integrateGL(math::ScaledIntegrand<math::ScalingCub>
				(math::ScalingCub(0, z), pzfunc), 0, 1, math::MAX_GL_ORDER) - I;
			if (deriv) {
				double pz;
				pzfunc.evalDeriv(z, &pz);
				*deriv = pz;
			}
		}
		virtual unsigned int numDerivs() const { return 1; }
};
		//Bubble sort (x,y) in ascending order of x.
void sort(std::vector<double>& x, std::vector<double>& y) {
	int n = x.size();
	for (int i = 0; i < n - 1; i++) {
		bool swap = false;
		for (int j = 0; j < n - i - 1; j++) {
			if (x[j] > x[j + 1]) {
				std::swap(x[j], x[j + 1]);
				std::swap(y[j], y[j + 1]);
				swap = true;
			}
		}
		if (!swap) break;
	}
}
// computes area inside the curve defined by  the points (x,y) - used to compute Jr at box/loop
// transition. On return biggest values of x and y in x2max & y2max
// and in ysh the y-value of the curve when x=Rsh
double Area(std::vector<double> x, std::vector<double> y,
	    double& x2max, double& ymv2) {
	std::vector<double> xn, yn;
	for (int i = 0; i < x.size(); i++) {//Assume 4-fold symmetry
		xn.push_back(fabs(x[i])); yn.push_back(fabs(y[i]));
	}
	sort(xn, yn);
	int n = yn.size();
	double I = yn[0] * (xn[1] - xn[0])
		   + yn[n - 1] * (xn[n - 1] - xn[n - 2]);
	for (int i = 1; i < n - 1; i++) 
		I +=  yn[i] * (xn[i + 1] - xn[i - 1]);
	I *= .5;
	x2max = xn[xn.size()-1];
	ymv2 = yn[yn.size()-1];
	return I;
}
// computes Actions J of the orbit with enegy E at the box/loop
// transition. On return vR has Rdot(Rsh)
actions::Actions BoxLoopTrAct(const potential::BasePotential& pot, double E) {
	
	double zmax = potential::z_max(pot, E);
			//Angle to z axis at origin
	double t = 0;
	double Rmax0 = potential::R_max(pot,E);
	double x0 = 1e-3*Rmax0;
	pzf pzt(pot,E,x0);
	double Phi1 = pot.value(coord::PosCyl(x0, 0, 0));
	double v = sqrt(2 * (E - Phi1));
	double vz = v * cos(t), vx = v * sin(t);
	coord::PosVelCyl xv0(x0, 0, 0, vx, vz, 0);
	std::vector<double> R, pR;
	orbit::makeSoS(xv0, pot, R, pR, 1000);
	double Rmax, pRm;
	double Jr = Area(R, pR, Rmax, pRm) / M_PI;
	double v1 = sqrt(2*(E-pot.value(coord::PosCyl(Rmax,0,0))));
	if(pRm/v1>0.8||R.size()<3||abs(Rmax)<10*x0){//From Rmax to Rmax0 approx motion as on x axis 
		pxf pxfunc(pot,E);
		Jr+=math::integrateGK(math::ScaledIntegrand<math::ScalingCub>
				      (math::ScalingCub(Rmax, Rmax0), pxfunc), 0, 1, 1e-8)/M_PI;
	}
	//Get Jfast=Jx+Jz=2*Jr+Jz from motion along the z axis
	double Jfast = 2 * math::integrateGK(math::ScaledIntegrand<math::ScalingCub>
					     (math::ScalingCub(0, zmax), pzt), 0, 1, 1e-8) / M_PI;
	double Jz = Jfast - 2 * Jr;
	if(Jz<0){
		Jz=0;
		Jr=.5*Jfast;
	}
	return actions::Actions(Jr, Jz, 0);
}

/// function to use in ODE integrator
class OrbitIntegratorMeridionalPlane: public math::IOdeSystem {
	public:
		OrbitIntegratorMeridionalPlane(const potential::BasePotential& p, double Lz) :
		    poten(p), Lz2(Lz*Lz) {};

    /** apply the equations of motion in R,z plane without tracking the azimuthal motion.
        Integration variables are: R, z, vR, vz, dR, dz, dvR, dvz
        R here can have a negative sign (this happens for Lz=0, when the orbit flips to x<0
        and crosses the z=0 plane at negative 'R', but we compute the potential derivatives at |R|,
        and multiply by sign(R) when necessary.
    */
		virtual void eval(const double /*t*/, const double x[], double dxdt[],double*) const
		{
			coord::GradCyl grad;
			coord::HessCyl hess;
			double signR = x[0]>=0 ? 1 : -1;
			coord::PosCyl pos(fabs(x[0]), x[1], 0);
			poten.eval(pos, NULL, &grad, &hess);
			double Lz2ovR4 = Lz2>0 ? Lz2/pow_2(pow_2(pos.R)) : 0;
			dxdt[0] = x[2];
			dxdt[1] = x[3];
			dxdt[2] = -(grad.dR - Lz2ovR4 * pos.R) * signR;
			dxdt[3] = - grad.dz;
			dxdt[4] = x[6];
			dxdt[5] = x[7];
			dxdt[6] = -(hess.dR2 + 3*Lz2ovR4) * x[4] - hess.dRdz * signR * x[5];
			dxdt[7] = - hess.dRdz * signR * x[4] - hess.dz2 * x[5];
			dxdt[8] = pow_2(x[3]) + pow_2(x[2]);//dJz = vz*dz + vR*dR 
		}

		virtual unsigned int size() const { return 9; }  // two coordinates and two velocities
	private:
		const potential::BasePotential& poten;
		const double Lz2;
};

/// function to use in locating the exact time of the x-y plane crossing
class FindCrossingPointZequal0: public math::IFunction {
	public:
		FindCrossingPointZequal0(const math::BaseOdeStepper& _solver) :
		    solver(_solver) {};
    /** used in root-finder to locate the root z(t)=0 */
		virtual void evalDeriv(const double timeOffset, double* val, double* der, double*) const
		{
			if(val)
				*val = solver.getSol(timeOffset, 1);  // z
			if(der)
				*der = solver.getSol(timeOffset, 3);  // vz
		}
		virtual unsigned int numDerivs() const { return 1; }
	private:
		const math::BaseOdeStepper& solver;
};

/// function to use in locating the exact time vz goes -ve
class FindCrossingPointVZequal0: public math::IFunction {
	public:
		FindCrossingPointVZequal0(const math::BaseOdeStepper& _solver) :
		    solver(_solver) {};
    /** used in root-finder to locate the root z(t)=0 */
		virtual void evalDeriv(const double timeOffset, double* val, double* der, double*) const
		{
			if(val)	*val = solver.getSol(timeOffset, 3);  // Vz
		}
		virtual unsigned int numDerivs() const { return 0; }
	private:
		const math::BaseOdeStepper& solver;
};

/** launch an orbit perpendicularly to x-y plane from radius R0 with vz>0,
    and record the radius at which it crosses this plane downward (vz<0).
    \param[in]  poten  is the potential;
    \param[in]  E  is the orbit energy;
    \param[in]  Lz  is the z-component of angular momentum;
    \param[in]  R0  is the radius of the starting point;
    \param[out] timeCross stores the time required to complete the half-oscillation in z;
    \param[out] traj stores the trajectory recorded at equal intervals of time;
    \param[out] Rcross stores the radius of the crossing point;
    \param[out] dRcrossdR0  stores the derivative dRcross/dR0, computed from the variational equation.
*/
void findCrossingPointR(
			const potential::BasePotential& poten, double E, double Lz, double R0,
			double& timeCross,
			std::vector<std::pair<coord::PosVelCyl, double> >& traj, double& Rcross,
			double& dRcrossdR0, double& Jz)
{
	double Phi;
	coord::GradCyl grad;
	poten.eval(coord::PosCyl(R0, 0, 0), &Phi, &grad);
    // initial vertical velocity
	double vz0 = sqrt(fmax( 2 * (E-Phi) - (Lz>0 ? pow_2(Lz/R0) : 0), 0));
    // initial R-component of the deviation vector
	double dR0 = 1.;
    // initial vz-component (assigned from the requirement that E=const)
	double dvz0= vz0>0 ? ((Lz>0 ? pow_2(Lz) / pow_3(R0) : 0) - grad.dR) / vz0 * dR0 : 0;
	double vars[9] = {R0, 0, 0, vz0, dR0, 0, 0, dvz0, 0};
    OrbitIntegratorMeridionalPlane odeSystem(poten, Lz);
	math::OdeStepperDOP853 stepper(odeSystem, ACCURACY_INTEGR);
	stepper.init(vars);
	bool finished = false;
	unsigned int numStepsODE = 0;
	// time at the beginning of the current timestep
    double timeBegin = 0;
	// time offset of the next point on the stored trajectory from the beginning of the completed timestep
    double timeOffsetTraj = 0;
	//Store the rising quarter of the orbit
	const double timeStepTraj = timeCross*0.5/(NUM_STEPS_TRAJ-1);
	double tcurr=0;
	traj.clear();
	double vR,vz,dR,dz;
	while(!finished) {
        double timeStep = stepper.doStep(INFINITY);
		if(timeStep <= 0 || numStepsODE >= MAX_NUM_STEPS_ODE) { // signal of error
			 FILTERMSG(utils::VL_WARNING, "estimateFocalDistanceShellOrbit",
                "Failed to compute orbit for E="+utils::toString(E,16)+
                ", Lz="+utils::toString(Lz,16)+", R="+utils::toString(R0,16));
			timeCross  = 0;
			Rcross     = R0;   // this would terminate the root-finder, but we have no better option..
			dRcrossdR0 = NAN;
			if(R0>1e-8) printf("After %d steps failed to compute shell orbit for E=%g, Lz=%g, R=%g\n",
			       numStepsODE,E,Lz,R0);
			return;
		} else {
			numStepsODE++;
			if(timeStepTraj!=INFINITY)
			{   // store first part of trajectory
				while(timeOffsetTraj <= timeStep && traj.size() < NUM_STEPS_TRAJ) {
		    // store R, z, vR, vz at equal intervals of time
					double R = stepper.getSol(timeOffsetTraj, 0);
					double z = stepper.getSol(timeOffsetTraj, 1);
					double vR = stepper.getSol(timeOffsetTraj, 2);
					double vz = stepper.getSol(timeOffsetTraj, 3);
					traj.push_back(std::make_pair(coord::PosVelCyl(fabs(R), z, 0, vR, vz, 0),tcurr));
                	timeOffsetTraj += timeStepTraj;
					tcurr+=timeStepTraj;
            }
            timeOffsetTraj -= timeStep;  // prepare for the next timestep
			}
			if(stepper.getSol(timeStep, 1) <= 0) {  // z<=0 - we're done
				finished = true;
				double timeOffsetCross = math::findRoot(FindCrossingPointZequal0(stepper),
					0, timeStep, ACCURACY_RSHELL);
				timeCross = timeBegin + timeOffsetCross;
				Rcross    = stepper.getSol(timeOffsetCross, 0);
                vR = stepper.getSol(timeOffsetCross, 2);
                vz = stepper.getSol(timeOffsetCross, 3);
            // components of the deviation vector (dR,dz) at the crossing point
            	dR = stepper.getSol(timeOffsetCross, 4);
            	dz = stepper.getSol(timeOffsetCross, 5);
				Jz=stepper.getSol(timeOffsetCross, 8)/M_PI;
				finished=true;
			}
        timeBegin += timeStep;
		}
	}
	dRcrossdR0= dR - dz * vR / vz;
	if(Rcross < 0) {  // this happens for Lz=0, when the orbit crosses the x axis at negative x
		Rcross     = -Rcross;
		dRcrossdR0 = -dRcrossdR0;
	}
}
}

void FindClosedOrbitRZplane::evalDeriv(const double R0, double* val, double* der, double*) const {
	// first two calls in root-finder are for the boundary points, we already know the answer
	if(R0==Rmin || R0==Rmax) {
		if(val) *val = R0==Rmin ? Rmax-Rmin : Rmin-Rmax;
		if(der) *der = NAN;
		return;
	}
	double Rcross, dRcrossdR=NAN;
	findCrossingPointR(poten, E, Jphi, R0, timeCross, traj, Rcross, dRcrossdR, Jz);
	if(val)
		*val = Rcross-R0;
	if(der)
		*der = dRcrossdR-1;
}

const PolarInterpolator PolarInterpolatorSph;

void isMonotone(std::vector<double> x){
	for(int i=1; i<x.size(); i++){
		if(x[i]<=x[i-1])
			printf("isn't Monontone %g %g\n",x[i-1],x[i]);
	}
}
ShellInterpolator::ShellInterpolator(const BasePotential& pot,
					 const std::string logfname){
	FILE* logfile = NULL;
	if(logfname.size()>0){
		logfile=fopen(logfname.c_str(), "w");
			//printf("I can't open logfile %s\n", logfname.c_str());
		//else printf("logfile %s opened\n",logfname.c_str());
	}
	const double invPhi0 = 1/pot.value(coord::PosCyl(0,0,0));
	std::vector<double> gridR = potential::createInterpolationGrid
				    (pot, ACCURACY_INTERP2);
	const int sizeE=gridR.size(), sizeXi = 25;//sizeE/2;
	std::vector<double> gridE(sizeE), gridEscaled(sizeE);
	std::vector<double> gridL(sizeE), gridLscaled(sizeE);
	const math::ScalingSemiInf scalingSemi;
	for (int i = 0; i < sizeE; i++) {
		gridE[i] = pot.value(coord::PosCyl(gridR[i],0,0));
		gridEscaled[i] = scaleE(gridE[i], invPhi0);
		gridL[i] = L_circ(pot, gridE[i]);
		gridLscaled[i] = scale(scalingSemi,gridL[i]);
	}
	std::vector<double> gridXi(sizeXi), gridXiscaled(sizeXi);
	math::ScalingCub scalingXi(0, 1);//fix Xi grid so it's densest near ends
	for(int i=0; i<sizeXi; i++){
		gridXiscaled[i] = math::unscale(scalingXi, i/(double)(sizeXi-1));
		gridXi[i] = math::unscale(scalingXi, gridXiscaled[i]);  // Lzrel = u(chi)
	}
	math::Matrix<double> grid2dDE(sizeE,sizeXi);  // focal distance
	math::Matrix<double> grid2dRE(sizeE,sizeXi);  // Rshell / Rcirc(E)
	math::Matrix<double> grid2dL(sizeE,sizeXi);  // L = Jz+Jphi values
	math::Matrix<double> grid2dDL(sizeE,sizeXi);  // focal distance
	math::Matrix<double> grid2dRL(sizeE,sizeXi);  // Rshell
	std::string errorMessage;  // store the error text in case of an exception in the openmp block
	#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for(int iE = 0; iE < sizeE-1; iE++) {//avoid almost free orbits
		/*int nth=0;
#ifdef _OPENMP
		nth = omp_get_thread_num();
#endif
		//*/
		double E  = gridE[iE];
		double Rc = R_circ(pot, E);
		double Jz, Lc = gridL[iE];
		std::vector<double> L_vals(sizeXi);
		std::vector<double> Xi_vals(sizeXi);
		std::vector<double> D_vals(sizeXi);
		std::vector<double> R_vals(sizeXi);
		for(int iXi=0; iXi < sizeXi-1; iXi++){//omit circ orbits
			try{
				double Jphi = Lc * gridXi[iXi];
				std::vector<coord::PosVelCyl> shell;
				double Rsh, FD0 = estimateFocalDistanceShellOrbit
					(pot, E, Jphi, &Rsh, &Jz, &shell);
				//printf("FD:%f\n",FD0);
				double FD = FD0;
				//if(iXi==0) writeShell(iE,Rsh,FD0,shell);
				grid2dDE(iE, iXi) = FD;
				grid2dRE(iE, iXi) = Rsh / Rc;
				double L = Jphi + Jz;
				L_vals[iXi] = L; Xi_vals[iXi] = Jphi / L;
				D_vals[iXi] = FD; R_vals[iXi] = Rsh;
			}
			catch(std::exception& ex) {
				printf("%s",ex.what());
			}
		}
	// values for circular (planar) orbits
		grid2dDE(iE, 0) = grid2dDE(iE, 1);//don't trust Jphi=0 result
		grid2dDE(iE,sizeXi-1) = grid2dDE(iE, sizeXi-2);
		grid2dRE(iE,sizeXi-1) = 1;// Rsh = Rc
		D_vals[sizeXi - 1] = D_vals[sizeXi - 2];
		R_vals[sizeXi - 1] = R_vals[sizeXi - 2];
		L_vals[sizeXi - 1] = Lc;
		Xi_vals[sizeXi - 1] = 1;
	//Now interpolate values L, D, R at given E onto
	//regular grid in Jz/L
		isMonotone(Xi_vals);
		math::LinearInterpolator interpL(Xi_vals, L_vals);
		math::LinearInterpolator interpD(Xi_vals, D_vals);
		math::LinearInterpolator interpR(Xi_vals, R_vals);
		for (int iXi = 0; iXi < sizeXi; iXi++) {
			interpL.evalDeriv(gridXi[iXi], &grid2dL(iE, iXi));
			interpD.evalDeriv(gridXi[iXi], &grid2dDL(iE, iXi));
			interpR.evalDeriv(gridXi[iXi], &grid2dRL(iE, iXi));
		}
	}
    // limiting case of E=0 - assume a Keplerian potential at large radii
	for(int iXi=0; iXi<sizeXi; iXi++) {
		grid2dDE(sizeE-1, iXi) = grid2dDE(sizeE-2, iXi);
		grid2dRE(sizeE-1, iXi) = 1;  // Rshell = Rcirc
		grid2dL (sizeE-1, iXi) = L_circ(pot, gridE[sizeE-1]);
		grid2dDL(sizeE-1, iXi) = grid2dDL(sizeE-2,iXi);
		grid2dRL(sizeE-1, iXi) = R_circ(pot, gridE[sizeE-1]);// Rshell = Rcirc
	}
	if(logfile){
		fprintf(logfile,"%d\n",sizeE);
		fprintf(logfile,"  E  Delta(E,0) Rsh(E,0) R(E,0)\n");      
		for(int i=0; i<sizeE; i++) fprintf(logfile,"%g %g %g %g\n",
			gridE[i],grid2dDE(i,0),grid2dRE(i,0),gridR[i]);
	}
	//grid2dD contains D on regular grid in Xi but irregular
	//values of L that are stored in grid2dL
	for (int iXi = 0; iXi < sizeXi; iXi++) {
		std::vector<double> L_vals(sizeE);
		std::vector<double> D_vals(sizeE);
		std::vector<double> R_vals(sizeE);
		for (int iE = 0; iE < sizeE; iE++) {
			L_vals[iE] = grid2dL(iE, iXi);
			D_vals[iE] = grid2dDL(iE, iXi);
			R_vals[iE] = grid2dRL(iE, iXi);
		}
		//We now have D and R at series of L values
		isMonotone(L_vals);
		math::LinearInterpolator DL(L_vals, D_vals);
		math::LinearInterpolator RL(L_vals, R_vals);
		for (int iE = 0; iE < sizeE; iE++) {
			DL.evalDeriv(gridL[iE], &grid2dDL(iE, iXi));
			RL.evalDeriv(gridL[iE], &grid2dRL(iE, iXi));
		}

	}
	if(!errorMessage.empty())
		throw std::runtime_error(errorMessage);
	std::vector<double> gridJr, gridJz, gridI3, gridFD, scaledJ;
	interpDE = math::LinearInterpolator2d(gridEscaled, gridXiscaled, grid2dDE);
	interpRE = math::LinearInterpolator2d(gridEscaled, gridXiscaled, grid2dRE);
	interpDL = math::LinearInterpolator2d(gridLscaled, gridXiscaled, grid2dDL);
	interpRL = math::LinearInterpolator2d(gridLscaled, gridXiscaled, grid2dRL);
	if(logfile) fclose(logfile);
}
PolarInterpolator::PolarInterpolator(const potential::BasePotential& pot){
					 //const PtrShellInterpolator PtrShellI) {
			//spherical potential only has loop orbits by conservation of angular momentum.
			// Also if potential infinite as centre also always loop orbits.
	double Phi0 = pot.value(coord::PosCar(0, 0, 0));//potential at centre
	std::vector<double> gridR = potential::createInterpolationGrid
				    (pot, ACCURACY_INTERP2);
	const int sizeE=gridR.size();
	std::vector<double> gridE(sizeE), gridEscaled(sizeE);
	for (int i = 0; i < sizeE; i++) {
		gridE[i] = pot.value(coord::PosCyl(gridR[i],0,0));
		gridEscaled[i] = scaleE(gridE[i], 1/Phi0);
	}
	math::ScalingSemiInf Sc;
	std::vector<double> gridJz(sizeE), gridJfScaled(sizeE);
	//gridI3, gridFD, gridUmin;
	if (potential::isSpherical(pot)|| std::isnan(Phi0) || std::isinf(Phi0)) {
		double fac=sizeE<=1?1:1/((double)(sizeE-1));
		for (int i = 0; i < sizeE; i++) {
			gridJz[i] = 0;
			/*gridFD.push_back(0);
			gridI3.push_back(0);//change this
			//*/
			gridJfScaled[i]=i*fac;
		}
	} else {
		std::vector<double> gridJr(sizeE);
		int N = 7;//number of points to be fitted
		bool fitted = false;//gets if fitted straight point in E,z1 plane
		bool interp = false;
		double E0, z0;
		double b;//z=b(E-E0)+z0
		int i = 0;
		std::vector<double> prevz(N),prevx(N);
		//double Rsh, vR, Umin, d2pu2du2;
		while(i<sizeE) {
			//Rsh = PtrShellI->getRsh(gridE[i], 0, 1/Phi0) * R_circ(pot, gridE[i]);
			//double FD = PtrShellI->getDelta(gridE[i], 0, 1/Phi0);
			if (gridE[i] / Phi0 > 0.2 && !interp) {
				actions::Actions Jcrit = BoxLoopTrAct(pot, gridE[i]);
				/*FDfinder FDf(gridE[i], Rsh, vR, FD, pot);
				if(i==5) testFDfinder(gridE[i], Rsh, vR, FD, pot);
				gridFD.push_back(FDf.bestFD(Umin, d2pu2du2));
				gridUmin.push_back(Umin);
				const coord::ProlSph coordsys(gridFD.back());
				double vzsq = 2*(gridE[i] - pot.value(coord::PosCyl(Rsh,0,0))) - pow_2(vR); 
				double vz = vzsq>0? sqrt(vzsq) : 0;
				coord::PosVelCyl point(Rsh,0,0,vR,vz,0);
				//gridI3.push_back(actions::getI3(pot, point, coordsys));
				//*/
				gridJr[i] = Jcrit.Jr; gridJz[i] = Jcrit.Jz;
				//printf("J:%f %f\n",Jcrit.Jr*2+Jcrit.Jz,Jcrit.Jz);
				if ((i > N+2) && (1-gridE[i]/Phi0) > 1e-4) {
					if (gridJz[i] < gridJz[i-1]) {
						interp = true;
					}
				}
			}
			else {
				if (!fitted) {
					interp = true;
					int index = i - 2;
					if (N > index)N = index;
					E0 = gridE[index-N+1];
					Ifind fI(.5 * M_PI * gridJz[index-N+1], pot, gridE[index-N+1]);
					double zmax0 = potential::z_max(pot, gridE[index]);
					std::vector<double> zmaxs(N-1);
					z0 = math::findRoot(fI, 0, zmax0, 1e-6);
					std::vector<double> x(N - 1), y(N - 1);
					for (int j = 0;j < N-1;j++) {
						int j1 = index + j - N + 2;
						Ifind fI(.5 * M_PI * gridJz[j1], pot, gridE[j1]);
						double zmax = potential::z_max(pot, gridE[j1]);
						x[j] = (gridE[j1] - E0);
						y[j] = math::findRoot(fI, 0, zmax, 1e-6)- z0;
						zmaxs[j]=zmax;
					}
					b=math::linearFitZero(x, y,NULL);
					for(int j=0;j<N-1;j++){
						int j1 = index + j - N + 2;
						pzf pzfunc(pot,gridE[j1]);
						double z=z0+b*(gridE[j1]-E0);
						if(z>zmaxs[j])z=zmaxs[j];
						if(z<0)z=0;
						double gridJf=2*gridJr[j1]+gridJz[j1];
						gridJz[j1]=2*math::integrateGK(pzfunc,0,z,1e-8)/M_PI;
						gridJr[j1]=.5*(gridJf-gridJz[j1]);
					}
					fitted = true;
					i--;
				}
				pzf pzfunc(pot,gridE[i]);
				double zmax = potential::z_max(pot,gridE[i]);
				double z1 = fmin(zmax,b * (gridE[i] - E0) + z0);
				if(z1<0)z1=0;
				gridJz[i] = 2 * math::integrateGK(pzfunc, 0, z1, 1e-8) / M_PI;
				gridJr[i] = math::integrateGK(pzfunc, z1, zmax,1e-8) / M_PI;
			}
			i++;
		}
		//int sizeFD = gridFD.size();
		for(int i=0; i<sizeE; i++){
			gridJfScaled[i]=scale(Sc,2*gridJr[i]+gridJz[i]);
			/*if(i>=sizeFD){
				gridFD.push_back(0);
				gridI3.push_back(0);
				gridUmin.push_back(0);
			}
            //*/
		}
	}
	/*interpI3 = math::LinearInterpolator(gridEscaled, gridI3);
	interpFD = math::LinearInterpolator(gridEscaled, gridFD);
	interpUmin=math::LinearInterpolator(gridEscaled, gridUmin);
    //*/
	interpJzJ=math::LinearInterpolator(gridJfScaled,gridJz);
	interpJzE=math::LinearInterpolator(gridEscaled,gridJz);
	//coeffsJz = math::fitPoly(15, gridJfScaled, gridJz);
}double estimateFocalDistanceShellOrbit(
					   const potential::BasePotential& poten, double E, double Jphi,
					   double* _Rshell, double* _Jz, std::vector<coord::PosVelCyl>* shell)
{
	double Rmin, Rmax, FD;
	findPlanarOrbitExtent(poten, E, Jphi, Rmin, Rmax);
	//double timeCross = 5*pow_2(Rmin)/Jphi, Jz;
	double Rbar=.5*(Rmin+Rmax), Jz;
	double timeCross = 5*pow_2(Rbar)/(2*(E-poten.value(coord::PosCyl(Rbar,0,0))));
	//if(Jphi<1e-4) printf("Rmin..%g %g %g %g\n",Rmin,Rmax,Jphi,timeCross);
	std::vector<std::pair<coord::PosVelCyl,double> > traj;
	FindClosedOrbitRZplane finder(poten, E, Jphi, Rmin, Rmax, timeCross, traj, Jz);
    // locate the radius of a shell orbit;  as a by-product, store the orbit in 'traj'
	double Rshell = math::findRoot(finder, Rmin, Rmax, ACCURACY_RSHELL);
	if(shell){//return shell orbit up to p_theta=0
		for(int i=0; i<traj.size(); i++)
			shell->push_back(traj[i].first);
	}
	
//#define OLD_METHOD 1
#ifdef OLD_METHOD
	/*if(traj.size() >= 2){
	// now find the best-fit value of delta for this orbit
		FD = fitFocalDistanceShellOrbit(*shell);
		printf("FD:%f\n",FD);
	} else {
	// something went wrong; use a backup solution
		if(!isFinite(Rshell))
			Rshell = 0.5 * (Rmin+Rmax);
		utils::msg(utils::VL_WARNING, FUNCNAME,
			   "Could not find a thin orbit for E="+utils::toString(E,16)+", Jphi="+utils::toString(Jphi,16)+
			   " - assuming Rthin="+utils::toString(Rshell,16));
		FD = actions::estimateFocalDistancePoints(poten, std::vector<coord::PosCyl>(1,
            coord::PosCyl(Rshell,Rshell*ACCURACY_RSHELL, 0)))
		
	}
			//*/
	FD = actions::estimateFocalDistancePoints(poten, std::vector<coord::PosCyl>(1,
            coord::PosCyl(Rshell, /* z=very small number */ Rshell*ACCURACY_RSHELL, 0)));
#else
	double Phi;
	coord::GradCyl grad;
	poten.eval(coord::PosCyl(Rshell,0,0), &Phi, &grad);
	double vphi = Jphi!=0 ? Jphi / Rshell : 0;
	FD = Rshell * sqrt( math::clip((2 * (E-Phi) - Rshell * grad.dR) / ( Rshell * grad.dR - vphi*vphi), 0., 1e6) );
#endif
	if(_Jz) *_Jz = Jz;
	if(_Rshell) *_Rshell = Rshell; 
	return FD;
	
}
}
