#include "potential_interpolators.h"
#include "potential_utils.h"
namespace potential{

/// accuracy parameter determining the spacing of the interpolation grid along the energy axis
static const double ACCURACY_INTERP2 = 1e-5;

namespace{
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
	if(pRm/v1>0.8){//From Rmax to Rmax0 approx motion as on x axis 
		pxf pxfunc(pot,E);
		Jr+=math::integrateGK(math::ScaledIntegrand<math::ScalingCub>
				      (math::ScalingCub(Rmax, Rmax0), pxfunc), 0, 1, 1e-8)/M_PI;
	}
	//Get Jfast=Jx+Jz=2*Jr+Jz from motion along the z axis
	double Jfast = 2 * math::integrateGK(math::ScaledIntegrand<math::ScalingCub>
					     (math::ScalingCub(0, zmax), pzt), 0, 1, 1e-8) / M_PI;
	double Jz = Jfast - 2 * Jr;
	if(Jz<0){
		Jz=Jfast;
		Jr=0;
	}
	return actions::Actions(Jr, Jz, 0);
}
const PolarInterpolator PolarInterpolatorSph;
PolarInterpolator::PolarInterpolator(const potential::BasePotential& pot){
					 //const PtrShellInterpolator PtrShellI) {
			//spherical potential only has loop orbits by conservation of angular momentum.
			// Also if potential infinite as centre also always loop orbits.
	double Phi0 = pot.value(coord::PosCar(0, 0, 0));//potential at centre
	if (potential::isSpherical(pot)|| std::isnan(Phi0) || std::isinf(Phi0)) {
		isSph=true;
		return;
	}else{
		std::vector<double> gridR = potential::createInterpolationGrid
						(pot, ACCURACY_INTERP2);
		const int sizeE=gridR.size();
		std::vector<double> gridE(sizeE), gridEscaled(sizeE);
		for (int i = 0; i < sizeE; i++) {
			gridE[i] = pot.value(coord::PosCyl(gridR[i],0,0));
			gridEscaled[i] = scaleE(gridE[i], 1/Phi0);
		}
		math::ScalingSemiInf Sc;
		std::vector<double> gridJr(sizeE), gridJz(sizeE), gridJfScaled(sizeE);
		//gridI3, gridFD, gridUmin;
		isSph=false;
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
						gridJz[j1]=2*math::integrateGK(pzfunc,0,z,1e-8)/M_PI;
						gridJr[j1]=math::integrateGK(pzfunc,z,zmaxs[j],1e-8)/M_PI;
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
		/*interpI3 = math::LinearInterpolator(gridEscaled, gridI3);
		interpFD = math::LinearInterpolator(gridEscaled, gridFD);
		interpUmin=math::LinearInterpolator(gridEscaled, gridUmin);
		//*/
		interpJzJ=math::LinearInterpolator(gridJfScaled,gridJz);
		interpJzE=math::LinearInterpolator(gridEscaled,gridJz);
		//coeffsJz = math::fitPoly(15, gridJfScaled, gridJz);
	}
}
}
