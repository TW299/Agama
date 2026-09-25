#include "actions_factory.h"
#include "actions_newisochrone.h"
#include "actions_spherical.h"
#include "actions_staeckel.h"
#include "actions_staeckel3d.h"
#include "actions_torusmapper.h"
#include "potential_analytic.h"
#include "potential_perfect_ellipsoid.h"
#include "potential_perfect_ellipsoid_triaxial.h"

namespace actions {

void eval(const potential::BasePotential& pot, const coord::PosVelCyl& point,
    Actions* act, Angles* ang, Frequencies* freq, double focalDistance)
{
    const potential::Isochrone* potIso = dynamic_cast<const potential::Isochrone*>(&pot);
    if(potIso) {
        evalIsochrone(potIso->totalMass(), potIso->getRadius(), point, act, ang, freq);
        return;
    }

    if(isSpherical(pot)) {
        evalSpherical(pot, point, act, ang, freq);
        return;
    }

    const potential::OblatePerfectEllipsoid* potPE =
        dynamic_cast<const potential::OblatePerfectEllipsoid*>(&pot);
    if(potPE) {
        evalAxisymStaeckel(*potPE, point, act, ang, freq);
        return;
    }

    evalAxisymFudge(pot, point, act, ang, freq, focalDistance);
}

PtrActionFinder createActionFinder(const potential::PtrPotential& pot, bool interpolate)
{
    const potential::Isochrone* potIso = dynamic_cast<const potential::Isochrone*>(pot.get());
    if(potIso)
        return PtrActionFinder(new ActionFinderIsochrone(potIso->totalMass(), potIso->getRadius()));

    if(isSpherical(*pot))
       return PtrActionFinder(new ActionFinderSpherical(*pot));

#if __cplusplus >= 201103L   // aliasing constructor for shared pointers only works in C++11
    const potential::OblatePerfectEllipsoid* potPE =
        dynamic_cast<const potential::OblatePerfectEllipsoid*>(pot.get());
    if(potPE)
        return PtrActionFinder(new ActionFinderAxisymStaeckel(
            potential::PtrOblatePerfectEllipsoid(pot, potPE)));
    const potential::PerfectEllipsoidTriaxial* potTri =
        dynamic_cast<const potential::PerfectEllipsoidTriaxial*>(pot.get());
    if(potTri)
        return PtrActionFinder(new ActionFinderTriaxialStaeckel(
            potential::PtrPerfectEllipsoidTriaxial(pot, potTri)));
#endif
             
    return PtrActionFinder(new ActionFinderAxisymFudge(pot, interpolate));
}

PtrActionMapper createActionMapper(const potential::PtrPotential& pot, double tol)
{
    const potential::Isochrone* potIso = dynamic_cast<const potential::Isochrone*>(pot.get());
    if(potIso)
        return PtrActionMapper(new ActionMapperIsochrone(potIso->totalMass(), potIso->getRadius()));

    if(isSpherical(*pot))
        return PtrActionMapper(new ActionMapperSpherical(*pot));

    if(tol==tol)  // non-default value for tol
        return PtrActionMapper(new ActionMapperTorus(*pot, tol));
    else
        return PtrActionMapper(new ActionMapperTorus(*pot /*, default_value_for_tol */));
}

}  // namespace actions
