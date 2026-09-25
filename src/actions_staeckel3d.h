/** \file    actions_staeckel2.h
    \brief   Action-angle finders using triaxial Staeckel potential approximation
    \author  Eugene Vasiliev (for actions_staeckel code which is slightly modified for triaxial case)
    \date    2015

Computation of actions and angles for non axisymmetric perfect ellipsoid

Most sub-steps are shared between the two methods; 
only the computation of integrals of motion, and the auxiliary function that enters 
the expression for canonical momentum, are specific to each method.

The implementation is inspired by the code written by Jason Sanders,
but virtually nothing of the original code remains.
*/
#pragma once
#include "actions_base.h"
#include "potential_perfect_ellipsoid_triaxial.h"
#include "smart.h"

namespace actions {

/// \name  ------- Stand-alone driver routines that compute actions for a single point -------
///@{

/** Evaluate (exactly up to integration errors) any combination of actions, angles and frequencies
    in the Staeckel potential of an oblate Perfect Ellipsoid.
    \param[in]  potential is the input Staeckel potential.
    \param[in]  point is the position/velocity point.
    \param[out] act   if not NULL, will contain computed actions (Jr=Jz=NAN if E>=0).
    \param[out] ang   if not NULL, will contain corresponding angles (NAN if E>=0).
    \param[out] freq  if not NULL, will contain corresponding frequencies (NAN if E>=0).
*/
void evalTriaxialStaeckel(
    const potential::PerfectEllipsoidTriaxial& potential, 
    const coord::PosVelCar& point,
    Actions* act=NULL,
    Angles* ang=NULL,
    Frequencies* freq=NULL);

class ActionFinderTriaxialStaeckel:public BaseActionFinder{
public:
    explicit ActionFinderTriaxialStaeckel(const potential::PtrPerfectEllipsoidTriaxial& potential):pot(potential){}
    
    virtual std::string name() const{return "TriaxialStaeckel";}

    virtual void eval(const coord::PosVelCyl& point,
        Actions* act=NULL, Angles* ang=NULL, Frequencies* freq=NULL) const{
        evalTriaxialStaeckel(*pot,coord::toPosVelCar(point),act,ang,freq);
    }

private:
    const potential::PtrPerfectEllipsoidTriaxial pot;
};
///@}
}  // namespace actions
