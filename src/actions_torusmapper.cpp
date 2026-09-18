#include "actions_torusmapper.h"
//From actions_torus.cpp in AGAMA. used to cache torus objects;
#if __cplusplus >= 201103L
// with C++11 use unordered map as it is faster
#include <unordered_map>
#include "math_random.h"

namespace{
struct ActionsHash {
    size_t operator() (const actions::Actions& a) const {
        return math::hash((const void*)(&a), 3);
    }
};

struct ActionsEqual {
    bool operator() (const actions::Actions& lhs, const actions::Actions& rhs) const {
        return lhs.Jr == rhs.Jr && lhs.Jz == rhs.Jz && lhs.Jphi == rhs.Jphi;
    }
};
typedef std::unordered_map<actions::Actions, shared_ptr<actions::Torus>, ActionsHash, ActionsEqual>
    TorusCache;
}
#else
// use ordinary map pre-C++11
#include <map>
namespace {
struct ActionsLess {
    bool operator() (const actions::Actions& lhs, const actions::Actions& rhs) const {
        if(lhs.Jr < rhs.Jr) return true;
        if(lhs.Jr > rhs.Jr) return false;
        if(lhs.Jz < rhs.Jz) return true;
        if(lhs.Jz > rhs.Jz) return false;
        return lhs.Jphi < rhs.Jphi;
    }
};

typedef std::map<Actions, shared_ptr<actions::Torus>, ActionsLess> TorusCache;
}
#endif
namespace actions{
class ActionMapperTorus::Impl {
public:
    Impl(const potential::BasePotential& pot,const double tol):TG(pot,tol){}
    Torus* getTorus(const Actions J){
        TorusCache::iterator it=cache.find(J);
        if(it != cache.end())
        return it->second.get();
        std::shared_ptr<Torus> ptrT(new Torus(TG.fitTorus(J)));
        cache.insert(std::make_pair(J, ptrT));
        return ptrT.get();
    }
private:
    const TorusGenerator TG;
    TorusCache cache; 
};
ActionMapperTorus::ActionMapperTorus(const potential::BasePotential& pot,const double tol):impl(new Impl(pot,tol)){}
coord::PosVelCyl ActionMapperTorus::map(const ActionAngles& actAng, Frequencies* freq) const{
    Torus* T=impl->getTorus(actAng);
    if(T){
        if(freq!=NULL){
            *freq=T->freqs;
        }
        coord::PosMomCyl dPH=T->from_true(actAng);
        return coord::toPosVelCyl(dPH);
    }
}
}