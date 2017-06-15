#pragma once

#include "physics/Physics.h"
#include <map>
#include "base/WrapTTree.h"
#include "utils/TriggerSimulation.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

/**
 * @brief A ppi0_2gamma class playground for learning stuff
 * Class uses basic kinematic cuts to find gp->pi0p (pi0->2g)
 *
 */

class scratch_lheijken_gppi0p: public Physics {
protected:

public:
    scratch_lheijken_gppi0p(const std::string& name, OptionsPtr opts);
    virtual ~scratch_lheijken_gppi0p() {}

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;

};

}
}
}
