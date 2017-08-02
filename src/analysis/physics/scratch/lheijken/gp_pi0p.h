#pragma once

#include "physics/Physics.h"
#include <map>
#include "base/WrapTTree.h"
#include "utils/TriggerSimulation.h"
#include "analysis/plot/PromptRandomHist.h"

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
    struct steps_t : WrapTTree {

        ADD_BRANCH_T(bool, isSignal)
        ADD_BRANCH_T(std::string, cut)
        ADD_BRANCH_T(double, promptrandom)

        void AddStep(bool _isSignal, const std::string& _cut, double fillweight = 1)
        {
            isSignal = _isSignal;
            cut = _cut;
            promptrandom = fillweight;
            Tree->Fill();
        }
    };

    steps_t steps;
    PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;
    TH1D* TrueGammaE = nullptr;
    TH1D* TrueGammaIM = nullptr;
    TH1D* RecMesIM = nullptr;

protected:
    static constexpr auto radtodeg = std_ext::radian_to_degree(1.0);

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
