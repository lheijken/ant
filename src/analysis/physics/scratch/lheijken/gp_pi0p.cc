#include "gp_pi0p.h"

#include "base/ParticleType.h"
#include "plot/HistogramFactory.h"
#include "utils/Combinatorics.h"
#include "utils/ParticleTools.h"
#include "base/ParticleTypeTree.h"
#include "utils/uncertainties/FitterSergey.h"
#include "expconfig/ExpConfig.h"
#include "base/Logger.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "base/WrapTFile.h"

#include <iostream>
#include <memory>

/**
 * Single pion photoproduction, decay to 2gamma
 */

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

scratch_lheijken_gppi0p::scratch_lheijken_gppi0p(const std::string& name, OptionsPtr opts):
    Physics(name, opts),
{

}


void scratch_lheijken_gppi0p::ProcessEvent(const TEvent& event, manager_t&)
{
    // If MC, track detection efficiency
    // *************************************
    if((event.Reconstructed().ID.isSet(ant::TID::Flags_t::MC)) && (signal)){

        auto pi0 = utils::ParticleTools::FindParticle(ParticleTypeDatabase::Pi0,event.MCTrue().ParticleTree);
        if(!pi0) LOG(ERROR) << "Did not find pi0 in MC True ... strange!!!";

        for (const auto &tc : event.Reconstructed().TaggerHits){

            promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(tc));
            detection_efficiency.TrackSignalEvent(*pi0, 0, tc, promptrandom);
        }
    }



    // Get full list of particles, create neutral/charged list
    // *******************************************************
    const auto& candidates = event.Reconstructed().Candidates;
    TCandidatePtrList neutral;
    TCandidatePtrList charged;

    for(const auto& cand : candidates.get_iter()) {
        if(cand->VetoEnergy == 0.0) neutral.emplace_back(cand);
        else charged.emplace_back(cand);
    }



    // ===================================================================
    // Basic event selection
    if (neutral.size() < 2) return; // neutral == 2, 3 is ok
    if (neutral.size() > 3) return; // neutral == 2, 3 is ok
    steps.AddStep(signal, "# Neutral == 2 or 3");
    steps.AddStep(signal, "# Neutral == "+to_string(neutral.size()));

    if (charged.size()  > 1) return; // charged == 0, 1 is ok
    steps.AddStep(signal, "# Charged == 0 or 1");
    steps.AddStep(signal, "# Charged == "+to_string(charged.size()));
    // ===================================================================



    // Create pion candidate
    // *********************
    TParticle Meson;
    double Meson_time = 0;
    for (const auto& photon : neutral)
    {
        Meson+=TParticle(ParticleTypeDatabase::Photon, photon);
        Meson_time+=photon->Time;
    }
    Meson_time = Meson_time / neutral.size(); // Avg the meson time


    // Loop over the Tagger hits
    // *************************
    int i = 0;
    for (const auto &tc : event.Reconstructed().TaggerHits)
    {
        auto missing = tc.GetPhotonBeam() + LorentzVec(vec3(0,0,0),ParticleTypeDatabase::Proton.Mass()) - Meson;

        // ===================================================================
        // Advanced event selection
        if(Meson.M() < 80)  continue;
        if(Meson.M() > 200) continue;
        steps.AddStep(signal, "IM = (80-200)", promptrandom.FillWeight());

        if(charged.size() > 0)
        {
            auto coplanarity  = abs(Meson.Phi() - charged.at(0)->Phi)*radtodeg;
            if (coplanarity < 150) continue;
            if (coplanarity > 210) continue;
            steps.AddStep(signal, "Cop = (150-210)", promptrandom.FillWeight());

        }

        if(missing.M() > 1300) continue;
        steps.AddStep(signal, "MM < 1300", promptrandom.FillWeight());

        // ===================================================================

        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(tc));

        detection_efficiency.AcceptEvent(Meson,Meson_time, tc,promptrandom);
        cross_section.AcceptEvent(Meson,Meson_time,tc,promptrandom);

        i++;
    }


}

void scratch_lheijken_gppi0p::Finish()
{
}

void scratch_lheijken_gppi0p::ShowResult()
{
    ant::canvas(GetName()+": Analysis cuts")
            << TTree_drawable(steps.Tree, "cut.c_str()>>cuts_signal","promptrandom*isSignal")
            << TTree_drawable(steps.Tree, "cut.c_str()>>cuts_background","promptrandom*(!isSignal)")
            << TTree_drawable(steps.Tree, "isSignal","promptrandom")
            << TTree_drawable(steps.Tree, "isSignal","promptrandom","signalcount","0 = bkg, 1 = sig","",BinSettings(2))
            << endc; // actually draws the canvas
}

AUTO_REGISTER_PHYSICS(scratch_lheijken_gppi0p)

