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
    promptrandom(ExpConfig::Setup::Get())
{
    CreateHistos();
    steps.CreateBranches(HistFac.makeTTree("analysis_cuts"));
}


void scratch_lheijken_gppi0p::ProcessEvent(const TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event);

    //-- Check the decay string for signal MC pattern
    bool MCsignal  = false;
    string decay;
    if(event.Reconstructed().ID.isSet(ant::TID::Flags_t::MC))
            decay = utils::ParticleTools::GetDecayString(event.MCTrue().ParticleTree);
    else    decay = "data" + ExpConfig::Setup::Get().GetName();
    if(decay == "(#gamma p) #rightarrow #pi^{0} [ #gamma #gamma ] p ") MCsignal = true;

    steps.AddStep(MCsignal, "All events");

    //-- MC true stuff
    if(event.Reconstructed().ID.isSet(ant::TID::Flags_t::MC)){
        //--- Fetches a list of all gammas in the MCTrue tree
        auto GammasTrue = utils::ParticleTools::FindParticles(ParticleTypeDatabase::Photon,event.MCTrue().ParticleTree);
        for(const auto& ph : GammasTrue) {
            hTrueGammaE->Fill(ph->Ek());
        }
        utils::ParticleTools::FillIMCombinations(hTrueGammaIM,2,GammasTrue);
    }

    //-- Get full list of particles, create neutral/charged list
    const auto& candidates = event.Reconstructed().Candidates;
    TCandidatePtrList neutral;
    TCandidatePtrList charged;
    for(const auto& cand : candidates.get_iter()) {
        if(cand->VetoEnergy == 0.0) neutral.emplace_back(cand);
        else charged.emplace_back(cand);
    }

    //-- Particle combinatorics selection
    if (neutral.size() < 2) return; // neutral == 2, 3 is ok
    if (neutral.size() > 3) return; // neutral == 2, 3 is ok
    steps.AddStep(MCsignal, "# Neutral == 2 or 3");
    steps.AddStep(MCsignal, "# Neutral == "+to_string(neutral.size()));
    if (charged.size()  > 1) return; // charged == 0, 1 is ok
    steps.AddStep(MCsignal, "# Charged == 0 or 1");
    steps.AddStep(MCsignal, "# Charged == "+to_string(charged.size()));


    //-- Create pion candidate
    TParticle Meson;
    double Meson_time = 0;
    for (const auto& photon : neutral) {
        Meson+=TParticle(ParticleTypeDatabase::Photon, photon);
        Meson_time+=photon->Time;
    }
    Meson_time = Meson_time / neutral.size(); // Avg the meson time
    //--- Have a look at it
    hRecMesIMN1->Fill(Meson.M());
    if (neutral.size()==2) hRecMesIMN2->Fill(Meson.M());
    if (neutral.size()==3) hRecMesIMN3->Fill(Meson.M());
    //--- Ignore the events with three clusters
    if (neutral.size()==3) return;


    //-- Loop over the Tagger hits
    int i=0;
    for (const auto &tc : event.Reconstructed().TaggerHits) {
        //--- Have a look at the missing mass of the meson
        auto missing = tc.GetPhotonBeam() + LorentzVec(vec3(0,0,0),ParticleTypeDatabase::Proton.Mass()) - Meson;
        hMMMeson->Fill(missing.M(),promptrandom.FillWeight());
        //--- Have a look at the coplanarity of the meson and the proton
        double coplanarity = -1;
        if(charged.size() > 0) {
            coplanarity  = abs(Meson.Phi() - charged.at(0)->Phi)*radtodeg;
            hCoplMesProt->Fill(coplanarity,promptrandom.FillWeight());
        }

        if(missing.M() > 1300) continue;
        steps.AddStep(MCsignal, "MM < 1300", promptrandom.FillWeight());
        if(coplanarity < 150 || coplanarity > 210) continue;
        steps.AddStep(MCsignal, "Cop = (150-210)", promptrandom.FillWeight());

        //--- Have a look at the mesons after the requirements on the proton
        hRecMesIMPCut->Fill(Meson.M(), promptrandom.FillWeight());
        i++;

        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(tc)); // - vad hander har?
    }
    if(i>0) hRecMesIMPCut2->Fill(Meson.M());


}

void scratch_lheijken_gppi0p::Finish()
{
}

void scratch_lheijken_gppi0p::ShowResult()
{
/*
    canvas(GetName())
             <<  hTrueGammaE
             << hTrueGammaIM
             << endc; // actually draws the canvas

    ant::canvas(GetName()+": Analysis cuts")
            << TTree_drawable(steps.Tree, "cut.c_str()>>cuts_signal","promptrandom*isSignal")
            << TTree_drawable(steps.Tree, "cut.c_str()>>cuts_background","promptrandom*(!isSignal)")
            << TTree_drawable(steps.Tree, "isSignal","promptrandom")
            << TTree_drawable(steps.Tree, "isSignal","promptrandom","signalcount","0 = bkg, 1 = sig","",BinSettings(2))
            << endc; // actually draws the canvas
*/
}

void scratch_lheijken_gppi0p::CreateHistos()
{
    auto hfTrueMC = new HistogramFactory("TrueMC", HistFac, "");
    auto hfPartKin = new HistogramFactory("PartKin", HistFac, "");

    const BinSettings ETGBins(100,0.,500.);
    const BinSettings ProBins(200,0.,2000.);
    const BinSettings PhiBins(100,0.,360.);

    hTrueGammaE   = hfTrueMC->makeTH1D("energies of true gammas","E_{#gamma}","",ETGBins,"TrueGammaE");
    hTrueGammaIM  = hfTrueMC->makeTH1D("IM of true gammas","m_{#gamma#gamma}","",ETGBins, "TrueGammasIM");
    hRecMesIMN1 = hfPartKin->makeTH1D("IM of rec meson 2+3 neutral","m_{#gamma#gamma}","",ETGBins, "RecMesIMN1");
    hRecMesIMN2 = hfPartKin->makeTH1D("IM of rec meson 2 neutral","m_{#gamma#gamma}","",ETGBins, "RecMesIMN2");
    hRecMesIMN3 = hfPartKin->makeTH1D("IM of rec meson 3 neutral","m_{#gamma#gamma}","",ETGBins, "RecMesIMN3");
    hRecMesIMPCut = hfPartKin->makeTH1D("IM of rec meson cut on proton","m_{#gamma#gamma}","",ETGBins, "RecMesIMPCut");
    hRecMesIMPCut2 = hfPartKin->makeTH1D("IM of rec meson cut on proton","m_{#gamma#gamma}","",ETGBins, "RecMesIMPCut2");
    hMMMeson = hfPartKin->makeTH1D("MM of rec meson","MM(#pi^{0})","",ProBins, "MMMeson");
    hCoplMesProt = hfPartKin->makeTH1D("Coplanarity of rec meson/proton","","abs(#phi_{#pi^{0}}-#phi_{p}) [deg]",PhiBins, "CoplMesProt");
}
AUTO_REGISTER_PHYSICS(scratch_lheijken_gppi0p)

