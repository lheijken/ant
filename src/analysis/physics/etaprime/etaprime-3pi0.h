#pragma once

#include "analysis/utils/fitter/TreeFitter.h"
#include "analysis/utils/uncertainties/FitterSergey.h"
#include "analysis/physics/Physics.h"
#include "analysis/plot/PromptRandomHist.h"
#include "utils/A2GeoAcceptance.h"
#include "analysis/utils/TriggerSimulation.h"


#include "TLorentzVector.h"

class TH1D;
class TTree;

namespace ant {
namespace analysis {
namespace physics {

class Etap3pi0 : public Physics {

protected:
    utils::TriggerSimulation triggersimu;
    //geometry
    ant::analysis::utils::A2SimpleGeometry geometry;

    // =======================   constants =====================================================

    const double IM_mean_etaprime = 895.0;
    const double IM_sigma_etaprime = 26.3;

    const double IM_mean_eta = 515.5;
    const double IM_sigma_eta = 19.4;

    const double IM_mean_pi0 = 126.0;
    const double IM_sigma_pi0 = 15.0;


    const double copl_opening_sigma = 7.80;

    ParticleTypeTree signal_tree;
    ParticleTypeTree reference_tree;
    ParticleTypeTree bkg_tree;

    struct settings_t
    {
        std::map<int,std::string> EventTypes= {{0,"signal"},        // etaprime -> pi0 pi0 pi0
                                               {1,"reference"},     // etaprime -> eta pi0 pi0
                                               {2,"background"},    // 3 pi0 photoproduction
                                               {-1,"other"}};
        const double EMBChi2Cut= std::numeric_limits<double>::infinity();
        const double Chi2CutEMB = 22.5;
        const double Chi2CutSig = 17.5;
        const double Chi2CutRef = 20;
        const double coplCut    = 15;
        const double etaprimeThreshold = 1445.6;
    };
    settings_t phSettings;

    //obsolete...
    const std::vector<std::vector<std::pair<unsigned,unsigned>>> combinations =
    {
        { {0, 1}, {2, 3}, {4, 5} },
        { {0, 1}, {2, 4}, {3, 5} },
        { {0, 1}, {2, 5}, {3, 4} },

        { {0, 2}, {1, 3}, {4, 5} },
        { {0, 2}, {1, 4}, {3, 5} },
        { {0, 2}, {1, 5}, {3, 4} },

        { {0, 3}, {1, 2}, {4, 5} },
        { {0, 3}, {1, 4}, {2, 5} },
        { {0, 3}, {1, 5}, {2, 4} },

        { {0, 4}, {1, 2}, {3, 5} },
        { {0, 4}, {1, 3}, {2, 5} },
        { {0, 4}, {1, 5}, {2, 3} },

        { {0, 5}, {1, 2}, {3, 4} },
        { {0, 5}, {1, 3}, {2, 4} },
        { {0, 5}, {1, 4}, {2, 3} }
    };

    //===================== KinFitting ========================================================

    std::shared_ptr<utils::UncertaintyModel> uncertModel = std::make_shared<utils::UncertaintyModels::FitterSergey>();

    utils::TreeFitter fitterSig;
    std::vector<utils::TreeFitter::tree_t> intermediatesTreeSig= std::vector<utils::TreeFitter::tree_t>(3);

    utils::TreeFitter fitterRef;
    std::vector<utils::TreeFitter::tree_t> intermediatesTreeRef= std::vector<utils::TreeFitter::tree_t>(3);


    utils::KinFitter kinFitterEMB;

    ant::analysis::PromptRandom::Switch promptrandom;





    //========================  Storage  ============================================================
    TTree* tree;

    struct branches {
        struct kinFitReturn_t
        {
            double          beamE= {};
            std::vector<TLorentzVector> intermediatesSig= std::vector<TLorentzVector>(3);
            std::vector<TLorentzVector> gammasSig= std::vector<TLorentzVector>(6);
            std::vector<TLorentzVector> intermediatesRef= std::vector<TLorentzVector>(3);
            std::vector<TLorentzVector> gammasRef= std::vector<TLorentzVector>(6);
            TLorentzVector  etaprimeCand= {};
            TLorentzVector  p= {};
        };

        kinFitReturn_t kinfitted= {};

        TLorentzVector etaprimeCand= {};

        TLorentzVector proton= {};
        double protonTime= std_ext::NaN;

        TLorentzVector trueProton= {};

        TLorentzVector MM= {};

        double coplanarity= {};
        double EsumCB= {};

        double      taggWeight= {};
        double      taggE= {};
        unsigned    taggCh= {};
        double      taggTime= {};

        double EMB_chi2= std::numeric_limits<double>::infinity();

        double chi2_ref= std::numeric_limits<double>::infinity();
        double prob_ref= {};
        int    iteration_ref= {};
        int    status_ref= {};
        double chi2_sig= std::numeric_limits<double>::infinity();
        double prob_sig= {};
        int    iteration_sig= {};
        int    status_sig= {};

        int type= -1;
        int truetype= -1;

        std::string decayString= {};

        void SetBranches(TTree* tree);
        void FillKinfitBeamProton(double beamE, const TParticlePtr& proton);

    };
    branches vars;

    //======================= histograms ===========================================================
    std::map<std::string,std::map<std::string,TH1*>> hists;
    void AddHist1D(const std::string& category, const std::string& hname,
                   const std::string& title,
                   const std::string& xlabel, const std::string& ylabel,
                   const BinSettings& bins);
    void AddHist2D(const std::string& category, const std::string& hname,
                   const std::string& title,
                   const std::string& xlabel, const std::string& ylabel,
                   const BinSettings& xbins, const BinSettings& ybins);

    // functions
//    void MakeSignal(const TParticleList& photonLeaves);
//    void MakeReference(const TParticleList& photonLeaves);
    bool MakeMCProton(const TEventData& mcdata, TParticlePtr& proton);

    double applyEnergyMomentumConservation(double EBeam, const TParticleList& photons, const TParticlePtr& proton);

public:
    Etap3pi0(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;

};


}}}
