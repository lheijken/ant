#include "EtapOmegaG.h"
#include "physics/Plotter.h"

#include "analysis/plot/CutTree.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/EPT.h"

#include "base/interval.h"
#include "base/Logger.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::plot;
using namespace std;

// define the structs containing the histograms
// and the cuts. for simple branch variables, that could
// be combined...

struct CommonHist_t {
    static OptionsPtr opts;

    using Tree_t = physics::EtapOmegaG::TreeCommon;
    using ProtonPhoton_t = physics::EtapOmegaG::ProtonPhotonTree_t;
    struct Fill_t {
        const Tree_t& Common;
        const ProtonPhoton_t& ProtonPhoton;
        const utils::MCWeighting::tree_t& MCWeighting;

        Fill_t(const Tree_t& common, const ProtonPhoton_t& protonphoton,
               const utils::MCWeighting::tree_t& mcWeighting) :
            Common(common), ProtonPhoton(protonphoton), MCWeighting(mcWeighting) {}

        double Weight() const {
            if(MCWeighting.Tree)
                return MCWeighting.MCWeight;
            return Common.TaggW;
        }

    };

    const BinSettings bins_FitProb{100, 0, 1};
    const BinSettings bins_LogFitProb{100, -17, -1};
    TH1D* h_CBSumE = nullptr;
    TH1D* h_CBSumVetoE = nullptr;
    TH1D* h_PIDSumE = nullptr;
    TH2D* h_CBSumVetoE_PIDSumE = nullptr;
    TH1D* h_MissingMass = nullptr;
    TH1D* h_DiscardedEk = nullptr;
    TH1D* h_nTouchesHole = nullptr;
    TH1D* h_MCMissedBkg = nullptr;

    TH2D* h_ProtonTOF = nullptr;
    TH2D* h_ProtonTOFFitted = nullptr;
    TH2D* h_ProtonVetoE = nullptr;
    TH2D* h_ProtonShortE = nullptr;

    const bool isLeaf;
    const bool includeProtonHists;
    const bool moreCutsLessPlots;


    CommonHist_t(HistogramFactory HistFac, cuttree::TreeInfo_t treeInfo) :
        isLeaf(treeInfo.nDaughters==0),
        includeProtonHists(opts->Get<bool>("IncludeProtonHists")),
        moreCutsLessPlots(opts->Get<bool>("MoreCutsLessPlots"))
    {
        if(moreCutsLessPlots)
            return;

        const AxisSettings axis_CBSumVetoE{"CBSumVetoE / MeV", BinSettings(100,0,4)};
        const AxisSettings axis_PIDSumE{"PIDSumE / MeV", BinSettings(50,0,10)};

        h_CBSumE = HistFac.makeTH1D("CB Sum E","E / MeV","",BinSettings(100,500,1600),"h_CBSumE");
        h_CBSumVetoE = HistFac.makeTH1D("CB Sum VetoE",axis_CBSumVetoE,"h_CBSumVetoE");
        h_PIDSumE = HistFac.makeTH1D("PID Sum E",axis_PIDSumE,"h_PIDSumE");
        h_CBSumVetoE_PIDSumE = HistFac.makeTH2D("PIDSumE vs. CBSumVetoE",axis_PIDSumE,axis_CBSumVetoE,"h_CBSumVetoE_PIDSumE");
        h_MissingMass = HistFac.makeTH1D("MissingMass","m / MeV","",BinSettings(200,600,1300),"h_MissingMass");
        h_DiscardedEk = HistFac.makeTH1D("DiscardedEk","E / MeV","",BinSettings(100,0,100),"h_DiscardedEk");
        h_nTouchesHole = HistFac.makeTH1D("nTouchesHole","nTouchesHole","",BinSettings(5),"h_nTouchesHole");
        h_MCMissedBkg = HistFac.makeTH1D("MCMissedBkg","","",BinSettings(15),"h_MCMissedBkg");

        if(!isLeaf)
            return;

        if(includeProtonHists) {
            BinSettings bins_protonE(100,0,600);
            h_ProtonTOF = HistFac.makeTH2D("ProtonTOF","t","E",
                                           BinSettings(50,-5,20),bins_protonE,"h_ProtonTOF");
            h_ProtonTOFFitted = HistFac.makeTH2D("ProtonTOFFitted","t","E fitted",
                                                 BinSettings(50,-5,20),bins_protonE,"h_ProtonTOFFitted");
            h_ProtonVetoE = HistFac.makeTH2D("ProtonVeto","E fitted","Veto E",
                                             bins_protonE,BinSettings(50,0,8),"h_ProtonVeto");
            h_ProtonShortE = HistFac.makeTH2D("ProtonShortE","E fitted","E short",
                                              bins_protonE,BinSettings(100,0,300),"h_ProtonShortE");
        }
    }


    void Fill(const Fill_t& f) const {
        if(moreCutsLessPlots)
            return;

        h_CBSumE->Fill(f.Common.CBSumE, f.Weight());
        h_CBSumVetoE->Fill(f.ProtonPhoton.CBSumVetoE, f.Weight());
        h_PIDSumE->Fill(f.Common.PIDSumE, f.Weight());
        h_CBSumVetoE_PIDSumE->Fill(f.Common.PIDSumE, f.ProtonPhoton.CBSumVetoE, f.Weight());
        h_MissingMass->Fill(f.ProtonPhoton.MissingMass, f.Weight());
        h_DiscardedEk->Fill(f.ProtonPhoton.DiscardedEk, f.Weight());
        h_nTouchesHole->Fill(f.ProtonPhoton.nTouchesHole, f.Weight());

        if(!f.Common.MCTrueMissed().empty())
            h_MCMissedBkg->Fill(f.Common.MCTrueMissed().c_str(), f.Weight());

        if(!isLeaf)
            return;

        if(includeProtonHists) {
            h_ProtonTOF->Fill(f.ProtonPhoton.Proton().Time, f.ProtonPhoton.Proton().Ek(), f.Weight());
            h_ProtonTOFFitted->Fill(f.ProtonPhoton.Proton().Time, f.ProtonPhoton.FittedProtonE, f.Weight());
            h_ProtonVetoE->Fill(f.ProtonPhoton.FittedProtonE, f.ProtonPhoton.Proton().VetoE, f.Weight());
            h_ProtonShortE->Fill(f.ProtonPhoton.FittedProtonE, f.ProtonPhoton.Proton().ShortE, f.Weight());
        }
    }

    std::vector<TH1*> GetHists() const {
        return {h_CBSumE, h_CBSumVetoE, h_PIDSumE, h_MissingMass,
                    h_DiscardedEk, h_nTouchesHole};
    }

    // Sig and Ref channel (may) share some cuts...
    static cuttree::Cuts_t<Fill_t> GetCuts() {
        using cuttree::MultiCut_t;
        cuttree::Cuts_t<Fill_t> cuts;
        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"DiscardedEk=0", [] (const Fill_t& f) { return f.ProtonPhoton.DiscardedEk == 0; } },
                              {"DiscardedEk<20", [] (const Fill_t& f) { return f.ProtonPhoton.DiscardedEk < 20; } },
                              {"DiscardedEk<50", [] (const Fill_t& f) { return f.ProtonPhoton.DiscardedEk < 50; } },
                          });
        return cuts;
    }
};

struct SigHist_t : CommonHist_t {
    using SharedTree_t = physics::EtapOmegaG::Sig_t::SharedTree_t;
    using Tree_t = physics::EtapOmegaG::Sig_t::Fit_t::BaseTree_t;

    struct Fill_t : CommonHist_t::Fill_t {
        const SharedTree_t& Shared;
        const Tree_t& Tree;
        Fill_t(const CommonHist_t::Tree_t& common,
               const SharedTree_t& shared,
               const Tree_t& tree,
               const utils::MCWeighting::tree_t& mcWeighting) :
            CommonHist_t::Fill_t(common, tree, mcWeighting),
            Shared(shared),
            Tree(tree)
        {}

    };

    TH1D* h_IM_4g = nullptr;        // EtaPrime IM
    TH2D* h_IM_4g_TaggCh = nullptr;

    TH1D* h_KinFitProb = nullptr;
    TH1D* h_AntiPi0FitProb = nullptr;
    TH1D* h_AntiEtaFitProb = nullptr;
    TH1D* h_TreeFitProb = nullptr;

    TH1D* h_AntiPi0ZVertex = nullptr;
    TH1D* h_AntiEtaZVertex = nullptr;
    TH1D* h_TreeZVertex = nullptr;

    TH2D* h_gNonPi0_CaloE_Theta = nullptr;
    TH1D* h_gNonPi0_TouchesHoles = nullptr;
    TH1D* h_gNonPi0_CBSumVetoE = nullptr;

    TH1D* h_Bachelor_E = nullptr;

    const BinSettings bins_IM_Etap {1020-910, 910, 1020};
    const BinSettings bins_IM_Omega{100, 700, 900};
    const BinSettings bins_ZVertex{100, -15, 15};

    SigHist_t(HistogramFactory HistFac, cuttree::TreeInfo_t treeInfo)
        : CommonHist_t(HistFac, treeInfo)
    {
        if(moreCutsLessPlots && !isLeaf)
            return;

        h_IM_4g = HistFac.makeTH1D("#eta' IM", "IM(#pi^{0}#gamma#gamma) / MeV","",bins_IM_Etap,"h_IM_4g");

        if(moreCutsLessPlots)
            return;

        auto ept = ExpConfig::Setup::GetDetector<expconfig::detector::EPT>();
        h_IM_4g_TaggCh = HistFac.makeTH2D("IM 4g vs. TaggCh","IM(#pi^{0}#gamma#gamma) / MeV","Tagger Channel",
                                          bins_IM_Etap, BinSettings(ept->GetNChannels()),
                                          "h_IM_4g_TaggCh",true);

        h_KinFitProb = HistFac.makeTH1D("KinFitProb","p","",bins_FitProb,"h_KinFitProb");
        h_AntiPi0FitProb = HistFac.makeTH1D("AntiPi0FitProb", "log_{10} p","",bins_LogFitProb,"h_AntiPi0FitProb");
        h_AntiEtaFitProb = HistFac.makeTH1D("AntiEtaFitProb", "log_{10} p","",bins_LogFitProb,"h_AntiEtaFitProb");
        h_TreeFitProb = HistFac.makeTH1D("TreeFitProb", "p","",bins_FitProb,"h_TreeFitProb");

        h_AntiPi0ZVertex = HistFac.makeTH1D("AntiPi0ZVertex", "z / cm","",bins_ZVertex,"h_AntiPi0ZVertex");
        h_AntiEtaZVertex = HistFac.makeTH1D("AntiEtaZVertex", "z / cm","",bins_ZVertex,"h_AntiEtaZVertex");
        h_TreeZVertex = HistFac.makeTH1D("TreeZVertex", "z / cm","",bins_ZVertex,"h_TreeZVertex");

        h_gNonPi0_CaloE_Theta = HistFac.makeTH2D("gNonPi0 E_{k} #theta","#theta / #circ","E_{k} / MeV",
                                                 BinSettings(180,0,180), BinSettings(50,0,400), "h_gNonPi0_CaloE_Theta");

        h_gNonPi0_TouchesHoles = HistFac.makeTH1D("gNonPi0_TouchesHole","nTouchesHole","",
                                                  BinSettings(3),"h_gNonPi0_TouchesHole");
        h_gNonPi0_CBSumVetoE = HistFac.makeTH1D("CB Veto Sum E Bachelor Photons","E / MeV","",
                                                BinSettings(100,0,2),"h_gNonPi0_CBSumVetoE");

        BinSettings bins_BachelorE(100,100,200);
        h_Bachelor_E = HistFac.makeTH1D("E_#gamma in #eta' frame","E_{#gamma} / MeV","",
                                        bins_BachelorE,"h_Bachelor_E");


    }

    void Fill(const Fill_t& f) const {
        CommonHist_t::Fill(f);
        const SharedTree_t& s = f.Shared;
        const Tree_t& tree = f.Tree;

        if(h_IM_4g)
            h_IM_4g->Fill(tree.IM_Pi0gg, f.Weight());

        if(moreCutsLessPlots)
            return;

        h_IM_4g_TaggCh->Fill(tree.IM_Pi0gg, f.Common.TaggCh, f.Weight());

        h_KinFitProb->Fill(s.KinFitProb, f.Weight());
        const auto get_log_prob = [this] (double prob) {
            const auto log_prob = std::log10(prob);
            return isfinite(log_prob) ? log_prob : bins_LogFitProb.Start();
        };
        h_AntiPi0FitProb->Fill(get_log_prob(s.AntiPi0FitProb), f.Weight());
        h_AntiEtaFitProb->Fill(get_log_prob(s.AntiEtaFitProb), f.Weight());
        h_TreeFitProb->Fill(tree.TreeFitProb, f.Weight());

        h_AntiPi0ZVertex->Fill(s.AntiPi0FitZVertex, f.Weight());
        h_AntiEtaZVertex->Fill(s.AntiEtaFitZVertex, f.Weight());
        h_TreeZVertex->Fill(tree.TreeFitZVertex, f.Weight());

        {
            const auto& theta0 = std_ext::radian_to_degree(tree.gNonPi0()[0].Theta());
            h_gNonPi0_CaloE_Theta->Fill(theta0, tree.gNonPi0()[0].Ek(), f.Weight());\
            const auto& theta1 = std_ext::radian_to_degree(tree.gNonPi0()[1].Theta());
            h_gNonPi0_CaloE_Theta->Fill(theta1, tree.gNonPi0()[1].Ek(), f.Weight());
        }

        h_gNonPi0_TouchesHoles->Fill(tree.gNonPi0()[0].TouchesHole+tree.gNonPi0()[1].TouchesHole, f.Weight());
        h_gNonPi0_CBSumVetoE->Fill(tree.gNonPi0()[0].VetoE+tree.gNonPi0()[1].VetoE, f.Weight());
    }

    std::vector<TH1*> GetHists() const {
        auto hists = CommonHist_t::GetHists();
        hists.insert(hists.end(), {
                         h_IM_4g, h_KinFitProb,
                         h_AntiPi0FitProb, h_AntiEtaFitProb, h_TreeFitProb,
                         h_AntiPi0ZVertex, h_AntiEtaZVertex, h_TreeZVertex,
                         h_gNonPi0_TouchesHoles, h_gNonPi0_CBSumVetoE,
                         h_Bachelor_E
                     });
        return hists;
    }

    static cuttree::Cuts_t<Fill_t> GetCuts() {
        const auto moreCutsLessPlots = opts->Get<bool>("MoreCutsLessPlots");

        using cuttree::MultiCut_t;
        auto cuts = cuttree::ConvertCuts<Fill_t, CommonHist_t::Fill_t>(CommonHist_t::GetCuts());

        // anti tree fits
        {
            const auto cut = [] (double logcutval, double prob) {
                return std::isnan(prob) || std::log10(prob)<logcutval;
            };

            // anti pi0
            MultiCut_t<Fill_t> cut_AntiPi0{
                {"AntiPi0FitProb<10^{-5}||nan",  [cut] (const Fill_t& f) { return cut(-5, f.Shared.AntiPi0FitProb); } },
                {"AntiPi0FitProb<10^{-7}||nan",  [cut] (const Fill_t& f) { return cut(-7, f.Shared.AntiPi0FitProb); } },
            };
            if(moreCutsLessPlots)
                cut_AntiPi0.emplace_back("AntiPi0FitProb<10^{-3}||nan",  [cut] (const Fill_t& f) { return cut(-3, f.Shared.AntiPi0FitProb); });
            cuts.emplace_back(cut_AntiPi0);

            // anti eta
            MultiCut_t<Fill_t> cut_AntiEta{
                {"AntiEtaFitProb<10^{-4}||nan", [cut] (const Fill_t& f) { return cut(-4, f.Shared.AntiEtaFitProb); } },
                {"AntiEtaFitProb<10^{-6}||nan", [cut] (const Fill_t& f) { return cut(-6, f.Shared.AntiEtaFitProb); } },
            };
            if(moreCutsLessPlots)
                cut_AntiEta.emplace_back("AntiEtaFitProb<10^{-2}||nan", [cut] (const Fill_t& f) { return cut(-2, f.Shared.AntiEtaFitProb); });
            cuts.emplace_back(cut_AntiEta);
        }

        // tree fit cut
        {
            MultiCut_t<Fill_t> cut_TreeFit{
                {"TreeFitProb>0.2", [] (const Fill_t& f) { return f.Tree.TreeFitProb>0.2; } },
                {"TreeFitProb>0.1", [] (const Fill_t& f) { return f.Tree.TreeFitProb>0.1; } },
            };
            if(moreCutsLessPlots)
                cut_TreeFit.emplace_back("TreeFitProb>0.05", [] (const Fill_t& f) { return f.Tree.TreeFitProb>0.05; });

            cuts.emplace_back(cut_TreeFit);
        }

        // kinematic bachelor photon cut
        {
            auto gNonPi0_cut_1 = [] (const Fill_t& f) {
                const auto cut = [] (const TSimpleParticle& p) {
                    const auto& theta = std_ext::radian_to_degree(p.Theta());
                    const auto& caloE = p.Ek();
                    return caloE > 230.0*(1.0-theta/160.0);
                };
                return cut(f.Tree.gNonPi0()[0]) && cut(f.Tree.gNonPi0()[1]);
            };

            const auto cut_simple = [] (const TSimpleParticle& p, double factor = 1.0) {
                const auto& theta = std_ext::radian_to_degree(p.Theta());
                const auto& caloE = p.Ek();
                if(theta<22) // decide if TAPS or CB
                    return caloE > factor*140;
                else
                    return caloE > factor*60;
            };

            auto gNonPi0_cut_2 = [cut_simple] (const Fill_t& f) {
                return cut_simple(f.Tree.gNonPi0()[0]) && cut_simple(f.Tree.gNonPi0()[1]);
            };

            auto gNonPi0_cut_3 = [cut_simple] (const Fill_t& f) {
                return cut_simple(f.Tree.gNonPi0()[0], 0.5) && cut_simple(f.Tree.gNonPi0()[1], 0.5);
            };

            MultiCut_t<Fill_t> cut_gNonPi0{
                {"gNonPi0_1", gNonPi0_cut_1},
                {"gNonPi0_2", gNonPi0_cut_2},
                {"gNonPi0_3", gNonPi0_cut_3},
            };

            cuts.emplace_back(cut_gNonPi0);
        }

        // PID cuts
        {
            auto cut_gNonPi0 = [] (const Fill_t& f, double cut) {
                auto& v = f.Tree.gNonPi0();
                return (v.front().VetoE + v.back().VetoE)<=cut;
            };

            MultiCut_t<Fill_t> pid_cut{
                {"CBSumVetoE_gNonPi0<0.2", [cut_gNonPi0] (const Fill_t& f) { return cut_gNonPi0(f, 0.2); }},
            };

            if(moreCutsLessPlots) {
                // concentrate on gNonPi0 VetoE
                pid_cut.emplace_back("CBSumVetoE_gNonPi0=0", [cut_gNonPi0] (const Fill_t& f) { return cut_gNonPi0(f, 0.0); });
                pid_cut.emplace_back("CBSumVetoE_gNonPi0<0.1", [cut_gNonPi0] (const Fill_t& f) { return cut_gNonPi0(f, 0.1); });
                pid_cut.emplace_back("CBSumVetoE_gNonPi0<0.4", [cut_gNonPi0] (const Fill_t& f) { return cut_gNonPi0(f, 0.4); });
            }
            else {
                // some standard test cuts
                pid_cut.emplace_back("NoCBSumVetoE", [] (const Fill_t&) { return true; });
                pid_cut.emplace_back("CBSumVetoE<0.2", [] (const Fill_t& f) { return f.ProtonPhoton.CBSumVetoE<0.2; });
            }


            cuts.emplace_back(pid_cut);
        }
        return cuts;
    }
};

struct SigPi0Hist_t : SigHist_t {
    using Tree_t = physics::EtapOmegaG::Sig_t::Pi0_t::Tree_t;

    TH2D* h_IM_3g_4g_high = nullptr;    // Omega IM vs. EtaPrime IM

    struct Fill_t : SigHist_t::Fill_t {
        const Tree_t& Pi0;
        Fill_t(const CommonHist_t::Tree_t& common,
               const SharedTree_t& shared,
               const Tree_t& pi0,
               const utils::MCWeighting::tree_t& mcWeighting) :
            SigHist_t::Fill_t(common, shared, pi0, mcWeighting),
            Pi0(pi0)
        {}
    };

    SigPi0Hist_t(HistogramFactory HistFac, cuttree::TreeInfo_t treeInfo) : SigHist_t(HistFac, treeInfo) {
        if(moreCutsLessPlots)
            return;

        h_IM_3g_4g_high = HistFac.makeTH2D("#omega vs. #eta' IM",
                                           "IM(#pi^{0}#gamma#gamma) / MeV",
                                           "IM(#pi^{0}#gamma) / MeV",
                                           bins_IM_Etap, bins_IM_Omega,"h_IM_3g_4g_high"
                                           );
    }

    void Fill(const Fill_t& f) const {
        SigHist_t::Fill(f);
        if(moreCutsLessPlots)
            return;
        const Tree_t& pi0 = f.Pi0;
        h_IM_3g_4g_high->Fill(pi0.IM_Pi0gg, pi0.IM_Pi0g()[1], f.Weight());
        h_Bachelor_E->Fill(pi0.Bachelor_E()[0], f.Weight());
    }

    static cuttree::Cuts_t<Fill_t> GetCuts() {
        using cuttree::MultiCut_t;
        auto cuts = cuttree::ConvertCuts<Fill_t, SigHist_t::Fill_t>(SigHist_t::GetCuts());
        const auto moreCutsLessPlots = opts->Get<bool>("MoreCutsLessPlots");
        if(moreCutsLessPlots) {
            cuts.emplace_back(MultiCut_t<Fill_t>{
                                  {"IM_Pi0g[1] 40", [] (const Fill_t& f) {
                                       const auto& window = ParticleTypeDatabase::Omega.GetWindow(40);
                                       return window.Contains(f.Pi0.IM_Pi0g()[1]);
                                   }},
                                  {"IM_Pi0g[1] 30", [] (const Fill_t& f) {
                                       const auto& window = ParticleTypeDatabase::Omega.GetWindow(30);
                                       return window.Contains(f.Pi0.IM_Pi0g()[1]);
                                   }},
                                  {"IM_Pi0g[1] 50", [] (const Fill_t& f) {
                                       const auto& window = ParticleTypeDatabase::Omega.GetWindow(50);
                                       return window.Contains(f.Pi0.IM_Pi0g()[1]);
                                   }},
                              });
        }
        else {
            cuts.emplace_back(MultiCut_t<Fill_t>{
                                  {"IM_Pi0g[1]", [] (const Fill_t& f) {
                                       const auto& window = ParticleTypeDatabase::Omega.GetWindow(40);
                                       return window.Contains(f.Pi0.IM_Pi0g()[1]);
                                   }},
                              });
        }
        return cuts;
    }

};

struct SigOmegaPi0Hist_t : SigHist_t {
    using Tree_t = physics::EtapOmegaG::Sig_t::OmegaPi0_t::Tree_t;
    struct Fill_t : SigHist_t::Fill_t {
        const Tree_t& OmegaPi0;
        Fill_t(const CommonHist_t::Tree_t& common,
               const SharedTree_t& shared,
               const Tree_t& omegapi0,
               const utils::MCWeighting::tree_t& mcWeighting) :
            SigHist_t::Fill_t(common, shared, omegapi0, mcWeighting),
            OmegaPi0(omegapi0)
        {}
    };




    SigOmegaPi0Hist_t(HistogramFactory HistFac, cuttree::TreeInfo_t treeInfo) : SigHist_t(HistFac, treeInfo) {

    }

    std::vector<TH1*> GetHists() const {
        return SigHist_t::GetHists();
    }

    void Fill(const Fill_t& f) const {
        if(moreCutsLessPlots)
            return;

        SigHist_t::Fill(f);
        const Tree_t& omegapi0 = f.OmegaPi0;
        h_Bachelor_E->Fill(omegapi0.Bachelor_E, f.Weight());
    }

    static cuttree::Cuts_t<Fill_t> GetCuts() {
        return cuttree::ConvertCuts<Fill_t, SigHist_t::Fill_t>(SigHist_t::GetCuts());
    }

};

struct RefHist_t : CommonHist_t {
    using Tree_t = physics::EtapOmegaG::Ref_t::Tree_t;

    struct Fill_t : CommonHist_t::Fill_t {
        const Tree_t& Tree;
        Fill_t(const CommonHist_t::Tree_t& common, const Tree_t& tree,
               const utils::MCWeighting::tree_t& mcWeighting) :
            CommonHist_t::Fill_t(common, tree, mcWeighting),
            Tree(tree)
        {}

    };


    TH1D* h_IM_2g = nullptr;
    TH2D* h_IM_2g_TaggCh = nullptr;

    TH1D* h_KinFitProb = nullptr;
    TH1D* h_TaggT = nullptr;

    RefHist_t(HistogramFactory HistFac, cuttree::TreeInfo_t treeInfo) : CommonHist_t(HistFac, treeInfo) {
        BinSettings bins_im(1050-875,875,1050);

        h_KinFitProb = HistFac.makeTH1D("KinFitProb","p","",bins_FitProb,"h_KinFitProb");
        h_IM_2g = HistFac.makeTH1D("IM 2g","IM / MeV","",bins_im,"h_IM_2g");

        auto ept = ExpConfig::Setup::GetDetector<expconfig::detector::EPT>();
        h_IM_2g_TaggCh = HistFac.makeTH2D("IM 2g vs. TaggCh","IM / MeV","Tagger Channel",
                                          bins_im, BinSettings(ept->GetNChannels()), "h_IM_2g_TaggCh", true);

        h_TaggT = HistFac.makeTH1D("Tagger Time","t / ns","",BinSettings(400,-50,50),"h_TaggT");
    }

    void Fill(const Fill_t& f) const {
        CommonHist_t::Fill(f);
        const Tree_t& tree = f.Tree;

        h_KinFitProb->Fill(tree.KinFitProb, f.Weight());

        h_IM_2g->Fill(tree.IM_2g, f.Weight());
        h_IM_2g_TaggCh->Fill(tree.IM_2g, f.Common.TaggCh, f.Weight());

        h_TaggT->Fill(f.Common.TaggT-f.Common.CBAvgTime, f.Weight()<0 ? -1.0 : 1.0);
    }

    std::vector<TH1*> GetHists() const {
        auto hists = CommonHist_t::GetHists();
        hists.insert(hists.end(), {h_KinFitProb, h_IM_2g});
        return hists;
    }

    static cuttree::Cuts_t<Fill_t> GetCuts() {
        using cuttree::MultiCut_t;
        auto cuts = cuttree::ConvertCuts<Fill_t, CommonHist_t::Fill_t>(CommonHist_t::GetCuts());

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"KinFitProb>0.02", [] (const Fill_t& f) { return f.Tree.KinFitProb>0.02; } },
                              {"KinFitProb>0.05", [] (const Fill_t& f) { return f.Tree.KinFitProb>0.05; } },
                              {"KinFitProb>0.1", [] (const Fill_t& f) { return f.Tree.KinFitProb>0.1; } },
                              {"KinFitProb>0.2", [] (const Fill_t& f) { return f.Tree.KinFitProb>0.2; } },
                          });

        cuts.emplace_back(MultiCut_t<Fill_t>{
                              {"PIDSumE=0", [] (const Fill_t& f) { return f.Common.PIDSumE==0; }},
                              {"PIDSumE<0.2", [] (const Fill_t& f) { return f.Common.PIDSumE<0.2; }},
                              {"CBSumVetoE=0", [] (const Fill_t& f) { return f.ProtonPhoton.CBSumVetoE==0; }},
                              {"CBSumVetoE<0.2", [] (const Fill_t& f) { return f.ProtonPhoton.CBSumVetoE<0.2; }},
                          });
        return cuts;
    }
};

OptionsPtr CommonHist_t::opts;

template<typename Hist_t>
struct MCTrue_Splitter : cuttree::StackedHists_t<Hist_t> {

    const bool moreCutsLessPlots;

    // Hist_t should have that type defined
    using Fill_t = typename Hist_t::Fill_t;

    MCTrue_Splitter(const HistogramFactory& histFac,
                    const cuttree::TreeInfo_t& treeInfo) :
        cuttree::StackedHists_t<Hist_t>(histFac, treeInfo),
        moreCutsLessPlots(CommonHist_t::opts->Get<bool>("MoreCutsLessPlots"))
    {
        using histstyle::Mod_t;
        this->GetHist(0, "Data", Mod_t::MakeDataPoints(kBlack));
        if(!moreCutsLessPlots) {
            this->GetHist(5, "D07", Mod_t::MakeDataPoints(kGray));
            this->GetHist(6, "D10", Mod_t::MakeDataPoints(kGray));
            this->GetHist(7, "D12", Mod_t::MakeDataPoints(kGray));
        }

        this->GetHist(1, "Sig",  Mod_t::MakeLine(kRed, 2.0));

        if(moreCutsLessPlots)
            return;

        this->GetHist(2, "Ref",  Mod_t::MakeLine(kRed, 2.0));
        // mctrue is never >=3 (and <9) in tree, use this to sum up all MC and all bkg MC
        // see also Fill()
        this->GetHist(3, "Sum_MC", Mod_t::MakeLine(kBlack, 2.0));
        this->GetHist(4, "Bkg_MC", Mod_t::MakeFill(kGray+2));
    }

    void Fill(const Fill_t& f) {

        const unsigned mctrue = f.Common.MCTrue;

        auto get_bkg_name = [] (unsigned mctrue) {
            const string& name = mctrue>=10 ?
                        physics::EtapOmegaG::ptreeBackgrounds[mctrue-10].Name
                    : "Other";
            return "Bkg_"+name;
        };

        using histstyle::Mod_t;
        const Hist_t& hist = mctrue<9 ? this->GetHist(mctrue) :
                                        this->GetHist(mctrue,
                                                      get_bkg_name(mctrue),
                                                      Mod_t::MakeLine(histstyle::color_t::GetLight(mctrue-9), 1, kGray+2)
                                                      );

        hist.Fill(f);

        if(moreCutsLessPlots)
            return;

        // handle MC_all and MC_bkg
        if(mctrue>0) {
            this->GetHist(3).Fill(f);
            if(mctrue >= 9)
                this->GetHist(4).Fill(f);
        }

        // handle D07/D10/D12
        if(mctrue == 0) {
            this->GetHist(4+f.Common.BeamTime).Fill(f);
        }
    }
};

struct EtapOmegaG_plot : Plotter {

    CommonHist_t::Tree_t treeCommon;
    utils::MCWeighting::tree_t treeMCWeighting;

    EtapOmegaG_plot(const string& tag, const string& name, const WrapTFileInput& input, OptionsPtr opts) :
        Plotter(name, input, opts)
    {
        /// \todo using a static field is actually quite ugly...
        CommonHist_t::opts = opts;

        init_tree(input, treeCommon, "EtapOmegaG/"+tag+"/Common");

        if(input.GetObject("EtapOmegaG/"+tag+"/"+utils::MCWeighting::treeName, treeMCWeighting.Tree)) {
            LOG(INFO) << "Found " << tag << " MCWeighting tree";
            treeMCWeighting.LinkBranches();
            check_entries(treeMCWeighting);
        }

        create_MCTrue_generated(tag, HistFac, input);
    }

    static void create_MCTrue_generated(const string& tag, const HistogramFactory& HistFac, const WrapTFileInput& input) {
        utils::MCWeighting::tree_t treeMCWeighting;
        physics::EtapOmegaG::TreeMCWeighting treeMCWeighting_extra;

        if(input.GetObject("EtapOmegaG/"+utils::MCWeighting::treeName, treeMCWeighting.Tree) &&
           input.GetObject("EtapOmegaG/"+utils::MCWeighting::treeName+"_extra", treeMCWeighting_extra.Tree)) {
            treeMCWeighting.LinkBranches();
            treeMCWeighting_extra.LinkBranches();
            if(treeMCWeighting.Tree->GetEntries() != treeMCWeighting.Tree->GetEntries()) {
                LOG(ERROR) << "Mismatch in trees for MCTrue generated hist";
                return;
            }
            LOG(INFO) << "Creating generated histogram from MCWeighting trees...";
            auto ept = ExpConfig::Setup::GetDetector<expconfig::detector::EPT>();

            auto h = HistFac.makeTH1D(tag+": Generated events",{"TaggCh",{ept->GetNChannels()}},"h_mctrue_generated",true);
            for(long long entry=0;entry<treeMCWeighting.Tree->GetEntries();entry++) {
                treeMCWeighting_extra.Tree->GetEntry(entry);
                if(tag == "Sig" && treeMCWeighting_extra.MCTrue != 1)
                    continue;
                if(tag == "Ref" && treeMCWeighting_extra.MCTrue != 2)
                    continue;
                treeMCWeighting.Tree->GetEntry(entry);
                h->Fill(treeMCWeighting_extra.TaggCh(), treeMCWeighting.MCWeight());
            }
        }

    }

    static void init_tree(const WrapTFileInput& input, WrapTTree& tree, const string& name) {
        if(!input.GetObject(name, tree.Tree))
            throw Exception("Cannot find tree "+name);
        tree.LinkBranches();
    }

    void check_entries(const WrapTTree& tree) {
        if(tree.Tree->GetEntries() == treeCommon.Tree->GetEntries())
            return;
        throw Exception(std_ext::formatter() << "Tree " << tree.Tree->GetName()
                        << " does not have expected entries=" << treeCommon.Tree->GetEntries()
                        << " but " << tree.Tree->GetEntries());
    }

    virtual long long GetNumEntries() const override
    {
        return treeCommon.Tree->GetEntries();
    }

    virtual void ProcessEntry(const long long entry) override
    {
        treeCommon.Tree->GetEntry(entry);
        if(treeMCWeighting.Tree)
            treeMCWeighting.Tree->GetEntry(entry);
    }

};

struct EtapOmegaG_plot_Ref : EtapOmegaG_plot {

    RefHist_t::Tree_t          treeRef;

    using MCRefHist_t = MCTrue_Splitter<RefHist_t>;
    cuttree::Tree_t<MCRefHist_t> cuttreeRef;

    EtapOmegaG_plot_Ref(const string& name, const WrapTFileInput& input, OptionsPtr opts) :
        EtapOmegaG_plot("Ref", name, input, opts)
    {
        init_tree(input, treeRef, "EtapOmegaG/Ref/Ref");
        check_entries(treeRef);

        cuttreeRef = cuttree::Make<MCRefHist_t>(HistFac);
    }

    virtual void ProcessEntry(const long long entry) override
    {
        EtapOmegaG_plot::ProcessEntry(entry);
        treeRef.Tree->GetEntry(entry);
        cuttree::Fill<MCRefHist_t>(cuttreeRef, {treeCommon, treeRef, treeMCWeighting});
    }
};

struct EtapOmegaG_plot_Sig : EtapOmegaG_plot {

    SigHist_t::SharedTree_t   treeSigShared;
    SigPi0Hist_t::Tree_t      treeSigPi0;
    SigOmegaPi0Hist_t::Tree_t treeSigOmegaPi0;

    using MCSigPi0Hist_t = MCTrue_Splitter<SigPi0Hist_t>;
    using MCSigOmegaPi0Hist_t = MCTrue_Splitter<SigOmegaPi0Hist_t>;
    cuttree::Tree_t<MCSigPi0Hist_t>      cuttreeSigPi0;
    cuttree::Tree_t<MCSigOmegaPi0Hist_t> cuttreeSigOmegaPi0;

    EtapOmegaG_plot_Sig(const string& name, const WrapTFileInput& input, OptionsPtr opts) :
        EtapOmegaG_plot("Sig", name, input, opts)
    {
        init_tree(input, treeSigShared,   "EtapOmegaG/Sig/Shared");
        init_tree(input, treeSigPi0,      "EtapOmegaG/Sig/Pi0");
        init_tree(input, treeSigOmegaPi0, "EtapOmegaG/Sig/OmegaPi0");

        check_entries(treeSigShared);
        check_entries(treeSigPi0);
        check_entries(treeSigOmegaPi0);

        cuttreeSigPi0 = cuttree::Make<MCSigPi0Hist_t>(HistogramFactory("SigPi0",HistFac,"SigPi0"));
        if(!CommonHist_t::opts->Get<bool>("MoreCutsLessPlots"))
            cuttreeSigOmegaPi0 = cuttree::Make<MCSigOmegaPi0Hist_t>(HistogramFactory("SigOmegaPi0",HistFac,"SigOmegaPi0"));
    }

    virtual void ProcessEntry(const long long entry) override
    {
        EtapOmegaG_plot::ProcessEntry(entry);
        treeSigShared.Tree->GetEntry(entry);
        treeSigPi0.Tree->GetEntry(entry);
        treeSigOmegaPi0.Tree->GetEntry(entry);
        cuttree::Fill<MCSigPi0Hist_t>(cuttreeSigPi0, {treeCommon, treeSigShared, treeSigPi0, treeMCWeighting});
        if(cuttreeSigOmegaPi0)
            cuttree::Fill<MCSigOmegaPi0Hist_t>(cuttreeSigOmegaPi0, {treeCommon, treeSigShared, treeSigOmegaPi0, treeMCWeighting});
    }
};

AUTO_REGISTER_PLOTTER(EtapOmegaG_plot_Ref)
AUTO_REGISTER_PLOTTER(EtapOmegaG_plot_Sig)
