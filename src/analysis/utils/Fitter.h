#pragma once

#include "base/vec/LorentzVec.h"
#include "base/ParticleTypeTree.h"
#include "tree/TParticle.h"
#include "analysis/utils/combinatorics.h"
#include "base/std_ext/math.h"

#include "APLCON.hpp"

#include <stdexcept>
#include <string>
#include <vector>
#include <memory>
#include <set>
#include <functional>

class TTree;
class TH1D;
class TF1;

namespace ant {

class WrapTFile;

namespace analysis {
namespace utils {

class Fitter {

public:

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    /**
     * @brief Uncertainties for E, theta, and phi
     */
    struct Uncertainties_t {
        double sigmaE     = {};
        double sigmaTheta = {};
        double sigmaPhi   = {};

        Uncertainties_t() = default;
        Uncertainties_t(const double E, const double Theta, const double Phi) : sigmaE(E), sigmaTheta(Theta), sigmaPhi(Phi) {}

        bool operator==(const Uncertainties_t& other) const noexcept {
            return sigmaE == other.sigmaE && sigmaTheta == other.sigmaTheta && sigmaPhi == other.sigmaPhi;
        }
    };

    /**
     * @brief Virtual base class for different Uncertainty Models for kion fitter.
     *        Derive and implement the GetSigmas() method.
     * @see UncertaintyModels::Constant
     * @see UncertaintyModels::MCExtracted
     */
    class UncertaintyModel {
    public:
        virtual ~UncertaintyModel();
        virtual Uncertainties_t GetSigmas(const TParticle& particle) const =0;
    };
    using UncertaintyModelPtr = std::shared_ptr<const Fitter::UncertaintyModel>;

    static const APLCON::Fit_Settings_t DefaultSettings;

    Fitter(const Fitter&) = delete;
    Fitter& operator=(const Fitter&) = delete;
    virtual ~Fitter();

    struct FitParticle
    {
        TParticlePtr Particle; // pointer to unfitted value

        struct Var_t {
            double Value = 0;
            double Sigma = 0;
            double Pull = 0;
        protected:
            friend struct FitParticle;
            void SetupBranches(TTree* tree, const std::string& prefix);
        };

        Var_t Ek;
        Var_t Theta;
        Var_t Phi;

        FitParticle(const std::string& name): Name(name) {}


    protected:
        friend class Fitter;
        friend class KinFitter;
        friend class TreeFitter;


        const std::string Name;
        void SetupBranches(TTree* tree, const std::string& prefix);
        static ant::LorentzVec GetVector(const std::vector<double>& EkThetaPhi, double m);

        TParticlePtr AsFitted(const ParticleTypeDatabase::Type& type);

        std::vector<double*> Addresses()
        {
            return { std::addressof(Ek.Value),
                     std::addressof(Theta.Value),
                     std::addressof(Phi.Value)};
        }
        std::vector<double*> Addresses_Sigma()
        {
            return { std::addressof(Ek.Sigma),
                     std::addressof(Theta.Sigma),
                     std::addressof(Phi.Sigma)};
        }

        std::vector<double*> Addresses_Pulls()
        {
            return { std::addressof(Ek.Pull),
                     std::addressof(Theta.Pull),
                     std::addressof(Phi.Pull)};
        }
    };

protected:

    Fitter(const std::string& fittername,
           const APLCON::Fit_Settings_t& settings, std::shared_ptr<const Fitter::UncertaintyModel>& uncertainty_model);

    Fitter(Fitter&&) = default;
    Fitter& operator=(Fitter&&) = default;

    std::shared_ptr<const Fitter::UncertaintyModel> uncertainty;
    std::unique_ptr<APLCON> aplcon;

    void LinkVariable(FitParticle& particle);

    void SetPhotonEkThetaPhi(FitParticle& photon, const TParticlePtr& p) const;

    static double fct_TaggerEGausSigma(double E);

private:
    static APLCON::Fit_Settings_t MakeDefaultSettings();

};


class KinFitter : public Fitter
{
public:

    KinFitter(const std::string& name,
              unsigned numGammas,
              UncertaintyModelPtr Uncertainty_model,
              const APLCON::Fit_Settings_t& settings = DefaultSettings
              );

    virtual ~KinFitter();

    KinFitter(const KinFitter&) = delete;
    KinFitter& operator=(const KinFitter&) = delete;
    KinFitter(KinFitter&&) = default;
    KinFitter& operator=(KinFitter&&) = default;

    void SetEgammaBeam(double ebeam);
    void SetProton(const TParticlePtr& proton);
    virtual void SetPhotons(const TParticleList& photons);

    TParticlePtr GetFittedProton() const;
    TParticleList GetFittedPhotons() const;
    double GetFittedBeamE() const;

    std::vector<FitParticle> GetFitParticles() const;

    void SetupBranches(TTree* tree, std::string branch_prefix="");

    APLCON::Result_t DoFit();

protected:

    struct PhotonBeamVector {
        double E_before = std_ext::NaN;
        double E     = std_ext::NaN;
        double Sigma = std_ext::NaN;

        const std::string Name;

        PhotonBeamVector(const std::string& name);

        std::vector<double*> Adresses()
        {
            return { std::addressof(E) };
        }

        std::vector<double*> Adresses_Sigma()
        {
            return { std::addressof(Sigma)};
        }

    };

    // it's pretty important that those things are pointers,
    // since the members are linked to APLCON in ctor!
    // A move/copy of those members may not happen, so we just
    // point to their fixed location in memory.
    std::vector<std::shared_ptr<FitParticle>> Photons;
    std::unique_ptr<FitParticle> Proton;
    std::unique_ptr<PhotonBeamVector> Beam;

private:
    double result_chi2ndof       =  0.0;
    int result_iterations        =  0;
    int result_status            = -1;
    double result_probability    =  0.0;
};

class TreeFitter : public KinFitter
{

public:
    struct nodesetup_t {
        double IM_Sigma;
        bool Excluded;
        nodesetup_t(double IM_sigma = 1.0, bool excluded = false) :
            IM_Sigma(IM_sigma),
            Excluded(excluded)
        {}

        // the getter lets the user decide how to setup the node
        struct getter : std::function<nodesetup_t(ParticleTypeTree)> {
            // use base class constructors
            using std::function<nodesetup_t(ParticleTypeTree)>::function;
            // provide default which gets the default nodesetup_t
            getter() :
                std::function<nodesetup_t(ParticleTypeTree)>([] (ParticleTypeTree) {return nodesetup_t{};})
            {}
        };
    };

    // construct TreeFitter with KinFit
    TreeFitter(const std::string& name,
               ParticleTypeTree ptree,
               unsigned kinFitGammas,
               UncertaintyModelPtr uncertainty_model,
               nodesetup_t::getter nodeSetup = {},
               const APLCON::Fit_Settings_t& settings = DefaultSettings
              );

    // construct TreeFitter without additional KinFit
    TreeFitter(const std::string& name,
               ParticleTypeTree ptree,
               UncertaintyModelPtr uncertainty_model,
               nodesetup_t::getter nodeSetup = {},
               const APLCON::Fit_Settings_t& settings = DefaultSettings
              ) : TreeFitter(name, ptree, 0, uncertainty_model, nodeSetup, settings)
    {}

    TreeFitter(const TreeFitter&) = delete;
    TreeFitter& operator=(const TreeFitter&) = delete;
    TreeFitter(TreeFitter&&) = default;
    TreeFitter& operator=(TreeFitter&&) = default;

    virtual ~TreeFitter();

    struct node_t {
        node_t(const ParticleTypeTree& ptree) : TypeTree(ptree) {}
        const ParticleTypeTree TypeTree;
        LorentzVec LVSum;
        std::shared_ptr<FitParticle> Leave;
        bool operator<(const node_t& rhs) const {
            return TypeTree->Get() < rhs.TypeTree->Get();
        }
    };

    using tree_t = Tree<node_t>::node_t;

    virtual void SetPhotons(const TParticleList& photons) override;

    tree_t GetTreeNode(const ParticleTypeDatabase::Type& type) const {
        tree_t treenode = nullptr;
        tree->Map_nodes([&treenode, &type] (tree_t t) {
            if(t->Get().TypeTree->Get() == type)
                treenode = t;
        });
        return treenode;
    }

    bool NextFit(APLCON::Result_t& fit_result);

    using current_comb_t = KofNvector<TParticlePtr>;
    const current_comb_t& GetCurrentCombination() const {
        return *current_comb_ptr;
    }



protected:

    // use while(NextFit()) {} instead
    // to run all fits
    using KinFitter::DoFit;

    static tree_t MakeTree(ParticleTypeTree ptree);

    const tree_t tree;
    using permutations_t = std::vector<std::vector<size_t>>;
    permutations_t permutations;
    permutations_t::const_iterator current_perm;

    // use unique_ptr since KofNvector does not have default ctor
    std::unique_ptr<current_comb_t> current_comb_ptr;

};




namespace UncertaintyModels {

/**
 * @brief Constant/static kin fitter uncertainty model. has a fixed value for {cb,taps}{photon,proton}{E,theta,phi}
 */
struct Constant : public Fitter::UncertaintyModel {
public:

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    Fitter::Uncertainties_t photon_cb;
    Fitter::Uncertainties_t photon_taps;
    Fitter::Uncertainties_t proton_cb;
    Fitter::Uncertainties_t proton_taps;

    Constant();
    virtual ~Constant();

    Fitter::Uncertainties_t GetSigmas(const TParticle &particle) const override;

    /**
     * @brief Create a new instance and return a shared pointer to it
     * @return
     */
    static std::shared_ptr<Constant> make();
};

/**
 * @brief Simple kin fitter uncertainty model. has a fixed value for {cb,taps}{photon,proton}{theta,phi}, Energries are relative values and get multipied with the particle energy on GetSigmas()
 */
struct ConstantRelativeE : public Constant {
public:

    ConstantRelativeE();
    virtual ~ConstantRelativeE();

    Fitter::Uncertainties_t GetSigmas(const TParticle &particle) const override;

    /**
     * @brief Create a new, empty instance and return a shared pointer to it
     * @return
     */
    static std::shared_ptr<ConstantRelativeE> make();

    /**
     * @brief Create a new instance filled with global values determined from MC and return a shared ptr to it
     * @return
     */
    static std::shared_ptr<ConstantRelativeE> makeMCLongTarget();
};

struct ConstantRelativeEpow : public Constant {
public:
    double Eexp_cb   = 1.0;
    double Eexp_taps = 1.0;

    ConstantRelativeEpow();
    virtual ~ConstantRelativeEpow();

    Fitter::Uncertainties_t GetSigmas(const TParticle &particle) const override;

    /**
     * @brief Create a new, empty instance and return a shared pointer to it
     * @return
     */
    static std::shared_ptr<ConstantRelativeEpow> make();

    /**
     * @brief Create a new instance filled with global values determined from MC and return a shared ptr to it
     * @return
     */
    static std::shared_ptr<ConstantRelativeEpow> makeMCLongTarget();
};

/**
 * @brief Kin fitter uncertainties, uses histograms. Energy depenent values for each detector element. Histograms can be loaded from root files in setup database.
 */
class MCExtracted : public Fitter::UncertaintyModel {
public:
    struct angular_sigma {
        using Hist = std::shared_ptr<TH1D>;
        Hist p0 = nullptr;
        Hist p1 = nullptr;
        Hist p2 = nullptr;

        double GetSigma(const unsigned element, const double E) const;

        static double f(const double x, const double p0, const double p1, const double p2) noexcept;
        static double f_root(const double* x, const double* p) noexcept;

        static TF1* GetTF1(const std::string& name="SigmaFit");

        void Load(ant::WrapTFile& f, const std::string& prefix, const int bins);
        Hist LoadHist(ant::WrapTFile& f, const std::string& name, const int bins);

        angular_sigma();
        ~angular_sigma();

    };

protected:

    angular_sigma cb_sigma_theta;
    angular_sigma cb_sigma_phi;
    angular_sigma taps_sigma_theta;
    angular_sigma taps_sigma_phi;

    Fitter::Uncertainties_t GetSigmasProton(const TParticle &proton) const;
    Fitter::Uncertainties_t GetSigmasPhoton(const TParticle &photon) const;

public:

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    MCExtracted();
    virtual ~MCExtracted();

    /**
     * @brief Load Sigmas from histograms in ROOT file
     * @param path Path to root file
     */
    void LoadSigmas(const std::string& path);

    Fitter::Uncertainties_t GetSigmas(const TParticle &particle) const override;

    /**
     * @brief Create an instance of this model and directly load sigmas from ROOT file of current setup
     * @return new instance
     */
    static std::shared_ptr<MCExtracted> makeAndLoad();

};

/**
 * @brief Optimized / adjusted to data uncertainty model:
 *
 * Uses fixed/constant values for protons: Energy always 0 -> unmeasured,
 * different values for theta and phi in CB and TAPS.
 *
 * Photons:
 * CB:
 *  simga E = A*E / (E[GeV])**b
 *
 *  sigma Theta = constant + A*sin(theta)
 *
 *  sigma Phi   = sPhi / sin(theta)
 *
 * TAPS:
 *   constant values for now.
 *
 */
struct Optimized : Fitter::UncertaintyModel {

    double cb_photon_theta_const = 0.0;
    double cb_photon_theta_Sin   = 0.0;
    double cb_photon_phi         = 0.0;

    double cb_photon_E_rel = 0.0;
    double cb_photon_E_exp = 0.0;
    double cb_photon_E_lin = 0.0;

    Fitter::Uncertainties_t cb_proton = {};


    double taps_photon_E_rel = 0.0;
    double taps_photon_E_exp = 0.0;
    double taps_photon_E_lin = 0.0;
    double taps_photon_theta = 0.0;
    double taps_photon_phi   = 0.0;

    Fitter::Uncertainties_t taps_proton = {};

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    Optimized();

    virtual ~Optimized();

    Fitter::Uncertainties_t GetSigmas(const TParticle& particle) const;
    static double dThetaSin(const double theta, const double offset, const double thetapart) noexcept;

    static double dE(const double E, const double rel, const double exp, const double reloffset) noexcept;

    std::string to_string_simple() const;

    /**
     * @brief serialize to JSON string
     * @return
     */
    std::string to_string() const;

    /**
     * @brief serialzie to JSON string without whitespace
     * @return
     */
    std::string to_string_short() const;

    /**
     * @brief load from a JSON string
     * @param data
     */
    void load_from_string(const std::string& data);

    void load_from_string_simple(const std::string& data);

    bool operator==(const Optimized& other) const noexcept;
    bool operator!=(const Optimized& other) const noexcept;

protected:
    const static std::string sepatator;

    void ReadToken(const std::string& token);
};

/**
 * @brief Values for the optimized model, obtainded from KinFitPi0 on data
 */
struct Optimized_Oli1 : Optimized {
    Optimized_Oli1();
};

}



}}} // namespace ant::analysis::utils
