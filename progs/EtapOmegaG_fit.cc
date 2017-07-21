#include "base/Logger.h"

#include "tclap/CmdLine.h"
#include "tclap/ValuesConstraintExtra.h"
#include "base/interval.h"
#include "base/WrapTFile.h"
#include "base/std_ext/string.h"
#include "base/std_ext/system.h"
#include "base/std_ext/memory.h"
#include "base/std_ext/math.h"
#include "base/std_ext/container.h"
#include "base/ParticleType.h"

#include "analysis/plot/RootDraw.h"
#include "analysis/plot/HistogramFactory.h"
#include "analysis/utils/ParticleTools.h"
#include "root-addons/analysis_codes/Math.h"
#include "expconfig/ExpConfig.h"
#include "base/Detector_t.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TVectorT.h"

#include "TSystem.h"
#include "TRint.h"
#include "TROOT.h"
#include "TFrame.h"
#include "TLegend.h"

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooChebychev.h"
#include "RooConstVar.h"
#include "RooDerivative.h"
#include "RooFFTConvPdf.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "RooHist.h"
#include "RooPlotable.h"

#include "APLCON.hpp"

#include <fstream>

auto debug = false;

using namespace ant;
using namespace std;
using namespace RooFit;

// use APLCON to calculate the total sum with error propagation
// Value, Sigma, Pull
struct N_t {
    explicit N_t(const RooRealVar& var) : Value(var.getVal()), Sigma(var.getError()) {}
    explicit N_t(double v=0, double s=0) : Value(v), Sigma(s) {}

    double Value;
    double Sigma;
    double Pull = std_ext::NaN;

    template<std::size_t N>
    std::tuple<double&> linkFitter() noexcept {
        // the following get<N> assumes this order of indices
        static_assert(APLCON::ValueIdx==0,"");
        static_assert(APLCON::SigmaIdx==1,"");
        static_assert(APLCON::PullIdx ==2,"");
        // the extra std::tie around std::get is for older compilers...
        return std::tie(std::get<N>(std::tie(Value, Sigma, Pull)));
    }

    friend ostream& operator<<(ostream& s, const N_t& o) {
        return s << o.Value << "+/-" << o.Sigma << "(" << 100.0*o.Sigma/o.Value << "%)";
    }

    static N_t fromIntegral(const TH1D& h) {
        /// \todo check if binwidth is taken into account correctly
        N_t n;
        n.Value = h.IntegralAndError(1, h.GetNbinsX(), n.Sigma, "");
        return n;
    }

    static N_t fromBin(const TH1D& h, int bin) {
        return N_t(h.GetBinContent(bin), h.GetBinError(bin));
    }


};

// helper structs to pass parameters for fitting around

struct fit_params_t {
    static constexpr interval<double> etap_region{920, 990};
    static constexpr int nSamplingBins{10000};
    static constexpr int interpOrder{4};

    static constexpr const char* p_N = "N_sig";
    static constexpr const char* p_sigma = "sigma";
    static constexpr const char* p_delta = "delta";
    static constexpr const char* p_argus_chi = "argus_chi";

    int TaggCh = -1;
    double Eg = std_ext::NaN;
    double start_Nsig = 1e3;

    TH1D* h_mc = nullptr;
    TH1D* h_data = nullptr;
};

struct fit_return_t : ant::root_drawable_traits {

    fit_params_t p;

    RooFitResult* fitresult = nullptr;
    double chi2ndf = std_ext::NaN;
    double peakpos = std_ext::NaN;
    double threshold = std_ext::NaN;

    int numParams() const {
        return fitresult->floatParsFinal().getSize();
    }

    double residualSignalIntegral() {
        auto& h = residual;
        return h->Integral(h->GetXaxis()->FindBin(p.etap_region.Start()),
                           h->GetXaxis()->FindBin(p.etap_region.Stop()))
                /h->getNominalBinWidth();
    }

    N_t getPar(const char* name) const {
        auto& pars = fitresult->floatParsFinal();
        return N_t(dynamic_cast<const RooRealVar&>(*pars.at(pars.index(name))));
    }

    N_t getPar_N() const { return getPar(fit_params_t::p_N); }
    N_t getPar_sigma() const { return getPar(fit_params_t::p_sigma); }
    N_t getPar_delta() const { return getPar(fit_params_t::p_delta); }
    N_t getPar_argus_chi() const { return getPar(fit_params_t::p_argus_chi); }


    RooPlot* fitplot = nullptr;
    RooHist*  h_data = nullptr;
    RooCurve* f_sum = nullptr;
    RooCurve* f_sig = nullptr;
    RooCurve* f_bkg = nullptr;
    RooHist* residual = nullptr;

    N_t N_effcorr;

    double ymax = 160;

    void Draw(const string& option) const override
    {
        auto nsig = getPar_N();
        auto sigma = getPar_sigma();
        auto shift = getPar_delta();

        auto lbl = new TPaveText();
        lbl->SetX1NDC(0.63);
        lbl->SetX2NDC(0.9);
        lbl->SetY1NDC(0.55);
        lbl->SetY2NDC(0.95);
        lbl->SetBorderSize(0);
        lbl->SetFillColor(kWhite);
        lbl->SetTextSize(0.04);
        if(isfinite(p.Eg))
            lbl->AddText(static_cast<string>(std_ext::formatter() << setprecision(1) << fixed << "E_{#gamma} = " << p.Eg << " MeV").c_str());
        lbl->AddText(static_cast<string>(std_ext::formatter() << setprecision(0) << fixed << "N = " << nsig.Value << "#pm " << nsig.Sigma).c_str());
        lbl->AddText(static_cast<string>(std_ext::formatter() << setprecision(0) << fixed << "N/#varepsilon = " << N_effcorr.Value << "#pm " << N_effcorr.Sigma).c_str());
        lbl->AddText(static_cast<string>(std_ext::formatter() << setprecision(2) << fixed << "#chi^{2}_{red} = " << chi2ndf).c_str());
        lbl->AddText(static_cast<string>(std_ext::formatter() << setprecision(1) << fixed << "#sigma = " << sigma.Value << " MeV").c_str());
        lbl->AddText(static_cast<string>(std_ext::formatter() << setprecision(1) << fixed << "#delta = " << shift.Value << " MeV").c_str());

        if(debug) {
            std_ext::formatter extra;
            extra << "Status: ";
            for(unsigned i=0;i<fitresult->numStatusHistory();i++) {
                auto code = fitresult->statusCodeHistory(i);
                if(code != 0)
                    lbl->SetTextColor(kRed);
                extra << code;
            }
            if(p.TaggCh>=0)
                extra << " TaggCh=" << p.TaggCh;
            lbl->AddText(static_cast<string>(extra).c_str());

            TVectorT<double> eigenvalues;
            const auto detCov = fitresult->covarianceMatrix().EigenVectors(eigenvalues);
            eigenvalues.Print();
            LOG(INFO) << eigenvalues.NonZeros();
            LOG(INFO) << eigenvalues.GetNrows();
//            lbl->AddText(static_cast<string>(std_ext::formatter() << setprecision(0) << fixed << "det(Cov) = " << detCov).c_str());

            if(fitresult->status())
                lbl->SetTextColor(kRed);

            // also sets error if Value=0 (division is nan)
            if(!(nsig.Sigma/nsig.Value<0.5))
                lbl->SetTextColor(kRed);

            if(detCov<=0)
                lbl->SetTextColor(kRed);

            if(chi2ndf > 2)
                lbl->SetTextColor(kRed);
        }

        fitplot->SetMinimum(0);
        fitplot->SetMaximum(ymax);
        fitplot->SetTitle("");
        fitplot->Draw(option.c_str());
        lbl->Draw();

    }

    friend ostream& operator<<(ostream& s, const fit_return_t& o) {
        const auto options = "v";
        o.fitresult->printStream(s,o.fitresult->defaultPrintContents(options),o.fitresult->defaultPrintStyle(options));
        return s;
    }
};

// helper functions (mostly templates which can't be lambdas...)

template<typename T, typename Transform>
N_t calcSum(const std::vector<T>& input, Transform transform) {
    std::vector<N_t> Ns(input.size());
    std::transform(input.begin(), input.end(), Ns.begin(), transform);
    N_t N_sum; // sigma=0 means unmeasured

    APLCON::Fit_Settings_t fit_settings;
    fit_settings.ConstraintAccuracy = 1e-2;
    APLCON::Fitter<std::vector<N_t>, N_t> fitter(fit_settings);
    fitter.DoFit(Ns, N_sum, [] (const vector<N_t>& Ns, const N_t& Nsum) {
        double sum = 0.0;
        for(auto& n : Ns)
            sum += n.Value;
        return Nsum.Value - sum;
    });
    return N_sum;
}

template<typename T>
struct draw_TGraph_t : ant::root_drawable_traits {
    T* graph;
    string xlabel;
    string ylabel;
    interval<double> yrange;

    explicit draw_TGraph_t(T* g, const string& xlabel_, const string& ylabel_ = "",
                           const interval<double>& yrange_ = {0,-1}) :
        graph(g), xlabel(xlabel_), ylabel(ylabel_), yrange(yrange_)
    {}

    void Draw(const string& opt) const override {
        graph->Draw(opt.c_str());
        graph->GetXaxis()->SetTitle(xlabel.c_str());
        graph->GetYaxis()->SetTitle(ylabel.c_str());
        if(yrange.IsSane()) {
            graph->SetMinimum(yrange.Start());
            graph->SetMaximum(yrange.Stop());
        }
        // necesary to immediately show changes to multigraph after drawing in canvas
        gPad->Modified();
        gPad->Update();
    }
};

template<typename T, typename... Args>
draw_TGraph_t<T> draw_TGraph(T* g, Args&&... args) {
    return draw_TGraph_t<T>(g, std::forward<Args>(args)...);
}

void calcNEffCorr(N_t N_fit, N_t N_mcreco, N_t N_mcgen, N_t& N_effcorr) {
    // determine efficiency corrected N_effcorr = N_fit/efficiency = N_fit * N_mcgen/N_mcreco;
    APLCON::Fit_Settings_t fit_settings;
    fit_settings.ConstraintAccuracy = 1e-2;
    APLCON::Fitter<N_t, N_t, N_t, N_t> fitter(fit_settings);
    fitter.DoFit(N_fit, N_mcreco, N_mcgen, N_effcorr,
                 [] (const N_t& N_fit, const N_t& N_mcreco, const N_t& N_mcgen, const N_t& N_effcorr) {
        return N_effcorr.Value - N_fit.Value * N_mcgen.Value / N_mcreco.Value;
    });

    if(debug) {
        LOG(INFO) << " N_fit=" << N_fit
                  << " N_effcorr=" << N_effcorr;
    }
}

void traverseCuts(TDirectory* dir, vector<vector<string>>& cuts) {
    auto keys = dir->GetListOfKeys();
    if(!keys)
        return;

    vector<string> dirnames;
    bool h_found = false;
    TIter nextk(keys);
    TKey* key;
    TKey* nextdir = nullptr;
    while((key = (TKey*)nextk()))
    {
        auto classPtr = TClass::GetClass(key->GetClassName());
        if(classPtr->InheritsFrom(TDirectory::Class())) {
            const string dirname(key->GetName());
            if(dirname == "h")
                h_found = true;
            else {
                nextdir = key;;
                dirnames.emplace_back(dirname);
            }
        }
    }

    if(h_found && !dirnames.empty()) {
        cuts.emplace_back(dirnames);
        if(nextdir) {
            traverseCuts(dynamic_cast<TDirectory*>(nextdir->ReadObj()), cuts);
        }
    }
}

vector<vector<string>> extractCuts(const string& prefix, const WrapTFileInput& input) {
    TDirectory* prefixDir = nullptr;
    if(!input.GetObject(prefix, prefixDir))
        throw runtime_error("Cannot find prefix dir " + prefix);
    vector<vector<string>> cuts;
    traverseCuts(prefixDir, cuts);
    return cuts;
}

string pickCutString(const vector<string>& cutchoice, const vector<vector<string>>& cuts) {
    // do some sanity checking for provided cutchoice
    for(auto& choice : cutchoice) {
        bool found = false;
        for(auto& cut : cuts) {
            if(std_ext::contains(cut, choice)) {
                found = true;
                break;
            }
        }
        if(!found)
            throw std::runtime_error(std_ext::formatter()
                                     << "Cut choice '" << choice << "' not found in " << cuts);
    }

    bool foundNothing = false;
    vector<string> pickedCuts;
    for(auto& cut : cuts) {
        auto it_cutchoice = std::find_if(cutchoice.begin(), cutchoice.end(), [cut] (const string& s) {
                       return std_ext::contains(cut, s);
                       });
        if(it_cutchoice != cutchoice.end())
        {
            if(foundNothing)
                throw std::runtime_error("The choice of cuts cannot be realized");
            pickedCuts.emplace_back(*it_cutchoice);
        }
        else {
            foundNothing = true;
        }
    }

    return std_ext::concatenate_string(pickedCuts, "/");
}

string formatCutString(string s) {
    std_ext::replace(s, "/",", ");
    return "\""+s+"\"";
}

void saveAllPads(ant::canvas& c, const std::string& prefix) {
    auto old_pad = gPad;
    c.cd();
    TIter next(gPad->GetListOfPrimitives());
    auto n = 0;
    while(TObject* o = next()) {
        if(TPad* pad = dynamic_cast<TPad*>(o)) {

            TCanvas* c = new TCanvas("tmp_canvas_","tmp_canvas_",600,600);
            c->Modified();
            c->Update();
            c->cd();
            auto p = dynamic_cast<TPad*>(pad->DrawClone());
            p->SetPad(0,0,1,1);

            const auto adjustPad = [] (TPad* p) {
                p->SetTopMargin(0.02);
                p->SetBottomMargin(0.12);
                p->SetRightMargin(0.05);
                p->SetLeftMargin(0.14);

                const auto adjustAxis = [] (TAxis* ax) {
                    ax->SetLabelSize(0.05);
                    ax->SetTitleSize(0.05);
                    ax->SetNdivisions(6,5,0,true);
                };

                TIter next(p->GetListOfPrimitives());

                while(TObject* o = next()) {
                    if(auto h = dynamic_cast<TH1*>(o)) {
                        adjustAxis(h->GetXaxis());
                        adjustAxis(h->GetYaxis());
                        h->GetYaxis()->SetTitleOffset(1.35);
                        h->SetTitle("");
                    }
                    else if(auto h = dynamic_cast<TGraphErrors*>(o)) {
                        adjustAxis(h->GetXaxis());
                        adjustAxis(h->GetYaxis());
                        h->GetYaxis()->SetTitleOffset(1.35);
                        h->SetTitle("");
                    }
                    else if(auto f = dynamic_cast<TFrame*>(o)) {
                        // delete some hidden frames...
                        delete f;
                    }
                }
            };
            adjustPad(p);

            const string filename = prefix+"_"+to_string(n)+".pdf";

            c->SaveAs(filename.c_str());
            delete c;
            n++;
        }
    }
    old_pad->cd();
}

// start reference routines

fit_return_t doReferenceFit(const fit_params_t& p) {

    fit_return_t r;
    r.p = p; // remember params

    const auto calcIMThresh = [] (double Eg) {
        const auto mp = ParticleTypeDatabase::Proton.Mass();
        return std::sqrt(std_ext::sqr(mp) + 2*mp*Eg) - mp;
    };
    r.threshold = calcIMThresh(p.Eg);

    // define observable and ranges
    RooRealVar x("IM","IM", p.h_data->GetXaxis()->GetXmin(), p.h_data->GetXaxis()->GetXmax(), "MeV");
    x.setBins(p.nSamplingBins);

    // for close-to-threshold tagger bins, it's important to handle
    // the fit range and the MC-data shift properly
    const auto bin_threshold = p.h_data->FindBin(r.threshold);
    const auto threshold_upedge = p.h_data->GetXaxis()->GetBinUpEdge(bin_threshold);
    x.setRange("full",x.getMin(), threshold_upedge);

    // load data to be fitted
    RooDataHist h_roo_data("h_roo_data","dataset",x,p.h_data);

    // build shifted mc lineshape
    auto x_shift_range = ParticleTypeDatabase::EtaPrime.GetWindow(30);
    if(x_shift_range.Stop()>threshold_upedge)
        x_shift_range.Stop() = threshold_upedge;
    const auto maxPosMC = p.h_mc->GetXaxis()->GetBinCenter(p.h_mc->GetMaximumBin()); // MC has enough statistics to make this number reliable
    x_shift_range -= maxPosMC;
    if(debug)
        LOG(INFO) << "delta range: " << x_shift_range;
    RooRealVar x_shift(fit_params_t::p_delta, "shift in IM", 0.0, x_shift_range.Start(), x_shift_range.Stop());
    RooProduct x_shift_invert("x_shift_invert","shifted IM",RooArgSet(x_shift, RooConst(-1.0)));
    RooAddition x_shifted("x_shifted","shifted IM",RooArgSet(x,x_shift_invert));
    RooDataHist h_roo_mc("h_roo_mc","MC lineshape", x, p.h_mc);
    RooHistPdf pdf_mc_lineshape("pdf_mc_lineshape","MC lineshape as PDF", x_shifted, x, h_roo_mc, fit_params_t::interpOrder);

    // build detector resolution smearing
    RooRealVar  sigma(fit_params_t::p_sigma,"detector resolution",  2.0, 0.03, 10.0);
    RooGaussian pdf_smearing("pdf_smearing","Single Gaussian", x, RooConst(0.0), sigma);

    // build signal as convolution, note that the gaussian must be the second PDF (see documentation)
    RooFFTConvPdf pdf_signal("pdf_signal","MC_lineshape (X) gauss", x, pdf_mc_lineshape, pdf_smearing);

    // build background with ARGUS function
    RooRealVar argus_cutoff("argus_cutoff","argus pos param", r.threshold);
    RooRealVar argus_shape(fit_params_t::p_argus_chi,"argus shape param #chi", -5, -25.0, 0.0);
    RooRealVar argus_p("argus_p","argus p param", 0.5);
    RooArgusBG pdf_background("pdf_background","bkg argus",x,argus_cutoff,argus_shape,argus_p);

    // build sum
    RooRealVar nsig(fit_params_t::p_N,"number signal events", p.start_Nsig, 0, 1e6);
    RooRealVar nbkg("N_bkg","number background events", 1e3, 0, 1e6);
    RooAddPdf pdf_sum("pdf_sum","total sum",RooArgList(pdf_signal,pdf_background),RooArgList(nsig,nbkg));

    // do some pre-fitting to obtain better starting values, make sure function is non-zero in range
//    x.setRange("nonzero",x.getMin(), threshold-5);
//    pdf_sum.chi2FitTo(h_roo_data, Range("full"), PrintLevel(-1)); // using Range(..., ...) does not work here (bug in RooFit, sigh)

    // do the actual maximum likelihood fit
    // use , Optimize(false), Strategy(2) for double gaussian...?!
    r.fitresult = pdf_sum.fitTo(h_roo_data, Extended(), SumW2Error(kTRUE), Range("full"), Save(), PrintLevel(debug ? 3 : -1));

    // draw output and remember pointer
    r.fitplot = x.frame();
    r.fitplot->SetName(("FitPlot_TaggCh="+to_string(p.TaggCh)).c_str());

    h_roo_data.plotOn(r.fitplot);
    r.h_data = dynamic_cast<RooHist*>(r.fitplot->findObject(0));

    // need to figure out chi2nds and stuff after plotting data and finally fitted pdf_sum
    // also the residHist must be created here (and rememebered for later use)
    pdf_sum.plotOn(r.fitplot, LineColor(kRed));
    //    pdf_sum.plotOn(frame, LineColor(kRed), VisualizeError(*fr));
    r.f_sum = dynamic_cast<RooCurve*>(r.fitplot->findObject(0));
    r.chi2ndf = r.fitplot->chiSquare(r.numParams());

    auto pdf_sum_tf = pdf_sum.asTF(RooArgList(x), RooArgList(*pdf_sum.getParameters(x)), RooArgSet(x));
    constexpr interval<double> peak_region(fit_params_t::etap_region);
    r.peakpos = pdf_sum_tf->GetMaximumX(peak_region.Start(), peak_region.Stop());
    r.residual = r.fitplot->residHist();

    pdf_sum.plotOn(r.fitplot, Components(pdf_background), LineColor(kBlue));
    r.f_bkg = dynamic_cast<RooCurve*>(r.fitplot->findObject(0));
    pdf_sum.plotOn(r.fitplot, Components(pdf_signal), LineColor(kGreen));
    r.f_sig = dynamic_cast<RooCurve*>(r.fitplot->findObject(0));
    return r;
}

N_t doReference(const WrapTFileInput& input,
                const WrapTFileInput& mctestinput,
                const unique_ptr<ofstream>& textout,
                const std::string& imgdir,
                vector<string> cutchoice,
                const interval<int>& taggChRange) {
    analysis::HistogramFactory::DirStackPush HistFacDir(analysis::HistogramFactory("Ref"));

    auto Tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();

    // add some default cuts
    cutchoice.emplace_back("DiscardedEk=0");
    cutchoice.emplace_back("KinFitProb>0.02");

    const string ref_prefix   = "EtapOmegaG_plot_Ref";
    const auto pickedCut = pickCutString(cutchoice, extractCuts(ref_prefix, input));
    LOG(INFO) << "Reference picked cut: " << pickedCut;
    const string ref_histpath = ref_prefix + (pickedCut.empty() ? "" : "/") + pickedCut;
    const string ref_histname = "h_IM_2g_TaggCh";

    TH2D* ref_data;
    TH2D* ref_mc;
    TH1D* ref_mctrue_generated;

    if(mctestinput.NumberOfFiles()>0)
    {
        const string histpath = ref_histpath+"/h/Sum_MC/"+ref_histname;
        if(!mctestinput.GetObject(histpath, ref_data)) {
            throw runtime_error("Cannot find " + histpath);
        }
    }
    else
    {
        const string histpath = ref_histpath+"/h/Data/"+ref_histname;
        if(!input.GetObject(histpath, ref_data)) {
            throw runtime_error("Cannot find " + histpath);
        }
    }
    {
        const string histpath = ref_histpath+"/h/Ref/"+ref_histname;
        if(!input.GetObject(histpath, ref_mc)) {
            throw runtime_error("Cannot find " + histpath);
        }
    }
    {
        const string histpath = ref_prefix+"/h_mctrue_generated";
        if(!input.GetObject(histpath, ref_mctrue_generated)) {
            throw runtime_error("Cannot find " + histpath);
        }
    }

    // start creating the overview (more will be added after fits)
    ant::canvas c_overview("Ref Overview");
    c_overview << drawoption("colz")
               << ref_mc << ref_data;

    std::vector<fit_return_t> fit_results;

    if(textout) {
        *textout << formatCutString(pickedCut) << '\n';
    }

    for(auto taggch=taggChRange.Stop();taggch>=taggChRange.Start();taggch--) {
        LOG(INFO) << "Fitting TaggCh=" << taggch;
        fit_params_t p;
        p.TaggCh = taggch;
        p.Eg = Tagger->GetPhotonEnergy(taggch);

        // fit MC lineshape to data
        const auto taggbin = taggch+1;
        p.h_mc   = ref_mc->ProjectionX("h_mc",taggbin,taggbin);
        p.h_data = ref_data->ProjectionX("h_data",taggbin,taggbin);

        // higher tagger channels have quite low number of signal events
        // provide better starting value in this case
        if(p.TaggCh>=38)
            p.start_Nsig = 1e2;

        auto r = doReferenceFit(p);

        calcNEffCorr(r.getPar_N(), // N from fit
                     N_t::fromIntegral(*p.h_mc), // N_mcreco
                     N_t::fromBin(*ref_mctrue_generated, taggbin), // N_mcgen
                     r.N_effcorr // output
                     );

        // save/plot values
        fit_results.emplace_back(r);

        if(debug) {
            LOG(INFO) << r;
        }

        if(textout) {
            const auto Eg_err = Tagger-> GetPhotonEnergyWidth(r.p.TaggCh)/2;
            *textout << p.Eg << " " << Eg_err << " "
                     << r.getPar_N().Value << " " << r.getPar_N().Sigma << " "
                     << r.N_effcorr.Value << " " << r.N_effcorr.Sigma
                     << '\n';
        }
    }

    if(textout) {
        *textout << '\n' << '\n';
    }

    // plot all single fits
    ant::canvas c_plots("Ref Plots");
    for(auto& r : fit_results)
        c_plots << r;
    c_plots << endc;

    if(!imgdir.empty())
        saveAllPads(c_plots, imgdir+"/Ref_TaggCh");


    // sum up the N_data and N_effcorr

    auto N_fit_sum = calcSum(fit_results, [] (const fit_return_t& r) {
        return r.getPar_N(); // N from fit
    });
    LOG(INFO) << "Sum of N_fit: " << N_fit_sum;

    auto N_effcorr_sum = calcSum(fit_results, [] (const fit_return_t& r) {
        return r.N_effcorr;
    });

    // plot summation over tagg channels, calc total chi2
    {
        /// \todo maybe interpolate the sums again? summing them up all time
        /// is rather slow, and we also need to copy the fit_results several times, sigh...
        const auto  makeTF1sum = [fit_results] (RooCurve* fit_return_t::* PtrToMember) -> TF1* {
            if(fit_results.empty())
                return nullptr;
            auto& r = fit_results.front();
            auto axis = r.fitplot->GetXaxis();
            auto f = new TF1("", [fit_results, PtrToMember] (double* x, double*) {
                double sum = 0;
                for(auto& r : fit_results) {
                    // individual thresholds of results are lower than global
                    if(x[0] < r.threshold)
                        sum += (r.*PtrToMember)->Eval(x[0]);
                }
                return sum;
            }, axis->GetXmin(), fit_results.back().threshold, 0);
            f->SetNpx(1000);
            f->SetLineColor((r.*PtrToMember)->GetLineColor());
            f->SetLineWidth((r.*PtrToMember)->GetLineWidth());
            return f;
        };

        auto& r = fit_results.front();

        auto h_proj_all_ch = ref_data->ProjectionX("proj_all_taggch",taggChRange.Start()+1,taggChRange.Stop()+1);
        h_proj_all_ch->GetYaxis()->SetTitle(r.fitplot->GetYaxis()->GetTitle());
        h_proj_all_ch->SetMarkerStyle(r.h_data->GetMarkerStyle());
        h_proj_all_ch->SetLineColor(r.h_data->GetLineColor());
        h_proj_all_ch->SetStats(0);

        auto f_sum = makeTF1sum(&fit_return_t::f_sum);

        const auto calcApproxNDF = [r, f_sum] () {
            double x_low, x_high;
            f_sum->GetRange(x_low, x_high);
            // the number of parameters for a single is used here,
            // don't know how to handle this 100% correct
            const auto binwidth = r.h_data->getNominalBinWidth();
            return (x_high-x_low)/binwidth - r.numParams();
        };

        const auto chi2 = h_proj_all_ch->Chisquare(f_sum,"R");
        const auto ndf = calcApproxNDF();
        auto lbl = new TPaveText();
        lbl->SetX1NDC(0.6);
        lbl->SetX2NDC(0.88);
        lbl->SetY1NDC(0.75);
        lbl->SetY2NDC(0.95);
        lbl->SetBorderSize(0);
        lbl->SetFillColor(kWhite);
        lbl->SetTextSize(0.04);
        lbl->AddText(static_cast<string>(std_ext::formatter() << setprecision(0) << fixed << "N = " << N_fit_sum.Value << "#pm " << N_fit_sum.Sigma).c_str());
        lbl->AddText(static_cast<string>(std_ext::formatter() << setprecision(0) << fixed << "N/#varepsilon = " << N_effcorr_sum.Value << "#pm " << N_effcorr_sum.Sigma).c_str());
        lbl->AddText(static_cast<string>(std_ext::formatter() << setprecision(0) << fixed << "#chi^{2}_{red} = " << chi2 << "/" << ndf << "#approx " << setprecision(2) << chi2/ndf).c_str());

        c_overview << drawoption("E1") << h_proj_all_ch
                   << samepad << makeTF1sum(&fit_return_t::f_sig)
                   << samepad << makeTF1sum(&fit_return_t::f_bkg)
                   << samepad << f_sum
                   << samepad << lbl;

    }

    auto makeTGraph = [] (const string& title = "", Color_t color = kBlack) {
        auto g = new TGraphErrors();
        g->SetTitle(title.c_str());
        g->SetFillColor(kWhite);
        g->SetLineColor(color);
        g->SetLineWidth(3);
        g->SetMarkerSize(0);
        return g;
    };

    auto setEgPoint = [Tagger] (const fit_return_t& r, TGraphErrors* g, const N_t& y) {
        auto n = g->GetN();
        g->SetPoint(n, r.p.Eg, y.Value);
        g->SetPointError(n, Tagger-> GetPhotonEnergyWidth(r.p.TaggCh)/2, y.Sigma);
    };

    // number of events and effcorr events as multigraph
    {
        auto g_N_fit = makeTGraph("N", kRed);
        auto g_N_effcorr = makeTGraph("N/#varepsilon", kBlue);

        for(const fit_return_t& r : fit_results) {
            setEgPoint(r, g_N_fit,  r.getPar_N());
            setEgPoint(r, g_N_effcorr, r.N_effcorr);
        }

        auto multigraph = new TMultiGraph("g_N_fit_effcorr","Ref: Number of events");
        multigraph->Add(g_N_fit);
        multigraph->Add(g_N_effcorr);

        c_overview << drawoption("AP") << padoption::Legend
                   << draw_TGraph(multigraph, "E_{#gamma} / MeV", "Events", interval<double>{0, 4400})
                   << endc;
        if(!imgdir.empty())
            saveAllPads(c_overview, imgdir+"/Ref_Overview");
    }

    // fit parameters vs photon energy
    {
        auto g_par_chi2  = makeTGraph();
        auto g_par_delta = makeTGraph();
        auto g_par_sigma = makeTGraph();
        auto g_par_argus_chi = makeTGraph();

        for(const auto& r : fit_results) {
            setEgPoint(r, g_par_chi2, N_t(r.chi2ndf));
            setEgPoint(r, g_par_delta, r.getPar_delta());
            setEgPoint(r, g_par_sigma, r.getPar_sigma());
            setEgPoint(r, g_par_argus_chi, r.getPar_argus_chi());
        }

        canvas c("Ref Fit Parameters");
        c << drawoption("AP")
          << draw_TGraph(g_par_chi2,      "E_{#gamma} / MeV", "#chi^{2}_{red}")
          << draw_TGraph(g_par_delta,     "E_{#gamma} / MeV", "#delta_{Data-MC} / MeV")
          << draw_TGraph(g_par_sigma,     "E_{#gamma} / MeV", "#sigma_{smearing} / MeV")
          << draw_TGraph(g_par_argus_chi, "E_{#gamma} / MeV", "#chi_{ARGUS}")
          << endc;
        if(!imgdir.empty())
            saveAllPads(c, imgdir+"/Ref_FitParams");
    }

    return N_effcorr_sum;
}

// start signal routines

const string sig_prefix   = "EtapOmegaG_plot_Sig";

N_t doSignal(const string& sig_histpath,
              const WrapTFileInput& input,
              const WrapTFileInput& mctestinput,
              N_t& N_fit, bool showcanvas)
{
    TH1D* sig_data;
    TH1D* sig_mc;
    TH1D* sig_mctrue_generated;
    {

        const string sig_histname = "h_IM_4g";

        if(mctestinput.NumberOfFiles()>0)
        {
            const string histpath = sig_histpath+"/h/Sum_MC/"+sig_histname;
            if(!mctestinput.GetObject(histpath, sig_data)) {
                throw runtime_error("Cannot find " + histpath);
            }
        }
        else
        {
            const string histpath = sig_histpath+"/h/Data/"+sig_histname;
            if(!input.GetObject(histpath, sig_data)) {
                throw runtime_error("Cannot find " + histpath);
            }
        }
        {
            const string histpath = sig_histpath+"/h/Sig/"+sig_histname;
            if(!input.GetObject(histpath, sig_mc)) {
                throw runtime_error("Cannot find " + histpath);
            }
        }
        {
            const string histpath = sig_prefix+"/h_mctrue_generated";
            if(!input.GetObject(histpath, sig_mctrue_generated)) {
                throw runtime_error("Cannot find " + histpath);
            }
        }

    }

    analysis::HistogramFactory::DirStackPush HistFacDir(analysis::HistogramFactory("Sig"));


    // start creating the overview (more will be added after fits)
    ant::canvas c_overview("Sig Overview");
    c_overview << drawoption("colz")
               << sig_mc << sig_data
               << sig_mctrue_generated;

    // define observable and ranges
    RooRealVar x("IM","IM", sig_data->GetXaxis()->GetXmin(), sig_data->GetXaxis()->GetXmax(), "MeV");
    x.setBins(10000);
    x.setRange("full",x.getMin(),x.getMax());

    // load data to be fitted
    RooDataHist h_roo_data("h_roo_data","dataset",x,sig_data);

    // build shifted mc lineshape
    RooRealVar x_shift(fit_params_t::p_delta, "shift in IM", 0.0, -10.0, 10.0);
    RooProduct x_shift_invert("x_shift_invert","shifted IM",RooArgSet(x_shift, RooConst(-1.0)));
    RooAddition x_shifted("x_shifted","shifted IM",RooArgSet(x,x_shift_invert));
    RooDataHist h_roo_mc("h_roo_mc","MC lineshape", x, sig_mc);
    RooHistPdf pdf_mc_lineshape("pdf_mc_lineshape","MC lineshape as PDF", x_shifted, x, h_roo_mc, fit_params_t::interpOrder);

    // build detector resolution smearing

    RooRealVar  sigma(fit_params_t::p_sigma,"detector resolution",  2.0, 0.01, 10.0);
    RooGaussian pdf_smearing("pdf_smearing","Single Gaussian", x, RooConst(0.0), sigma);

    // build signal as convolution, note that the gaussian must be the second PDF (see documentation)
    RooFFTConvPdf pdf_signal("pdf_signal","MC_lineshape (X) gauss", x, pdf_mc_lineshape, pdf_smearing);

    // build background with generalized ARGUS
    RooRealVar argus_cutoff("argus_cutoff","argus pos param", 1025);
    RooRealVar argus_shape(fit_params_t::p_argus_chi,"argus shape param #chi", -15, -25.0, 0.0);
    RooRealVar argus_p("argus_p","argus p param", 1.5, 0.5, 4);
    RooArgusBG pdf_background("pdf_background","bkg argus",x,argus_cutoff,argus_shape,argus_p);

    // build sum
    RooRealVar nsig(fit_params_t::p_N,"number signal events", 1e3, 0, 1e4);
    RooRealVar nbkg("N_bkg","number background events", 1e4, 0, 1e6);
    RooAddPdf pdf_sum("pdf_sum","total sum",RooArgList(pdf_signal,pdf_background),RooArgList(nsig,nbkg));

    // do some pre-fitting to obtain better starting values, make sure function is non-zero in range
    //    x.setRange("nonzero",x.getMin(), threshold-5);
    //    pdf_sum.chi2FitTo(h_roo_data, Range("nonzero"), PrintLevel(-1)); // using Range(..., ...) does not work here (bug in RooFit, sigh)

    // do the actual maximum likelihood fit
    // use , Optimize(false), Strategy(2) for double gaussian...?!

    fit_return_t r;
    r.ymax = 300;

    r.fitresult = pdf_sum.fitTo(h_roo_data, Extended(), SumW2Error(kTRUE), Range("full"), Save(), PrintLevel(debug ? 3 : -1));

    if(showcanvas)
        r.fitresult->Print();

    // draw output and remember pointer
    r.fitplot = x.frame();

    h_roo_data.plotOn(r.fitplot);
    r.h_data = dynamic_cast<RooHist*>(r.fitplot->findObject(0));

    // need to figure out chi2nds and stuff after plotting data and finally fitted pdf_sum
    // also the residHist must be created here (and rememebered for later use)
    pdf_sum.plotOn(r.fitplot, LineColor(kRed));
    //    pdf_sum.plotOn(frame, LineColor(kRed), VisualizeError(*fr));
    r.f_sum = dynamic_cast<RooCurve*>(r.fitplot->findObject(0));
    r.chi2ndf = r.fitplot->chiSquare(r.numParams());

    auto pdf_sum_tf = pdf_sum.asTF(RooArgList(x), RooArgList(*pdf_sum.getParameters(x)), RooArgSet(x));
    constexpr interval<double> peak_region(fit_params_t::etap_region);
    r.peakpos = pdf_sum_tf->GetMaximumX(peak_region.Start(), peak_region.Stop());
    r.residual = r.fitplot->residHist();

    pdf_sum.plotOn(r.fitplot, Components(pdf_background), LineColor(kBlue));
    r.f_bkg = dynamic_cast<RooCurve*>(r.fitplot->findObject(0));
    pdf_sum.plotOn(r.fitplot, Components(pdf_signal), LineColor(kGreen));
    r.f_sig = dynamic_cast<RooCurve*>(r.fitplot->findObject(0));

    // do efficiency correction
    // (simple here, as integrated over all tagger channels)
    N_fit = N_t(nsig);
    calcNEffCorr(N_fit,
                 N_t::fromIntegral(*sig_mc),
                 N_t::fromIntegral(*sig_mctrue_generated),
                 r.N_effcorr
                 );

    // r is completely filled now, then we can plot it
    c_overview << r;
    if(showcanvas)
        c_overview << endc;

    // subsequent analysis just needs effcorr number of events
    return r.N_effcorr;
}

N_t calcBranchingRatio(N_t N_sig_events, N_t N_etap) {

    N_t BR_pi0_2g(99.823/100.0,0.034/100.0); // branching ratio pi0->2g is about 100 % (PDG)
    N_t BR_omega_pi0g(8.28/100.0,0.28/100.0); // branching ratio omega->pi0 g is about 8.3 % (PDG)

    N_t BR_etap_omega_g(0,0);
    {
        APLCON::Fit_Settings_t fit_settings;
        fit_settings.ConstraintAccuracy = 1e-2;
        APLCON::Fitter<N_t, N_t, N_t, N_t, N_t> fitter(fit_settings);
        fitter.DoFit(N_sig_events, N_etap, BR_pi0_2g, BR_omega_pi0g, BR_etap_omega_g, [] (
                     const N_t& N_sig_events, const N_t& N_etap, const N_t& BR_pi0_2g, const N_t& BR_omega_pi0g, const N_t& BR_etap_omega_g) {
            return BR_etap_omega_g.Value - N_sig_events.Value/BR_pi0_2g.Value/BR_omega_pi0g.Value/N_etap.Value;
        });
    }

    return BR_etap_omega_g;
}

struct TCLAPInterval : interval<int> {
    using interval::interval;
    using ValueCategory = TCLAP::ValueLike;
};

int main(int argc, char** argv) {
    SetupLogger();

    // parse command line
    TCLAP::CmdLine cmd("EtapOmegaG_fit", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");
    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b","batch","Run in batch mode (no ROOT shell afterwards)",false);
    auto cmd_debug = cmd.add<TCLAP::MultiSwitchArg>("","debug","Enable debug mode",false);

    auto cmd_input = cmd.add<TCLAP::ValueArg<string>>("i","input","ROOT input file",true,"","rootfile");
    auto cmd_mctestinput = cmd.add<TCLAP::ValueArg<string>>("","mctestinput","Input for MC in/out test",false,"","rootfile");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o","output","Output file",false,"","filename");
    TCLAP::ValuesConstraintExtra<decltype(ExpConfig::Setup::GetNames())> allowedsetupnames(ExpConfig::Setup::GetNames());
    auto cmd_setup  = cmd.add<TCLAP::ValueArg<string>>("s","setup","Choose setup by name",true,"", &allowedsetupnames);

    auto cmd_skipref = cmd.add<TCLAP::MultiSwitchArg>("","skipref","Skip analysis of reference channel",false);
    auto cmd_skipsig = cmd.add<TCLAP::MultiSwitchArg>("","skipsig","Skip analysis of signal channel",false);
    auto cmd_taggerrange = cmd.add<TCLAP::ValueArg<TCLAPInterval>>("","taggerrange","tagger range for reference, ex. 4-8",false,TCLAPInterval{0,40},"channels");
    auto cmd_cut = cmd.add<TCLAP::MultiArg<string>>("c","cut","Select cuts instead of default provided", false, "");
    auto cmd_textout = cmd.add<TCLAP::ValueArg<string>>("","textout","Dump numbers to file as text (gnuplot compatible)",false,"", "");
    auto cmd_imgdir = cmd.add<TCLAP::ValueArg<string>>("","imgdir","Output folder for SaveMultiImages calls",false,"", "");

    cmd.parse(argc, argv);

    // verbosity management
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }
    debug = cmd_debug->isSet();
    RooMsgService::instance().setGlobalKillBelow(debug ? RooFit::WARNING : RooFit::ERROR);
    if(!debug)
        RooMsgService::instance().setSilentMode(true);

    auto textout_stream = cmd_textout->isSet() ? std_ext::make_unique<ofstream>(cmd_textout->getValue()) : nullptr;
    if(textout_stream && !*textout_stream) {
        LOG(ERROR) << "Cannot open text output " << cmd_textout->getValue();
        return EXIT_FAILURE;
    }

    // open input file, config setup, other stuff..
    // keep inputs open, otherwise hist pointers become invalid

    ExpConfig::Setup::SetByName(cmd_setup->getValue());

    WrapTFileInput input(cmd_input->getValue());
    WrapTFileInput mctestinput;
    if(cmd_mctestinput->isSet())
        mctestinput.OpenFile(cmd_mctestinput->getValue());

    const bool skipRef = cmd_skipref->isSet();
    const bool skipSig = cmd_skipsig->isSet();

    const auto taggChRange = cmd_taggerrange->getValue();
    if(!taggChRange.IsSane()) {
        LOG(ERROR) << "Provided tagger channel range " << taggChRange << " not sane.";
        return EXIT_FAILURE;
    }
    if(cmd_taggerrange->isSet()) {
        LOG(WARNING) << "Using non-default tagger channel range may not yield correct results (use for debugging only)";
    }

    // create TRint as RooFit internally creates functions/histograms,
    // prevents this stupid gStyle=0 related error, sigh...
    argc=0; // prevent TRint to parse any cmdline
    TRint app("EtapOmegaG_fit",&argc,argv,nullptr,0,true);
    if(cmd_batchmode->isSet())
        gROOT->SetBatch(true);

    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        // cd into masterFile upon creation
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(), true);
    }

    // get total number of reference events (effcorr), thus produced eta primes...
    N_t N_etap(0,0);
    if(skipRef) {
        // do not stay silent if reference is skipped
        LOG(WARNING) << "Skipping reference channel analysis, using pre-calculated value";
        N_etap.Value = 5.12072e+06;
        N_etap.Sigma = 191409;
    }
    else
    {
        N_t BR_etap_2g(2.20/100.0,0.08/100.0); // branching ratio eta'->2g is about 2.2 % (PDG)
        auto N_ref_events = doReference(input, mctestinput, textout_stream,
                                        cmd_imgdir->getValue(),
                                        cmd_cut->getValue(), taggChRange);
        LOG(INFO) << "Number of eta' -> 2g events (effcorr): " << N_ref_events;
        APLCON::Fit_Settings_t fit_settings;
        fit_settings.ConstraintAccuracy = 1e-2;
        APLCON::Fitter<N_t, N_t, N_t> fitter(fit_settings);
        fitter.DoFit(N_ref_events, BR_etap_2g, N_etap, [] (
                     const N_t& N_ref_events, const N_t& BR_etap_2g, const N_t& N_etap) {
            return N_etap.Value - N_ref_events.Value/BR_etap_2g.Value;
        });
    }
    LOG(INFO) << "Number of tagged eta' in EPT 2014 beamtime: " << N_etap;


    // get total number of signal events
    if(!skipSig) {
        const string sig_path = sig_prefix+"/SigPi0";

        vector<string> sig_histpaths;

        if(textout_stream) {
            const auto cuts = extractCuts(sig_path, input);
            LOG(INFO) << "Will run signal fits on all possible cuts made from " << cuts;

            // transform cuts into state using ranges
            using it_t = decltype(cuts.front().cbegin());
            struct state_item_t {
                it_t begin;
                it_t end;
                it_t curr;
                state_item_t(it_t begin_, it_t end_) : begin(begin_), end(end_), curr(begin_) {}
            };

            vector<state_item_t> state;
            for(auto& cut : cuts)
                state.emplace_back(cut.cbegin(), cut.cend());

            bool running = !state.empty();
            while(running) {

                // add the current state to histpaths
                sig_histpaths.emplace_back(sig_prefix+"/SigPi0");
                for(const auto& s : state) {
                    sig_histpaths.back() += "/" + *s.curr;
                }

                // go to next state
                auto it = state.begin();
                while(running) {
                    ++(it->curr);
                    if(it->curr != it->end)
                        break;

                    it->curr = it->begin;
                    ++it;
                    running = it != state.end();
                }
            }
        }
        else {
            // no textout stream, then run on pre-selected cut string
            sig_histpaths.emplace_back(
                        sig_path+"/DiscardedEk=0"
                                 "/AntiPi0FitProb<10^{-5}||nan"
                                 "/AntiEtaFitProb<10^{-4}||nan"
                                 "/TreeFitProb>0.1"
                                 "/gNonPi0_2"
                                 "/CBSumVetoE_gNonPi0<0.2"
                                 "/IM_Pi0g[1]");
        }

        for(const auto& sig_histpath : sig_histpaths) {
            LOG(INFO) << "Cut selection: " << sig_histpath;

            N_t N_fit;
            const auto N_sig_events = doSignal(sig_histpath, input, mctestinput, N_fit, sig_histpaths.size()==1);

            const auto BR_etap_omega_g = calcBranchingRatio(N_sig_events, N_etap);
            if(BR_etap_omega_g.Sigma/BR_etap_omega_g.Value>0.5)
                LOG(WARNING) << "Very large error: " << sig_histpath;

            if(textout_stream) {
                *textout_stream << "# " << sig_histpath << endl
                                << BR_etap_omega_g.Value << " " << BR_etap_omega_g.Sigma << " "
                                << N_fit.Value << " " << N_fit.Sigma << endl;
            }
            else {
                LOG(INFO) << "Number of eta' -> omega g events          : " << N_fit;
                LOG(INFO) << "Number of eta' -> omega g events (effcorr): " << N_sig_events;
                LOG(INFO) << "BR(eta' -> omega g)          : " << BR_etap_omega_g;
                const N_t BR_etap_omega_g_expected(2.75/100.0,0.23/100.0); // branching ratio eta'->omega g is about 2.8 % (PDG)
                LOG(INFO) << "BR(eta' -> omega g) expected : " << BR_etap_omega_g_expected;
            }
        }
    }

    // run TRint
    if(!cmd_batchmode->isSet()) {
        if(!std_ext::system::isInteractive()) {
            LOG(INFO) << "No TTY attached. Not starting ROOT shell.";
        }
        else {
            if(masterFile)
                LOG(INFO) << "Close ROOT properly to write data to disk.";

            app.Run(kTRUE); // really important to return...
            if(masterFile)
                LOG(INFO) << "Writing output file...";
            masterFile = nullptr;   // and to destroy the master WrapTFile before TRint is destroyed
        }
    }

    return EXIT_SUCCESS;
}
