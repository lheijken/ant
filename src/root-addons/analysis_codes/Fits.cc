#include "Fits.h"

#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "root-addons/analysis_codes/Math.h"
#include "base/ParticleType.h"
#include "TStyle.h"
#include "TLatex.h"
#include "base/std_ext/string.h"
#include "TGraph.h"
#include "analysis/plot/RootDraw.h"
#include <TMultiGraph.h>
#include "base/ParticleType.h"
#include "base/math_functions/CrystalBall.h"
#include "TGraphErrors.h"

using namespace ant;
using namespace std;

Fits::FitResult Fits::FitEtaCalib(TH1* h, const double r_min, const double r_max) {
    return FitPeakCrystalBallPol6(h, ParticleTypeDatabase::Eta.Mass(), 21.0, r_min, r_max);
}

Fits::FitResult Fits::FitPi0Calib(TH1* h, const double r_min, const double r_max) {
    return FitPeakCrystalBallPol4(h, ParticleTypeDatabase::Pi0.Mass(), 11.0, r_min, r_max);
}

Fits::FitResult Fits::FitPi0Calib0(TH1* h, const double r_min, const double r_max) {
    return FitPeakCrystalBallPol0(h, ParticleTypeDatabase::Pi0.Mass(), 11.0, r_min, r_max);
}
Fits::FitResult Fits::FitPi0CalibGaussian(TH1* h, const double r_min, const double r_max) {
    return FitPeakPol0(h, ParticleTypeDatabase::Pi0.Mass(), 11.0, r_min, r_max);
}



Fits::FitResult Fits::FitPeakPol0(TH1* h, const double mass, const double expected_width, const double r_min, const double r_max) {

    const int    npx   = 500;

    TF1* sig = new TF1("sig", "gaus", r_min, r_max);
    sig->SetLineColor(kGreen);
    sig->SetNpx(npx);

    // height
    sig->SetParameter(0, 0.5 * h->GetMaximum());

    // peak position
    sig->SetParameter(1, mass);

    // width
    sig->SetParameter(2, expected_width);



    TF1* bg = new TF1("bg", "pol0", r_min, r_max);
    bg->SetLineColor(kBlue);

    bg->SetParameter(0,0);
    bg->SetParName(0, "BG p_{0}");
    bg->SetParameter(1,0);
    bg->SetParName(1, "BG p_{1}");
    bg->SetParameter(2,0);
    bg->SetParName(2, "BG p_{2}");
    bg->SetParameter(3,0);
    bg->SetParName(3, "BG p_{3}");
    bg->SetParameter(4,0);
    bg->SetParName(4, "BG p_{4}");

    const auto peak_range = interval<double>::CenterWidth(mass,2*2*mass);

    TFSum::FitRanged(h, bg, r_min, peak_range.Start(), peak_range.Stop(), r_max);


    TFSum* sum = new TFSum("sum", sig, bg, r_min, r_max);
    sum->SetNpx(npx);


    //TCanvas* c = new TCanvas();
    //c->SetTitle(Form("Fit to %s", h->GetName()));
    h->SetStats(true);
    gStyle->SetOptFit(1);
    //h->Draw();
    h->Fit(sum->Function(), "REM0NB");

    sum->SyncToFcts();

    sum->Draw();


    const double total_area = sum->Function()->Integral(r_min, r_max);
    const double bg_area    =  bg->Integral(r_min, r_max);
    const double sig_area   = total_area - bg_area;
    const double peak_pos  = sig->GetParameter(1);
    const double peak_width = sig->GetParameter(2);
    const double relwidth   = peak_width / peak_pos;
    const double mean_error = sig->GetParError(1);

    cout << "Mass offset = " << peak_pos - mass << " MeV\n";
    cout << "sigma/pos   = " <<  relwidth << "\n";
    cout << "Sig         = " << sig_area << endl;

    const string text = std_ext::formatter()
            << "N = " << sig_area << " "
            << "#sigma/#mu = " << relwidth;
    auto l = new TLatex(peak_pos, sum->Function()->Eval(peak_pos), text.c_str());
    l->Draw();


    return FitResult(peak_pos, sig_area, sig->GetParameter(2), (sum->Function()->GetChisquare()) / (sum->Function()->GetNDF()), sum->Function(), bg, sig,mean_error);

}




Fits::FitResult Fits::FitPeakPol4(TH1* h, const double mass, const double expected_width, const double r_min, const double r_max) {

    const int    npx   = 500;

    TF1* sig = new TF1("sig", "gaus", r_min, r_max);
    sig->SetLineColor(kGreen);
    sig->SetNpx(npx);

    // height
    sig->SetParameter(0, 0.5 * h->GetMaximum());

    // peak position
    sig->SetParameter(1, mass);

    // width
    sig->SetParameter(2, expected_width);



    TF1* bg = new TF1("bg", "pol4", r_min, r_max);
    bg->SetLineColor(kBlue);

    bg->SetParameter(0,0);
    bg->SetParName(0, "BG p_{0}");
    bg->SetParameter(1,0);
    bg->SetParName(1, "BG p_{1}");
    bg->SetParameter(2,0);
    bg->SetParName(2, "BG p_{2}");
    bg->SetParameter(3,0);
    bg->SetParName(3, "BG p_{3}");
    bg->SetParameter(4,0);
    bg->SetParName(4, "BG p_{4}");

    const auto peak_range = interval<double>::CenterWidth(mass,2*2*mass);

    TFSum::FitRanged(h, bg, r_min, peak_range.Start(), peak_range.Stop(), r_max);


    TFSum* sum = new TFSum("sum", sig, bg, r_min, r_max);
    sum->SetNpx(npx);


    //TCanvas* c = new TCanvas();
    //c->SetTitle(Form("Fit to %s", h->GetName()));
    h->SetStats(true);
    gStyle->SetOptFit(1);
    //h->Draw();
    h->Fit(sum->Function(), "REM0NB");

    sum->SyncToFcts();

    sum->Draw();


    const double total_area = sum->Function()->Integral(r_min, r_max);
    const double bg_area    =  bg->Integral(r_min, r_max);
    const double sig_area   = total_area - bg_area;
    const double peak_pos  = sig->GetParameter(1);
    const double peak_width = sig->GetParameter(2);
    const double relwidth   = peak_width / peak_pos;
    const double mean_error = sig->GetParameter(3);


    cout << "Mass offset = " << peak_pos - mass << " MeV\n";
    cout << "sigma/pos   = " <<  relwidth << "\n";
    cout << "Sig         = " << sig_area << endl;

    const string text = std_ext::formatter()
            << "N = " << sig_area << " "
            << "#sigma/#mu = " << relwidth;
    auto l = new TLatex(peak_pos, sum->Function()->Eval(peak_pos), text.c_str());
    l->Draw();


    return FitResult(peak_pos, sig_area, sig->GetParameter(2), (sum->Function()->GetChisquare()) / (sum->Function()->GetNDF() ), sum->Function(), bg, sig,mean_error);

}


Fits::FitResult Fits::FitPeakCrystalBallPol6(TH1* h, const double mass, const double expected_width, const double r_min, const double r_max) {

    const int    npx   = 500;

    TF1* sig = ant::math::CrystalBall::GetTF1();
    //CrystalBall-Function f(x ; alpha , n , sigma, x- ) ;
    //alpha=Positon of the transition between normal destribution and power law
    //n= free parameter of the power law
    //sigma= standard deviation
    //x-= expected value
    sig->SetRange(r_min,r_max);

    sig->SetLineColor(kGreen);
    sig->SetNpx(npx);

    //alpha
    sig->SetParameter(0, 50);

    // n
    sig->SetParameter(1, 6);

    //sigma
    sig->SetParameter(2, expected_width);

    //x-
    sig->SetParameter(3, mass);

    //hight
    sig->SetParameter(4, 0.5 * h->GetMaximum());


    TF1* bg = new TF1("bg", "pol6", r_min, r_max);
    bg->SetLineColor(kBlue);

    bg->SetParameter(0,0);
    bg->SetParName(0, "BG p_{0}");
    bg->SetParameter(1,0);
    bg->SetParName(1, "BG p_{1}");
    bg->SetParameter(2,0);
    bg->SetParName(2, "BG p_{2}");
    bg->SetParameter(3,0);
    bg->SetParName(3, "BG p_{3}");
    bg->SetParameter(4,0);
    bg->SetParName(4, "BG p_{4}");
    bg->SetParameter(5,0);
    bg->SetParName(5, "BG p_{5}");
    bg->SetParameter(6,0);
    bg->SetParName(6, "BG p_{6}");

    const auto peak_range = interval<double>::CenterWidth(mass,2*2*mass);

    TFSum::FitRanged(h, bg, r_min, peak_range.Start(), peak_range.Stop(), r_max);


    TFSum* sum = new TFSum("sum", sig, bg, r_min, r_max);
    sum->SetNpx(npx);


    //TCanvas* c = new TCanvas();
    //c->SetTitle(Form("Fit to %s", h->GetName()));
    h->SetStats(true);
    gStyle->SetOptFit(1);
    //h->Draw();
    h->Fit(sum->Function(), "REM0NB");

    sum->SyncToFcts();

    sum->Draw();


    const double total_area = sum->Function()->Integral(r_min, r_max);
    const double bg_area    =  bg->Integral(r_min, r_max);
    const double sig_area   = total_area - bg_area;
    const double peak_pos  = sig->GetParameter(3);
    const double peak_width = sig->GetParameter(2);
    const double relwidth   = peak_width / peak_pos;
    const double mean_error = sig->GetParError(3);

    cout << "Mass offset = " << peak_pos - mass << " MeV\n";
    cout << "sigma/pos   = " <<  relwidth << "\n";
    cout << "Sig         = " << sig_area << endl;

    const string text = std_ext::formatter()
            << "N = " << sig_area << " "
            << "#sigma/#mu = " << relwidth;
    auto l = new TLatex(peak_pos, sum->Function()->Eval(peak_pos), text.c_str());
    l->Draw();


    return FitResult(peak_pos, sig_area, sig->GetParameter(2), (sum->Function()->GetChisquare()) / (sum->Function()->GetNDF() ), sum->Function(), bg, sig, mean_error);

}
Fits::FitResult Fits::FitPeakCrystalBallPol4(TH1* h, const double mass, const double expected_width, const double r_min, const double r_max) {

    const int    npx   = 500;

    TF1* sig = ant::math::CrystalBall::GetTF1();
    //CrystalBall-Function f(x ; alpha , n , sigma, x- ) ;
    //alpha=Positon of the transition between normal destribution and power law
    //n= free parameter of the power law
    //sigma= standard deviation
    //x-= expected value
    sig->SetRange(r_min,r_max);

    sig->SetLineColor(kGreen);
    sig->SetNpx(npx);

    //alpha
    sig->SetParameter(0, 3);

    // n
    sig->SetParameter(1, 6);

    //sigma
    sig->SetParameter(2, expected_width);

    //x-
    sig->SetParameter(3, mass);

    //hight
    sig->SetParameter(4, 0.5 * h->GetMaximum());


    TF1* bg = new TF1("bg", "pol0", r_min, r_max);
    bg->SetLineColor(kBlue);

    bg->SetParameter(0,0);
    bg->SetParName(0, "BG p_{0}");
    bg->SetParameter(1,0);
    bg->SetParName(1, "BG p_{1}");
    bg->SetParameter(2,0);
    bg->SetParName(2, "BG p_{2}");
    bg->SetParameter(3,0);
    bg->SetParName(3, "BG p_{3}");
    bg->SetParameter(4,0);
    bg->SetParName(4, "BG p_{4}");


    const auto peak_range = interval<double>::CenterWidth(mass,2*2*mass);

    TFSum::FitRanged(h, bg, r_min, peak_range.Start(), peak_range.Stop(), r_max);


    TFSum* sum = new TFSum("sum", sig, bg, r_min, r_max);
    sum->SetNpx(npx);


    //TCanvas* c = new TCanvas();
    //c->SetTitle(Form("Fit to %s", h->GetName()));
    h->SetStats(true);
    gStyle->SetOptFit(1);
    //h->Draw();
    h->Fit(sum->Function(), "REM0NB");

    sum->SyncToFcts();

    sum->Draw();


    const double total_area = sum->Function()->Integral(r_min, r_max);
    const double bg_area    =  bg->Integral(r_min, r_max);
    const double sig_area   = total_area - bg_area;
    const double peak_pos  = sig->GetParameter(3);
    const double peak_width = sig->GetParameter(2);
    const double relwidth   = peak_width / peak_pos;
    const double mean_error = sig->GetParError(3);

    cout << "Mass offset = " << peak_pos - mass << " MeV\n";
    cout << "sigma/pos   = " <<  relwidth << "\n";
    cout << "Sig         = " << sig_area << endl;

    const string text = std_ext::formatter()
            << "N = " << sig_area << " "
            << "#sigma/#mu = " << relwidth;
    auto l = new TLatex(peak_pos, sum->Function()->Eval(peak_pos), text.c_str());
    l->Draw();


    return FitResult(peak_pos, sig_area, sig->GetParameter(2), (sum->Function()->GetChisquare()) / (sum->Function()->GetNDF() ), sum->Function(), bg, sig, mean_error);

}

Fits::FitResult Fits::FitPeakCrystalBallPol0(TH1* h, const double mass, const double expected_width, const double r_min, const double r_max) {

    const int    npx   = 500;

    TF1* sig = ant::math::CrystalBall::GetTF1();
    //CrystalBall-Function f(x ; alpha , n , sigma, x- ) ;
    //alpha=Positon of the transition between normal destribution and power law
    //n= free parameter of the power law
    //sigma= standard deviation
    //x-= expected value
    sig->SetRange(r_min,r_max);

    sig->SetLineColor(kGreen);
    sig->SetNpx(npx);

    //alpha
    sig->SetParameter(0, 50);

    // n
    sig->SetParameter(1, 6);

    //sigma
    sig->SetParameter(2, expected_width);

    //x-
    sig->SetParameter(3, mass);

    //hight
    sig->SetParameter(4, 0.5 * h->GetMaximum());


    TF1* bg = new TF1("bg", "pol0", r_min, r_max);
    bg->SetLineColor(kBlue);

    bg->SetParameter(0,0);
//    bg->SetParName(0, "BG p_{0}");



    const auto peak_range = interval<double>::CenterWidth(mass,expected_width * 4);

    TFSum::FitRanged(h, bg, r_min, peak_range.Start(), peak_range.Stop(), r_max);


    TFSum* sum = new TFSum("sum", sig, bg, r_min, r_max);
    sum->SetNpx(npx);


    //TCanvas* c = new TCanvas();
    //c->SetTitle(Form("Fit to %s", h->GetName()));
    h->SetStats(true);
    gStyle->SetOptFit(1);
    //h->Draw();
    h->Fit(sum->Function(), "REM0NB");

    sum->SyncToFcts();

    sum->Draw();


    const double total_area = sum->Function()->Integral(r_min, r_max);
    const double bg_area    =  bg->Integral(r_min, r_max);
    const double sig_area   = total_area - bg_area;
    const double peak_pos  = sig->GetParameter(3);
    const double peak_width = sig->GetParameter(2);
    const double relwidth   = peak_width / peak_pos;
    const double mean_error = sig->GetParError(3);

    cout << "Mass offset = " << peak_pos - mass << " MeV\n";
    cout << "sigma/pos   = " <<  relwidth << "\n";
    cout << "Sig         = " << sig_area << endl;

    const string text = std_ext::formatter()
            << "N = " << sig_area << " "
            << "#sigma/#mu = " << relwidth;
    auto l = new TLatex(peak_pos, sum->Function()->Eval(peak_pos), text.c_str());
    l->Draw();


    return FitResult(peak_pos, sig_area, sig->GetParameter(2), (sum->Function()->GetChisquare()) / (sum->Function()->GetNDF() ), sum->Function(), bg, sig,mean_error);

}

Fits::FitResult Fits::FitPeakPol6(TH1* h, const double mass, const double expected_width, const double r_min, const double r_max) {

    const int    npx   = 500;

    TF1* sig = new TF1("sig", "gaus", r_min, r_max);
    sig->SetLineColor(kGreen);
    sig->SetNpx(npx);

    // height
    sig->SetParameter(0, 0.5 * h->GetMaximum());

    // peak position
    sig->SetParameter(1, mass);

    // width
    sig->SetParameter(2, expected_width);



    TF1* bg = new TF1("bg", "pol6", r_min, r_max);
    bg->SetLineColor(kBlue);

    bg->SetParameter(0,0);
    bg->SetParName(0, "BG p_{0}");
    bg->SetParameter(1,0);
    bg->SetParName(1, "BG p_{1}");
    bg->SetParameter(2,0);
    bg->SetParName(2, "BG p_{2}");
    bg->SetParameter(3,0);
    bg->SetParName(3, "BG p_{3}");
    bg->SetParameter(4,0);
    bg->SetParName(4, "BG p_{4}");
    bg->SetParameter(5,0);
    bg->SetParName(5, "BG p_{5}");
    bg->SetParameter(6,0);
    bg->SetParName(6, "BG p_{6}");

    const auto peak_range = interval<double>::CenterWidth(mass,2*2*mass);

    TFSum::FitRanged(h, bg, r_min, peak_range.Start(), peak_range.Stop(), r_max);


    TFSum* sum = new TFSum("sum", sig, bg, r_min, r_max);
    sum->SetNpx(npx);


    //TCanvas* c = new TCanvas();
    //c->SetTitle(Form("Fit to %s", h->GetName()));
    h->SetStats(true);
    gStyle->SetOptFit(1);
    //h->Draw();
    h->Fit(sum->Function(), "REM0NB");

    sum->SyncToFcts();

    sum->Draw();


    const double total_area = sum->Function()->Integral(r_min, r_max);
    const double bg_area    =  bg->Integral(r_min, r_max);
    const double sig_area   = total_area - bg_area;
    const double peak_pos  = sig->GetParameter(1);
    const double peak_width = sig->GetParameter(2);
    const double relwidth   = peak_width / peak_pos;
    const double mean_error = sig->GetParError(3);

    cout << "Mass offset = " << peak_pos - mass << " MeV\n";
    cout << "sigma/pos   = " <<  relwidth << "\n";
    cout << "Sig         = " << sig_area << endl;

    const string text = std_ext::formatter()
            << "N = " << sig_area << " "
            << "#sigma/#mu = " << relwidth;
    auto l = new TLatex(peak_pos, sum->Function()->Eval(peak_pos), text.c_str());
    l->Draw();


    return FitResult(peak_pos, sig_area, sig->GetParameter(2), 0, sum->Function(), bg, sig, mean_error);

}

void FixZeroBins(TH1* h) {
    for(int i=1;i<=h->GetNbinsX(); ++i) {
        if(h->GetBinContent(i) == 0.0) {
            h->SetBinError(i,1.0);
        }
    }
}

void Fits::FitSlicesPi0(TH2 *h2)
{

    const vector<interval<double>> ranges = {
//        Ranges needed for beamtime data
//        {80, 180},
//        {100, 180},
//        {100, 170},
//        {100, 170},
//        {100, 180},
//        {100, 180},
//        {100,160},
//        {140, 200},
//        {10, 230},
//        {10, 230},
//        {10, 280},
//        {10, 300},
//        {10, 300}

//        Ranges for MC data
        {80, 190},
        {80, 190},
        {80, 190},
        {80, 190},
        {80, 190},
        {80, 190},
        {80, 190},
        {80, 190},
        {80, 190},
        {80, 190},
        {80, 190},
        {80, 190},
        {80, 190}


    };


    gStyle->SetOptStat(0);
    double minEnergy=125;
    double maxEnergy=450;
    TGraph* g1 = new TGraph();
    TGraphErrors* g1_rel = new TGraphErrors();
    int k=0;

    canvas fits(string("Fits for ")+h2->GetTitle());

    for(int i=1; i>0; ++i) {
        TH1* b = h2->ProjectionX(Form("x%d",i),i,i+1);
        double e = h2->GetYaxis()->GetBinCenter(i);
        double elow = h2->GetYaxis()->GetBinLowEdge(i);
        double eup = h2->GetYaxis()->GetBinUpEdge(i);

        if (e <= maxEnergy && e > minEnergy)
        {
            fits << b;

            FixZeroBins(b);
//            auto result = FitPi0Calib(b,ranges.at(k).Start(),ranges.at(k).Stop());  //Fit with crystal ball function
            auto result = FitPi0CalibGaussian(b,ranges.at(k).Start(),ranges.at(k).Stop()); //Fit now with regular Gaussian Function


            fits << samepad << result.bkg << samepad << result.sum << samepad <<result.sig;

            g1->SetPoint(k,e,result.pos);
            g1_rel ->SetPoint(k,e,(result.pos/ParticleTypeDatabase::Pi0.Mass()-1) * 100);
            g1_rel->SetPointError(k,0,(result.position_error/ParticleTypeDatabase::Pi0.Mass()) * 100);
            const string title = std_ext::formatter() << "Energy from " << elow <<" to "<<eup<< " MeV "<<" Chi^2/Ndf: "<<result.chi2dof;
            b->SetTitle(title.c_str());
            k++;
        }

        if(e > maxEnergy){
            break;
        }
    }
    fits << g1 << g1_rel << endc;

    g1->SetTitle("Position of the pi0 peak in different energy intervals of 25 MeV");
    g1->GetXaxis()->SetTitle("Energy of the photons [MeV]");
    g1->GetYaxis()->SetTitle("Position of pi0 peak [MeV]");

    g1_rel->SetTitle("Position of the pi0 peak in different energy intervals of 25 MeV");
    g1_rel->GetXaxis()->SetTitle("Energy of the photons [MeV]");
    g1_rel->GetYaxis()->SetTitle("Deviation from the 135 MeV peak [%] ");
    g1_rel->SetMarkerStyle(21);
    g1_rel->SetMarkerSize(1.5);

}



void Fits::FitSlicesZVertex(TH3 *h3)
{
    gStyle->SetOptStat(0);
    double minEnergy=125;
    double maxEnergy=450;

    int k=0;
    canvas fits(string("Fits for ")+h3->GetTitle());

    auto c_canvas = new TCanvas();
    auto mg = new TMultiGraph();

    const auto mkGrapth = [] (int f) {
        auto c = new TGraphErrors();
        c->SetLineColor(f);
        c->SetFillStyle(0);
        return c;
    };

    int farbe=1;

    for(int z=1;z<= h3->GetNbinsZ();z++)
    {
   double zVer = h3->GetZaxis()->GetBinCenter(z);

        TGraph* g1     = mkGrapth(farbe);
        TGraphErrors* g1_rel = mkGrapth(farbe);
        g1_rel->SetMarkerStyle(20);
        farbe++ ;
        if(farbe==10)
        {
            farbe=1;
        }

        for(int i=1; i>0; ++i) {
            TH1* b = h3->ProjectionX(Form("x%d",i),i,i+1,z,z+1);
            double e = h3->GetYaxis()->GetBinCenter(i);
            double elow = h3->GetYaxis()->GetBinLowEdge(i);
            double eup = h3->GetYaxis()->GetBinUpEdge(i);

            if (e < maxEnergy && e > minEnergy)
            {
                fits << b;


                auto result = FitPi0Calib0(b);
                fits << samepad << result.bkg << samepad << result.sum << samepad <<result.sig;
                g1->SetPoint(k,e,result.pos);
                g1_rel ->SetPoint(k,e,(result.pos/ParticleTypeDatabase::Pi0.Mass()-1) * 100);
                g1_rel->SetPointError(k,0,(result.position_error/ParticleTypeDatabase::Pi0.Mass()) * 100);
                g1_rel->SetMarkerStyle(21);
                const string title = std_ext::formatter() << "Energy from " << elow <<" to "<<eup<< " MeV"<<" Chi^2/Ndf: "<<result.chi2dof;
                b->SetTitle(title.c_str());
                k++;
            }
            else {
                k = 0;
            }
            if(e > maxEnergy){
                break;
            }
        }
        fits << endr;


        const string deviationtitle = std_ext::formatter() << "ZVertexinterval from " << zVer- 0.5<< "cm to "<<zVer + 0.5<< "cm";

        c_canvas->cd();
        g1_rel->Draw();
        g1_rel->SetTitle(deviationtitle.c_str());
        g1_rel->GetXaxis()->SetTitle("Energy of the photons [MeV]");
        g1_rel->GetYaxis()->SetTitle("Deviation from the 135 MeV peak [%] ");
        g1_rel->SetMarkerStyle(20);
        g1_rel->SetMarkerSize(1);

        mg->Add(g1_rel);

    }
    fits << endc;

    auto mg_c = new TCanvas();
    mg->Draw("aLP");
    mg->GetYaxis()->SetTitle("Deviation from the 135 MeV peak [%] ");
    mg->GetXaxis()->SetTitle("Energy of the Photons [MeV]");
    mg_c->SetTitle("Z-Vertex Dependence of the Pi0 Position");


    mg_c->BuildLegend();
}


void Fits::FitSlicesEta(TH2 *h2)
{
    double minEnergy=200;
    double maxEnergy=600;
    TGraph* g1 = new TGraph();
    TGraph* g1_rel = new TGraph();
    int k=0;

    canvas fits(string("Fits for ")+h2->GetTitle());

    for(int i=1; i>0; ++i) {
        TH1* b = h2->ProjectionX(Form("x%d",i),i,i+1);
        double e = h2->GetYaxis()->GetBinCenter(i);
        double elow = h2->GetYaxis()->GetBinLowEdge(i);
        double eup = h2->GetYaxis()->GetBinUpEdge(i);

        if (e <= maxEnergy && e > minEnergy)
        {
            fits << b;

            auto result = FitEtaCalib(b);
            fits << samepad << result.bkg << samepad << result.sum << samepad <<result.sig;
            g1->SetPoint(k,e,result.pos);
            g1_rel ->SetPoint(k,e,(result.pos/ParticleTypeDatabase::Eta.Mass()-1) * 100);
            const string title = std_ext::formatter() << "Energy from " << elow <<" to "<<eup<< "  MeV"<<"Chi^2/Ndf: "<<result.chi2dof;
            b->SetTitle(title.c_str());
            k++;
        }
        else {
            k = 0;
        }
        if(e > maxEnergy){
            break;
        }
    }
    fits << endc;
    new TGraph();
    g1->Draw();
    g1->SetTitle("Position of the Eta peak in different energy intervals of 25 MeV");
    g1->GetXaxis()->SetTitle("Energy of the photons [MeV]");
    g1->GetYaxis()->SetTitle("Position of Eta peak [MeV]");

    new TGraph();
    g1_rel->Draw();
    g1_rel->SetTitle("Position of the Eta peak in different energy intervals of 25 MeV");
    g1_rel->GetXaxis()->SetTitle("Energy of the photons [MeV]");
    g1_rel->GetYaxis()->SetTitle("Deviation from the 548 MeV peak [%] ");
}

