#include <iostream>
#include "TMath.h"
#include "TMinuit.h"
#include <string>

void fit_peak_v50_study()
{
    std::string workdir("/data/HPGe/USERS/kkw/DAQ/ANA125");

    char runnumber3[256];
    char runnumber6[256];
    
    int runnum = 624;

    sprintf(runnumber6, "%06d", runnum);
    sprintf(runnumber3, "%03d", runnum);

    int binnum = 8000;

    int const N = 1; // peaks
    double peak[N] = {1553.8};

    double respa[3] = {0.000158, 0.547238, 0.006808};

    TF1 *gfit = new TF1("gfit", "gaus");
    TF1 *agfit = new TF1("agfit", "[0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2])");
    agfit->SetParNames("Area", "Mean", "Sigma"); 

    TF1 *p1fit = new TF1("p1fit", "pol1");
    TF1 *fit = new TF1("fit", "([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+[3]+[4]*x");
    fit->SetParNames("Area", "Mean", "Sigma", "Intercept", "Slope");

    TCanvas *c1 = new TCanvas("c1", "fitting", 1200, 800);

    fit->SetLineColor(kGreen);
    p1fit->SetLineColor(kGreen);
    p1fit->SetLineStyle(7);

    TFile *hf = new TFile(Form("%s/result/RUN%s/Bin%i/his_%s.rrs_on.root", workdir.c_str(), runnumber3, binnum, runnumber6));
    TH1D *his_temp = (TH1D*)hf->Get("his_tot3");

    his_temp->Draw();
    his_temp->SetName("his_temp");
    his_temp->SetLineColor(1);

    double sigma_cal;

    for (int i = 0; i < N; i++) {
        c1->cd();
        double energy = peak[i];
        int xbin = his_temp->FindBin(peak[i]);

        printf("\n=========================\n");
        printf("%.2f keV bin : %i\n", peak[i], xbin);

        sigma_cal = peak[i] * (respa[0] + (respa[1] / peak[i]) + (respa[2] / sqrt(peak[i])));
        printf("Cal. Sigma : %f (3 sigma : %f)\n", sigma_cal, 3 * sigma_cal);

        his_temp->GetXaxis()->SetRange(his_temp->FindBin(peak[i] - (3 * sigma_cal)), his_temp->FindBin(peak[i] + (3 * sigma_cal)));

        int maxb = his_temp->GetMaximumBin();
        double maxx = his_temp->GetBinCenter(maxb);
        printf("maxb(bin) : %i\n", maxb);
        printf("maxx(kev) : %f\n", maxx);

        his_temp->GetXaxis()->SetRange(0, binnum);
        gfit->SetRange(maxx - 3 * sigma_cal, maxx + 3 * sigma_cal);
        agfit->SetRange(maxx - 3 * sigma_cal, maxx + 3 * sigma_cal);
        fit->SetRange(maxx - 6 * sigma_cal, maxx + 6 * sigma_cal);

        his_temp->Fit("p1fit", "RQ+");
        double fixbg[2] = {p1fit->GetParameter(0), p1fit->GetParameter(1)};

        his_temp->Fit("gfit", "RQ0+");
        double gfitpa[3] = {gfit->GetParameter(0), gfit->GetParameter(1), gfit->GetParameter(2)};

        agfit->SetParameters(2.5 * gfitpa[1] * gfitpa[2], peak[i], sigma_cal);
        his_temp->Fit("agfit", "R0Q+");

        fit->SetParameter(0, agfit->GetParameter(0));
        fit->SetParameter(1, maxx);
        fit->SetParLimits(1, maxx - 3 * sigma_cal, maxx + 3 * sigma_cal);
        fit->SetParameter(2, sigma_cal);
        fit->SetParameter(3, fixbg[0]);
        fit->SetParameter(4, fixbg[1]);

        his_temp->Fit("fit", "R+");

        double fitpa[5] = {fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2), fit->GetParameter(3), fit->GetParameter(4)};
        double fitchi2 = fit->GetChisquare() / fit->GetNDF();

        printf("\nfit chi^2/NDF = %f\n", fitchi2);
    }

    his_temp->SetTitle("source peaks Fitting");
    his_temp->SetXTitle("Energy[keV]");
    his_temp->SetYTitle("Counts/kev");

    c1->cd();
    his_temp->GetXaxis()->SetRangeUser(1535, 1575);

    his_temp->Draw();

    TFile rf(Form("%s/result/RUN%s/Bin%i/res_%s.rrs_on.root", workdir.c_str(), runnumber3, binnum, runnumber6), "RECREATE");
    rf.cd();
    c1->Write();
    rf.Close();
}
