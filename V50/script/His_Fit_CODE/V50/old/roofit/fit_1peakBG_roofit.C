
#include <iostream>
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <string>
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooCmdArg.h"
#include "RooChi2Var.h"

void fit_1peakBG_roofit()
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

    TFile *hf = new TFile(Form("%s/result/RUN%s/Bin%i/his_%s.rrs_on.root", workdir.c_str(), runnumber3, binnum, runnumber6));
    TH1D *his_temp = (TH1D*)hf->Get("his_tot3");

    if (!his_temp) {
        std::cerr << "Histogram not found in file!" << std::endl;
        return;
    }

    TCanvas *cbg = new TCanvas("cbg", "bg fitting", 1000, 800);
    TCanvas *cfit = new TCanvas("cfit", "peak fitting", 1000, 800);

    for (int i = 0; i < N; i++) {

        double energy = peak[i];
        double sigma_cal = energy * (respa[0] + (respa[1] / energy) + (respa[2] / sqrt(energy)));

        int bin_low = his_temp->FindBin(energy - 3 * sigma_cal);
        int bin_high = his_temp->FindBin(energy + 3 * sigma_cal);

// Define background fitting range excluding peak region (energy +/- 5 keV)
RooRealVar x_bkg("x_bkg", "Energy [keV]", energy - 50, energy + 50);

// Define ranges to exclude the peak region
//x_bkg.setRange("lowRange", energy - 50, energy - 5); // Below the peak
//x_bkg.setRange("highRange", energy + 5, energy + 50); // Above the peak
x_bkg.setRange("Range", energy -50, energy + 50); // Above the peak


RooRealVar p0("p0", "Constant", 0, -1e3, 1e3);
RooRealVar p1("p1", "Linear Coefficient", 0, -10, 10); // Added linear coefficient
RooPolynomial bkg("bkg", "Background", x_bkg, RooArgList(p0, p1)); // Use p0 and p1 for background

RooDataHist data_bkg("data_bkg", "Background Data", RooArgList(x_bkg), his_temp);
// Perform the background-only fit in the combined range (lowRange + highRange)
//bkg.fitTo(data_bkg, RooFit::Range("lowRange,highRange")); // Use separate range names
bkg.fitTo(data_bkg, RooFit::Range("Range")); // Use separate range names

	// Plot background fit results
	RooPlot *frame_bkg = x_bkg.frame();
	data_bkg.plotOn(frame_bkg);
	bkg.plotOn(frame_bkg);
	frame_bkg->SetTitle("Background Fit");
	frame_bkg->GetXaxis()->SetTitle("Energy [keV]");
	frame_bkg->GetYaxis()->SetTitle("Counts / keV");

	cbg->cd();
	frame_bkg->Draw();

	// Calculate and print Reduced Chi2 directly from frame
	double reducedChi2_bkg = frame_bkg->chiSquare();
	std::cout << "Reduced Chi2: " << reducedChi2_bkg << std::endl;


	// Get and fix the optimized p0 value
	double fixedP0 = p0.getVal();
	double fixedP1 = p1.getVal();
	p0.setConstant(true); // Fix p0 for the combined fit
	p1.setConstant(true); // Fix p1 for the combined fit

	// Define the fitting range for the signal + background fit
	double x_low = energy - 5;
	double x_high = energy + 5;
	RooRealVar x("x", "Energy [keV]", x_low, x_high);

	// Define the Gaussian signal component
	RooRealVar mean("mean", "Mean", energy, x_low, x_high);
	RooRealVar sigma("sigma", "Sigma", sigma_cal, 0.1 * sigma_cal, 2 * sigma_cal);
	RooRealVar gauss_area("gauss_area", "Gaussian Area", his_temp->Integral(bin_low, bin_high), 0, 1.5 * his_temp->Integral());
	RooGaussian gauss("gauss", "Gaussian", x, mean, sigma);

	// Define the background component with fixed p0
	RooPolynomial bkg_fixed("bkg_fixed", "Background", x, RooArgList(p0, p1));
	RooRealVar bkg_area("bkg_area", "Background Area", 0.1 * his_temp->Integral(bin_low, bin_high), 0, 1.5 * his_temp->Integral());

	// Combine signal and background models
	RooAddPdf model("model", "Signal + Background", RooArgList(gauss, bkg_fixed), RooArgList(gauss_area, bkg_area));

	// Perform the combined fit
	RooDataHist data("data", "Data", RooArgList(x), his_temp);
	model.fitTo(data);

	// Plot the results
	RooPlot *frame = x.frame();
	data.plotOn(frame);
	model.plotOn(frame);
	model.plotOn(frame, RooFit::Components("bkg_fixed"), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));
	model.plotOn(frame, RooFit::Components("gauss"), RooFit::LineStyle(kDotted), RooFit::LineColor(kBlue));

	frame->SetTitle("Source Peaks Fitting");
	frame->GetXaxis()->SetTitle("Energy [keV]");
	frame->GetYaxis()->SetTitle("Counts / keV");

	// Calculate and print Reduced Chi2 directly from frame
	double reducedChi2 = frame->chiSquare();
	std::cout << "Reduced Chi2: " << reducedChi2 << std::endl;
	std::cout << "BG Reduced Chi2: " << reducedChi2_bkg << std::endl;
	cout <<"P1 pa0 "<<fixedP0<<endl;
	cout <<"P1 pa1 "<<fixedP1<<endl;


	cfit->cd();
	frame->Draw();

    }
}
