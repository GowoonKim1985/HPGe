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
#include "RooMsgService.h"


void fit_V50_roofit(){

  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    std::string workdir("/data/HPGe/USERS/kkw/DAQ/ANA125");

    char runnumber3[256];
    char runnumber6[256];
    
    int runnum = 624;

    sprintf(runnumber6, "%06d", runnum);
    sprintf(runnumber3, "%03d", runnum);

    int binnum = 8000;

    int const N = 6; // peaks
    //    double peak[N] = {1553.8};
    double peak[N]={768.36, 772.26, 783.3, 785.4, 785.96, 794.95};
    double respa[3] = {0.000158, 0.547238, 0.006808};

    TFile *hf = new TFile(Form("%s/result/RUN%s/Bin%i/his_%s.rrs_on.root", workdir.c_str(), runnumber3, binnum, runnumber6));
    TH1D *his_temp = (TH1D*)hf->Get("his_tot3");

    if (!his_temp) {
        std::cerr << "Histogram not found in file!" << std::endl;
        return;
    }

    TCanvas *cbg = new TCanvas("cbg", "bg fitting", 1000, 800);
    TCanvas *cfit = new TCanvas("cfit", "peak fitting", 1000, 800);

    int bin_low[N], bin_high[N];
    double energy[N], sigma_cal[N];
    for(int i=0; i<N; i++){
        energy[i] = peak[i];
        sigma_cal[i] = energy[i] * (respa[0] + (respa[1] / energy[i]) + (respa[2] / sqrt(energy[i])));
	bin_low[i] = his_temp->FindBin(energy[i] - 3 * sigma_cal[i]);
        bin_high[i] = his_temp->FindBin(energy[i] + 3 * sigma_cal[i]);
    }

    cbg->cd();
    his_temp->GetXaxis()->SetRangeUser(energy[0]-50, energy[N-1]+50);
	his_temp->Draw("hist");
    cfit->cd();
    his_temp->GetXaxis()->SetRangeUser(energy[0]-5, energy[N-1]+5);
	his_temp->Draw("hist");


	// Define background fitting range excluding peak region (energy +/- 5 keV)
	RooRealVar x_bkg("x_bkg", "Energy [keV]", energy[0] - 50, energy[N-1] + 50);

	// Define ranges to exclude the peak region
	//	x_bkg.setRange("BGRange1", energy[0] - 50, energy[0] - 5); // Below the peak
	//	x_bkg.setRange("BGRange2", energy[0] + 5, energy[1] - 5); // Above the peak
	//	x_bkg.setRange("BGRange3", energy[1] +5 , energy[1] + 50); // Below the peak
	x_bkg.setRange("BGRange", energy[0] - 50, energy[N-1] + 50); // Below the peak


	RooRealVar p0("p0", "Constant", 0, -1e3, 1e3);
	RooRealVar p1("p1", "Linear Coefficient", 0, -10, 10); // Added linear coefficient
	//	RooRealVar p0("p0", "Constant", 30, 0, 60);
	//	RooRealVar p1("p1", "Linear Coefficient", 0, -1, 1); // Added linear coefficient
	RooPolynomial bkg("bkg", "Background", x_bkg, RooArgList(p0, p1)); // Use p0 and p1 for background
	RooDataHist data_bkg("data_bkg", "Background Data", RooArgList(x_bkg), his_temp);
	//	bkg.fitTo(data_bkg, RooFit::Range("BGRange1,BGRange2,BGRange3")); // Use separate range names
	bkg.fitTo(data_bkg, RooFit::Range("BGRange"));

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

	RooRealVar x("x", "Energy [keV]", energy[0]-5, energy[N-1]+5);

	RooRealVar* mean[N];
	RooRealVar* sigma[N];
	RooRealVar* gauss_area[N];
	RooGaussian* gauss[N];

	RooArgList gaussList;
	RooArgList gaussAreaList;

	for(int i=0; i<N; i++){
	  mean[i] = new RooRealVar(Form("mean%i",i+1), Form("Mean%i",i+1), energy[i], energy[i]-3*sigma_cal[i], energy[i]+3*sigma_cal[i]);
	  sigma[i] = new RooRealVar(Form("sigma%i",i+1), Form("Sigma%i",i+1), sigma_cal[i], 0.1 * sigma_cal[i], 2 * sigma_cal[i]);
	  gauss_area[i] = new RooRealVar(Form("gauss_area%i",i+1), Form("Gaussian Area%i",i+1), his_temp->Integral(bin_low[i], bin_high[i]), 0, 1.5 * his_temp->Integral());
	  gauss[i] = new RooGaussian(Form("gauss%i",i+1), Form("Gaussian%i",i), x, *mean[i], *sigma[i]);

	  gaussList.add(*gauss[i]);
	  gaussAreaList.add(*gauss_area[i]);
}

	// Define the background component with fixed p0
	RooPolynomial bkg_fixed("bkg_fixed", "Background", x, RooArgList(p0, p1));
	RooRealVar bkg_area("bkg_area", "Background Area", 0.1 * his_temp->Integral(bin_low[0], bin_high[N-1]), 0, 1.5 * his_temp->Integral());

	  gaussList.add(bkg_fixed);
	  gaussAreaList.add(bkg_area);

	  RooAddPdf model("model", "M-Signal + Background", gaussList, gaussAreaList);

	// Perform the combined fit
	RooDataHist data("data", "Data", RooArgList(x), his_temp);
	model.fitTo(data);

	// Plot the results
	RooPlot *frame = x.frame();
	data.plotOn(frame, RooFit::Invisible());

	model.plotOn(frame);
	model.plotOn(frame, RooFit::Components("bkg_fixed"), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));

	for(int i=0; i<N; i++){
	model.plotOn(frame,
		     RooFit::Components(Form("gauss%i",i+1)),
		     RooFit::LineStyle(kDotted), 
		     RooFit::LineColor(kBlue));
	}

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
	frame->SetMinimum(energy[0]-5);
	frame->SetMinimum(energy[N-1]+5);
	frame->Draw("SAME");
}
