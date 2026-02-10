#include <iostream>
#include "TMath.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooPolynomial.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"

void fit_peak_multi_study2() {
    // Define observable (energy in keV)
    RooRealVar x("x", "Energy (keV)", 0, 4000);

    // General resolution parameters for sigma calculation
    double respa[3] = {0.000158, 0.547238, 0.006808};

    // Peak positions and initializations
    const int N = 6;
    double peak[N] = {768.36, 772.26, 783.3, 785.4, 785.96, 794.95};
    double peakfix[N] = {767.403, 771.402, 782.723, 784.877, 785.451, 794.67};

    // File and histogram handling
    std::string workdir("/data/HPGe/USERS/kkw/DAQ/ANA125");
    char hisfile[256];
    char runnumber3[256];
    char runnumber6[256];

    int runnum = 624;
    sprintf(runnumber6, "%06d", runnum);
    sprintf(runnumber3, "%03d", runnum);

    int binnum = 8000; // Example default value, can be changed dynamically
    sprintf(hisfile, "%s/result/RUN%s/Bin%i/cal_v1/his_%s.root", workdir.c_str(), runnumber3, binnum, runnumber6);

    TFile* hf = new TFile(hisfile);
    TH1D* his_temp = (TH1D*)hf->Get("his_tot3");

    // Import histogram into RooFit
    RooDataHist data("data", "Dataset from histogram", x, RooFit::Import(*his_temp));

    // Variables for peaks (mean, sigma, and area)
    RooArgList pdfList;
    RooArgList coeffList;

    for (int i = 0; i < N; i++) {
        // Calculate sigma for each peak
        double sigma_cal = peak[i] * (respa[0] + respa[1] / peak[i] + respa[2] / sqrt(peak[i]));

        // Create mean, sigma, and area parameters
        RooRealVar* mean = new RooRealVar(Form("mean%d", i), Form("Mean of peak %d", i), peak[i], peak[i] - 1.5, peak[i] + 1.5);
        RooRealVar* sigma = new RooRealVar(Form("sigma%d", i), Form("Sigma of peak %d", i), sigma_cal, 0.8 * sigma_cal, 3.0 * sigma_cal);
        RooRealVar* area = new RooRealVar(Form("area%d", i), Form("Area of peak %d", i), 1000, 0, 1e6);

        // Create Gaussian PDF for each peak
        RooGaussian* gauss = new RooGaussian(Form("gauss%d", i), Form("Gaussian for peak %d", i), x, *mean, *sigma);

        // Add to PDF and coefficient lists
        pdfList.add(*gauss);
        coeffList.add(*area);
    }

    // Background model (polynomial)
    RooRealVar slope("slope", "Slope of background", 0, -1, 1);
    RooPolynomial background("background", "Background", x, RooArgList(slope));

    // Add background to the PDF list and coefficients
    pdfList.add(background);
    coeffList.add(RooFit::RooConst(1.0)); // Coefficient for background

    // Combine all peaks and background into a single model
    RooAddPdf model("model", "Combined model", pdfList, coeffList);

    // Fit model to data
    RooFitResult* fitResult = model.fitTo(data, RooFit::Save());

    // Extract fit results and output to console
    for (int i = 0; i < N; i++) {
        std::cout << "Peak " << i+1 << ": "
                  << "Mean = " << model.getParameters(data)->getRealValue(Form("mean%d", i)) << ", "
                  << "Sigma = " << model.getParameters(data)->getRealValue(Form("sigma%d", i)) << ", "
                  << "Area = " << model.getParameters(data)->getRealValue(Form("area%d", i)) << std::endl;
    }

    // Print Chi^2 / NDF
    RooPlot* frame = x.frame();
    data.plotOn(frame);
    model.plotOn(frame);
    double chi2 = frame->chiSquare();
    std::cout << "Chi^2 / NDF = " << chi2 << std::endl;
}
