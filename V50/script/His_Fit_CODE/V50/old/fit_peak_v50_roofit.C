#include "TMath.h"
#include "TMinuit.h"
#include <string>
#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooPolynomial.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <iostream>

void fit_peak_v50_roofit()
{   
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // ERROR 이상만 출력

    std::string workdir("/data/HPGe/USERS/kkw/DAQ/ANA125");

    char hisfile[256];
    char rawname[256];
    char resfile[256]; // fit result

    char runnumber3[256];
    char runnumber6[256];
    char bintag[256];
    char peakfile[258];

    int runnum = 624;

    sprintf(runnumber6, "%06d", runnum);
    sprintf(runnumber3, "%03d", runnum);

    int binnum = 8000; // 1bin = 1keV
    int bw = binnum / binnum; // bin width, keV

    int const N = 1; // peaks
    double peak[N] = {1553.8};

    // General resolution
    double respa[3];
    respa[0] = 0.000158;
    respa[1] = 0.547238;
    respa[2] = 0.006808;

    sprintf(hisfile, "%s/result/RUN%s/Bin%i/his_%s.rrs_on.root", workdir.c_str(), runnumber3, binnum, runnumber6);
    sprintf(resfile, "%s/result/RUN%s/Bin%i/res_%s.rrs_on.root", workdir.c_str(), runnumber3, binnum, runnumber6);

    TFile* hf = new TFile(hisfile);
    TH1D* his_temp = (TH1D*)hf->Get("his_tot3"); // 히스토그램 로드
    his_temp->SetName("his_temp");
    his_temp->SetLineColor(1);

    // 히스토그램 정규화
    his_temp->Scale(1.0 / his_temp->GetBinWidth(1));

    for (int i = 0; i < N; i++) {
        double energy = peak[i];
        double sigma_cal = peak[i] * (respa[0] + (respa[1] / peak[i]) + (respa[2] / (sqrt(peak[i]))));
        std::cout << "Cal. Sigma : " << sigma_cal << " (3 sigma : " << 3 * sigma_cal << ")\n";
        std::cout << "energy : " << energy << "\n";

        RooRealVar x("x", "Energy (keV)", energy - 20, energy + 20); // fitting range 설정
        RooDataHist data("data", "Dataset", x, his_temp); // RooFit 데이터

    double initial_area = his_temp->Integral(his_temp->FindBin(energy - 30), his_temp->FindBin(energy + 30));
    cout<<"initial_area "<<initial_area<<endl;

        RooRealVar mean1("mean1", "Mean of Gaussian", energy, energy - 5, energy + 5);
	//	RooRealVar area1("area1", "Area of Gaussian 1", initial_area, 0, 1.5 * initial_area);
        RooRealVar sigma1("sigma1", "Sigma of Gaussian 1", 0.5, 0.1, 1.5);
	RooGaussian gauss1("gauss1", "Gaussian Signa1", x, mean1, sigma1);
	
	RooRealVar area1("area1", "Area of Gaussian Signal", initial_area/2, 0, initial_area); // 면적
	RooAddPdf model("model", "Signal Model", RooArgList(gauss1), RooArgList(area1));



	//    gauss1.fitTo(data, RooFit::Range(energy - 5, energy + 5), RooFit::Strategy(2)); 

    RooFitResult* fitResult = model.fitTo(data, RooFit::Save(), RooFit::Range(energy - 5, energy + 5));

    // 면적값 및 에러 출력
    double fitted_area = area1.getVal();       // 피팅된 면적값
    double fitted_area_error = area1.getError(); // 면적값의 에러

    std::cout << "Fitted Gaussian Area: " << fitted_area << " ± " << fitted_area_error << std::endl;
    std::cout << "Gaussian Mean: " << mean1.getVal() << " ± " << mean1.getError() << std::endl;
    std::cout << "Gaussian Sigma: " << sigma1.getVal() << " ± " << sigma1.getError() << std::endl;


    TCanvas *c_gauss = new TCanvas("c_gauss", "Gaussian Fit Result", 800, 600); // 가우시
    RooPlot *frame_gauss = x.frame(); 
    TCanvas *c = new TCanvas("c", "Gaussian Fit Result", 800, 600); // 캔버스 생성
    RooPlot *frame = x.frame(); // 플롯 범위 설정
    data.plotOn(frame); // 히스토그램 데이터 플롯
    model.plotOn(frame); // 전체 모델 플롯
    model.plotOn(frame, RooFit::Components("gauss1"), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed)); // 가우시안 신호 플롯
    frame->SetTitle("Gaussian Fit with RooFit");
    frame->Draw(); // 플롯 그리기



    //    std::cout << "Gaussian 1 Mean (fitted): " << mean1.getVal() << " ± " << mean1.getError() << std::endl;
    }

}


