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

void fit_prac() {
    // histogram setting
    const int bins = 8000;
    const double xMin = 0;
    const double xMax = 4000;
    TH1D* h = new TH1D("h", "Gamma Energy Spectrum", bins, xMin, xMax);

    
    for (int i = 0; i < 1e6; ++i) {
        h->Fill(gRandom->Gaus(530, 5));
        h->Fill(gRandom->Gaus(560, 7));
        h->Fill(gRandom->Uniform(xMin, xMax));
    }

    // RooFit 변수 설정
    RooRealVar x("x", "Energy (keV)", 500, 600);

    // 히스토그램을 RooDataHist로 변환
    RooDataHist data("data", "Dataset", x, h);

    // 피팅 변수 정의
    RooRealVar mean1("mean1", "Mean of Gaussian 1", 530, 520, 540);
    RooRealVar sigma1("sigma1", "Sigma of Gaussian 1", 5, 1, 10);
    RooRealVar area1("area1", "Area of Gaussian 1", 1000, 0, 1e6);

    RooRealVar mean2("mean2", "Mean of Gaussian 2", 560, 550, 570);
    RooRealVar sigma2("sigma2", "Sigma of Gaussian 2", 7, 1, 15);
    RooRealVar area2("area2", "Area of Gaussian 2", 1000, 0, 1e6);

    RooRealVar a("a", "Linear coefficient", 0, -10, 10);
    RooRealVar b("b", "Constant coefficient", 0, -1e4, 1e4);

    // 가우시안과 백그라운드 정의
    RooGaussian gauss1("gauss1", "Gaussian 1", x, mean1, sigma1);
    RooGaussian gauss2("gauss2", "Gaussian 2", x, mean2, sigma2);
    RooPolynomial background("background", "Background", x, RooArgList(a, b));

    // 전체 피팅 함수 정의
    RooAddPdf model("model", "Total Fit Model",
                    RooArgList(gauss1, gauss2, background),
                    RooArgList(area1, area2));

    // 피팅 수행
    model.fitTo(data, RooFit::Range(500, 600));

    // 결과 플롯
    TCanvas* c = new TCanvas("c", "Fit Result", 800, 600);
    RooPlot* frame = x.frame();
    data.plotOn(frame);
    model.plotOn(frame);
    model.plotOn(frame, RooFit::Components("gauss1"), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));
    model.plotOn(frame, RooFit::Components("gauss2"), RooFit::LineStyle(kDashed), RooFit::LineColor(kBlue));
    model.plotOn(frame, RooFit::Components("background"), RooFit::LineStyle(kDashed), RooFit::LineColor(kGreen));
    frame->Draw();

    // 가우시안 면적 계산 및 출력
    double binWidth = (xMax - xMin) / bins; // 히스토그램의 bin width 계산
    double realArea1 = area1.getVal() / binWidth;
    double realArea2 = area2.getVal() / binWidth;

    double realArea1Error = area1.getError() / binWidth;
    double realArea2Error = area2.getError() / binWidth;

    std::cout << "Gaussian 1 Area (real): " << realArea1 << " \u00b1 " << realArea1Error << "\n";
    std::cout << "Gaussian 2 Area (real): " << realArea2 << " \u00b1 " << realArea2Error << "\n";
}
