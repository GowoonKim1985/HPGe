#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TROOT.h"

// RooFit headers
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooConstVar.h"

using namespace RooFit;

void step3_likelihood_fitting() {
    gStyle->SetOptStat(0);
    
    std::cout << "=== RooFit Likelihood Fitting Analysis ===" << std::endl;
    std::cout << "Fitting Be-7 peak (478 keV) with background peaks" << std::endl;
    
    // Load pseudo data
    std::cout << "\n1. Loading pseudo data..." << std::endl;
    TFile* data_file = new TFile("pseudo_data.root", "READ");
    if (!data_file || !data_file->IsOpen()) {
        std::cerr << "Error: Cannot open pseudo_data.root!" << std::endl;
        return;
    }
    
    TH1D* h_data = (TH1D*)data_file->Get("h_pseudo_data");
    TTree* event_tree = (TTree*)data_file->Get("event_tree");
    
    if (!h_data || !event_tree) {
        std::cerr << "Error: Cannot find required data objects!" << std::endl;
        return;
    }
    
    std::cout << "Data loaded successfully:" << std::endl;
    std::cout << "  Histogram bins: " << h_data->GetNbinsX() << std::endl;
    std::cout << "  Tree entries: " << event_tree->GetEntries() << std::endl;
    
    // Load resolution function
    std::cout << "\n2. Loading resolution function..." << std::endl;
    TFile* res_file = new TFile("resolution_results.root", "READ");
    TF1* resolution_func = nullptr;
    
    if (res_file && res_file->IsOpen()) {
        resolution_func = (TF1*)res_file->Get("resolution_function");
    }
    
    if (!resolution_func) {
        std::cout << "Warning: Using default resolution parameters..." << std::endl;
        resolution_func = new TF1("resolution_func_default", "[0]/sqrt(x) + [1] + [2]*x", 300, 2700);
        resolution_func->SetParameters(20.0, 0.5, 0.0);
    }
    
    // Calculate resolutions for each peak
    std::vector<double> peak_energies = {463.0, 478.0, 481.0, 486.0};
    std::vector<double> peak_sigmas;
    
    std::cout << "Peak resolutions:" << std::endl;
    for (double energy : peak_energies) {
        double resolution = resolution_func->Eval(energy); // %
        double fwhm = (resolution / 100.0) * energy; // keV
        double sigma = fwhm / 2.35; // keV
        peak_sigmas.push_back(sigma);
        std::cout << "  " << energy << " keV: #sigma = " << sigma << " keV" << std::endl;
    }
    
    // Define fitting range
    double fit_min = 456.0;
    double fit_max = 500.0;
    
    std::cout << "\n3. Setting up RooFit model..." << std::endl;
    std::cout << "Fitting range: " << fit_min << " - " << fit_max << " keV" << std::endl;
    
    // === BINNED LIKELIHOOD FITTING ===
    std::cout << "\n4. Binned likelihood fitting..." << std::endl;
    
    // Create RooDataHist from histogram (restrict to fitting range)
    int bin_min = h_data->FindBin(fit_min);
    int bin_max = h_data->FindBin(fit_max);
    int n_bins_fit = bin_max - bin_min + 1;  // 피팅에 사용된 빈 개수 계산
    
    std::cout << "  Fitting range bins: " << n_bins_fit << std::endl;
    
    // Create reduced histogram for fitting range
    TH1D* h_fit = new TH1D("h_fit", "Data for fitting", 
                          bin_max - bin_min + 1,
                          h_data->GetBinLowEdge(bin_min),
                          h_data->GetBinLowEdge(bin_max + 1));
    
    for (int i = bin_min; i <= bin_max; ++i) {
        int new_bin = i - bin_min + 1;
        h_fit->SetBinContent(new_bin, h_data->GetBinContent(i));
        h_fit->SetBinError(new_bin, h_data->GetBinError(i));
    }
    
    // RooFit Variables for BINNED fit
    RooRealVar energy_b("energy_b", "Energy (keV)", fit_min, fit_max);
    
    // Peak positions (1 sigma 범위로 변동 가능)
    RooRealVar mean_463_b("mean_463_b", "463 keV peak", 463.0, 463.0-peak_sigmas[0], 463.0+peak_sigmas[0]);
    RooRealVar mean_478_b("mean_478_b", "478 keV peak (Be-7)", 478.0, 478.0-peak_sigmas[1], 478.0+peak_sigmas[1]);
    RooRealVar mean_481_b("mean_481_b", "481 keV peak", 481.0, 481.0-peak_sigmas[2], 481.0+peak_sigmas[2]);
    RooRealVar mean_486_b("mean_486_b", "486 keV peak", 486.0, 486.0-peak_sigmas[3], 486.0+peak_sigmas[3]);
    
    // Peak widths (0.5sigma ~ 1.5sigma 범위로 변동 가능)
    RooRealVar sigma_463_b("sigma_463_b", "463 keV width", peak_sigmas[0], 0.5*peak_sigmas[0], 1.5*peak_sigmas[0]);
    RooRealVar sigma_478_b("sigma_478_b", "478 keV width", peak_sigmas[1], 0.5*peak_sigmas[1], 1.5*peak_sigmas[1]);
    RooRealVar sigma_481_b("sigma_481_b", "481 keV width", peak_sigmas[2], 0.5*peak_sigmas[2], 1.5*peak_sigmas[2]);
    RooRealVar sigma_486_b("sigma_486_b", "486 keV width", peak_sigmas[3], 0.5*peak_sigmas[3], 1.5*peak_sigmas[3]);
    
    // Peak amplitudes (free parameters) - 상한을 3배로 증가
    RooRealVar n_463_b("n_463_b", "463 keV amplitude", 10, 0, 3000);
    RooRealVar n_478_b("n_478_b", "478 keV amplitude (Be-7)", 50, 0, 1500);
    RooRealVar n_481_b("n_481_b", "481 keV amplitude", 10, 0, 3000);
    RooRealVar n_486_b("n_486_b", "486 keV amplitude", 10, 0, 3000);
    
    // Linear background parameters
    RooRealVar n_bkg_b("n_bkg_b", "Background amplitude", 100, 0, 10000);
    RooRealVar a0_b("a0_b", "Background constant", 100, 0, 10000);
    RooRealVar a1_b("a1_b", "Background slope", -1, -10, 10);
    
    // PDF Components for BINNED
    RooGaussian gauss_463_b("gauss_463_b", "463 keV peak", energy_b, mean_463_b, sigma_463_b);
    RooGaussian gauss_478_b("gauss_478_b", "478 keV peak (Be-7)", energy_b, mean_478_b, sigma_478_b);
    RooGaussian gauss_481_b("gauss_481_b", "481 keV peak", energy_b, mean_481_b, sigma_481_b);
    RooGaussian gauss_486_b("gauss_486_b", "486 keV peak", energy_b, mean_486_b, sigma_486_b);
    
    // Linear background
    RooPolynomial background_b("background_b", "Linear background", energy_b, RooArgList(a0_b, a1_b));
    
    // Total model for BINNED
    RooAddPdf model_b("model_b", "4 Gaussians + Linear background",
                     RooArgList(gauss_463_b, gauss_478_b, gauss_481_b, gauss_486_b, background_b),
                     RooArgList(n_463_b, n_478_b, n_481_b, n_486_b, n_bkg_b));
    
    RooDataHist data_hist("data_hist", "Binned data", energy_b, h_fit);
    
    std::cout << "  Data points in fit range: " << data_hist.sumEntries() << std::endl;
    
    // Perform binned fit
    RooFitResult* fit_result_binned = model_b.fitTo(data_hist, Save(true), PrintLevel(-1));
    
    std::cout << "  Binned fit completed!" << std::endl;
    std::cout << "  Fit status: " << fit_result_binned->status() << std::endl;
    std::cout << "  -log(L): " << fit_result_binned->minNll() << std::endl;
    
    // Store binned fit results
    double be7_binned = n_478_b.getVal();
    double be7_err_binned = n_478_b.getError();
    double nll_binned = fit_result_binned->minNll();
    
    // === UNBINNED LIKELIHOOD FITTING ===
    std::cout << "\n5. Unbinned likelihood fitting..." << std::endl;
    
    // 빈드 데이터에서 언빈드 데이터 생성 (빈센터 x 빈컨텐츠)
    std::vector<double> unbinned_data;
    for (int i = bin_min; i <= bin_max; ++i) {
        double bin_center = h_data->GetBinCenter(i);
        int bin_content = static_cast<int>(h_data->GetBinContent(i));
        
        // 빈 컨텐츠만큼 빈센터 에너지값 반복
        for (int j = 0; j < bin_content; ++j) {
            unbinned_data.push_back(bin_center);
        }
    }
    
    std::cout << "  Created unbinned dataset with " << unbinned_data.size() << " events" << std::endl;
    
    // RooFit Variables for UNBINNED fit (독립적인 변수들)
    RooRealVar energy_u("energy_u", "Energy (keV)", fit_min, fit_max);
    
    // Peak positions (1 sigma 범위로 변동 가능)
    RooRealVar mean_463_u("mean_463_u", "463 keV peak", 463.0, 463.0-peak_sigmas[0], 463.0+peak_sigmas[0]);
    RooRealVar mean_478_u("mean_478_u", "478 keV peak (Be-7)", 478.0, 478.0-peak_sigmas[1], 478.0+peak_sigmas[1]);
    RooRealVar mean_481_u("mean_481_u", "481 keV peak", 481.0, 481.0-peak_sigmas[2], 481.0+peak_sigmas[2]);
    RooRealVar mean_486_u("mean_486_u", "486 keV peak", 486.0, 486.0-peak_sigmas[3], 486.0+peak_sigmas[3]);
    
    // Peak widths (0.5sigma ~ 1.5sigma 범위로 변동 가능)
    RooRealVar sigma_463_u("sigma_463_u", "463 keV width", peak_sigmas[0], 0.5*peak_sigmas[0], 1.5*peak_sigmas[0]);
    RooRealVar sigma_478_u("sigma_478_u", "478 keV width", peak_sigmas[1], 0.5*peak_sigmas[1], 1.5*peak_sigmas[1]);
    RooRealVar sigma_481_u("sigma_481_u", "481 keV width", peak_sigmas[2], 0.5*peak_sigmas[2], 1.5*peak_sigmas[2]);
    RooRealVar sigma_486_u("sigma_486_u", "486 keV width", peak_sigmas[3], 0.5*peak_sigmas[3], 1.5*peak_sigmas[3]);
    
    // Peak fractions (free parameters for unbinned)
    RooRealVar frac_463_u("frac_463_u", "463 keV fraction", 0.1, 0.0, 1.0);
    RooRealVar frac_478_u("frac_478_u", "478 keV fraction (Be-7)", 0.2, 0.0, 1.0);
    RooRealVar frac_481_u("frac_481_u", "481 keV fraction", 0.1, 0.0, 1.0);
    RooRealVar frac_486_u("frac_486_u", "486 keV fraction", 0.1, 0.0, 1.0);
    
    // PDF Components for UNBINNED
    RooGaussian gauss_463_u("gauss_463_u", "463 keV peak", energy_u, mean_463_u, sigma_463_u);
    RooGaussian gauss_478_u("gauss_478_u", "478 keV peak (Be-7)", energy_u, mean_478_u, sigma_478_u);
    RooGaussian gauss_481_u("gauss_481_u", "481 keV peak", energy_u, mean_481_u, sigma_481_u);
    RooGaussian gauss_486_u("gauss_486_u", "486 keV peak", energy_u, mean_486_u, sigma_486_u);
    
    // Linear background parameters for unbinned
    RooRealVar a0_u("a0_u", "Background constant", 100, 0, 10000);
    RooRealVar a1_u("a1_u", "Background slope", -1, -10, 10);
    RooPolynomial background_u("background_u", "Linear background", energy_u, RooArgList(a0_u, a1_u));
    
    // Total model for UNBINNED
    RooAddPdf model_u("model_u", "4 Gaussians + Linear background",
                     RooArgList(gauss_463_u, gauss_478_u, gauss_481_u, gauss_486_u, background_u),
                     RooArgList(frac_463_u, frac_478_u, frac_481_u, frac_486_u));
    
    // RooDataSet 생성
    RooDataSet data_unbinned("data_unbinned", "Unbinned data", RooArgSet(energy_u));
    
    for(double e : unbinned_data) {
        energy_u.setVal(e);
        data_unbinned.add(RooArgSet(energy_u));
    }
    
    std::cout << "  Data points in fit range: " << data_unbinned.sumEntries() << std::endl;
    
    // Perform unbinned fit
    RooFitResult* fit_result_unbinned = model_u.fitTo(data_unbinned, Save(true), PrintLevel(-1));
    
    std::cout << "  Unbinned fit completed!" << std::endl;
    std::cout << "  Fit status: " << fit_result_unbinned->status() << std::endl;
    std::cout << "  -log(L): " << fit_result_unbinned->minNll() << std::endl;
    
    // Store unbinned fit results (fraction * total events)
    double total_events = data_unbinned.sumEntries();
    double be7_unbinned = frac_478_u.getVal() * total_events;
    double be7_err_unbinned = frac_478_u.getError() * total_events;
    double nll_unbinned = fit_result_unbinned->minNll();
    
    // === RESULTS COMPARISON ===
    std::cout << "\n6. Fit results comparison:" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Parameter        Binned        Unbinned" << std::endl;
    std::cout << "========================================" << std::endl;
    
    std::cout << Form("Be-7 (478 keV)   %.1f#pm%.1f     %.1f#pm%.1f", 
                     be7_binned, be7_err_binned, be7_unbinned, be7_err_unbinned) << std::endl;
    
    // Print fitted peak positions and widths
    std::cout << "\nFitted peak positions (Binned):" << std::endl;
    std::cout << Form("  463 keV: %.2f #pm %.2f keV", mean_463_b.getVal(), mean_463_b.getError()) << std::endl;
    std::cout << Form("  478 keV: %.2f #pm %.2f keV", mean_478_b.getVal(), mean_478_b.getError()) << std::endl;
    std::cout << Form("  481 keV: %.2f #pm %.2f keV", mean_481_b.getVal(), mean_481_b.getError()) << std::endl;
    std::cout << Form("  486 keV: %.2f #pm %.2f keV", mean_486_b.getVal(), mean_486_b.getError()) << std::endl;
    
    std::cout << "\nFitted peak widths (Binned):" << std::endl;
    std::cout << Form("  463 keV: %.2f #pm %.2f keV", sigma_463_b.getVal(), sigma_463_b.getError()) << std::endl;
    std::cout << Form("  478 keV: %.2f #pm %.2f keV", sigma_478_b.getVal(), sigma_478_b.getError()) << std::endl;
    std::cout << Form("  481 keV: %.2f #pm %.2f keV", sigma_481_b.getVal(), sigma_481_b.getError()) << std::endl;
    std::cout << Form("  486 keV: %.2f #pm %.2f keV", sigma_486_b.getVal(), sigma_486_b.getError()) << std::endl;
    
    // === PLOTTING ===
    std::cout << "\n7. Creating fit plots..." << std::endl;
    
    TCanvas* c_fits = new TCanvas("c_fits", "Likelihood Fits Comparison", 1400, 700);
    c_fits->Divide(2, 1);
    
    // Plot 1: Binned fit
    c_fits->cd(1);
    RooPlot* frame_binned = energy_b.frame(Title("Binned Likelihood Fit"));
    data_hist.plotOn(frame_binned, MarkerStyle(20), MarkerSize(0.8));
    
    model_b.plotOn(frame_binned, LineColor(kRed), LineWidth(2));
    
    // Plot individual components
    model_b.plotOn(frame_binned, Components(gauss_478_b), LineColor(kGreen+2), LineStyle(2));
    model_b.plotOn(frame_binned, Components(RooArgSet(gauss_463_b, gauss_481_b, gauss_486_b)), 
                LineColor(kMagenta), LineStyle(2));
    model_b.plotOn(frame_binned, Components(background_b), LineColor(kBlue), LineStyle(3));
    
    frame_binned->Draw();
    
    // Add legend
    TLegend* leg1 = new TLegend(0.6, 0.6, 0.9, 0.9);
    leg1->AddEntry((TObject*)0, "Binned Fit", "");
    leg1->AddEntry((TObject*)0, Form("Be-7: %.1f#pm%.1f", be7_binned, be7_err_binned), "");
    leg1->AddEntry((TObject*)0, Form("-log(L): %.1f", nll_binned), "");
    leg1->SetTextSize(0.03);
    leg1->Draw();
    
    // Plot 2: Unbinned fit
    c_fits->cd(2);
    RooPlot* frame_unbinned = energy_u.frame(Title("Unbinned Likelihood Fit"));
    
    // 데이터 플롯 (빈드와 같은 빈 개수 사용)
    data_unbinned.plotOn(frame_unbinned, Binning(n_bins_fit), MarkerStyle(20), MarkerSize(0.8)); 
    
    model_u.plotOn(frame_unbinned, LineColor(kRed), LineWidth(2));
    
    // 각 컴포넌트별 플롯
    model_u.plotOn(frame_unbinned, Components(gauss_478_u), LineColor(kGreen+2), LineStyle(2));
    model_u.plotOn(frame_unbinned, Components(RooArgSet(gauss_463_u, gauss_481_u, gauss_486_u)), 
                LineColor(kMagenta), LineStyle(2));
    model_u.plotOn(frame_unbinned, Components(background_u), LineColor(kBlue), LineStyle(3));
    
    // 전체 모델 다시 그리기 (맨 위에 표시)
    model_u.plotOn(frame_unbinned, LineColor(kRed), LineWidth(2));
    
    frame_unbinned->Draw();
    
    // Add legend
    TLegend* leg2 = new TLegend(0.6, 0.6, 0.9, 0.9);
    leg2->AddEntry((TObject*)0, "Unbinned Fit", "");
    leg2->AddEntry((TObject*)0, Form("Be-7: %.1f#pm%.1f", be7_unbinned, be7_err_unbinned), "");
    leg2->AddEntry((TObject*)0, Form("-log(L): %.1f", nll_unbinned), "");
    leg2->SetTextSize(0.03);
    leg2->Draw();
    
    c_fits->Update();
    
    // === SAVE RESULTS ===
    std::cout << "\n8. Saving results..." << std::endl;
    
    TFile* outfile = new TFile("likelihood_results.root", "RECREATE");
    if (outfile->IsOpen()) {
        // Save fit results
        fit_result_binned->Write("fit_result_binned");
        fit_result_unbinned->Write("fit_result_unbinned");
        
        // Save plots
        c_fits->Write();
        frame_binned->Write("frame_binned");
        frame_unbinned->Write("frame_unbinned");
        
        // Save data
        h_fit->Write("h_data_fit_range");
        
        // Save summary
        TPaveText* summary = new TPaveText(0, 0, 1, 1);
        summary->SetName("fit_summary");
        summary->AddText("=== Likelihood Fitting Results ===");
        summary->AddText(Form("Fit range: %.0f - %.0f keV", fit_min, fit_max));
        summary->AddText("");
        summary->AddText("Be-7 Peak (478 keV) Results:");
        summary->AddText(Form("Binned:   %.1f #pm %.1f counts", be7_binned, be7_err_binned));
        summary->AddText(Form("Unbinned: %.1f #pm %.1f counts", be7_unbinned, be7_err_unbinned));
        summary->AddText("");
        summary->AddText(Form("Binned -log(L):   %.2f", nll_binned));
        summary->AddText(Form("Unbinned -log(L): %.2f", nll_unbinned));
        summary->AddText("");
        summary->AddText("Fitted Peak Positions (Binned):");
        summary->AddText(Form("463: %.2f#pm%.2f, 478: %.2f#pm%.2f", 
                             mean_463_b.getVal(), mean_463_b.getError(),
                             mean_478_b.getVal(), mean_478_b.getError()));
        summary->AddText(Form("481: %.2f#pm%.2f, 486: %.2f#pm%.2f", 
                             mean_481_b.getVal(), mean_481_b.getError(),
                             mean_486_b.getVal(), mean_486_b.getError()));
        summary->Write();
        
        std::cout << "Results saved to 'likelihood_results.root'" << std::endl;
        outfile->Close();
    }
    
    std::cout << "\n=== Likelihood Fitting Analysis Completed! ===" << std::endl;
    std::cout << "Be-7 detection results:" << std::endl;
    std::cout << "  Binned likelihood:   " << be7_binned << " #pm " << be7_err_binned << " counts" << std::endl;
    std::cout << "  Unbinned likelihood: " << be7_unbinned << " #pm " << be7_err_unbinned << " counts" << std::endl;

    // 추가: 피팅된 가우시안 함수들의 평균값과 시그마 출력
    std::cout << "\n=== Fitted Gaussian Parameters ===" << std::endl;
    std::cout << "BINNED FIT RESULTS:" << std::endl;
    std::cout << "Peak Positions (Mean values):" << std::endl;
    std::cout << Form("  463 keV peak: %.3f ± %.3f keV", mean_463_b.getVal(), mean_463_b.getError()) << std::endl;
    std::cout << Form("  478 keV peak: %.3f ± %.3f keV", mean_478_b.getVal(), mean_478_b.getError()) << std::endl;
    std::cout << Form("  481 keV peak: %.3f ± %.3f keV", mean_481_b.getVal(), mean_481_b.getError()) << std::endl;
    std::cout << Form("  486 keV peak: %.3f ± %.3f keV", mean_486_b.getVal(), mean_486_b.getError()) << std::endl;
    
    std::cout << "\nPeak Widths (Sigma values):" << std::endl;
    std::cout << Form("  463 keV peak: %.3f ± %.3f keV", sigma_463_b.getVal(), sigma_463_b.getError()) << std::endl;
    std::cout << Form("  478 keV peak: %.3f ± %.3f keV", sigma_478_b.getVal(), sigma_478_b.getError()) << std::endl;
    std::cout << Form("  481 keV peak: %.3f ± %.3f keV", sigma_481_b.getVal(), sigma_481_b.getError()) << std::endl;
    std::cout << Form("  486 keV peak: %.3f ± %.3f keV", sigma_486_b.getVal(), sigma_486_b.getError()) << std::endl;
    
    std::cout << "\nUNBINNED FIT RESULTS:" << std::endl;
    std::cout << "Peak Positions (Mean values):" << std::endl;
    std::cout << Form("  463 keV peak: %.3f ± %.3f keV", mean_463_u.getVal(), mean_463_u.getError()) << std::endl;
    std::cout << Form("  478 keV peak: %.3f ± %.3f keV", mean_478_u.getVal(), mean_478_u.getError()) << std::endl;
    std::cout << Form("  481 keV peak: %.3f ± %.3f keV", mean_481_u.getVal(), mean_481_u.getError()) << std::endl;
    std::cout << Form("  486 keV peak: %.3f ± %.3f keV", mean_486_u.getVal(), mean_486_u.getError()) << std::endl;
    
    std::cout << "\nPeak Widths (Sigma values):" << std::endl;
    std::cout << Form("  463 keV peak: %.3f ± %.3f keV", sigma_463_u.getVal(), sigma_463_u.getError()) << std::endl;
    std::cout << Form("  478 keV peak: %.3f ± %.3f keV", sigma_478_u.getVal(), sigma_478_u.getError()) << std::endl;
    std::cout << Form("  481 keV peak: %.3f ± %.3f keV", sigma_481_u.getVal(), sigma_481_u.getError()) << std::endl;
    std::cout << Form("  486 keV peak: %.3f ± %.3f keV", sigma_486_u.getVal(), sigma_486_u.getError()) << std::endl;
    
    std::cout << "\n=== Initial vs Fitted Values Comparison ===" << std::endl;
    std::cout << "Peak Positions (Initial → Fitted Binned):" << std::endl;
    std::cout << Form("  463 keV: 463.000 → %.3f (Δ = %.3f)", mean_463_b.getVal(), mean_463_b.getVal() - 463.0) << std::endl;
    std::cout << Form("  478 keV: 478.000 → %.3f (Δ = %.3f)", mean_478_b.getVal(), mean_478_b.getVal() - 478.0) << std::endl;
    std::cout << Form("  481 keV: 481.000 → %.3f (Δ = %.3f)", mean_481_b.getVal(), mean_481_b.getVal() - 481.0) << std::endl;
    std::cout << Form("  486 keV: 486.000 → %.3f (Δ = %.3f)", mean_486_b.getVal(), mean_486_b.getVal() - 486.0) << std::endl;
    
    std::cout << "\nPeak Widths (Initial → Fitted Binned):" << std::endl;
    std::cout << Form("  463 keV: %.3f → %.3f (Δ = %.3f)", peak_sigmas[0], sigma_463_b.getVal(), sigma_463_b.getVal() - peak_sigmas[0]) << std::endl;
    std::cout << Form("  478 keV: %.3f → %.3f (Δ = %.3f)", peak_sigmas[1], sigma_478_b.getVal(), sigma_478_b.getVal() - peak_sigmas[1]) << std::endl;
    std::cout << Form("  481 keV: %.3f → %.3f (Δ = %.3f)", peak_sigmas[2], sigma_481_b.getVal(), sigma_481_b.getVal() - peak_sigmas[2]) << std::endl;
    std::cout << Form("  486 keV: %.3f → %.3f (Δ = %.3f)", peak_sigmas[3], sigma_486_b.getVal(), sigma_486_b.getVal() - peak_sigmas[3]) << std::endl;
    
    
    // Wait for user
    std::cout << "\nPress Enter to close..." << std::endl;
    c_fits->WaitPrimitive();
    
    // Cleanup
    delete fit_result_binned;
    delete fit_result_unbinned;
}

