#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TROOT.h"

void step1_eneResol() {
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    
    std::cout << "=== Energy Resolution Analysis ===" << std::endl;
    
    // Load calibration results
    std::cout << "Loading calibration results..." << std::endl;
    TFile* infile = new TFile("calibration_results.root", "READ");
    if (!infile || !infile->IsOpen()) {
        std::cerr << "Error: Cannot open calibration_results.root!" << std::endl;
        std::cerr << "Please run the calibration code first." << std::endl;
        return;
    }
    
    // Get calibrated energy histogram
    TH1D* h_energy = (TH1D*)infile->Get("h_energy_norm");
    if (!h_energy) {
        std::cerr << "Error: Cannot find h_energy_norm histogram!" << std::endl;
        infile->Close();
        return;
    }
    
    // Clone histogram to avoid issues
    TH1D* h_work = (TH1D*)h_energy->Clone("h_work");
    
    // Get calibration function
    TF1* calib_func = (TF1*)infile->Get("calibration_function");
    if (!calib_func) {
        std::cerr << "Error: Cannot find calibration_function!" << std::endl;
        infile->Close();
        return;
    }
    
    std::cout << "Calibration function loaded: " << calib_func->GetExpFormula() << std::endl;
    std::cout << "Slope = " << calib_func->GetParameter(1) << " keV/ch" << std::endl;
    std::cout << "Intercept = " << calib_func->GetParameter(0) << " keV" << std::endl;
    
    // Check histogram range
    std::cout << "Histogram range: " << h_work->GetXaxis()->GetXmin() 
              << " to " << h_work->GetXaxis()->GetXmax() << " keV" << std::endl;
    std::cout << "Number of bins: " << h_work->GetNbinsX() << std::endl;
    
    // Define known peak energies for resolution analysis
    std::vector<double> peak_energies = {2614.5, 1764.5, 1460.8, 1238.1, 609.3, 352.0}; // keV
    std::vector<std::string> peak_names = {"2615 keV (Tl-208)", "1765 keV (Bi-214)", 
                                          "1461 keV (K-40)", "1238 keV (unidentified)",
                                          "609 keV (Bi-214)", "352 keV (Pb-214)"};
    
    std::cout << "\nFitting peaks for energy resolution analysis..." << std::endl;
    std::cout << "Target energies: ";
    for (double energy : peak_energies) {
        std::cout << energy << " ";
    }
    std::cout << "keV" << std::endl;
    std::cout << "Fit range: +- 10 keV around each peak" << std::endl;
    
    // Vectors to store results
    std::vector<double> fitted_energies;
    std::vector<double> fitted_sigmas;
    std::vector<double> fitted_fwhms;
    std::vector<double> resolutions;
    std::vector<double> energy_errors;
    std::vector<double> sigma_errors;
    std::vector<double> resolution_errors;  // 해상도 에러 벡터 추가
    
    // Create canvas for peak fits
    TCanvas* c_peaks = new TCanvas("c_peaks", "Peak Fits for Resolution", 1200, 800);
    c_peaks->Divide(3, 2); // 3x2 grid for 6 peaks
    
    // Fit each peak with Gaussian + Linear background
    for (int i = 0; i < peak_energies.size(); ++i) {
        double target_energy = peak_energies[i];
        std::cout << "\n--- Fitting " << peak_names[i] << " ---" << std::endl;
        
        // Check if target energy is within histogram range
        if (target_energy < h_work->GetXaxis()->GetXmin() || 
            target_energy > h_work->GetXaxis()->GetXmax()) {
            std::cout << "Warning: " << target_energy << " keV is outside histogram range. Skipping..." << std::endl;
            continue;
        }
        
        // Find peak region (±10 keV around target for searching)
        double search_range = 10.0;
        double search_min = std::max(target_energy - search_range, h_work->GetXaxis()->GetXmin());
        double search_max = std::min(target_energy + search_range, h_work->GetXaxis()->GetXmax());
        
        int bin_left = h_work->FindBin(search_min);
        int bin_right = h_work->FindBin(search_max);
        
        std::cout << "Searching for peak in range: " << search_min << " to " << search_max << " keV" << std::endl;
        std::cout << "Bin range: " << bin_left << " to " << bin_right << std::endl;
        
        // Find actual peak position in the region
        double max_content = 0;
        int max_bin = h_work->FindBin(target_energy);
        for (int bin = bin_left; bin <= bin_right; ++bin) {
            double content = h_work->GetBinContent(bin);
            if (content > max_content) {
                max_content = content;
                max_bin = bin;
            }
        }
        
        if (max_content <= 0) {
            std::cout << "Warning: No peak found for " << peak_names[i] << ". Skipping..." << std::endl;
            continue;
        }
        
        double peak_pos = h_work->GetBinCenter(max_bin);
        std::cout << "Peak found at bin " << max_bin << ", energy = " << peak_pos << " keV, counts = " << max_content << std::endl;
        
        // Set fitting range: ±10 keV around found peak
        double fit_range = 10.0;
        double fit_min = std::max(peak_pos - fit_range, h_work->GetXaxis()->GetXmin());
        double fit_max = std::min(peak_pos + fit_range, h_work->GetXaxis()->GetXmax());
        
        std::cout << "Fit range: " << fit_min << " to " << fit_max << " keV (±" << fit_range << " keV)" << std::endl;
        
        // Estimate background level at the edges
        int bin_min = h_work->FindBin(fit_min);
        int bin_max = h_work->FindBin(fit_max);
        double bg_left = 0, bg_right = 0;
        int bg_bins = 3; // Average over 3 bins at each edge (smaller range now)
        
        for (int j = 0; j < bg_bins; ++j) {
            if (bin_min + j <= h_work->GetNbinsX()) {
                bg_left += h_work->GetBinContent(bin_min + j);
            }
            if (bin_max - j >= 1) {
                bg_right += h_work->GetBinContent(bin_max - j);
            }
        }
        bg_left /= bg_bins;
        bg_right /= bg_bins;
        double avg_bg = (bg_left + bg_right) / 2.0;
        
        std::cout << "Background estimate: left=" << bg_left << ", right=" << bg_right << ", avg=" << avg_bg << std::endl;
        
        // Create Gaussian + Linear background fit function
        // f(x) = [0]*exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*x
        TF1* gaus_bg_fit = new TF1(Form("gaus_bg_fit_%d", i), 
                                   "[0]*exp(-0.5*((x-[1])/[2])**2) + [3] + [4]*x", 
                                   fit_min, fit_max);
        
        // Set initial parameters
        double initial_sigma = (target_energy < 600) ? 0.5 : 3.0; // keV, smaller for low energy
        double net_amplitude = max_content - avg_bg;
        
        gaus_bg_fit->SetParameter(0, net_amplitude);           // Gaussian amplitude
        gaus_bg_fit->SetParameter(1, peak_pos);                // Gaussian mean
        gaus_bg_fit->SetParameter(2, initial_sigma);           // Gaussian sigma
        gaus_bg_fit->SetParameter(3, avg_bg);                  // Background constant
        gaus_bg_fit->SetParameter(4, 0.0);                     // Background slope
        
        // Set parameter names
        gaus_bg_fit->SetParNames("Amplitude", "Mean", "Sigma", "BG_Const", "BG_Slope");
        
        // Set reasonable parameter limits
        gaus_bg_fit->SetParLimits(0, 0, 100 * max_content);     // Amplitude > 0
        gaus_bg_fit->SetParLimits(1, fit_min, fit_max);        // Mean within fit range
        gaus_bg_fit->SetParLimits(2, .5, 10.0);               // Sigma reasonable range (smaller for ±20 keV range)
        gaus_bg_fit->SetParLimits(3, 0, 2 * avg_bg);           // Background positive
        
        std::cout << "Initial parameters:" << std::endl;
        std::cout << "  Amplitude: " << net_amplitude << std::endl;
        std::cout << "  Mean: " << peak_pos << std::endl;
        std::cout << "  Sigma: " << initial_sigma << std::endl;
        std::cout << "  BG_Const: " << avg_bg << std::endl;
        std::cout << "  BG_Slope: 0.0" << std::endl;
        
        // Perform fit
        Int_t fit_status = h_work->Fit(gaus_bg_fit, "RQN", "", fit_min, fit_max);
        
        std::cout << "Fit status: " << fit_status << std::endl;
        
        if (fit_status == 0) { // Successful fit
            double fitted_mean = gaus_bg_fit->GetParameter(1);
            double fitted_sigma = gaus_bg_fit->GetParameter(2);
            double fitted_fwhm = 2.35 * fitted_sigma;
            double resolution = (fitted_fwhm / fitted_mean) * 100.0; // %
            
            double mean_error = gaus_bg_fit->GetParError(1);
            double sigma_error = gaus_bg_fit->GetParError(2);
            
            // Check if fit results are reasonable
            if (fitted_sigma > 0 && fitted_sigma < 50 && 
                fitted_mean > fit_min && fitted_mean < fit_max &&
                resolution > 0 && resolution < 50) {
                
                // 해상도 에러 계산 (에러 전파)
                // R = (FWHM/E) * 100 = (2.35*σ/E) * 100
                // δR = R * sqrt((δσ/σ)² + (δE/E)²)
                double rel_sigma_error = sigma_error / fitted_sigma;
                double rel_energy_error = mean_error / fitted_mean;
                double resolution_error = resolution * sqrt(rel_sigma_error*rel_sigma_error + 
                                                           rel_energy_error*rel_energy_error);
                
                // Store results
                fitted_energies.push_back(fitted_mean);
                fitted_sigmas.push_back(fitted_sigma);
                fitted_fwhms.push_back(fitted_fwhm);
                resolutions.push_back(resolution);
                energy_errors.push_back(mean_error);
                sigma_errors.push_back(sigma_error);
                resolution_errors.push_back(resolution_error);
                
                std::cout << "✓ Peak fitted successfully:" << std::endl;
                std::cout << "  Mean: " << fitted_mean << " ± " << mean_error << " keV" << std::endl;
                std::cout << "  Sigma: " << fitted_sigma << " ± " << sigma_error << " keV" << std::endl;
                std::cout << "  FWHM: " << fitted_fwhm << " keV" << std::endl;
                std::cout << "  Resolution: " << resolution << " ± " << resolution_error << " %" << std::endl;
                std::cout << "  Chi2/NDF: " << gaus_bg_fit->GetChisquare() << "/" << gaus_bg_fit->GetNDF() << std::endl;
                std::cout << "  BG_Const: " << gaus_bg_fit->GetParameter(3) << " ± " << gaus_bg_fit->GetParError(3) << std::endl;
                std::cout << "  BG_Slope: " << gaus_bg_fit->GetParameter(4) << " ± " << gaus_bg_fit->GetParError(4) << std::endl;
                
                // Draw fit
                c_peaks->cd(fitted_energies.size());
                
                // Create a copy of histogram for this plot
                TH1D* h_plot = (TH1D*)h_work->Clone(Form("h_plot_%d", i));
                h_plot->GetXaxis()->SetRangeUser(fit_min - 10, fit_max + 10);
                h_plot->SetTitle(Form("%s Fit;Energy (keV);Counts", peak_names[i].c_str()));
                h_plot->Draw();
                
                // Draw total fit
                gaus_bg_fit->SetLineColor(kRed);
                gaus_bg_fit->SetLineWidth(2);
                gaus_bg_fit->Draw("same");
                
                // Draw background component
                TF1* bg_component = new TF1(Form("bg_%d", i), "[0] + [1]*x", fit_min, fit_max);
                bg_component->SetParameter(0, gaus_bg_fit->GetParameter(3));
                bg_component->SetParameter(1, gaus_bg_fit->GetParameter(4));
                bg_component->SetLineColor(kBlue);
                bg_component->SetLineStyle(2);
                bg_component->SetLineWidth(2);
                bg_component->Draw("same");
                
                // Draw Gaussian component
                TF1* gaus_component = new TF1(Form("gaus_%d", i), 
                                              "[0]*exp(-0.5*((x-[1])/[2])**2)", fit_min, fit_max);
                gaus_component->SetParameter(0, gaus_bg_fit->GetParameter(0));
                gaus_component->SetParameter(1, gaus_bg_fit->GetParameter(1));
                gaus_component->SetParameter(2, gaus_bg_fit->GetParameter(2));
                gaus_component->SetLineColor(kGreen);
                gaus_component->SetLineStyle(3);
                gaus_component->SetLineWidth(2);
                gaus_component->Draw("same");
                
                // Add legend
                TLegend* leg = new TLegend(0.6, 0.75, 0.9, 0.9);
                leg->AddEntry(gaus_bg_fit, "Total Fit", "l");
                leg->AddEntry(gaus_component, "Gaussian", "l");
                leg->AddEntry(bg_component, "Background", "l");
                leg->SetTextSize(0.03);
                leg->Draw();
                
                // Add fit info
                TPaveText* fit_info = new TPaveText(0.05, 0.65, 0.45, 0.9, "NDC");
                fit_info->AddText(Form("Mean: %.1f #pm %.1f keV", fitted_mean, mean_error));
                fit_info->AddText(Form("Sigma: %.1f #pm %.1f keV", fitted_sigma, sigma_error));
                fit_info->AddText(Form("FWHM: %.1f keV", fitted_fwhm));
                fit_info->AddText(Form("Resolution: %.1f #pm %.1f %%", resolution, resolution_error));
                fit_info->AddText(Form("#chi^{2}/NDF: %.1f/%d", gaus_bg_fit->GetChisquare(), gaus_bg_fit->GetNDF()));
                fit_info->SetFillColor(kWhite);
                fit_info->SetBorderSize(1);
                fit_info->SetTextSize(0.025);
                fit_info->Draw();
                
            } else {
                std::cout << "✗ Fit results unreasonable. Skipping..." << std::endl;
                std::cout << "  Fitted sigma: " << fitted_sigma << " keV" << std::endl;
                std::cout << "  Fitted mean: " << fitted_mean << " keV" << std::endl;
                std::cout << "  Resolution: " << resolution << " %" << std::endl;
            }
        } else {
            std::cout << "✗ Fit failed for " << peak_names[i] << std::endl;
        }
        
        //delete gaus_bg_fit;
    }
    
    c_peaks->Update();
    
    if (fitted_energies.size() < 2) {
        std::cerr << "Error: Not enough successful fits (" << fitted_energies.size() 
                  << ") for resolution analysis!" << std::endl;
        std::cerr << "Need at least 2 peaks for resolution function fitting." << std::endl;
        infile->Close();
        return;
    }
    
    // Create resolution vs energy graph with error bars
    std::cout << "\n=== Energy Resolution Analysis ===" << std::endl;
    std::cout << "Successfully fitted " << fitted_energies.size() << " peaks:" << std::endl;
    
    TGraphErrors* gr_resolution = new TGraphErrors();
    gr_resolution->SetName("gr_resolution");
    gr_resolution->SetTitle("Energy Resolution vs Energy;Energy (keV);Resolution (%)");
    gr_resolution->SetMarkerStyle(20);
    gr_resolution->SetMarkerSize(1.2);
    gr_resolution->SetMarkerColor(kBlue);
    gr_resolution->SetLineColor(kBlue);
    
    for (int i = 0; i < fitted_energies.size(); ++i) {
        gr_resolution->SetPoint(i, fitted_energies[i], resolutions[i]);
        gr_resolution->SetPointError(i, energy_errors[i], resolution_errors[i]);
        std::cout << "  E = " << fitted_energies[i] << " ± " << energy_errors[i] 
                  << " keV, Resolution = " << resolutions[i] << " ± " << resolution_errors[i] << " %" << std::endl;
    }
    
    // Fit resolution function: R(E) = a/√E + b + c*E
    TF1* resolution_fit = new TF1("resolution_fit", "[0]/sqrt(x) + [1] + [2]*x", 300, 2700);
    resolution_fit->SetParameters(20.0, 0.5, 0.0); // Initial guess
    resolution_fit->SetParNames("a", "b", "c");
    resolution_fit->SetLineColor(kRed);
    resolution_fit->SetLineWidth(2);
    
    std::cout << "\nFitting resolution function R(E) = a/√E + b + c×E..." << std::endl;
    Int_t res_fit_status = gr_resolution->Fit(resolution_fit, "RQN");
    
    // Create resolution canvas
    TCanvas* c_resolution = new TCanvas("c_resolution", "Energy Resolution Analysis", 800, 600);
    gr_resolution->Draw("AP");
    
    if (res_fit_status == 0) {
        resolution_fit->Draw("same");
        
        // Add fit results
        TPaveText* res_info = new TPaveText(0.55, 0.65, 0.9, 0.9, "NDC");
        res_info->AddText("Resolution Function:");
        res_info->AddText("R(E) = a/#sqrt{E} + b + c#timesE");
        res_info->AddText(Form("a = %.2f #pm %.2f", resolution_fit->GetParameter(0), resolution_fit->GetParError(0)));
        res_info->AddText(Form("b = %.3f #pm %.3f", resolution_fit->GetParameter(1), resolution_fit->GetParError(1)));
        res_info->AddText(Form("c = %.2e #pm %.2e", resolution_fit->GetParameter(2), resolution_fit->GetParError(2)));
        res_info->AddText(Form("#chi^{2}/NDF = %.2f/%d", resolution_fit->GetChisquare(), resolution_fit->GetNDF()));
        res_info->SetFillColor(kWhite);
        res_info->SetBorderSize(1);
        res_info->SetTextSize(0.03);
        res_info->Draw();
        
        // Calculate resolution at 478 keV
        double resolution_478 = resolution_fit->Eval(478.0);
        double sigma_478 = (resolution_478 / 100.0) * 478.0 / 2.35; // Convert % to keV, then to sigma
        
        std::cout << "\n=== Resolution at 478 keV ===" << std::endl;
        std::cout << "Resolution at 478 keV: " << resolution_478 << " %" << std::endl;
        std::cout << "FWHM at 478 keV: " << (resolution_478 / 100.0) * 478.0 << " keV" << std::endl;
        std::cout << "Sigma at 478 keV: " << sigma_478 << " keV" << std::endl;
    } else {
        std::cout << "Warning: Resolution function fit failed!" << std::endl;
    }
    
    c_resolution->Update();

    std::cout << "\n=== Resolution Peak Identification ===" << std::endl;
    std::cout << "Zoom and inspect the spectrum as needed." << std::endl;
    std::cout << "When ready, close the canvas window or press 'q' in canvas to continue..." << std::endl;

    c_resolution->WaitPrimitive();
    
    // Save results
    std::cout << "\nDo you want to save the resolution analysis results? (y/n): ";
    char save_choice;
    std::cin >> save_choice;
    
    if (save_choice == 'y' || save_choice == 'Y') {
        TFile* outfile = new TFile("resolution_results.root", "RECREATE");
        
        if (outfile->IsOpen()) {
            std::cout << "Saving resolution results to 'resolution_results.root'..." << std::endl;
            
            // Save resolution function if fit was successful
            if (res_fit_status == 0) {
                resolution_fit->SetName("resolution_function");
                resolution_fit->SetTitle("Energy Resolution Function: R(E) = a/sqrt(E) + b + c*E");
                resolution_fit->Write();
            }
            
            // Save resolution graph
            gr_resolution->Write();
            
            // Save canvases
            c_peaks->Write();
            c_resolution->Write();
            
            // Save calibrated histogram for reference
            h_work->Write();
            
            //outfile->Close();
            std::cout << "Resolution analysis results saved successfully!" << std::endl;
            std::cout << "Saved objects:" << std::endl;
            if (res_fit_status == 0) {
                std::cout << "  - resolution_function (TF1)" << std::endl;
            }
            std::cout << "  - gr_resolution (TGraphErrors)" << std::endl;
            std::cout << "  - c_peaks (TCanvas)" << std::endl;
            std::cout << "  - c_resolution (TCanvas)" << std::endl;
            std::cout << "  - h_work (TH1D)" << std::endl;
        } else {
            std::cerr << "Error: Could not create output file!" << std::endl;
        }
        
        //delete outfile;
    }
    
    //infile->Close();
    std::cout << "\nEnergy resolution analysis completed!" << std::endl;
}

