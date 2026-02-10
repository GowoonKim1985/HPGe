#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <random>
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TROOT.h"

void step2_mkPseudo() {
    gStyle->SetOptStat(1111);
    
    std::cout << "=== Pseudo Data Generation ===" << std::endl;
    std::cout << "Creating 6-day measurement simulation with Be-7 peak (478 keV)" << std::endl;
    
    // Load calibration results for background
    std::cout << "\n1. Loading background data (count/day)..." << std::endl;
    TFile* calib_file = new TFile("calibration_results.root", "READ");
    if (!calib_file || !calib_file->IsOpen()) {
        std::cerr << "Error: Cannot open calibration_results.root!" << std::endl;
        std::cerr << "Please run the calibration code first." << std::endl;
        return;
    }
    
    TH1D* h_bg_norm = (TH1D*)calib_file->Get("h_energy_norm");
    if (!h_bg_norm) {
        std::cerr << "Error: Cannot find h_energy_norm histogram!" << std::endl;
        return;
    }
    
    // Clone for safety
    TH1D* h_background = (TH1D*)h_bg_norm->Clone("h_background");
    
    std::cout << "Background histogram loaded:" << std::endl;
    std::cout << "  Range: " << h_background->GetXaxis()->GetXmin() 
              << " to " << h_background->GetXaxis()->GetXmax() << " keV" << std::endl;
    std::cout << "  Bins: " << h_background->GetNbinsX() << std::endl;
    std::cout << "  Total rate: " << h_background->Integral() << " counts/day" << std::endl;
    
    // Load resolution function
    std::cout << "\n2. Loading resolution function..." << std::endl;
    TFile* res_file = new TFile("resolution_results.root", "READ");
    if (!res_file || !res_file->IsOpen()) {
        std::cerr << "Error: Cannot open resolution_results.root!" << std::endl;
        std::cerr << "Please run the resolution analysis first." << std::endl;
        return;
    }
    
    TF1* resolution_func = (TF1*)res_file->Get("resolution_function");
    if (!resolution_func) {
        std::cerr << "Error: Cannot find resolution_function!" << std::endl;
        std::cerr << "Using default resolution parameters..." << std::endl;
        // Create default resolution function if not found
        resolution_func = new TF1("resolution_func_default", "[0]/sqrt(x) + [1] + [2]*x", 300, 2700);
        resolution_func->SetParameters(20.0, 0.5, 0.0); // Typical values
    }
    
    // Calculate resolution at 478 keV
    double be7_energy = 478.0; // keV
    double resolution_478 = resolution_func->Eval(be7_energy); // %
    double fwhm_478 = (resolution_478 / 100.0) * be7_energy; // keV
    double sigma_478 = fwhm_478 / 2.35; // keV
    
    std::cout << "Resolution at 478 keV:" << std::endl;
    std::cout << "  Resolution: " << resolution_478 << " %" << std::endl;
    std::cout << "  FWHM: " << fwhm_478 << " keV" << std::endl;
    std::cout << "  Sigma: " << sigma_478 << " keV" << std::endl;
    
    // Simulation parameters
    int measurement_days = 6;
    //double be7_total_counts = 26.1; // Expected Be-7 counts in 6 days
    double be7_total_counts = 100.1; // Expected Be-7 counts in 6 days
    
    std::cout << "\n3. Simulation parameters:" << std::endl;
    std::cout << "  Measurement period: " << measurement_days << " days" << std::endl;
    std::cout << "  Expected Be-7 counts: " << be7_total_counts << std::endl;
    std::cout << "  Be-7 peak position: " << be7_energy << " keV" << std::endl;
    std::cout << "  Be-7 peak sigma: " << sigma_478 << " keV" << std::endl;
    
    // Create pseudo data histogram (same binning as background)
    TH1D* h_pseudo = (TH1D*)h_background->Clone("h_pseudo_data");
    h_pseudo->Reset(); // Clear contents
    h_pseudo->SetTitle("Pseudo Data (6-day measurement);Energy (keV);Counts");
    
    // Initialize random number generator
    TRandom3* rnd = new TRandom3(0); // Use time-based seed
    
    std::cout << "\n4. Generating background events..." << std::endl;
    
    // Convert background histogram to PDF and sample from it
    double total_bg_rate = h_background->Integral(); // counts/day
    double total_bg_6days = total_bg_rate * measurement_days; // expected counts in 6 days
    
    // Sample actual number of background events from Poisson distribution
    int actual_bg_events = rnd->Poisson(total_bg_6days);
    
    std::cout << "  Expected background: " << total_bg_6days << " counts" << std::endl;
    std::cout << "  Sampled background: " << actual_bg_events << " counts" << std::endl;
    
    // Generate background events
    for (int i = 0; i < actual_bg_events; ++i) {
        double energy = h_background->GetRandom(); // Sample from PDF
        h_pseudo->Fill(energy);
        
        if (i % 10000 == 0 && i > 0) {
            std::cout << "    Generated " << i << " background events..." << std::endl;
        }
    }
    
    std::cout << "  Background generation completed!" << std::endl;
    
    std::cout << "\n5. Generating Be-7 peak events..." << std::endl;
    
    // Sample actual number of Be-7 events from Poisson distribution
    int actual_be7_events = rnd->Poisson(be7_total_counts);
    
    std::cout << "  Expected Be-7 events: " << be7_total_counts << std::endl;
    std::cout << "  Sampled Be-7 events: " << actual_be7_events << std::endl;
    
    // Generate Be-7 peak events (Gaussian around 478 keV)
    for (int i = 0; i < actual_be7_events; ++i) {
        double energy = rnd->Gaus(be7_energy, sigma_478);
        
        // Only fill if within histogram range
        if (energy >= h_pseudo->GetXaxis()->GetXmin() && 
            energy <= h_pseudo->GetXaxis()->GetXmax()) {
            h_pseudo->Fill(energy);
        }
    }
    
    std::cout << "  Be-7 peak generation completed!" << std::endl;
    
    // Summary statistics
    double total_counts = h_pseudo->Integral();
    double be7_region_min = be7_energy - 3*sigma_478;
    double be7_region_max = be7_energy + 3*sigma_478;
    int bin_min = h_pseudo->FindBin(be7_region_min);
    int bin_max = h_pseudo->FindBin(be7_region_max);
    double be7_region_counts = h_pseudo->Integral(bin_min, bin_max);
    
    std::cout << "\n=== Pseudo Data Summary ===" << std::endl;
    std::cout << "Total counts: " << total_counts << std::endl;
    std::cout << "Background counts: " << actual_bg_events << std::endl;
    std::cout << "Be-7 peak counts: " << actual_be7_events << std::endl;
    std::cout << "Counts in Be-7 region (±3σ): " << be7_region_counts << std::endl;
    
    // Create comparison plots
    std::cout << "\n6. Creating comparison plots..." << std::endl;
    
    TCanvas* c_comparison = new TCanvas("c_comparison", "Pseudo Data vs Background", 1200, 600);
    c_comparison->Divide(2, 1);
    
    // Prepare 6-day background for comparison
    TH1D* h_bg_6days = (TH1D*)h_background->Clone("h_bg_6days");
    h_bg_6days->Scale(measurement_days);
    h_bg_6days->SetTitle("Pseudo Data vs 6-day Background;Energy (keV);Counts");
    h_bg_6days->SetLineColor(kBlue);
    h_bg_6days->SetFillColor(kBlue);
    h_bg_6days->SetFillStyle(3004);
    
    // Set pseudo data style
    h_pseudo->SetLineColor(kRed);
    h_pseudo->SetMarkerColor(kRed);
    h_pseudo->SetMarkerStyle(20);
    h_pseudo->SetMarkerSize(0.5);
    
    // Plot 1: Full energy range
    c_comparison->cd(1);
    gPad->SetLogy(); // Log scale for better visibility
    h_bg_6days->Draw("HIST");
    h_pseudo->Draw("E SAME");
    
    TLegend* leg1 = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg1->AddEntry(h_bg_6days, "6-day Background", "f");
    leg1->AddEntry(h_pseudo, "Pseudo Data", "pe");
    leg1->Draw();
    
    // Plot 2: Be-7 region (user range keV)
    c_comparison->cd(2);
    TH1D* h_bg_zoom = (TH1D*)h_bg_6days->Clone("h_bg_zoom");
    TH1D* h_pseudo_zoom = (TH1D*)h_pseudo->Clone("h_pseudo_zoom");
   

    double zoom_low = 440;
    double zoom_high = 500; 
    h_bg_zoom->GetXaxis()->SetRangeUser(zoom_low, zoom_high);
    h_pseudo_zoom->GetXaxis()->SetRangeUser(zoom_low, zoom_high);
    h_bg_zoom->SetTitle("Be-7 Region (450-520 keV);Energy (keV);Counts");
    
    h_bg_zoom->Draw("HIST");
    h_pseudo_zoom->Draw("SAME");
    
    // Add vertical line at 478 keV
    TLine* be7_line = new TLine(478, h_bg_zoom->GetMinimum(), 478, h_bg_zoom->GetMaximum()*1.2);
    be7_line->SetLineColor(kGreen+2);
    be7_line->SetLineWidth(2);
    be7_line->SetLineStyle(2);
    be7_line->Draw();
    
    TLegend* leg2 = new TLegend(0.15, 0.7, 0.45, 0.9);
    leg2->AddEntry(h_bg_zoom, "6-day Background", "f");
    leg2->AddEntry(h_pseudo_zoom, "Pseudo Data", "pe");
    leg2->AddEntry(be7_line, "Be-7 (478 keV)", "l");
    leg2->Draw();
    
    c_comparison->Update();
    
    // Save results
    std::cout << "\n7. Saving pseudo data..." << std::endl;
    
    TFile* outfile = new TFile("pseudo_data.root", "RECREATE");
    if (outfile->IsOpen()) {
        // Save pseudo data histogram
        h_pseudo->Write();
        
        // Save background for comparison
        h_background->Write();
        h_bg_6days->Write();
        
        // Save canvas
        c_comparison->Write();
        
        // === NEW: Save data for likelihood analysis ===
        std::cout << "  Preparing data for likelihood analysis..." << std::endl;
        
        // 1. Binned likelihood: Histogram for 440-500 keV region
        double analysis_low = 440.0;
        double analysis_high = 500.0;
        
        // Find bin range
        int bin_low = h_pseudo->FindBin(analysis_low);
        int bin_high = h_pseudo->FindBin(analysis_high);
        
        // Create analysis histogram
        TH1D* h_analysis = new TH1D("h_analysis", 
                                   "Pseudo Data for Analysis (440-500 keV);Energy (keV);Counts",
                                   bin_high - bin_low + 1,
                                   h_pseudo->GetBinLowEdge(bin_low),
                                   h_pseudo->GetBinLowEdge(bin_high + 1));
        
        // Fill analysis histogram
        for (int i = bin_low; i <= bin_high; ++i) {
            double bin_center = h_pseudo->GetBinCenter(i);
            double bin_content = h_pseudo->GetBinContent(i);
            int analysis_bin = h_analysis->FindBin(bin_center);
            h_analysis->SetBinContent(analysis_bin, bin_content);
            h_analysis->SetBinError(analysis_bin, sqrt(bin_content)); // Poisson error
        }
        
        h_analysis->Write();
        
        // Also save background histogram for the same region
        TH1D* h_bg_analysis = new TH1D("h_bg_analysis",
                                      "Background for Analysis (440-500 keV);Energy (keV);Counts/day",
                                      bin_high - bin_low + 1,
                                      h_background->GetBinLowEdge(bin_low),
                                      h_background->GetBinLowEdge(bin_high + 1));
        
        for (int i = bin_low; i <= bin_high; ++i) {
            double bin_center = h_background->GetBinCenter(i);
            double bin_content = h_background->GetBinContent(i);
            int analysis_bin = h_bg_analysis->FindBin(bin_center);
            h_bg_analysis->SetBinContent(analysis_bin, bin_content);
        }
        
        h_bg_analysis->Write();
        
        // 2. Unbinned likelihood: Tree with individual events
        std::cout << "  Creating event tree for unbinned analysis..." << std::endl;
        
        TTree* event_tree = new TTree("event_tree", "Individual Events for Unbinned Analysis");
        
        // Tree variables
        Double_t energy;
        Int_t count;
        
        event_tree->Branch("energy", &energy, "energy/D");
        event_tree->Branch("count", &count, "count/I");
        
        // Fill tree with bin centers and counts
        int total_events = 0;
        for (int i = bin_low; i <= bin_high; ++i) {
            energy = h_pseudo->GetBinCenter(i);
            count = (Int_t)h_pseudo->GetBinContent(i);
            
            if (count > 0) {
                event_tree->Fill();
                total_events += count;
            }
        }
        
        event_tree->Write();
        
        std::cout << "  Analysis data prepared:" << std::endl;
        std::cout << "    Energy range: " << analysis_low << " - " << analysis_high << " keV" << std::endl;
        std::cout << "    Number of bins: " << h_analysis->GetNbinsX() << std::endl;
        std::cout << "    Total events in range: " << total_events << std::endl;
        std::cout << "    Tree entries: " << event_tree->GetEntries() << std::endl;
        
        // Save generation parameters
        TPaveText* params = new TPaveText(0, 0, 1, 1);
        params->SetName("generation_parameters");
        params->AddText(Form("Measurement days: %d", measurement_days));
        params->AddText(Form("Be-7 energy: %.1f keV", be7_energy));
        params->AddText(Form("Be-7 expected counts: %.1f", be7_total_counts));
        params->AddText(Form("Be-7 actual counts: %d", actual_be7_events));
        params->AddText(Form("Be-7 sigma: %.2f keV", sigma_478));
        params->AddText(Form("Background expected: %.1f", total_bg_6days));
        params->AddText(Form("Background actual: %d", actual_bg_events));
        params->AddText(Form("Total counts: %.0f", total_counts));
        params->AddText(Form("Analysis range: %.0f-%.0f keV", analysis_low, analysis_high));
        params->AddText(Form("Events in analysis range: %d", total_events));
        params->Write();
        
        std::cout << "Pseudo data saved to 'pseudo_data.root'" << std::endl;
        std::cout << "Saved objects:" << std::endl;
        std::cout << "  - h_pseudo_data (TH1D): Full pseudo data" << std::endl;
        std::cout << "  - h_analysis (TH1D): Analysis region (440-500 keV) for binned likelihood" << std::endl;
        std::cout << "  - h_bg_analysis (TH1D): Background in analysis region" << std::endl;
        std::cout << "  - event_tree (TTree): Individual events for unbinned likelihood" << std::endl;
        std::cout << "  - h_background (TH1D): Original background" << std::endl;
        std::cout << "  - h_bg_6days (TH1D): 6-day background" << std::endl;
        std::cout << "  - c_comparison (TCanvas): Comparison plots" << std::endl;
        std::cout << "  - generation_parameters (TPaveText): Parameters" << std::endl;
        
        outfile->Close();
    } else {
        std::cerr << "Error: Could not create pseudo_data.root!" << std::endl;
    }
} 
