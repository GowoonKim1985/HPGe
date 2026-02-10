void step0_DataRead_calibration() {
  //    const char* filename = "/home/dhh/Linux_workDir/Yemilab_HPGe_Study/2025_Apr_Yemi_Sens_BG/20250515-0530_data/2025May_BGdata3_15.TKA";
    const char* filename = "/home/kkw/study/unbinned/2025May_BGdata3_15.TKA";

    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return;
    }

    std::string line;
    long liveTime = 0, realTime = 0;

    if (std::getline(infile, line)) liveTime = std::stol(line);
    if (std::getline(infile, line)) realTime = std::stol(line);

    std::vector<int> adcCounts;
    while (std::getline(infile, line)) {
        int val = std::stoi(line);
        adcCounts.push_back(val);
    }

    infile.close();

    const int nBins = 8192;
    double liveTime_day = liveTime / 86400.0;

    std::cout << "Live Time: " << liveTime << " sec (" << liveTime_day << " day)" << std::endl;
    std::cout << "Real Time: " << realTime << " sec" << std::endl;
    std::cout << "Data points read: " << adcCounts.size() << std::endl;
    std::cout << "Total histogram bins: " << nBins << std::endl;

    TH1D* h_raw = new TH1D("h_raw", "ADC Spectrum (Raw);Channel;Counts", nBins, 0.5, nBins + 0.5);
    TH1D* h_norm = new TH1D("h_norm", "ADC Spectrum (Normalized);Channel;Counts/day", nBins, 0.5, nBins + 0.5);

    // Fill histograms
    h_raw->SetBinContent(1, 0);
    h_raw->SetBinContent(2, 0);
    h_norm->SetBinContent(1, 0);
    h_norm->SetBinContent(2, 0);

    for (int i = 0; i < adcCounts.size(); ++i) {
        int binNumber = i + 3;
        if (binNumber <= nBins) {
            h_raw->SetBinContent(binNumber, adcCounts[i]);
            h_norm->SetBinContent(binNumber, adcCounts[i] / liveTime_day);
        }
    }

    // Draw histograms
    TCanvas* c1 = new TCanvas("c1", "ADC Spectra", 800, 800);
    c1->ToggleEventStatus();
    c1->Divide(1, 2);
    c1->cd(1);
    h_raw->Draw();
    c1->cd(2);
    h_norm->Draw();
    c1->Update();

    std::cout << "\n=== Calibration Peak Identification ===" << std::endl;
    std::cout << "Zoom and inspect the spectrum as needed." << std::endl;
    std::cout << "When ready, close the canvas window or press 'q' in canvas to continue..." << std::endl;
    
    c1->WaitPrimitive();

    std::cout << "Enter the ADC channel for 2614.5 keV peak: ";
    double adc_2614;
    std::cin >> adc_2614;

    // Reopen canvas for results
    if (c1->GetCanvasID() == -1) {
        c1 = new TCanvas("c1", "ADC Spectra with Calibration", 800, 800);
        c1->ToggleEventStatus();
        c1->Divide(1, 2);
        c1->cd(1);
        h_raw->Draw();
        c1->cd(2);
        h_norm->Draw();
    }

    double slope = 2614.5 / adc_2614;

    std::cout << "\n=== Predicted ADC channels for other peaks ===" << std::endl;
    std::cout << "Linear calibration: Energy = " << slope << " * ADC" << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;

    std::vector<double> energies = {2614.5, 1764.5, 1460.8, 1238.1, 609.3, 511.0, 352.0};
    std::vector<std::string> sources = {"Tl-208", "Bi-214", "K-40", "Bi-214", "Bi-214", "annihilation", "Pb-214"};
    std::vector<double> predicted_adcs;

    for (int i = 0; i < energies.size(); ++i) {
        double predicted_adc = energies[i] / slope;
        predicted_adcs.push_back(predicted_adc);
        std::cout << std::fixed << std::setprecision(1);
        std::cout << energies[i] << " keV (" << sources[i] << "): ~" 
                  << (int)predicted_adc << " ADC" << std::endl;
    }

    // Add vertical lines
    c1->cd(2);
    double ymax = h_norm->GetMaximum();
    
    for (int i = 0; i < predicted_adcs.size(); ++i) {
        TLine* line = new TLine(predicted_adcs[i], 0, predicted_adcs[i], ymax);
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        line->SetLineStyle(2);
        line->Draw("same");
        
        TText* text = new TText(predicted_adcs[i], ymax * 0.9, Form("%.1f keV", energies[i]));
        text->SetTextColor(kRed);
        text->SetTextSize(0.03);
        text->SetTextAngle(90);
        text->Draw("same");
    }
    
    c1->Update();
    std::cout << "\nRed dashed lines show predicted peak positions." << std::endl;
    
    // Wait before proceeding to fitting
    std::cout << "\nProceeding to Gaussian fitting..." << std::endl;
    std::cout << "Close this canvas or press 'q' to continue to fitting..." << std::endl;
    c1->WaitPrimitive();

    // === GAUSSIAN FITTING SECTION ===
    std::cout << "\n=== Gaussian Fitting for Peak Calibration ===" << std::endl;
    
    // Create new canvas for fitting results
    TCanvas* c2 = new TCanvas("c2", "Peak Fitting Results", 1200, 800);
    c2->ToggleEventStatus();
    c2->Divide(4, 2);
    
    std::vector<double> fitted_means;
    std::vector<double> fitted_sigmas;
    std::vector<double> fitted_amplitudes;
    
    for (int i = 0; i < predicted_adcs.size(); ++i) {
        c2->cd(i + 1);
        
        double center = predicted_adcs[i];
        double fit_min = center - 5.0;
        double fit_max = center + 5.0;
        
        // Create histogram for this peak region
        TH1D* h_peak = (TH1D*)h_norm->Clone(Form("h_peak_%d", i));
        h_peak->SetTitle(Form("%.1f keV Peak Fit;ADC Channel;Counts/day", energies[i]));
        h_peak->GetXaxis()->SetRangeUser(fit_min - 2, fit_max + 2);
        
        // Define Gaussian function
        TF1* gaus_fit = new TF1(Form("gaus_%d", i), "gaus", fit_min, fit_max);
        gaus_fit->SetParameters(h_peak->GetBinContent(h_peak->FindBin(center)), center, 1.0);
        gaus_fit->SetParNames("Amplitude", "Mean", "Sigma");
        
        // Perform fit
        h_peak->Fit(gaus_fit, "RQ");
        
        // Store fitted parameters
        fitted_means.push_back(gaus_fit->GetParameter(1));
        fitted_sigmas.push_back(gaus_fit->GetParameter(2));
        fitted_amplitudes.push_back(gaus_fit->GetParameter(0));
        
        // Draw
        h_peak->Draw();
        gaus_fit->SetLineColor(kRed);
        gaus_fit->SetLineWidth(2);
        gaus_fit->Draw("same");
        
        // Add text with fit results
        TPaveText* pt = new TPaveText(0.6, 0.7, 0.9, 0.9, "NDC");
        pt->AddText(Form("%.1f keV", energies[i]));
        pt->AddText(Form("Mean: %.2f", gaus_fit->GetParameter(1)));
        pt->AddText(Form("Sigma: %.2f", gaus_fit->GetParameter(2)));
        pt->AddText(Form("#chi^{2}/NDF: %.2f", gaus_fit->GetChisquare()/gaus_fit->GetNDF()));
        pt->SetFillColor(kWhite);
        pt->SetBorderSize(1);
        pt->Draw();
    }
    
    // Add empty pad for summary
    c2->cd(8);
    TPaveText* summary = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
    summary->AddText("Fitted Peak Positions (ADC):");
    summary->AddText("------------------------");
    for (int i = 0; i < energies.size(); ++i) {
        summary->AddText(Form("%.1f keV: %.2f ADC", energies[i], fitted_means[i]));
    }
    summary->SetFillColor(kYellow);
    summary->SetBorderSize(1);
    summary->Draw();
    
    c2->Update();
    
    // Print results to console
    std::cout << "\n=== Fitted Peak Positions ===" << std::endl;
    std::cout << "Energy (keV) | Predicted ADC | Fitted ADC | Sigma | Chi2/NDF" << std::endl;
    std::cout << "-------------|---------------|------------|-------|----------" << std::endl;
    
    for (int i = 0; i < energies.size(); ++i) {
        TF1* fit_func = (TF1*)gROOT->FindObject(Form("gaus_%d", i));
        double chi2ndf = (fit_func) ? fit_func->GetChisquare()/fit_func->GetNDF() : -1;
        
        std::cout << std::fixed << std::setprecision(1) << std::setw(12) << energies[i] 
                  << " | " << std::setw(13) << predicted_adcs[i]
                  << " | " << std::setprecision(2) << std::setw(10) << fitted_means[i]
                  << " | " << std::setw(5) << fitted_sigmas[i]
                  << " | " << std::setw(8) << chi2ndf << std::endl;
    }
    
    std::cout << "\nFitting complete! Proceeding to energy calibration..." << std::endl;
    std::cout << "Close this canvas or press 'q' to continue..." << std::endl;
    c2->WaitPrimitive();

    // === ENERGY CALIBRATION SECTION ===
    std::cout << "\n=== Energy Calibration Function ===" << std::endl;
    
    // Create TGraph for ADC vs Energy
    TGraph* gr_calib = new TGraph(fitted_means.size());
    for (int i = 0; i < fitted_means.size(); ++i) {
        gr_calib->SetPoint(i, fitted_means[i], energies[i]);
    }
    
    // Linear fit: Energy = slope * ADC + intercept
    TF1* linear_fit = new TF1("linear_fit", "pol1", 0, 8192);
    gr_calib->Fit(linear_fit, "Q");
    
    double calib_slope = linear_fit->GetParameter(1);
    double calib_intercept = linear_fit->GetParameter(0);
    // R² 수동 계산
    double sum_y = 0, sum_y2 = 0;
    for (int i = 0; i < energies.size(); ++i) {
        sum_y += energies[i];
        sum_y2 += energies[i] * energies[i];
    }
    double mean_y = sum_y / energies.size();
    double ss_tot = sum_y2 - energies.size() * mean_y * mean_y;
    double ss_res = linear_fit->GetChisquare();
    double r_squared = 1.0 - (ss_res / ss_tot);

    
    // Draw calibration graph
    TCanvas* c3 = new TCanvas("c3", "Energy Calibration", 800, 600);
    c3->ToggleEventStatus();
    gr_calib->SetTitle("Energy Calibration;ADC Channel;Energy (keV)");
    gr_calib->SetMarkerStyle(20);
    gr_calib->SetMarkerSize(1.2);
    gr_calib->SetMarkerColor(kBlue);
    gr_calib->Draw("AP");
    
    linear_fit->SetLineColor(kRed);
    linear_fit->SetLineWidth(2);
    linear_fit->Draw("same");
    
    // Add calibration info
    TPaveText* calib_info = new TPaveText(0.15, 0.75, 0.5, 0.9, "NDC");
    calib_info->AddText(Form("Energy = %.4f #times ADC + %.2f", calib_slope, calib_intercept));
    calib_info->AddText(Form("R^{2} = %.6f", r_squared));
    calib_info->AddText(Form("#chi^{2}/NDF = %.2f", linear_fit->GetChisquare()/linear_fit->GetNDF()));
    calib_info->SetFillColor(kWhite);
    calib_info->SetBorderSize(1);
    calib_info->Draw();
    
    c3->Update();
    
    // Print calibration results
    std::cout << "Calibration Function: Energy = " << calib_slope << " * ADC + " << calib_intercept << std::endl;
    std::cout << "R² = " << r_squared << std::endl;
    std::cout << "χ²/NDF = " << linear_fit->GetChisquare()/linear_fit->GetNDF() << std::endl;
    
    std::cout << "\nCalibration accuracy check:" << std::endl;
    std::cout << "Energy (keV) | Fitted ADC | Calibrated Energy | Difference" << std::endl;
    std::cout << "-------------|------------|-------------------|------------" << std::endl;
    
    for (int i = 0; i < energies.size(); ++i) {
        double calibrated_energy = calib_slope * fitted_means[i] + calib_intercept;
        double difference = calibrated_energy - energies[i];
        std::cout << std::fixed << std::setprecision(1) << std::setw(12) << energies[i]
                  << " | " << std::setprecision(2) << std::setw(10) << fitted_means[i]
                  << " | " << std::setprecision(1) << std::setw(17) << calibrated_energy
                  << " | " << std::setw(10) << difference << std::endl;
    }

    std::cout << "\nProceeding to create energy-calibrated histogram..." << std::endl;
    std::cout << "Close this canvas or press 'q' to continue..." << std::endl;
    c3->WaitPrimitive();

    // === CREATE ENERGY-CALIBRATED HISTOGRAM ===
    std::cout << "\n=== Creating Energy-Calibrated Histogram ===" << std::endl;
    
    // Calculate energy range
    double min_energy = calib_slope * 1 + calib_intercept;
    double max_energy = calib_slope * nBins + calib_intercept;
    
    // Create energy histogram with same number of bins
    TH1D* h_energy_raw = new TH1D("h_energy_raw", "Energy Spectrum (Raw);Energy (keV);Counts", 
                                  nBins, min_energy, max_energy);
    TH1D* h_energy_norm = new TH1D("h_energy_norm", "Energy Spectrum (Normalized);Energy (keV);Counts/day", 
                                   nBins, min_energy, max_energy);
    
    // Fill energy histograms by copying contents from ADC histograms
    for (int i = 1; i <= nBins; ++i) {
        double adc_content_raw = h_raw->GetBinContent(i);
        double adc_content_norm = h_norm->GetBinContent(i);
        
        h_energy_raw->SetBinContent(i, adc_content_raw);
        h_energy_norm->SetBinContent(i, adc_content_norm);
    }
    
    // Create final comparison canvas
    TCanvas* c4 = new TCanvas("c4", "Energy Spectrum and Calibration Differences", 1200, 600);
    c4->ToggleEventStatus();
    c4->Divide(2, 1);
    
    // Pad 1: Energy spectrum with theoretical peak positions
    c4->cd(1);
    h_energy_norm->SetTitle("Energy Spectrum (Normalized)");
    h_energy_norm->Draw();
    
    // Add theoretical peak positions to energy spectrum
    double ymax_energy = h_energy_norm->GetMaximum();
    for (int i = 0; i < energies.size(); ++i) {
        TLine* line = new TLine(energies[i], 0, energies[i], ymax_energy);
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        line->SetLineStyle(2);
        line->Draw("same");
        
        TText* text = new TText(energies[i], ymax_energy * 0.9, Form("%.1f", energies[i]));
        text->SetTextColor(kRed);
        text->SetTextSize(0.03);
        text->SetTextAngle(90);
        text->Draw("same");
    }
    
    // Pad 2: Calibration differences graph
    c4->cd(2);
    
    // Create difference graph (Calibrated - Theoretical)
    TGraph* gr_difference = new TGraph(energies.size());
    double max_difference = 0;
    for (int i = 0; i < energies.size(); ++i) {
        double calibrated_energy = calib_slope * fitted_means[i] + calib_intercept;
        double difference = calibrated_energy - energies[i];
        gr_difference->SetPoint(i, energies[i], difference);
        if (TMath::Abs(difference) > max_difference) max_difference = TMath::Abs(difference);
    }
    
    gr_difference->SetTitle("Calibration Differences;Theoretical Energy (keV);Calibrated - Theoretical (keV)");
    gr_difference->SetMarkerStyle(20);
    gr_difference->SetMarkerSize(1.2);
    gr_difference->SetMarkerColor(kBlue);
    gr_difference->SetLineColor(kBlue);
    gr_difference->GetYaxis()->SetRangeUser(-1,1);
    gr_difference->Draw("APL");
    
    // Add zero line
    double min_x, max_x, min_y, max_y;
    gr_difference->GetPoint(0, min_x, min_y);
    gr_difference->GetPoint(energies.size()-1, max_x, max_y);
    TLine* zero_line = new TLine(min_x, 0, max_x, 0);
    zero_line->SetLineColor(kRed);
    zero_line->SetLineStyle(2);
    zero_line->SetLineWidth(2);
    zero_line->Draw("same");
    
    // Add difference statistics
    TPaveText* difference_info = new TPaveText(0.15, 0.75, 0.5, 0.9, "NDC");
    difference_info->AddText(Form("Max |Difference|: %.2f keV", max_difference));
    
    // Calculate RMS difference
    double rms_difference = 0;
    for (int i = 0; i < energies.size(); ++i) {
        double calibrated_energy = calib_slope * fitted_means[i] + calib_intercept;
        double difference = calibrated_energy - energies[i];
        rms_difference += difference * difference;
    }
    rms_difference = TMath::Sqrt(rms_difference / energies.size());
    
    difference_info->AddText(Form("RMS Difference: %.2f keV", rms_difference));
    difference_info->SetFillColor(kWhite);
    difference_info->SetBorderSize(1);
    difference_info->Draw();
    
    c4->Update();
    
    std::cout << "\n=== Energy Calibration Complete ===" << std::endl;
    std::cout << "Final calibration function: Energy = " << calib_slope << " * ADC + " << calib_intercept << std::endl;
    std::cout << "Energy range: " << min_energy << " to " << max_energy << " keV" << std::endl;
    std::cout << "Energy bin width: " << (max_energy - min_energy) / nBins << " keV/bin" << std::endl;
    
    std::cout << "\nRed dashed lines in energy spectrum show theoretical peak positions." << std::endl;
    std::cout << "\nEnergy calibration analysis complete!" << std::endl;
    
    // === SAVE CALIBRATION RESULTS ===
    std::cout << "\nDo you want to save the calibration function and energy histogram to ROOT file? (y/n): ";
    char save_choice;
    std::cin >> save_choice;
    
    if (save_choice == 'y' || save_choice == 'Y') {
        // Create output ROOT file
        TFile* outfile = new TFile("calibration_results.root", "RECREATE");
        
        if (outfile->IsOpen()) {
            std::cout << "Saving calibration results to 'calibration_results.root'..." << std::endl;
            
            // Save calibration function
            linear_fit->SetName("calibration_function");
            linear_fit->SetTitle("Energy Calibration Function: Energy = slope*ADC + intercept");
            linear_fit->Write();
            
            // Save energy-calibrated histograms
            h_energy_raw->Write();
            h_energy_norm->Write();
            
            // Save calibration graph
            gr_calib->SetName("calibration_graph");
            gr_calib->SetTitle("Energy Calibration Graph");
            gr_calib->Write();
            
            // Save original ADC histograms for reference
            h_raw->Write();
            h_norm->Write();
            
            outfile->Close();
            std::cout << "Calibration results saved successfully!" << std::endl;
            std::cout << "Saved objects:" << std::endl;
            std::cout << "  - calibration_function (TF1)" << std::endl;
            std::cout << "  - h_energy_raw (TH1D)" << std::endl;
            std::cout << "  - h_energy_norm (TH1D)" << std::endl;
            std::cout << "  - calibration_graph (TGraph)" << std::endl;
            std::cout << "  - h_raw (TH1D)" << std::endl;
            std::cout << "  - h_norm (TH1D)" << std::endl;
        } else {
            std::cerr << "Error: Could not create output file!" << std::endl;
        }
        
        delete outfile;
    } else {
        std::cout << "Calibration results not saved." << std::endl;
    }
    
    std::cout << "\nProgram completed." << std::endl;
}

