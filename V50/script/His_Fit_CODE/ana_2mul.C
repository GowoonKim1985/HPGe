{
	gROOT->SetStyle("Plain");

	char s_filename[256];
	int s_runnum;
	char s_runnumber[256];
	cout << "sample run num : ";
	scanf("%i", &s_runnum);
	sprintf(s_runnumber, "%06d", s_runnum);

	char b_filename[256];
	int b_runnum;
	char b_runnumber[256];
	cout << "BKG run num : ";
	scanf("%i", &b_runnum);
	sprintf(b_runnumber, "%06d", b_runnum);

	double bn = 2; // 1bin = 0.5 keV
	double time_day = 0, b_time_day;
	double e4 = 583.187, e5 = 2614.511; // Tl208
	double sig_e4 = 1.17122, sig_e5 = 1.62939; // RUN13678, Tl208
//	double e4 = 898.042, e5 = 1836.063; // Y88
//	double sig_e4 = 1.0208667, sig_e5 = 1.3530827; // Y88
	double fpeak_e4 = sig_e4*3, fpeak_e5 = sig_e5*3;

	sprintf(s_filename,"../result/RUN%d/Gana_%s_2mul.root", s_runnum, s_runnumber);
	sprintf(b_filename,"../result/RUN%d/Gana_%s_2mul.root", b_runnum, b_runnumber);

//	if(runnum == 90) { time_day = 504122./(60.*60.*24.); }
	if(s_runnum == 136) { time_day = 3384917./(60.*60.*24.); }
	if(s_runnum == 137) { time_day = 1883220./(60.*60.*24.); }
	if(s_runnum == 138) { time_day = 1640511./(60.*60.*24.); }
	if(s_runnum == 13678) { time_day = 6908648./(60.*60.*24.); }
	if(s_runnum == 1368) { time_day = 5025428./(60.*60.*24.); }


	if(b_runnum == 151) { b_time_day = 2358370./(60.*60.*24.); }
	if(b_runnum == 155) { b_time_day = 3324515.25/(60.*60.*24.); }
	if(b_runnum == 1515) { b_time_day = 5682885.25/(60.*60.*24.); }

	TFile *file1 = new TFile(s_filename);
	file1->ls();

	//tot_day->Draw("COLZ");
	//tot->Draw("COLZ");

	//c1->SetGrid();

	int e4_bin = mul2->GetXaxis()->FindBin(e4);
	int minb_e4 = mul2->GetXaxis()->FindBin(e4-fpeak_e4);
	int maxb_e4 = mul2->GetXaxis()->FindBin(e4+fpeak_e4);

	int e5_bin = mul2->GetXaxis()->FindBin(e5);
	int minb_e5 = mul2->GetXaxis()->FindBin(e5-fpeak_e5);
	int maxb_e5 = mul2->GetXaxis()->FindBin(e5+fpeak_e5);

	//	int e6_bin = (int)e6+1;

	printf("\n------------------------------------\n");

	// Tl208

	cout << endl;
	cout << "E1 = " << e4 << " keV, " << "E2 = " << e5 << " keV" << endl;
	cout << endl;
	cout << "min x of E1 keV = " << e4 - fpeak_e4 << " keV" << endl;
	cout << "max x of E1 keV = " << e4 + fpeak_e4 << " keV" << endl;
	cout << "bin of min x of E1 keV = " << minb_e4 << endl;
	cout << "bin of E1 keV = " << e4_bin << endl;
	cout << "bin of max x of E1 keV = " << maxb_e4 << endl;
	cout << endl;
	cout << "min x of E2 keV = " << e5 - fpeak_e5 << " keV" << endl;
	cout << "max x of E2 keV = " << e5 + fpeak_e5 << " keV" << endl;
	cout << "bin of min x of E2 keV = " << minb_e5 << endl;
	cout << "bin of E2 keV = " << e5_bin << endl;
	cout << "bin of max x of E2 keV = " << maxb_e5 << endl;

	//EC P1, 894.979 ~ 901.105 keV, 1832 ~ 1840.12 keV
	double P1 = mul2->Integral(minb_e5, maxb_e5, minb_e4, maxb_e4);
	double P1_cpd = P1/time_day;
	double P1_dcpd = sqrt(P1)/time_day;

	printf("\n------------------------------------\n");
	printf("\nMeasrued time : %0.2f days \n", time_day);
	printf("\nP1 count : %0.2f \n", P1);
	printf("P1 cpd : %0.2f \n", P1_cpd);
	printf("P1 dcpd : %0.2f \n \n", P1_dcpd);


	//Y88 P2, 0 ~ 1832 keV, 894.979 ~ 901.105 keV
	double P2_1 = mul2->Integral(minb_e4, maxb_e4, 0, minb_e5);
	double P2_2 = mul2->Integral(0, minb_e5, minb_e4, maxb_e4);
	double P2_3 = mul2->Integral(minb_e4, maxb_e4, minb_e4, maxb_e4);
	double P2 = P2_1 + P2_2 - P2_3;
	double P2_cpd = P2/time_day;
	double P2_dcpd = sqrt(P2)/time_day;

	printf("583 keV full deposit");
	printf("P2 count : %0.2f \n", P2);
	printf("P2 cpd : %0.2f \n", P2_cpd);
	printf("P2 dcpd : %0.2f \n \n", P2_dcpd);


	//Y88 P3, 0 ~ 594.979 keV, 1832 ~ 1840.12 keV
	double P3 = mul2->Integral(minb_e5, maxb_e5, 0, minb_e4);
	double P3_cpd = P3/time_day;
	double P3_dcpd = sqrt(P3)/time_day;

	printf("2614 keV full deposit");
	printf("P3 count : %0.2f \n", P3);
	printf("P3 cpd : %0.2f \n", P3_cpd);
	printf("P3 dcpd : %0.2f \n \n", P3_dcpd);
/*
	if(runnum == 90) {

		P1_tl = 2.*P1_tl;
		P2_tl = 2.*P2_tl;
		P3_tl = 2.*P3_tl;

		P1_cpd = 2.*P1_cpd;
		P2_cpd = 2.*P2_cpd;
		P3_cpd = 2.*P3_cpd;

		P1_dcpd = sqrt(2.)*P1_dcpd;
		P2_dcpd = sqrt(2.)*P2_dcpd;
		P3_dcpd = sqrt(2.)*P3_dcpd;

		printf("\nMeasured with 7 detectors\n\n"); 

		printf("Y88 P1 count (normalized) : %0.2f \n", P1_tl);
		printf("Y88 P1 cpd (normalized) : %0.2f \n", P1_cpd);
		printf("Y88 P1 dcpd (normalized) : %0.2f \n \n", P1_dcpd);

		printf("Y88 P2 count (normalized) : %0.2f \n", P2_tl);
		printf("Y88 P2 cpd (normalized) : %0.2f \n", P2_cpd);
		printf("Y88 P2 dcpd (normalized) : %0.2f \n \n", P2_dcpd);

		printf("Y88 P3 count (normalized) : %0.2f \n", P3_tl);
		printf("Y88 P3 cpd (normalized) : %0.2f \n", P3_cpd);
		printf("Y88 P3 dcpd (normalized) : %0.2f \n \n", P3_dcpd);
	}
*/
	printf("------------------------------------\n");

	TFile *file2 = new TFile(b_filename);
	file2->ls();

	int b_e4_bin = mul2->GetXaxis()->FindBin(e4);
	int b_minb_e4 = mul2->GetXaxis()->FindBin(e4-fpeak_e4);
	int b_maxb_e4 = mul2->GetXaxis()->FindBin(e4+fpeak_e4);

	int b_e5_bin = mul2->GetXaxis()->FindBin(e5);
	int b_minb_e5 = mul2->GetXaxis()->FindBin(e5-fpeak_e5);
	int b_maxb_e5 = mul2->GetXaxis()->FindBin(e5+fpeak_e5);

	//	int e6_bin = (int)e6+1;

	printf("\n------------------------------------\n");


	cout << "BKG" << endl;
/*	cout << "min x of 898 keV = " << e4 - fpeak_e4 << " keV" << endl;
	cout << "max x of 898 keV = " << e4 + fpeak_e4 << " keV" << endl;
	cout << "bin of min x of 898 keV = " << b_minb_e4 << endl;
	cout << "bin of 898 keV = " << b_e4_bin << endl;
	cout << "bin of max x of 898 keV = " << b_maxb_e4 << endl;

	cout << "min x of 1836 keV = " << e5 - fpeak_e5 << " keV" << endl;
	cout << "max x of 1836 keV = " << e5 + fpeak_e5 << " keV" << endl;
	cout << "bin of min x of 1836 keV = " << b_minb_e5 << endl;
	cout << "bin of 1836 keV = " << b_e5_bin << endl;
	cout << "bin of max x of 1836 keV = " << b_maxb_e5 << endl;
*/

	//EC P1, 894.979 ~ 901.105 keV, 1832 ~ 1840.12 keV
	double b_P1 = mul2->Integral(b_minb_e5, b_maxb_e5, b_minb_e4, b_maxb_e4);
	double b_P1_cpd = b_P1/b_time_day;
	double b_P1_dcpd = sqrt(b_P1)/b_time_day;

//	printf("\n------------------------------------\n");
	printf("\nMeasrued time : %0.2f days \n", b_time_day);
	printf("\nP1 count : %0.2f \n", b_P1);
	printf("P1 cpd : %0.2f \n", b_P1_cpd);
	printf("P1 dcpd : %0.2f \n \n", b_P1_dcpd);


	//Y88 P2, 0 ~ 1832 keV, 894.979 ~ 901.105 keV
	double b_P2_1 = mul2->Integral(b_minb_e4, b_maxb_e4, 0, b_minb_e5);
	double b_P2_2 = mul2->Integral(0, b_minb_e5, b_minb_e4, b_maxb_e4);
	double b_P2_3 = mul2->Integral(b_minb_e4, b_maxb_e4, b_minb_e4, b_maxb_e4);
	double b_P2 = b_P2_1 + b_P2_2 - b_P2_3;
	double b_P2_cpd = b_P2/b_time_day;
	double b_P2_dcpd = sqrt(b_P2)/b_time_day;

	printf("583 keV full deposit");
	printf("P2 count : %0.2f \n", b_P2);
	printf("P2 cpd : %0.2f \n", b_P2_cpd);
	printf("P2 dcpd : %0.2f \n \n", b_P2_dcpd);

	double b_P3 = mul2->Integral(b_minb_e5, b_maxb_e5, 0, b_minb_e4);
	double b_P3_cpd = b_P3/b_time_day;
	double b_P3_dcpd = sqrt(b_P3)/b_time_day;

	printf("2614 keV full deposit");
	printf("P3 count : %0.2f \n", b_P3);
	printf("P3 cpd : %0.2f \n", b_P3_cpd);
	printf("P3 dcpd : %0.2f \n \n", b_P3_dcpd);


	///////////////// TGraph with cpd, dcpd

	int a=3;

	double en[3] = {1, 2, 3};
	double den[3] = {};

	double cpd[3] = {P1_cpd, P2_cpd, P3_cpd};
	double dcpd[3] = {P1_dcpd, P2_dcpd, P3_dcpd};

	double b_cpd[3] = {b_P1_cpd, b_P2_cpd, b_P3_cpd};
	double b_dcpd[3] = {b_P1_dcpd, b_P2_dcpd, b_P3_dcpd};

	TGraph *gr = new TGraphErrors(a, en, cpd, den, dcpd);
	TGraph *b_gr = new TGraphErrors(a, en, b_cpd, den, b_dcpd);
	TMultiGraph *mg = new TMultiGraph();

	gr->SetMarkerStyle(21);
	gr->SetMarkerSize(1.5);
	gr->SetMarkerColor(kRed);
	gr->SetLineColor(kRed);
	gr->GetYaxis()->SetTitle("counts/day");
	gr->GetXaxis()->SetTitle("");

	b_gr->SetMarkerStyle(23);
	b_gr->SetMarkerSize(1.5);
	b_gr->SetMarkerColor(kBlack);
	b_gr->SetLineColor(kBlack);
	b_gr->GetYaxis()->SetTitle("counts/day");
	b_gr->GetXaxis()->SetTitle("");

	TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1200);
	
	c1->cd(1);
	c1->SetGrid();
	mg->Add(gr);
	mg->Add(b_gr);
	mg->Draw("ap");

}
