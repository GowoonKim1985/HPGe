#include <iostream>
#include "TMath.h"
#include "TMinuit.h"

void ana_double_peak()
{

	char iso[256];
	char iso_dat[256];
/*	cout << "tho2? all? : ";
	scanf("%s", iso);
	sprintf(iso_dat, "/home/kkw/DAQ/ANA25/script/dat/%s_peak.dat", iso);
*/
	sprintf(iso_dat, "/home/kkw/DAQ/ANA25/script/dat/tho2_double_peak.dat");

	ifstream peakdata1(iso_dat); 
	string line_n1;
	string line_n2;

	int line_max1 = 0, line_max2 = 0;
	double file, evt;
	char rtemp1[256]; //nuclear info.
	double rtemp2; //peak energy

	while (peakdata1.good()) {
		getline(peakdata1, line_n1);
		peakdata1 >> rtemp1 >> rtemp2;
		++line_max1;
	}

	cout << "line max : " << line_max1 << endl;

	char nu[line_max1-1][256];
	double peak[line_max1-1];

	double readfile, readevt;
	ifstream peakdata2(iso_dat);

	while (peakdata2.good()) {
		if ( line_max2 < (line_max1-1) ) {
			getline(peakdata2, line_n2);
			peakdata2 >> rtemp1 >> rtemp2;
			sprintf(nu[line_max2], "%s", rtemp1);
			peak[line_max2] = rtemp2;
			cout << "line " << line_max2 << " : " << nu[line_max2] << " , " << peak[line_max2] << endl;
		}
		else { break; }
		++line_max2;
	}

	char hisfile[256];
	char rawname[256];
	char res2file[256]; //fit his
	char res3file[256]; //fit result
/*
	int runnum;
	char runnumber[256];
	cout << "run number : ";
	scanf("%i", &runnum);
	sprintf(runnumber, "%06d", runnum);

	char evt_case[256];
	cout << "event case (1mul / 2mul / tot) : ";
	scanf("%s", evt_case);
*/
//	int bin=8000;
	int bin=4000;
	//	cout << "binning (4000, 1bin=1keV / 8000, 1bin=0.5keV) : ";
	//	scanf("%d", &bin);
/*
	int det_num;
	cout << "the number of detector : ";
	scanf("%i", &det_num);
*/


	sprintf(hisfile, "/home/kkw/DAQ/ANA25/histogram/SP02/FWHM310_1kev/phis_tho2.root");
	sprintf(res2file, "/home/kkw/DAQ/ANA25/histogram/SP02/FIT310_1kev/fit_double.C");
	sprintf(res3file, "/home/kkw/DAQ/ANA25/histogram/SP02/FIT310_1kev/fit_double.root");
/*
	double time_day=0;

	if(runnum == 84) {time_day = 1974210./(60.*60.*24.);}
	if(runnum == 136) {time_day = 3384917./(60.*60.*24.);}
	if(runnum == 137) {time_day = 1883220./(60.*60.*24.);}
	if(runnum == 138) {time_day = 1640511./(60.*60.*24.);}
	if(runnum == 151) {time_day = 2358370./(60.*60.*24.);}
	if(runnum == 155) {time_day = 3324515.25/(60.*60.*24.);}
	if(runnum == 13678) {time_day = 6908648./(60.*60.*24.);}
	if(runnum == 1368) {time_day = 5025428./(60.*60.*24.);}
	if(runnum == 1515) {time_day = 5682885.25/(60.*60.*24.);}
	if(runnum == 21490) {time_day = 8463510./(60*60*24); }

	cout << "measured day = " << time_day << " days" << endl;
*/


	TF1 * gfit = new TF1("gfit", "gaus");
	TF1 * agfit = new TF1("agfit", "[0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2])");
	agfit->SetParNames("Area", "Mean", "Sigma"); //area = sqrt(2pi)*constant*sigma

	TF1 * efit = new TF1("efit", "expo");
	TF1 * e2fit = new TF1("e2fit", "expo");
	TF1 * p0fit = new TF1("p0fit", "pol0");
	TF1 * p1fit = new TF1("p1fit", "pol1");

	TF1 * m3fit = new TF1("m3fit", "[0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]) + [3]/(sqrt(2*TMath::Pi())*[5])*exp(-((x-[4])*(x-[4]))/2/[5]/[5]) + [6]");
	m3fit->SetParNames("Area1", "Mean1", "Sigma1", "Area2", "Mean2", "Sigma2", "Flat");

	TF1 * m4fit = new TF1("m4fit", "[0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]) + [3]/(sqrt(2*TMath::Pi())*[5])*exp(-((x-[4])*(x-[4]))/2/[5]/[5]) + exp([6]+[7]*x)");
	m4fit->SetParNames("Area1", "Mean1", "Sigma1", "Area2", "Mean2", "Sigma2", "expo1", "expo2");

	TF1 * m5fit = new TF1("m5fit", "[0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]) + [3]/(sqrt(2*TMath::Pi())*[5])*exp(-((x-[4])*(x-[4]))/2/[5]/[5]) + [6]+[7]*x");
	m5fit->SetParNames("Area1", "Mean1", "Sigma1", "Area2", "Mean2", "Sigma2", "Intercept", "Slope");

//	TCanvas * c1 = new TCanvas("c1", "CC1 Spectrum w/ fwhm cut", 800, 800);
	TCanvas * c1 = new TCanvas("c1", "CC1 Spectrum w/ fwhm cut", 1200, 800);

	agfit->SetLineColor(9);
	gfit->SetLineColor(3);
	efit->SetLineColor(5);
	e2fit->SetLineColor(6);
	p1fit->SetLineColor(2);
	m3fit->SetLineColor(kBlue);
	m4fit->SetLineColor(kRed);
	m5fit->SetLineColor(kGreen);


	int maxb, maxx, maxh;
	double bg1, bg2, bg;
	double cnt1, cnt2, cnt3, cnt1_2, cnt3_2;
	double pa1[line_max1-1][3];//gfit 
	//double pa2[line_max1-1][4];//m1fit
	//double pa2er[line_max1-1][4];//m1fit error
	//	double pa3[line_max1-1][7];//integral point1,2,3,4 & sum 1,2,3
	//	double pa4[line_max1-1][5];//m2fit
	//	double pa4er[line_max1-1][5];//m1fit error
	double cal[line_max1-1][2]; //integral fit func - m1fit,m2fit
	double epa[6]; //expo fit 
	int xbin[line_max1-1];
	double energy;

	double agfitpa1, agfitpa2, agfitpa3, agfiter1, agfiter2, agfiter3;

	double gfitpa1, gfitpa2, gfitpa3, gfiter1, gfiter2, gfiter3;

	double p1fitpa1, p1fitpa2, p1fiter1, p1fiter2;

	double m3fitpa1, m3fitpa2, m3fitpa3, m3fitpa4, m3fitpa5, m3fitpa6, m3fitpa7;
	double m3fiter1, m3fiter2, m3fiter3, m3fiter4, m3fiter5, m3fiter6, m3fiter7;
	double m3pcnt1, m3pcnter1, m3cpd1, m3dcpd1;
	double m3pcnt2, m3pcnter2, m3cpd2, m3dcpd2;

	double m4fitpa1, m4fitpa2, m4fitpa3, m4fitpa4, m4fitpa5, m4fitpa6, m4fitpa7, m4fitpa8;
	double m4fiter1, m4fiter2, m4fiter3, m4fiter4, m4fiter5, m4fiter6, m4fiter7, m4fiter8;
	double m4pcnt1, m4pcnter1, m4cpd1, m4dcpd1;
	double m4pcnt2, m4pcnter2, m4cpd2, m4dcpd2;

	double m5fitpa1, m5fitpa2, m5fitpa3, m5fitpa4, m5fitpa5, m5fitpa6, m5fitpa7, m5fitpa8;
	double m5fiter1, m5fiter2, m5fiter3, m5fiter4, m5fiter5, m5fiter6, m5fiter7, m5fiter8;
	double m5pcnt1, m5pcnter1, m5cpd1, m5dcpd1;
	double m5pcnt2, m5pcnter2, m5cpd2, m5dcpd2;

	double gfitchi2, agfitchi2, p1fitchi2, m3fitchi2, m4fitchi2, m5fitchi2;



	TString m1stat;

	TTree *gfitpara = new TTree("gfitpara", "gaus fitting result");
//	gfitpara -> Branch("run", &runnum, "run/I");
	gfitpara -> Branch("peak", &energy, "energy/D");
	gfitpara -> Branch("const", &gfitpa1, "gfitpa1/D");
	gfitpara -> Branch("mean", &gfitpa2, "gfitpa2/D");
	gfitpara -> Branch("chi2", &gfitchi2, "gfitchi2/D");
	gfitpara -> Branch("const_er", &gfiter1, "gfiter1/D");
	gfitpara -> Branch("mean_er", &gfiter2, "gfiter2/D");
	gfitpara -> Branch("sigma_er", &gfiter3, "gfiter3/D");

	TTree *p1fitpara = new TTree("p1fitpara", "pol1 fitting result");
//	p1fitpara -> Branch("run", &runnum, "run/I");
	p1fitpara -> Branch("peak", &energy, "energy/D");
	p1fitpara -> Branch("intercept", &p1fitpa1, "p1fitpa1/D");
	p1fitpara -> Branch("slope", &p1fitpa2, "p1fitpa2/D");
	p1fitpara -> Branch("chi2", &p1fitchi2, "p1fitchi2/D");
	p1fitpara -> Branch("intercept_er", &p1fiter1, "p1fiter1/D");
	p1fitpara -> Branch("slope_er", &p1fiter2, "p1fiter2/D");

	TTree *agfitpara = new TTree("agfitpara", "area-gaus fitting result");
//	agfitpara -> Branch("run", &runnum, "run/I");
	agfitpara -> Branch("peak", &energy, "energy/D");
	agfitpara -> Branch("area", &agfitpa1, "gfitpa1/D");
	agfitpara -> Branch("mean", &agfitpa2, "gfitpa2/D");
	agfitpara -> Branch("chi2", &agfitchi2, "gfitchi2/D");
	agfitpara -> Branch("area_er", &agfiter1, "gfiter1/D");
	agfitpara -> Branch("mean_er", &agfiter2, "gfiter2/D");
	agfitpara -> Branch("sigma_er", &agfiter3, "gfiter3/D");

	TTree *m3fitpara = new TTree("m3fitpara", "gaus+gaus+pol0 fitting result");
//	m3fitpara -> Branch("run", &runnum, "run/I");
	m3fitpara -> Branch("peak", &energy, "energy/D");
	m3fitpara -> Branch("area1", &m3fitpa1, "m3fitpa1/D");
	m3fitpara -> Branch("mean1", &m3fitpa2, "m3fitpa2/D");
	m3fitpara -> Branch("sigma1", &m3fitpa3, "m3fitpa3/D");
	m3fitpara -> Branch("area2", &m3fitpa4, "m3fitpa4/D");
	m3fitpara -> Branch("mean2", &m3fitpa5, "m3fitpa5/D");
	m3fitpara -> Branch("sigma2", &m3fitpa6, "m3fitpa6/D");
	m3fitpara -> Branch("chi2", &m3fitchi2, "m3fitchi2/D");
	m3fitpara -> Branch("pol0", &m3fitpa7, "m3fitpa7/D");
	m3fitpara -> Branch("area1_er", &m3fiter1, "m3fiter1/D");
	m3fitpara -> Branch("mean1_er", &m3fiter2, "m3fiter2/D");
	m3fitpara -> Branch("sigma1_er", &m3fiter3, "m3fiter3/D");
	m3fitpara -> Branch("area2_er", &m3fiter4, "m3fiter4/D");
	m3fitpara -> Branch("mean2_er", &m3fiter5, "m3fiter5/D");
	m3fitpara -> Branch("sigma2_er", &m3fiter6, "m3fiter6/D");
	m3fitpara -> Branch("pol0_er", &m3fiter7, "m3fiter7/D");
	m3fitpara -> Branch("pcount1", &m3pcnt1, "m3pcnt1/D");
	m3fitpara -> Branch("pcount2", &m3pcnt2, "m3pcnt2/D");
	m3fitpara -> Branch("pcount1_er", &m3pcnter1, "m3pcnter1/D");
	m3fitpara -> Branch("pcount2_er", &m3pcnter2, "m3pcnter2/D");
//	m3fitpara -> Branch("cpd1", &m3cpd1, "m3cpd1/D");
//	m3fitpara -> Branch("cpd2", &m3cpd2, "m3cpd2/D");
//	m3fitpara -> Branch("dcpd1", &m3dcpd1, "m3dcpd1/D");
//	m3fitpara -> Branch("dcpd2", &m3dcpd2, "m3dcpd2/D");

	TTree *m4fitpara = new TTree("m4fitpara", "gaus+gaus+exp fitting result");
//	m4fitpara -> Branch("run", &runnum, "run/I");
	m4fitpara -> Branch("peak", &energy, "energy/D");
	m4fitpara -> Branch("area1", &m4fitpa1, "m4fitpa1/D");
	m4fitpara -> Branch("mean1", &m4fitpa2, "m4fitpa2/D");
	m4fitpara -> Branch("sigma1", &m4fitpa3, "m4fitpa3/D");
	m4fitpara -> Branch("area2", &m4fitpa4, "m4fitpa4/D");
	m4fitpara -> Branch("mean2", &m4fitpa5, "m4fitpa5/D");
	m4fitpara -> Branch("sigma2", &m4fitpa6, "m4fitpa6/D");
	m4fitpara -> Branch("chi2", &m4fitchi2, "m4fitchi2/D");
	m4fitpara -> Branch("exp1", &m4fitpa7, "m4fitpa7/D");
	m4fitpara -> Branch("exp2", &m4fitpa8, "m4fitpa8/D");
	m4fitpara -> Branch("area1_er", &m4fiter1, "m4fiter1/D");
	m4fitpara -> Branch("mean1_er", &m4fiter2, "m4fiter2/D");
	m4fitpara -> Branch("sigma1_er", &m4fiter3, "m4fiter3/D");
	m4fitpara -> Branch("area2_er", &m4fiter4, "m4fiter4/D");
	m4fitpara -> Branch("mean2_er", &m4fiter5, "m4fiter5/D");
	m4fitpara -> Branch("sigma2_er", &m4fiter6, "m4fiter6/D");
	m4fitpara -> Branch("exp1_er", &m4fiter7, "m4fiter7/D");
	m4fitpara -> Branch("exp2_er", &m4fiter8, "m4fiter8/D");
	m4fitpara -> Branch("pcount1", &m4pcnt1, "m4pcnt1/D");
	m4fitpara -> Branch("pcount2", &m4pcnt2, "m4pcnt2/D");
	m4fitpara -> Branch("pcount1_er", &m4pcnter1, "m4pcnter1/D");
	m4fitpara -> Branch("pcount2_er", &m4pcnter2, "m4pcnter2/D");
//	m4fitpara -> Branch("cpd1", &m4cpd1, "m4cpd1/D");
//	m4fitpara -> Branch("cpd2", &m4cpd2, "m4cpd2/D");
//	m4fitpara -> Branch("dcpd1", &m4dcpd1, "m4dcpd1/D");
//	m4fitpara -> Branch("dcpd2", &m4dcpd2, "m4dcpd2/D");

	TTree *m5fitpara = new TTree("m5fitpara", "gaus+gaus+pol1 fitting result");
//	m5fitpara -> Branch("run", &runnum, "run/I");
	m5fitpara -> Branch("peak", &energy, "energy/D");
	m5fitpara -> Branch("area1", &m5fitpa1, "m5fitpa1/D");
	m5fitpara -> Branch("mean1", &m5fitpa2, "m5fitpa2/D");
	m5fitpara -> Branch("sigma1", &m5fitpa3, "m5fitpa3/D");
	m5fitpara -> Branch("area2", &m5fitpa4, "m5fitpa4/D");
	m5fitpara -> Branch("mean2", &m5fitpa5, "m5fitpa5/D");
	m5fitpara -> Branch("sigma2", &m5fitpa6, "m5fitpa6/D");
	m5fitpara -> Branch("chi2", &m5fitchi2, "m5fitchi2/D");
	m5fitpara -> Branch("pol0", &m5fitpa7, "m5fitpa7/D");
	m5fitpara -> Branch("pol1", &m5fitpa8, "m5fitpa8/D");
	m5fitpara -> Branch("area1_er", &m5fiter1, "m5fiter1/D");
	m5fitpara -> Branch("mean1_er", &m5fiter2, "m5fiter2/D");
	m5fitpara -> Branch("sigma1_er", &m5fiter3, "m5fiter3/D");
	m5fitpara -> Branch("area2_er", &m5fiter4, "m5fiter4/D");
	m5fitpara -> Branch("mean2_er", &m5fiter5, "m5fiter5/D");
	m5fitpara -> Branch("sigma2_er", &m5fiter6, "m5fiter6/D");
	m5fitpara -> Branch("pol0_er", &m5fiter7, "m5fiter7/D");
	m5fitpara -> Branch("pol1_er", &m5fiter8, "m5fiter8/D");
	m5fitpara -> Branch("pcount1", &m5pcnt1, "m5pcnt1/D");
	m5fitpara -> Branch("pcount2", &m5pcnt2, "m5pcnt2/D");
	m5fitpara -> Branch("pcount1_er", &m5pcnter1, "m5pcnter1/D");
	m5fitpara -> Branch("pcount2_er", &m5pcnter2, "m5pcnter2/D");
//	m5fitpara -> Branch("cpd1", &m5cpd1, "m5cpd1/D");
//	m5fitpara -> Branch("cpd2", &m5cpd2, "m5cpd2/D");
//	m5fitpara -> Branch("dcpd1", &m5dcpd1, "m5dcpd1/D");
//	m5fitpara -> Branch("dcpd2", &m5dcpd2, "m5dcpd2/D");

	int mul = bin/4000;

	TFile * hf = new TFile(hisfile);
	TH1D * his_temp = new TH1D("his_temp", "", bin, 0, 4000);
	char data_his[256];
//	sprintf(data_his, "%s_cnt", evt_case);
	his_temp = (TH1D*)hf->Get("fhis_nctot");
	his_temp -> SetLineColor(1);

	c1 -> cd();
/*
	efit -> SetRange(500, 2200);
	his_temp -> Fit("efit", "RQ0+");
	efit -> GetParameters(&epa[0]);
	efit -> SetRange(2500, 2800);
	his_temp -> Fit("efit", "RQ0+");
	efit -> GetParameters(&epa[2]);
	efit -> SetRange(3000, 4000);
	his_temp -> Fit("efit", "RQ0+");
	efit -> GetParameters(&epa[4]);
*/
	int odd_check = 1;
	double fixpa[3];//expo fit par
	double gpa1[line_max1]; //areas
	double gpa2[line_max1]; // sigmas

	for(int i=0; i<(line_max1-1); i++) {

		energy = peak[i];
		xbin[i] = his_temp->GetXaxis()->FindBin(peak[i]);
		printf("\n");
		printf("\n=========================\n");
		printf("%.2f keV bin : %i	\n", peak[i], xbin[i]);


		his_temp -> GetXaxis() -> SetRange(xbin[i]-5*mul, xbin[i]+5*mul);
		maxb = his_temp->GetMaximumBin();
		maxx = maxb/mul;
		cout << "peak[i] : "<< peak[i] <<endl;
		cout << "maxb(bin) : " << maxb << endl;
		cout << "maxx(kev) : " << maxx << endl;
		maxh = his_temp->GetMaximum();

		his_temp->GetXaxis()->SetRange(0, 4000*mul);

		int range_val = 20;

		gfit->SetRange((maxx-range_val),(maxx+range_val));
		agfit->SetRange((maxx-range_val),(maxx+range_val));
		p1fit->SetRange((maxx-range_val),(maxx+range_val));
		m3fit->SetRange((maxx-range_val),(maxx+range_val));
		m4fit->SetRange((maxx-range_val),(maxx+range_val));
		m5fit->SetRange((maxx-range_val),(maxx+range_val));
	
		efit->SetRange((maxx-150),(maxx+150));	
		his_temp->Fit("efit","RQ0+");
		fixpa[1] = efit->GetParameter(0);
		fixpa[2] = efit->GetParameter(1);

		//gaus fit
	//	gfit->FixParameter(1,peak[i]);
		his_temp->Fit("gfit","RQ0+");
		gfitpa1 = gfit->GetParameter(0);
		gfitpa2 = gfit->GetParameter(1);
		gfitpa3 = gfit->GetParameter(2);
		gfiter1 = gfit->GetParError(0);
		gfiter2 = gfit->GetParError(1);
		gfiter3 = gfit->GetParError(2);
		
		gpa1[i] = gfitpa1;
		gpa2[i] = gfitpa3;

		gfitchi2 = (gfit->GetChisquare())/(gfit->GetNDF());

		bg1 = 0; bg2 = 0; 
		for(int j1=0; j1<(5*mul); j1++) {
			bg1 = bg1 + (his_temp->GetBinContent((maxb)-(gfitpa3*3+5)*mul+j1));
			bg2 = bg2 + (his_temp->GetBinContent((maxb)+(gfitpa3*3+5)*mul-j1));
			//		cout<<maxx-10+j<<"kev: "<<bg1<<" / "<<maxx+10-j<<"kev : "<<bg2<<endl;
		}
		bg=(bg1+bg2)/(10*mul);


		//pol1 fit
		his_temp->Fit("p1fit","RQ0+");
		p1fitpa1 = p1fit->GetParameter(0);
		p1fitpa2 = p1fit->GetParameter(1);
		p1fiter1 = p1fit->GetParError(0);
		p1fiter2 = p1fit->GetParError(1);

		p1fitchi2 = (p1fit->GetChisquare())/(p1fit->GetNDF());


		//area gaus fit

		agfit->SetParameters(2.5*gfitpa2*gfitpa3, gfitpa2, gfitpa3);
		//agfit->SetParameters(2.5*gfitpa2*gfitpa3,maxx,2);
		his_temp->Fit("agfit","RQ0+");
		agfitpa1 = agfit->GetParameter(0);
		agfitpa2 = agfit->GetParameter(1);
		agfitpa3 = agfit->GetParameter(2);
		agfiter1 = agfit->GetParError(0);
		agfiter2 = agfit->GetParError(1);
		agfiter3 = agfit->GetParError(2);
		agfitchi2 = (agfit->GetChisquare())/(agfit->GetNDF());


	if(i==odd_check){	

cout<<"///////m3 fit///////"<<endl;
cout<<endl;

		odd_check = odd_check + 2;

		//m3fit
		//m3fit -> SetParameters(gfitpa1, peak[i], gfitpa3, gfitpa1, peak[i]+10, gfitpa3, bg);
		//m3fit -> SetParameters(gpa1[i-1], peak[i-1], gpa2[i-1], gpa1[i], peak[i], gpa2[i], bg);
		//m3fit -> SetParameters(gpa1[i-1], peak[i-1], gpa2[i-1], gpa1[i], peak[i], gpa2[i], bg);
		m3fit -> SetParameters(gfitpa1, peak[i-1], gfitpa3, gfitpa1, peak[i], gfitpa3, bg);
		
		his_temp -> Fit("m3fit", "R+"); 
		m3fitpa1 = m3fit->GetParameter(0);
		m3fitpa2 = m3fit->GetParameter(1);
		m3fitpa3 = m3fit->GetParameter(2);
		m3fitpa4 = m3fit->GetParameter(3);
		m3fitpa5 = m3fit->GetParameter(4);
		m3fitpa6 = m3fit->GetParameter(5);
		m3fitpa7 = m3fit->GetParameter(6);
		m3fiter1 = m3fit->GetParError(0);	
		m3fiter2 = m3fit->GetParError(1);	
		m3fiter3 = m3fit->GetParError(2);	
		m3fiter4 = m3fit->GetParError(3);	
		m3fiter5 = m3fit->GetParError(4);	
		m3fiter6 = m3fit->GetParError(5);	
		m3fiter7 = m3fit->GetParError(6);	
		m3fitchi2 = (m3fit->GetChisquare())/(m3fit->GetNDF());
/*
		m3pcnt1 = m3fitpa1*mul;
		m3pcnter1 = sqrt(m3fitpa1*mul);
		m3cpd1 = m3pcnt1/time_day;
		m3dcpd1 = m3pcnter1/time_day;

		m3pcnt2 = m3fitpa4*mul;
		m3pcnter2 = sqrt(m3fitpa4*mul);
		m3cpd2 = m3pcnt2/time_day;
		m3dcpd2 = m3pcnter2/time_day;
*/
		cout << "[m3fit]" << endl;
		cout << "peak1 : " << m3fitpa2 << endl;
		cout << "m3fit tot. counts (peak1) :" << m3pcnt1 << " +/- " << m3pcnter1 << endl;
//		cout << "m3fit cpd (peak1) :" << m3cpd1 << " +/- " << m3dcpd1 << endl;
		cout << "6 sigma : " << m3fitpa3*6 << endl;
		cout << "peak2 : " << m3fitpa5 << endl;
		cout << "m3fit tot. counts (peak2) :" << m3pcnt2 << " +/- " << m3pcnter2 << endl;
//		cout << "m3fit cpd (peak2) :" << m3cpd2 << " +/- " << m3dcpd2 << endl;
		cout << "6 sigma : " << m3fitpa6*6 << endl;
		cout << "chi2 : " << m3fitchi2 << endl;

		m3fitpara->Fill();

cout<<endl;
cout<<"///////m4 fit///////"<<endl;
cout<<endl;
		//m4fit
		//m4fit -> SetParameters(gfitpa1, peak[i], gfitpa3, gfitpa1, peak[i+1], gfitpa3, fixpa[1], fixpa[2]);
		//m4fit -> SetParameters(gfitpa1, peak[i], gfitpa3, gfitpa1, peak[i-1], gfitpa3, fixpa[1], fixpa[2]);
		//m4fit -> SetParameters(gpa1[i-1], peak[i-1], gpa2[i-1], gpa1[i], peak[i], gpa2[i], fixpa[1], fixpa[2]);
		m4fit -> SetParameters(gfitpa1, peak[i-1], gfitpa3, gfitpa1, peak[i], gfitpa3, fixpa[1], fixpa[2]);
/*
		if (peak[i]<2200) { 
			m4fit -> SetParameters(gfitpa1, peak[i]-10, gfitpa3, gfitpa1, peak[i], gfitpa3, epa[0], epa[1]); }
		else if (peak[i]>2500 && peak[i]<2800) { 
			m4fit -> SetParameters(gfitpa1, peak[i], gfitpa3, gfitpa1, peak[i]+5, gfitpa3, epa[2], epa[3]); }
		else if (peak[i]>3000) { 
			m4fit -> SetParameters(gfitpa1, peak[i], gfitpa3, gfitpa1, peak[i]+5, gfitpa3, epa[4], epa[5]); }
		//m4fit -> SetParameters(gfitpa1, peak[i]-5, gfitpa3, gfitpa1, peak[i]+10, gfitpa3, bg); // RUN1515, 238 keV
*/

		his_temp -> Fit("m4fit", "R+"); 
		m4fitpa1 = m4fit->GetParameter(0);
		m4fitpa2 = m4fit->GetParameter(1);
		m4fitpa3 = m4fit->GetParameter(2);
		m4fitpa4 = m4fit->GetParameter(3);
		m4fitpa5 = m4fit->GetParameter(4);
		m4fitpa6 = m4fit->GetParameter(5);
		m4fitpa7 = m4fit->GetParameter(6);
		m4fitpa8 = m4fit->GetParameter(7);
		m4fiter1 = m4fit->GetParError(0);	
		m4fiter2 = m4fit->GetParError(1);	
		m4fiter3 = m4fit->GetParError(2);	
		m4fiter4 = m4fit->GetParError(3);	
		m4fiter5 = m4fit->GetParError(4);	
		m4fiter6 = m4fit->GetParError(5);	
		m4fiter7 = m4fit->GetParError(6);	
		m4fiter8 = m4fit->GetParError(7);	
		m4fitchi2 = (m4fit->GetChisquare())/(m4fit->GetNDF());
/*
		m4pcnt1 = m4fitpa1*mul;
		m4pcnter1 = sqrt(m4fitpa1*mul);
		m4cpd1 = m4pcnt1/time_day;
		m4dcpd1 = m4pcnter1/time_day;

		m4pcnt2 = m4fitpa4*mul;
		m4pcnter2 = sqrt(m4fitpa4*mul);
		m4cpd2 = m4pcnt2/time_day;
		m4dcpd2 = m4pcnter2/time_day;
*/
		cout << "[m4fit]" << endl;
		cout << "peak1 : " << m4fitpa2 << endl;
		cout << "m4fit tot. counts (peak1) :" << m4pcnt1 << " +/- " << m4pcnter1 << endl;
//		cout << "m4fit cpd (peak1) :" << m4cpd1 << " +/- " << m4dcpd1 << endl;
		cout << "6 sigma : " << m4fitpa3*6 << endl;
		cout << "peak2 : " << m4fitpa5 << endl;
		cout << "m4fit tot. counts (peak2) :" << m4pcnt2 << " +/- " << m4pcnter2 << endl;
//		cout << "m4fit cpd (peak2) :" << m4cpd2 << " +/- " << m4dcpd2 << endl;
		cout << "6 sigma : " << m4fitpa6*6 << endl;
		cout << "chi2 : " << m4fitchi2 << endl;

		m4fitpara->Fill();


cout<<endl;
cout<<"///////m5 fit///////"<<endl;
cout<<endl;
		//m5fit
		//m5fit -> SetParameters(gfitpa1, peak[i], gfitpa3, gfitpa1, peak[i+1], gfitpa3, p1fitpa1, p1fitpa2);
		//m5fit -> SetParameters(gfitpa1, peak[i], gfitpa3, gfitpa1, peak[i-1], gfitpa3, p1fitpa1, p1fitpa2);
		//m5fit -> SetParameters(gpa1[i-1], peak[i-1], gpa2[i-1], gpa1[i], peak[i], gpa2[i], p1fitpa1, p1fitpa2);
		m5fit -> SetParameters(gfitpa1, peak[i-1], gfitpa3, gfitpa1, peak[i], gfitpa3, p1fitpa1, p1fitpa2);
//		m5fit -> SetParameters(gpa1[i], peak[i], gpa2[i], gpa1[i-1], peak[i-1], gpa2[i-1], p1fitpa1, p1fitpa2);


		his_temp -> Fit("m5fit", "R+"); 
		m5fitpa1 = m5fit->GetParameter(0);
		m5fitpa2 = m5fit->GetParameter(1);
		m5fitpa3 = m5fit->GetParameter(2);
		m5fitpa4 = m5fit->GetParameter(3);
		m5fitpa5 = m5fit->GetParameter(4);
		m5fitpa6 = m5fit->GetParameter(5);
		m5fitpa7 = m5fit->GetParameter(6);
		m5fitpa8 = m5fit->GetParameter(7);
		m5fiter1 = m5fit->GetParError(0);	
		m5fiter2 = m5fit->GetParError(1);	
		m5fiter3 = m5fit->GetParError(2);	
		m5fiter4 = m5fit->GetParError(3);	
		m5fiter5 = m5fit->GetParError(4);	
		m5fiter6 = m5fit->GetParError(5);	
		m5fiter7 = m5fit->GetParError(6);	
		m5fiter8 = m5fit->GetParError(7);	
		m5fitchi2 = (m5fit->GetChisquare())/(m5fit->GetNDF());
/*
		m5pcnt1 = m5fitpa1*mul;
		m5pcnter1 = sqrt(m5fitpa1*mul);
		m5cpd1 = m5pcnt1/time_day;
		m5dcpd1 = m5pcnter1/time_day;

		m5pcnt2 = m5fitpa4*mul;
		m5pcnter1 = sqrt(m5fitpa1*mul);
		m5cpd1 = m5pcnt1/time_day;
		m5dcpd1 = m5pcnter1/time_day;

		m5pcnt2 = m5fitpa4*mul;
		m5pcnter2 = sqrt(m5fitpa4*mul);
		m5cpd2 = m5pcnt2/time_day;
		m5dcpd2 = m5pcnter2/time_day;
*/
		cout << "[m5fit]" << endl;
		cout << "peak1 : " << m5fitpa2 << endl;
		cout << "m5fit tot. counts (peak1) :" << m5pcnt1 << " +/- " << m5pcnter1 << endl;
//		cout << "m5fit cpd (peak1) :" << m5cpd1 << " +/- " << m5dcpd1 << endl;
		cout << "6 sigma : " << m5fitpa3*6 << endl;
		cout << "peak2 : " << m5fitpa5 << endl;
		cout << "m5fit tot. counts (peak2) :" << m5pcnt2 << " +/- " << m5pcnter2 << endl;
//		cout << "m5fit cpd (peak2) :" << m5cpd2 << " +/- " << m5dcpd2 << endl;
		cout << "6 sigma : " << m5fitpa6*6 << endl;
//		cout << "chi2 : " << m5fitchi2 << endl;

		m5fitpara->Fill();

		}

	}
/*
	if (runnum == 136 || runnum == 137 || runnum == 138 || runnum == 13678 || runnum == 1368 || runnum ==21490) { his_temp->SetTitle("Mo-100 powder data"); }
	//	if (runnum == 151 || runnum == 155 || runnum == 1515 || runnum == 84) { his_temp->SetTitle("BKG"); }
	else { his_temp->SetTitle("BKG"); }
*/

	his_temp->SetTitle("ThO2 powder");
	his_temp->SetXTitle("Energy[keV]");
	his_temp->SetYTitle("Counts/kev");

	char his_name[256];
	sprintf(his_name, "double_peak_his");
	his_temp->SetName(his_name);

	c1->Modified();
	c1->Update();
	c1->SaveAs(res2file);

	TFile r3f (res3file,"RECREATE");
	r3f.cd();
	gfitpara->Write();
	m3fitpara->Write();
	m4fitpara->Write();
	m5fitpara->Write();
	r3f.Close();
	gfitpara->Reset();
	m3fitpara->Reset();
	m4fitpara->Reset();
	m5fitpara->Reset();
}
