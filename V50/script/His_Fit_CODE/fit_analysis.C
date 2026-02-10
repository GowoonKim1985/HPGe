#include <iostream>
#include "TMath.h"
#include "TMinuit.h"

void fit_peak() {

	int runnum;
	char runnumber3[256];
	char runnumber6[256];
	
        cout<< "run number : ";
        scanf("%i",&runnum);
        sprintf(runnumber6,"%06d",runnum);
        sprintf(runnumber3,"%03d",runnum);

        int bin;
        cout<< "x-bin number (default:4000(1ADC)) : ";
        scanf("%i",&bin);

        char binnumber[256];
        sprintf(binnumber,"%d",bin);

	char iso[256];
//	sprintf(iso, "./peak/peak_R%s.dat",runnumber3);
	sprintf(iso, "./peak/peak_BG.dat");

	ifstream peakdata1(iso);
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
	ifstream peakdata2(iso);

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
	char res2file[16][256]; //fit his
	char res3file[256]; //fit result

	//int runnum;
	//char runnumber[256];
	//cout << "run number : ";
	//scanf("%i", &runnum);
	//sprintf(runnumber, "%06d", runnum);

	//char evt_case[256];
	int evt_case;
	cout << "event case (0(=tot) / 1(1mul) / 2(2mul) : ";
	scanf("%i", &evt_case);

//	int bin;
//	cout << "binning (4000, 1bin=1keV / 8000, 1bin=0.5keV) : ";
//	scanf("%d", &bin);
	
	//	int det_num;
	//	cout << "the number of detector : ";
	//	scanf("%i", &det_num);


	if(evt_case==0){
		//sprintf(hisfile, "/home/kkw/DAQ/ANA500/result/RUN%d/Bin%i/his_%s.root", runnum, bin, runnumber);
		//sprintf(res2file, "/home/kkw/DAQ/ANA500/result/RUN%d/Bin%i/Tot_spectrum.C", runnum, bin);
		//sprintf(res3file, "/home/kkw/DAQ/ANA500/result/RUN%d/Bin%i/res_%s.root", runnum, bin, runnumber);
		sprintf(hisfile, "/data/HPGe/USERS/kkw/DAQ/ANA500/result/RUN%d/Bin%i/his_%s.root", runnum, bin, runnumber6);

		for(int res1 = 0; res1 <16; res1++){	
		sprintf(res2file[res1], "/data/HPGe/USERS/kkw/DAQ/ANA500/result/RUN%d/Bin%i/fit_his%i.C", runnum, bin, res1);
		}
		sprintf(res3file, "/data/HPGe/USERS/kkw/DAQ/ANA500/result/RUN%d/Bin%i/res_%s.root", runnum, bin, runnumber6);
	}
	else{
		sprintf(hisfile, "/home/kkw/DAQ/ANA500/result/RUN%d/his_%s_M%i.root", runnum, runnumber6, evt_case);
		for(int res2 = 0; res2 <16; res2++){	
		sprintf(res2file[res2], "/data/HPGe/USERS/kkw/DAQ/ANA500/result/RUN%d/Bin%i/fit_his%i.C", runnum, bin, res2);
		}
//		sprintf(res2file, "/home/kkw/DAQ/ANA500/result/RUN%d/M%i_spectrum.C", runnum, evt_case);
		sprintf(res3file, "/home/kkw/DAQ/ANA500/result/RUN%d/res_%s_M%i.root", runnum, runnumber6, evt_case);
	}


	TF1 * gfit = new TF1("gfit", "gaus");
	TF1 * agfit = new TF1("agfit", "[0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2])");
	agfit->SetParNames("Area", "Mean", "Sigma"); //area = sqrt(2pi)*constant*sigma

	TF1 * efit = new TF1("efit", "expo");
	//	TF1 * e2fit = new TF1("e2fit", "expo");
	TF1 * p0fit = new TF1("p0fit", "pol0");
	TF1 * p1fit = new TF1("p1fit", "pol1");

	TF1 * m1fit = new TF1("m1fit", "([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+[3]");
	m1fit->SetParNames("Area", "Mean", "Sigma", "Flat");

	TF1 * m2fit = new TF1("m2fit", "([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+exp([3]+[4]*x)");
	m2fit->SetParNames("Area", "Mean", "Sigma", "expo1", "expo2");

	TF1 * m3fit = new TF1("m3fit", "([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+[3]+[4]*x");
	m3fit->SetParNames("Area", "Mean", "Sigma", "Intercept", "Slope");

	agfit->SetLineColor(9);
	gfit->SetLineColor(3);
	efit->SetLineColor(5);
	//	e2fit->SetLineColor(6);
	p1fit->SetLineColor(2);
	m1fit->SetLineColor(kBlue);
	m2fit->SetLineColor(kRed);
	m3fit->SetLineColor(kGreen);

	int maxb, maxx, maxh, channel;
	double bg1, bg2, bg;
	double cnt1, cnt2, cnt3, cnt1_2, cnt3_2;
	double pa1[line_max1-1][3];//gfit 
	double cal[line_max1-1][2]; //integral fit func - m1fit,m2fit
	double epa_fix[2]; //expo fit 
	int xbin[line_max1-1];
	double energy;

	double agfitpa1, agfitpa2, agfitpa3, agfiter1, agfiter2, agfiter3;

	double gfitpa1, gfitpa2, gfitpa3, gfiter1, gfiter2, gfiter3;

	double p1fitpa1, p1fitpa2, p1fiter1, p1fiter2;

	double m1fitpa1, m1fitpa2, m1fitpa3, m1fitpa4;
	double m1fiter1, m1fiter2, m1fiter3, m1fiter4;
	double m1fitmean[16][line_max1-1], m1fitsig[16][line_max1-1], m1resolution[16][line_max1-1];
	double m1fitmeaner[16][line_max1-1], m1fitsiger[16][line_max1-1], m1resolutioner[16][line_max1-1];
	double m1pcnt, m1pcnter, m1pcnter1, m1pcnter2;  
	double m1res, m1reser;

	double m2fitpa1, m2fitpa2, m2fitpa3, m2fitpa4, m2fitpa5;
	double m2fiter1, m2fiter2, m2fiter3, m2fiter4, m2fiter5;
	double m2fitmean[16][line_max1-1], m2fitsig[16][line_max1-1], m2resolution[16][line_max1-1];
	double m2fitmeaner[16][line_max1-1], m2fitsiger[16][line_max1-1], m2resolutioner[16][line_max1-1];
	double m2pcnt, m2pcnter, m2pcnter1, m2pcnter2; 
	double m2res, m2reser;

	double m3fitpa1, m3fitpa2, m3fitpa3, m3fitpa4, m3fitpa5;
	double m3fiter1, m3fiter2, m3fiter3, m3fiter4, m3fiter5;
	double m3fitmean[16][line_max1-1], m3fitsig[16][line_max1-1], m3resolution[16][line_max1-1];
	double m3fitmeaner[16][line_max1-1], m3fitsiger[16][line_max1-1], m3resolutioner[16][line_max1-1];
	double m3pcnt, m3pcnter, m3pcnter1, m3pcnter2; 
	double m3res, m3reser;

	double gfitchi2, agfitchi2, p1fitchi2, m1fitchi2, m2fitchi2, m3fitchi2;

	double temp_m1fitmean[line_max1-1], temp_m1resol[line_max1-1], temp_m1fitmeaner[line_max1-1], temp_m1resoler[line_max1-1];
	double temp_m2fitmean[line_max1-1], temp_m2resol[line_max1-1], temp_m2fitmeaner[line_max1-1], temp_m2resoler[line_max1-1];
	double temp_m3fitmean[line_max1-1], temp_m3resol[line_max1-1], temp_m3fitmeaner[line_max1-1], temp_m3resoler[line_max1-1];

	int best_fit;
	double best_mean[16][line_max1-1], best_meaner[16][line_max1-1], best_res[16][line_max1-1], best_reser[16][line_max1-1], best_sig[16][line_max1-1];
	double final_mean[line_max1-1], final_meaner[line_max1-1], final_res[line_max1-1], final_reser[line_max1-1];
	double fin_res, fin_reser;
	double respa1[16], respa2[16], respa3[16];
	TString m1stat;

	TTree *gfitpara = new TTree("gfitpara", "gaus fitting result");
	gfitpara -> Branch("run", &runnum, "run/I");
	gfitpara -> Branch("channel", &channel, "channel/I");
	gfitpara -> Branch("peak", &energy, "energy/D");
	gfitpara -> Branch("const", &gfitpa1, "gfitpa1/D");
	gfitpara -> Branch("mean", &gfitpa2, "gfitpa2/D");
	gfitpara -> Branch("chi2", &gfitchi2, "gfitchi2/D");
	gfitpara -> Branch("const_er", &gfiter1, "gfiter1/D");
	gfitpara -> Branch("mean_er", &gfiter2, "gfiter2/D");
	gfitpara -> Branch("sigma_er", &gfiter3, "gfiter3/D");

	TTree *p1fitpara = new TTree("p1fitpara", "pol1 fitting result");
	p1fitpara -> Branch("run", &runnum, "run/I");
	p1fitpara -> Branch("channel", &channel, "channel/I");
	p1fitpara -> Branch("peak", &energy, "energy/D");
	p1fitpara -> Branch("intercept", &p1fitpa1, "p1fitpa1/D");
	p1fitpara -> Branch("slope", &p1fitpa2, "p1fitpa2/D");
	p1fitpara -> Branch("chi2", &p1fitchi2, "p1fitchi2/D");
	p1fitpara -> Branch("intercept_er", &p1fiter1, "p1fiter1/D");
	p1fitpara -> Branch("slope_er", &p1fiter2, "p1fiter2/D");

	TTree *agfitpara = new TTree("agfitpara", "area-gaus fitting result");
	agfitpara -> Branch("run", &runnum, "run/I");
	agfitpara -> Branch("channel", &channel, "channel/I");
	agfitpara -> Branch("peak", &energy, "energy/D");
	agfitpara -> Branch("area", &agfitpa1, "gfitpa1/D");
	agfitpara -> Branch("mean", &agfitpa2, "gfitpa2/D");
	agfitpara -> Branch("chi2", &agfitchi2, "gfitchi2/D");
	agfitpara -> Branch("area_er", &agfiter1, "gfiter1/D");
	agfitpara -> Branch("mean_er", &agfiter2, "gfiter2/D");
	agfitpara -> Branch("sigma_er", &agfiter3, "gfiter3/D");

	TTree *m1fitpara = new TTree("m1fitpara", "gaus+pol0 fitting result");
	m1fitpara -> Branch("run", &runnum, "run/I");
	m1fitpara -> Branch("channel", &channel, "channel/I");
	m1fitpara -> Branch("peak", &energy, "energy/D");
	m1fitpara -> Branch("area", &m1fitpa1, "m1fitpa1/D");
	m1fitpara -> Branch("mean", &m1fitpa2, "m1fitpa2/D");
	m1fitpara -> Branch("sigma", &m1fitpa3, "m1fitpa3/D");
	m1fitpara -> Branch("chi2", &m1fitchi2, "m1fitchi2/D");
	m1fitpara -> Branch("pol0", &m1fitpa4, "m1fitpa4/D");
	m1fitpara -> Branch("area_er", &m1fiter1, "m1fiter1/D");
	m1fitpara -> Branch("mean_er", &m1fiter2, "m1fiter2/D");
	m1fitpara -> Branch("sigma_er", &m1fiter3, "m1fiter3/D");
	m1fitpara -> Branch("pol0_er", &m1fiter4, "m1fitper4/D");
	m1fitpara -> Branch("pcount", &m1pcnt, "m1pcnt/D");
	m1fitpara -> Branch("pcount_er", &m1pcnter, "m1pcnter/D");
	m1fitpara -> Branch("resolution", &m1res, "m1resolution/D");
	m1fitpara -> Branch("resolution_er", &m1reser, "m1resolutioner/D");

	TTree *m2fitpara = new TTree("m2fitpara", "gaus+expo fitting result");
	m2fitpara -> Branch("run", &runnum, "run/I");
	m2fitpara -> Branch("channel", &channel, "channel/I");
	m2fitpara -> Branch("peak", &energy, "energy/D");
	m2fitpara -> Branch("area", &m2fitpa1, "m2fitpa1/D");
	m2fitpara -> Branch("mean", &m2fitpa2, "m2fitpa2/D");
	m2fitpara -> Branch("sigma", &m2fitpa3, "m2fitpa3/D");
	m2fitpara -> Branch("chi2", &m2fitchi2, "m2fitchi2/D");
	m2fitpara -> Branch("exp1", &m2fitpa4, "m2fitpa4/D");
	m2fitpara -> Branch("exp2", &m2fitpa5, "m2fitpa5/D");
	m2fitpara -> Branch("area_er", &m2fiter1, "m2fiter1/D");
	m2fitpara -> Branch("mean_er", &m2fiter2, "m2fiter2/D");
	m2fitpara -> Branch("sigma_er", &m2fiter3, "m2fiter3/D");
	m2fitpara -> Branch("exp1_er", &m2fiter4, "m2fitper4/D");
	m2fitpara -> Branch("exp2_er", &m2fiter5, "m2fitper5/D");
	m2fitpara -> Branch("pcount", &m2pcnt, "m2pcnt/D");
	m2fitpara -> Branch("pcount_er", &m2pcnter, "m2pcnter/D");
	m2fitpara -> Branch("resolution", &m2res, "m2resolution/D");
	m2fitpara -> Branch("resolution_er", &m2reser, "m2resolutioner/D");

	TTree *m3fitpara = new TTree("m3fitpara", "gaus+pol1 fitting result");
	m3fitpara -> Branch("run", &runnum, "run/I");
	m3fitpara -> Branch("channel", &channel, "channel/I");
	m3fitpara -> Branch("peak", &energy, "energy/D");
	m3fitpara -> Branch("area", &m3fitpa1, "m3fitpa1/D");
	m3fitpara -> Branch("mean", &m3fitpa2, "m3fitpa2/D");
	m3fitpara -> Branch("sigma", &m3fitpa3, "m3fitpa3/D");
	m3fitpara -> Branch("chi2", &m3fitchi2, "m3fitchi2/D");
	m3fitpara -> Branch("pol0", &m3fitpa4, "m3fitpa4/D");
	m3fitpara -> Branch("pol1", &m3fitpa5, "m3fitpa5/D");
	m3fitpara -> Branch("area_er", &m3fiter1, "m3fiter1/D");
	m3fitpara -> Branch("mean_er", &m3fiter2, "m3fiter2/D");
	m3fitpara -> Branch("sigma_er", &m3fiter3, "m3fiter3/D");
	m3fitpara -> Branch("pol0_er", &m3fiter4, "m3fitper4/D");
	m3fitpara -> Branch("pol1_er", &m3fiter5, "m3fitper5/D");
	m3fitpara -> Branch("pcount", &m3pcnt, "m3pcnt/D");
	m3fitpara -> Branch("pcount_er", &m3pcnter, "m3pcnter/D");
	m3fitpara -> Branch("resolution", &m3res, "m3resolution/D");
	m3fitpara -> Branch("resolution_er", &m3reser, "m3resolutioner/D");

	TTree *finalpara = new TTree("finalpara", "resolution checked");
	finalpara -> Branch("fitnum", &best_fit, "best_fit/I");
	finalpara -> Branch("run", &runnum, "run/I");
	finalpara -> Branch("channel", &channel, "channel/I");
	finalpara -> Branch("peak", &energy, "energy/D");
	finalpara -> Branch("resolution", &fin_res, "resolution/D");
	finalpara -> Branch("resolution_er", &fin_reser, "resolutioner/D");

	int mul = bin/4000;
	int agint;
	TFile * hf = new TFile(hisfile);
	TH1D * his_temp = new TH1D("his_temp", "", bin, 0, 4000);
	char data_his[256];

//	TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800);
	//	TCanvas *gc1 = new TCanvas("gc1", "gc1", 1200, 800);
	//	TCanvas *gc2 = new TCanvas("gc2", "gc2", 1200, 800);
	//	TCanvas *gc3 = new TCanvas("gc3", "gc3", 1200, 800);
	TCanvas *final_gc = new TCanvas("final_gc", "final_gc", 1200, 800);

	char canname[16][256];
	TCanvas *can[16];
		
	for(int ci=0;ci<16;ci++){
	sprintf(canname[ci],"can%i",ci);
	can[ci] = new TCanvas(canname[ci], canname[ci], 800, 600);
	}

//	c1->Divide(4, 4);
	//	gc1->Divide(4, 4);
	//	gc2->Divide(4, 4);
	//	gc3->Divide(4, 4);
	final_gc->Divide(4, 4);

	TGraphErrors *gr1[16], *gr2[16], *gr3[16], *final_gr[16];

	char result_hisname[256];
	TF1 * res_fit = new TF1("res_fit", "[0] + [1]*1./x + [2]/sqrt(x)", 250, 2500);

	for(int j=1; j<=16; j++) {
		if((j>=1&&j<=7)||(j>=9&&j<=15)){
			sprintf(data_his, "his%i", j);
			his_temp = (TH1D*)hf->Get(data_his);
			his_temp -> SetLineColor(1);

			channel = j;
			cout<<""<<endl;
			cout<<""<<endl;
			cout<<"[Ch : "<<j<<"]"<<endl;
			cout<<""<<endl;
			can[j]->cd();
	//		c1 -> cd(j);
			/*		efit -> SetRange(500, 2200);
					his_temp -> Fit("efit", "RQ0+");
					efit -> GetParameters(&epa[0]);
					efit -> SetRange(2500, 2800);
					his_temp -> Fit("efit", "RQ0+");
					efit -> GetParameters(&epa[2]);
					efit -> SetRange(3000, 4000);
					his_temp -> Fit("efit", "RQ0+");
					efit -> GetParameters(&epa[4]);
					*/
			for(int i=0; i<(line_max1-1); i++) {
				energy = peak[i];
				xbin[i] = his_temp->GetXaxis()->FindBin(peak[i]);
				printf("\n");
				printf("\n=========================\n");
				printf("%.2f keV bin : %i	\n", peak[i], xbin[i]);

				his_temp -> GetXaxis() -> SetRange(xbin[i]-5*mul, xbin[i]+5*mul);
				maxb = his_temp->GetMaximumBin();
				maxx = maxb/mul;
				cout << "maxb(bin) : " << maxb << endl;
				cout << "maxx(kev) : " << maxx << endl;
				maxh = his_temp->GetMaximum();

				his_temp->GetXaxis()->SetRange(0, 4000*mul);

				int range_val = 10;

				gfit->SetRange((maxx-range_val),(maxx+range_val));
				agfit->SetRange((maxx-range_val),(maxx+range_val));
				p1fit->SetRange((maxx-range_val),(maxx+range_val));
				m1fit->SetRange((maxx-range_val),(maxx+range_val));
				m2fit->SetRange((maxx-range_val),(maxx+range_val));
				m3fit->SetRange((maxx-range_val),(maxx+range_val));
				efit->SetRange((maxx-30),(maxx+30));

				//efit to get expo parameter
				epa_fix[0] = efit->GetParameter(0);
				epa_fix[1] = efit->GetParameter(1);

				//gaus fit
				his_temp->Fit("gfit","RQ0+");
				gfitpa1 = gfit->GetParameter(0);
				gfitpa2 = gfit->GetParameter(1);
				gfitpa3 = gfit->GetParameter(2);
				gfiter1 = gfit->GetParError(0);
				gfiter2 = gfit->GetParError(1);
				gfiter3 = gfit->GetParError(2);

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
				//agint = his_temp->Integral(maxx-int(3*mul*gfitpa3)-1,maxx-int(3*mul*gfitpa3)+1);
			//	cout<<"peak integral : "<<agint<<endl;
				//agfit->SetParameters(2.5*gfitpa2*gfitpa3, gfitpa2, gfitpa3);
				agfit->SetParameters(2.5*gfitpa2*gfitpa3, maxx, gfitpa3);
			//	agfit->SetParameters(agint, maxx, gfitpa3);
				his_temp->Fit("agfit","RQ0+");
				agfitpa1 = agfit->GetParameter(0);
				agfitpa2 = agfit->GetParameter(1);
				agfitpa3 = agfit->GetParameter(2);
				agfiter1 = agfit->GetParError(0);
				agfiter2 = agfit->GetParError(1);
				agfiter3 = agfit->GetParError(2);
				agfitchi2 = (agfit->GetChisquare())/(agfit->GetNDF());

				cout<<""<<endl;
				cout<<"[m1 fit]"<<endl;
				cout<<""<<endl;
	

				//m1fit->SetParameters(agfitpa1, peak[i], agfitpa3, bg);
        		        m1fit->ReleaseParameter(0);
        		        m1fit->ReleaseParameter(1);
        		        m1fit->ReleaseParameter(2);
        		        m1fit->ReleaseParameter(3);
				m1fit->SetParameters(agfitpa1, maxx, agfitpa3, bg);
				his_temp->Fit("m1fit", "RQ0+");
				m1fitpa2 = m1fit->GetParameter(1);
				if (abs(m1fitpa2-maxx)>5){m1fit->FixParameter(1,maxx);}
				his_temp->Fit("m1fit", "R+");

				m1fitpa1 = m1fit->GetParameter(0);
				m1fitpa2 = m1fit->GetParameter(1);
				m1fitpa3 = m1fit->GetParameter(2);
				m1fitpa4 = m1fit->GetParameter(3);
				m1fiter1 = m1fit->GetParError(0);
				m1fiter2 = m1fit->GetParError(1);
				m1fiter3 = m1fit->GetParError(2);
				m1fiter4 = m1fit->GetParError(3);
				m1fitchi2 = (m1fit->GetChisquare())/(m1fit->GetNDF());

				printf("BG : %.2f	\n", bg);
					if(m1fitpa3<0) {
						cout << "!!!!m1fit sigma is minus -> will be get ABS!!!" << endl;
					}
				m1pcnt = m1fitpa1;
				m1pcnter = m1fiter1;
				cout << "m1fit int.error : " << m1pcnter << endl;

				//m2fit
				/*
				   if (peak[i]<2200) {
				   m2fit->SetParameters(gfitpa1, peak[i], gfitpa3, epa[0], epa[1]);
				   }
				   else if (peak[i]>2500 && peak[i]<2700) {
				   m2fit->SetParameters(gfitpa1, peak[i], gfitpa3, epa[2], epa[3]);
				   }
				   else if (peak[i]>3000) {
				   m2fit->SetParameters(gfitpa1, peak[i], gfitpa3, epa[4], epa[5]);
				   }
				   */

				cout<<""<<endl;
				cout<<"[m2 fit]"<<endl;
				cout<<""<<endl;
			
				//m2fit->SetParameters(agfitpa1, peak[i], agfitpa3, epa_fix[0], epa_fix[1]);
        		        m2fit->ReleaseParameter(0);
        		        m2fit->ReleaseParameter(1);
        		        m2fit->ReleaseParameter(2);
        		        m2fit->ReleaseParameter(3);
        		        m2fit->ReleaseParameter(4);
				m2fit->SetParameters(agfitpa1, maxx, agfitpa3, epa_fix[0], epa_fix[1]);
				his_temp->Fit("m2fit", "RQ0+");
				m2fitpa2 = m2fit->GetParameter(1);
				if (abs(m2fitpa2-maxx)>5){m2fit->FixParameter(1,maxx);}
				his_temp->Fit("m2fit", "R+");
				m2fitpa1 = m2fit->GetParameter(0);
				m2fitpa2 = m2fit->GetParameter(1);
				m2fitpa3 = m2fit->GetParameter(2);
				m2fitpa4 = m2fit->GetParameter(3);
				m2fitpa5 = m2fit->GetParameter(4);
				m2fiter1 = m2fit->GetParError(0);
				m2fiter2 = m2fit->GetParError(1);
				m2fiter3 = m2fit->GetParError(2);
				m2fiter4 = m2fit->GetParError(3);
				m2fiter5 = m2fit->GetParError(4);
				m2fitchi2 = (m2fit->GetChisquare())/(m2fit->GetNDF());

					if(m2fitpa3 < 0) { cout << "!!!!m2fit sigma is minus -> will be get ABS!!!" << endl; }
				m2pcnt = m2fitpa1;
				m2pcnter = m2fiter1;


				//m3fit
				cout<<""<<endl;
				cout<<"[m3 fit]"<<endl;
				cout<<""<<endl;

        		        m3fit->ReleaseParameter(0);
        		        m3fit->ReleaseParameter(1);
        		        m3fit->ReleaseParameter(2);
        		        m3fit->ReleaseParameter(3);
        		        m3fit->ReleaseParameter(4);

				m3fit->SetParameters(agfitpa1, maxx, agfitpa3, p1fitpa1, p1fitpa2);
				his_temp->Fit("m3fit", "RQ0+");
				cout<<"test"<<endl;
				m3fitpa2 = m3fit->GetParameter(1);
					cout<<"mean : "<<m3fitpa2<<endl;
				if (abs(m3fitpa2-maxx)>5){
					cout<<"max height : "<<maxx<<endl;
					m3fit->FixParameter(1,maxx);}
				his_temp->Fit("m3fit", "R+");
				m3fitpa1 = m3fit->GetParameter(0);
				m3fitpa2 = m3fit->GetParameter(1);
				m3fitpa3 = m3fit->GetParameter(2);
				m3fitpa4 = m3fit->GetParameter(3);
				m3fitpa5 = m3fit->GetParameter(4);
				m3fiter1 = m3fit->GetParError(0);
				m3fiter2 = m3fit->GetParError(1);
				m3fiter3 = m3fit->GetParError(2);
				m3fiter4 = m3fit->GetParError(3);
				m3fiter5 = m3fit->GetParError(4);
				m3fitchi2 = (m3fit->GetChisquare())/(m3fit->GetNDF());

					if(m3fitpa3 < 0) { cout << "!!!!m3fit sigma is minus -> will be get ABS!!!" << endl; }
				m3pcnt = m3fitpa1;
				m3pcnter = m3fiter1;

				cout << "[m1fit]" << endl;
				cout <<  "m1fit tot. counts:" << m1pcnt << " +/- " << m1pcnter << endl;
				cout << "6 sigma : " << m1fitpa3*6 << endl;
				cout << "m1fit chi2 : " << m1fitchi2 << endl;
				cout << "" << endl;
				cout << "[m2fit]" << endl;
				cout <<  "m2fit tot. counts:" << m2pcnt << " +/- " << m2pcnter << endl;
				cout << "6 sigma : " << m2fitpa3*6 << endl;
				cout << "m2fit chi2 : " << m2fitchi2 << endl;
				cout << "" << endl;
				cout << "[m3fit]" << endl;
				cout <<  "m3fit tot. counts:" << m3pcnt << " +/- " << m3pcnter << endl;
				cout << "6 sigma : " << m3fitpa3*6 << endl;
				cout << "m3fit chi2 : " << m3fitchi2 << endl;
				cout << "" << endl;

				m1fitmean[j][i] = m1fitpa2;
				m1fitsig[j][i] = m1fitpa3; 
				m1fitmeaner[j][i] = m1fiter2;
				m1fitsiger[j][i] = m1fiter3;

				m2fitmean[j][i] = m2fitpa2;
				m2fitsig[j][i] = m2fitpa3; 
				m2fitmeaner[j][i] = m2fiter2;
				m2fitsiger[j][i] = m2fiter3;

				m3fitmean[j][i] = m3fitpa2;
				m3fitsig[j][i] = m3fitpa3; 
				m3fitmeaner[j][i] = m3fiter2;
				m3fitsiger[j][i] = m3fiter3;

				m1resolution[j][i] = m1fitsig[j][i]/m1fitmean[j][i];
				m2resolution[j][i] = m2fitsig[j][i]/m2fitmean[j][i];
				m3resolution[j][i] = m3fitsig[j][i]/m3fitmean[j][i];

				m1resolutioner[j][i] = m1resolution[j][i] * sqrt(pow(m1fitmeaner[j][i]/m1fitmean[j][i],2) + pow(m1fitsiger[j][i]/m1fitsig[j][i],2));
				m2resolutioner[j][i] = m2resolution[j][i] * sqrt(pow(m2fitmeaner[j][i]/m2fitmean[j][i],2) + pow(m2fitsiger[j][i]/m2fitsig[j][i],2));
				m3resolutioner[j][i] = m3resolution[j][i] * sqrt(pow(m3fitmeaner[j][i]/m3fitmean[j][i],2) + pow(m3fitsiger[j][i]/m3fitsig[j][i],2));

				m1res = m1resolution[j][i];
				m2res = m2resolution[j][i];
				m3res = m3resolution[j][i];

				m1reser = m1resolutioner[j][i];
				m2reser = m2resolutioner[j][i];
				m3reser = m3resolutioner[j][i];

				// To check the best fitting result

					if( abs(1-m1fitchi2) <= abs(1-m2fitchi2) ) { // m1 is better than m2
						if ( abs(1-m1fitchi2) <= abs(1-m3fitchi2) ){best_fit = 1;} // m1 is the best fit
						if ( abs(1-m1fitchi2) > abs(1-m3fitchi2) ) {best_fit = 3;} // m3 is the best fit
					}

					if( abs(1-m1fitchi2) > abs(1-m2fitchi2) ) { // m2 is better than m1
						if( abs(1-m2fitchi2) <= abs(1-m3fitchi2) ) {best_fit = 2;} // m2 is the best fit
						if( abs(1-m2fitchi2) > abs(1-m3fitchi2) ) {best_fit = 3;} //m3 is the best fit
					}
					
				if(j==7){best_fit=3;}//temp	
				
					if( best_fit == 1 ) {
						best_mean[j][i] = m1fitmean[j][i];
						best_meaner[j][i] = m1fitmeaner[j][i];
						best_res[j][i] = m1resolution[j][i];
						best_reser[j][i] = m1resolutioner[j][i];
						best_sig[j][i] = m1fitsig[j][i];
						cout<<"best fit : m1 fit"<<endl;
					}
	
					if( best_fit == 2 ) {
						best_mean[j][i] = m2fitmean[j][i];
						best_meaner[j][i] = m2fitmeaner[j][i];
						best_res[j][i] = m2resolution[j][i];
						best_reser[j][i] = m2resolutioner[j][i];
						best_sig[j][i] = m2fitsig[j][i];
						cout<<"best fit : m2 fit"<<endl;
					}

					if( best_fit == 3 ) {
						best_mean[j][i] = m3fitmean[j][i];
						best_meaner[j][i] = m3fitmeaner[j][i];
						best_res[j][i] = m3resolution[j][i];
						best_reser[j][i] = m3resolutioner[j][i];
						best_sig[j][i] = m3fitsig[j][i];
						cout<<"best fit : m3 fit"<<endl;
					}
cout<<"best resultion ="<<best_res[j][i]<<endl;
cout<<"best resultion error ="<<best_reser[j][i]<<endl;
					temp_m1fitmean[i] = m1fitmean[j][i];
					temp_m1resol[i] = m1resolution[j][i];
					temp_m1fitmeaner[i] = m1fitmeaner[j][i];
					temp_m1resoler[i] = m1resolutioner[j][i];

					temp_m2fitmean[i] = m2fitmean[j][i];
					temp_m2resol[i] = m2resolution[j][i];
					temp_m2fitmeaner[i] = m2fitmeaner[j][i];
					temp_m2resoler[i] = m2resolutioner[j][i];

					temp_m3fitmean[i] = m3fitmean[j][i];
					temp_m3resol[i] = m3resolution[j][i];
					temp_m3fitmeaner[i] = m3fitmeaner[j][i];
					temp_m3resoler[i] = m3resolutioner[j][i];

					fin_res = best_res[j][i];
					fin_reser = best_reser[j][i];

					final_mean[i] = best_mean[j][i];
					final_meaner[i] = best_meaner[j][i];
					final_res[i] = best_res[j][i];
					final_reser[i] = best_reser[j][i];

				//			gr1[j] = new TGraphErrors(line_max1-1, temp_m1fitmean, temp_m1resol, temp_m1fitmeaner, temp_m1resoler);
				//			gr2[j] = new TGraphErrors(line_max1-1, temp_m2fitmean, temp_m2resol, temp_m2fitmeaner, temp_m2resoler);
				//			gr3[j] = new TGraphErrors(line_max1-1, temp_m3fitmean, temp_m3resol, temp_m3fitmeaner, temp_m3resoler);
					final_gr[j] = new TGraphErrors(line_max1-1, final_mean, final_res, final_meaner, final_reser);

					m1fitpara->Fill();
					m2fitpara->Fill();
					m3fitpara->Fill();
					finalpara->Fill();
				}
			final_gc->cd(j);
			//final_gc[j]->cd();
			sprintf(result_hisname, "best fitting result ch%i", j);
			final_gr[j]->SetTitle(result_hisname);
			final_gr[j]->Draw("ALP");
			cout<<"[Ch : "<<j<<"]"<<endl;
			cout<<""<<endl;
			final_gr[j]->Fit("res_fit", "RE");
			respa1[j] = res_fit->GetParameter(0);
			respa2[j] = res_fit->GetParameter(1);
			respa3[j] = res_fit->GetParameter(2);
			cout<<""<<endl;
		}	

	}

        for(int rc=1;rc<16;rc++){
        printf("%.5f * %.5f * 1332.49\n",best_mean[rc][2], best_sig[rc][2]);
	}

			cout<<"[0] + [1]*1./x + [2]/sqrt(x)"<<endl;
		for(int l=1; l<=16 ;l++){
			if(l!=8){
			printf("ch %i * %.9f * %.9f * %.9f	\n",l,respa1[l],respa2[l],respa3[l]);
			}
		}
/*

		for(int l=1; l<=16 ;l++){
			cout<<""<<endl;
			cout<<"chn : "<<l<<endl;
			for(int m=0; m<5 ; m++){
			cout<<""<<endl;
			cout<< "mean "<<m<<" peak : "<< best_mean[l][m]<<endl; 
			cout<< "mean_er "<<m<<" peak : "<< best_meaner[l][m]<<endl; 
			cout<< "res "<<m<<" peak : "<< best_res[l][m]<<endl; 
			cout<< "res_er "<<m<<" peak : "<< best_reser[l][m]<<endl; 
			cout<<""<<endl;
			}
		}
*/

//	for(int j=1; j<16; j++){
//		if((j>=1&&j<=7)||(j>=9&&j<=15)){
/*			final_gc->cd(j);
			sprintf(result_hisname, "best fitting result ch%i", j);
			final_gr[j]->SetTitle(result_hisname);
			final_gr[j]->Draw("ALP");
			cout<<"[Ch : "<<j<<"]"<<endl;
			cout<<""<<endl;
			final_gr[j]->Fit("res_fit", "RE");
			cout<<""<<endl;
*/
//		}
//	}




	//	char result_hisname[256];
	//	TF1 * res_fit = new TF1("res_fit", "[0] + [1]*1./x + [2]/sqrt(x)", 300, 2500);
	/*
	   for(int j=1; j<16; j++) {
	   final_gc->cd(j);
	   sprintf(result_hisname, "best fitting result ch%i", j);
	   final_gr[j]->SetTitle(result_hisname);
	   final_gr[j]->Draw("ALP");

	   final_gr[j]->Fit("res_fit", "RE");
	   }
	   */
	//	his_temp->SetTitle("BKG");

	his_temp->SetXTitle("Energy[keV]");
	his_temp->SetYTitle("Counts/keV");

	char his_name[256];
	sprintf(his_name, "RUN%d_%s_single_peak_fitting", runnum, data_his); // ch4 check
	his_temp->SetName(his_name);

	//c1->Modified();
	//c1->Update();
	//c1->SaveAs(res2file);
//	can[j]->Modified();
//	can[j]->Update();
	for(int si=0;si<16;si++){
	can[si]->SaveAs(res2file[si]);
	}

	TFile r3f (res3file,"RECREATE");
	r3f.cd();
	gfitpara->Write();
	p1fitpara->Write();
	m1fitpara->Write();
	m2fitpara->Write();
	m3fitpara->Write();
	finalpara->Write();
	r3f.Close();
	gfitpara->Reset();
	p1fitpara->Reset();
	m1fitpara->Reset();
	m2fitpara->Reset();
	m3fitpara->Reset();
	finalpara->Reset();
}
