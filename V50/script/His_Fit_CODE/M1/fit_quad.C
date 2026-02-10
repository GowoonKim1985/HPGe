#include <iostream>
#include "TMath.h"
#include "TMinuit.h"
/*
Double_t fit_agaus (Double_t *x1, Double_t *par1){

Double_t gpart1 = par1[0]/(sqrt(2*TMath::Pi())*par1[2]);
Double_t gpart2 = exp(-((x1[0]-par1[1])*(x1[0]-par1[1]))/2/par1[2]/par1[2]);
Double_t gpart = gpart1*gpart2;

return gpart;
}

Double_t fit_quad (Double_t *x,  Double_t *par){

Double_t qpart1 = par[0]/(sqrt(2*TMath::Pi())*par[2])*exp(-((x[0]-par[1])*(x[0]-par[1]))/2/par[2]/par[2]);
Double_t qpart2 = par[3]/(sqrt(2*TMath::Pi())*par[5])*exp(-((x[1]-par[4])*(x[0]-par[4]))/2/par[5]/par[5]);
Double_t qpart3 = par[6]/(sqrt(2*TMath::Pi())*par[8])*exp(-((x[2]-par[7])*(x[0]-par[7]))/2/par[8]/par[8]);
Double_t qpart4 = par[9]/(sqrt(2*TMath::Pi())*par[11])*exp(-((x[3]-par[10])*(x[0]-par[10]))/2/par[11]/par[11]);

Double_t qpart = qpart1+qpart2+qpart3+qpart4+par[12];

return qpart;
}
*/
void fit_quad()
{

	ifstream peakdata1("/home/kkw/DAQ/ANA500/script/His_Fit_CODE/peak_quad.dat");
	string line_n1;
	string line_n2;
	int line_max1 = 0;
	int line_max2 = 0;
	double file, evt;
	char rtemp1[256];//nuclear info.
	double rtemp2;//peak energy
	while (peakdata1.good()){
		getline(peakdata1,line_n1);
		peakdata1>>rtemp1>>rtemp2;
		++line_max1;
	}
	cout<<"line max : "<<line_max1<<endl;

	char nu[line_max1-1][256];
	double peak[line_max1-1];

	double readfile, readevt;
	ifstream peakdata2("/home/kkw/DAQ/ANA500/script/His_Fit_CODE/peak_quad.dat");
	while (peakdata2.good()){
		if(line_max2<(line_max1-1)){
			getline(peakdata2, line_n2);
			peakdata2>>rtemp1>>rtemp2;
			sprintf(nu[line_max2],"%s",rtemp1);
			peak[line_max2]=rtemp2;
			cout<<"line "<<line_max2<<" : "<<nu[line_max2]<<" , "<<peak[line_max2]<<endl;
		}
		else{break;}
		++line_max2;
	}


	char hisfile[256];
	char rawname[256];
	char res2file[256];//fit his
	char res3file[256];//fit result

	//cout<<"file name(with out .root) :";
	//scanf("%s",rawname);
	//int hmax=17;
	//	int hmax=1;
	//	char rawhis[hmax][256];
	//	char rawhis[256];

	int runnum;
	char runnumber3[256];
	char runnumber6[256];
	char bintag[256];
	cout<< "run number : ";
	scanf("%i",&runnum);
	sprintf(runnumber6,"%06d",runnum);
	sprintf(runnumber3,"%03d",runnum);


	TF1 * gfit = new TF1("gfit", "gaus");
	TF1 * agfit = new TF1("agfit", "[0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2])");
	//TF1 * agfit = new TF1("agfit", fit_gaus, 0, 4000, 3);
	agfit->SetParNames("Area", "Mean", "Sigma"); //area = sqrt(2pi)*constant*sigma

	TF1 * efit = new TF1("efit", "expo");
	TF1 * e2fit = new TF1("e2fit", "expo");
	TF1 * p0fit = new TF1("p0fit", "pol0");
	TF1 * p1fit = new TF1("p1fit", "pol1");

	//TF1 * m9fit = new TF1("m9fit",fit_quad, 0, 4000, 13);
	TF1 * m9fit = new TF1("m9fit","[0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]) + [3]/(sqrt(2*TMath::Pi())*[5])*exp(-((x-[4])*(x-[4]))/2/[5]/[5]) + [6]/(sqrt(2*TMath::Pi())*[8])*exp(-((x-[7])*(x-[7]))/2/[8]/[8]) + [9]/(sqrt(2*TMath::Pi())*[11])*exp(-((x-[10])*(x-[10]))/2/[11]/[11]) + [12]");

	//m9fit->SetParNames("Area1", "Mean1", "Sigma1", "Area2", "Mean2", "Sigma2", "Area3", "Mean3", "Sigma3","Area4", "Mean4", "Sigma4", "Flat");

	TF1 * m10fit = new TF1("m10fit", "[0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]) + [3]/(sqrt(2*TMath::Pi())*[5])*exp(-((x-[4])*(x-[4]))/2/[5]/[5])+ [6]/(sqrt(2*TMath::Pi())*[8])*exp(-((x-[7])*(x-[7]))/2/[8]/[8]) + [9]/(sqrt(2*TMath::Pi())*[11])*exp(-((x-[10])*(x-[10]))/2/[11]/[11]) + exp([12]+[13]*x)");
//	m10fit->SetParNames("Area1", "Mean1", "Sigma1", "Area2", "Mean2", "Sigma2", "Area3", "Mean3", "Sigma3", "Area4", "Mean4", "Sigma4", "expo1", "expo2");



	TF1 * m11fit = new TF1("m11fit", "[0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]) + [3]/(sqrt(2*TMath::Pi())*[5])*exp(-((x-[4])*(x-[4]))/2/[5]/[5]) + [6]/(sqrt(2*TMath::Pi())*[8])*exp(-((x-[7])*(x-[7]))/2/[8]/[8]) + [9]/(sqrt(2*TMath::Pi())*[11])*exp(-((x-[10])*(x-[10]))/2/[11]/[11]) + [12]+[13]*x");
//	m11fit->SetParNames("Area1", "Mean1", "Sigma1", "Area2", "Mean2", "Sigma2", "Area3", "Mean3", "Sigma3", "Area4", "Mean4", "Sigma4","Intercept", "Slope");


	agfit->SetLineColor(9);
	gfit->SetLineColor(3);
	efit->SetLineColor(6);
	e2fit->SetLineColor(6);
	p1fit->SetLineColor(2);
	m9fit->SetLineColor(kBlue);
	m10fit->SetLineColor(kRed);
	m11fit->SetLineColor(kGreen);

	p0fit->SetLineStyle(7);
	efit->SetLineStyle(7);
	p1fit->SetLineStyle(7);

	TCanvas * c1 = new TCanvas("c1","fitting",1200,800);
	TCanvas * m9fit_peak = new TCanvas("Gaus+P0 fitting","Gaus+P0 fitting",600,400);
	TCanvas * m10fit_peak = new TCanvas("Gaus+Exp fitting","Gaus+Exp fitting",600,400);
	TCanvas * m11fit_peak = new TCanvas("Gaus+P1 fitting","Gaus+P1 fitting",600,400);

//	m9fit_peak->Divide(2,2);	
//	m10fit_peak->Divide(2,2);
//	m11fit_peak->Divide(2,2);



	int maxb, maxx, maxh;
	int maxb1, maxb2;
	double bg1, bg2, bg;
	double cnt1, cnt2, cnt3, cnt1_2, cnt3_2;
	double pa1[line_max1-1][3];//gfit 

	double cal[line_max1-1][2]; //integral fit func - m1fit,m2fit
	double epa[6]; //expo fit 
	int xbin[line_max1-1];
	double energy[4];

	double agfitpa1[4], agfitpa2[4], agfitpa3[4], agfiter1[4], agfiter2[4], agfiter3[4];

	double gfitpa1[4], gfitpa2[4], gfitpa3[4], gfiter1[4], gfiter2[4], gfiter3[4];

	double p1fitpa1, p1fitpa2, p1fiter1, p1fiter2;

	double m9fitpa1, m9fitpa2, m9fitpa3, m9fitpa4, m9fitpa5, m9fitpa6, m9fitpa7, m9fitpa8, m9fitpa9, m9fitpa10, m9fitpa11, m9fitpa12, m9fitpa13;
	double m9fiter1, m9fiter2, m9fiter3, m9fiter4, m9fiter5, m9fiter6, m9fiter7, m9fiter8, m9fiter9, m9fiter10, m9fiter11, m9fiter12, m9fiter13;
	double m9pcnt1, m9pcnter1, m9cpd1, m9dcpd1;
	double m9pcnt2, m9pcnter2, m9cpd2, m9dcpd2;

	double m10fitpa1, m10fitpa2, m10fitpa3, m10fitpa4, m10fitpa5, m10fitpa6, m10fitpa7, m10fitpa8, m10fitpa9, m10fitpa10, m10fitpa11, m10fitpa12, m10fitpa13, m10fitpa14;
	double m10fiter1, m10fiter2, m10fiter3, m10fiter4, m10fiter5, m10fiter6, m10fiter7, m10fiter8, m10fiter9, m10fiter10, m10fiter11, m10fiter12, m10fiter13, m10fiter14;
	double m10pcnt1, m10pcnter1, m10cpd1, m10dcpd1;
	double m10pcnt2, m10pcnter2, m10cpd2, m10dcpd2;

	double m11fitpa1, m11fitpa2, m11fitpa3, m11fitpa4, m11fitpa5, m11fitpa6, m11fitpa7, m11fitpa8, m11fitpa9, m11fitpa10, m11fitpa11, m11fitpa12, m11fitpa13, m11fitpa14;
	double m11fiter1, m11fiter2, m11fiter3, m11fiter4, m11fiter5, m11fiter6, m11fiter7, m11fiter8, m11fiter9, m11fiter10, m11fiter11, m11fiter12, m11fiter13, m11fiter14;
	double m11pcnt1, m11pcnter1, m11cpd1, m11dcpd1;
	double m11pcnt2, m11pcnter2, m11cpd2, m11dcpd2;

	double gfitchi2[4], agfitchi2[4], p1fitchi2, m9fitchi2, m10fitchi2, m11fitchi2;
	double fixpa[3];
	double m1ymax, m2ymax;

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

	TTree *m9fitpara = new TTree("m9fitpara", "gaus+pol0 fitting result");
//	m9fitpara -> Branch("run", &runnum, "run/I");
	m9fitpara -> Branch("peak", &energy, "energy/D");
	m9fitpara -> Branch("area1", &m9fitpa1, "m9fitpa1/D");
	m9fitpara -> Branch("mean1", &m9fitpa2, "m9fitpa2/D");
	m9fitpara -> Branch("sigma1", &m9fitpa3, "m9fitpa3/D");
	m9fitpara -> Branch("area2", &m9fitpa4, "m9fitpa4/D");
	m9fitpara -> Branch("mean2", &m9fitpa5, "m9fitpa5/D");
	m9fitpara -> Branch("sigma2", &m9fitpa6, "m9fitpa6/D");
	m9fitpara -> Branch("area3", &m9fitpa7, "m9fitpa7/D");
	m9fitpara -> Branch("mean3", &m9fitpa8, "m9fitpa8/D");
	m9fitpara -> Branch("sigma3", &m9fitpa9, "m9fitpa9/D");
	m9fitpara -> Branch("area4", &m9fitpa10, "m9fitpa10/D");
	m9fitpara -> Branch("mean4", &m9fitpa11, "m9fitpa11/D");
	m9fitpara -> Branch("sigma4", &m9fitpa12, "m9fitpa12/D");
	m9fitpara -> Branch("pol0", &m9fitpa13, "m9fitpa13/D");
	m9fitpara -> Branch("chi2", &m9fitchi2, "m9fitchi2/D");

	m9fitpara -> Branch("area1_er", &m9fiter1, "m9fiter1/D");
	m9fitpara -> Branch("mean1_er", &m9fiter2, "m9fiter2/D");
	m9fitpara -> Branch("sigma1_er", &m9fiter3, "m9fiter3/D");
	m9fitpara -> Branch("area2_er", &m9fiter4, "m9fiter4/D");
	m9fitpara -> Branch("mean2_er", &m9fiter5, "m9fiter5/D");
	m9fitpara -> Branch("sigma2_er", &m9fiter6, "m9fiter6/D");
	m9fitpara -> Branch("area3_er", &m9fiter7, "m9fiter7/D");
	m9fitpara -> Branch("mean3_er", &m9fiter8, "m9fiter8/D");
	m9fitpara -> Branch("sigma3_er", &m9fiter9, "m9fiter9/D");
	m9fitpara -> Branch("area4_er", &m9fiter10, "m9fiter10/D");
	m9fitpara -> Branch("mean4_er", &m9fiter11, "m9fiter11/D");
	m9fitpara -> Branch("sigma4_er", &m9fiter12, "m9fiter12/D");
	m9fitpara -> Branch("pol0_er", &m9fiter13, "m9fiter13/D");
	m9fitpara -> Branch("pcount1", &m9pcnt1, "m9pcnt1/D");
	m9fitpara -> Branch("pcount2", &m9pcnt2, "m9pcnt2/D");
	m9fitpara -> Branch("pcount1_er", &m9pcnter1, "m9pcnter1/D");
	m9fitpara -> Branch("pcount2_er", &m9pcnter2, "m9pcnter2/D");


	TTree *m10fitpara = new TTree("m10fitpara", "gaus+exp fitting result");
	m10fitpara -> Branch("peak", &energy, "energy/D");
	m10fitpara -> Branch("area1", &m10fitpa1, "m10fitpa1/D");
	m10fitpara -> Branch("mean1", &m10fitpa2, "m10fitpa2/D");
	m10fitpara -> Branch("sigma1", &m10fitpa3, "m10fitpa3/D");
	m10fitpara -> Branch("area2", &m10fitpa4, "m10fitpa4/D");
	m10fitpara -> Branch("mean2", &m10fitpa5, "m10fitpa5/D");
	m10fitpara -> Branch("sigma2", &m10fitpa6, "m10fitpa6/D");
	m10fitpara -> Branch("area3", &m10fitpa7, "m10fitpa7/D");
	m10fitpara -> Branch("mean3", &m10fitpa8, "m10fitpa8/D");
	m10fitpara -> Branch("sigma3", &m10fitpa9, "m10fitpa9/D");
	m10fitpara -> Branch("area4", &m10fitpa10, "m10fitpa10/D");
	m10fitpara -> Branch("mean4", &m10fitpa11, "m10fitpa11/D");
	m10fitpara -> Branch("sigma4", &m10fitpa12, "m10fitpa12/D");
	m10fitpara -> Branch("exp1", &m10fitpa13, "m10fitpa13/D");
	m10fitpara -> Branch("exp2", &m10fitpa14, "m10fitpa14/D");
	m10fitpara -> Branch("chi2", &m10fitchi2, "m10fitchi2/D");

	m10fitpara -> Branch("area1_er", &m10fiter1, "m10fiter1/D");
	m10fitpara -> Branch("mean1_er", &m10fiter2, "m10fiter2/D");
	m10fitpara -> Branch("sigma1_er", &m10fiter3, "m10fiter3/D");
	m10fitpara -> Branch("area2_er", &m10fiter4, "m10fiter4/D");
	m10fitpara -> Branch("mean2_er", &m10fiter5, "m10fiter5/D");
	m10fitpara -> Branch("sigma2_er", &m10fiter6, "m10fiter6/D");
	m10fitpara -> Branch("area3_er", &m10fiter7, "m10fiter7/D");
	m10fitpara -> Branch("mean3_er", &m10fiter8, "m10fiter8/D");
	m10fitpara -> Branch("sigma3_er", &m10fiter9, "m10fiter9/D");
	m10fitpara -> Branch("area4_er", &m10fiter10, "m10fiter10/D");
	m10fitpara -> Branch("mean4_er", &m10fiter11, "m10fiter11/D");
	m10fitpara -> Branch("sigma4_er", &m10fiter12, "m10fiter12/D");
	m10fitpara -> Branch("exp1_er", &m10fiter13, "m10fiter13/D");
	m10fitpara -> Branch("exp2_er", &m10fiter14, "m10fiter14/D");
	m10fitpara -> Branch("pcount1", &m10pcnt1, "m10pcnt1/D");
	m10fitpara -> Branch("pcount2", &m10pcnt2, "m10pcnt2/D");
	m10fitpara -> Branch("pcount1_er", &m10pcnter1, "m10pcnter1/D");
	m10fitpara -> Branch("pcount2_er", &m10pcnter2, "m10pcnter2/D");


	TTree *m11fitpara = new TTree("m11fitpara", "gaus+gaus+pol1 fitting result");
//	m11fitpara -> Branch("run", &runnum, "run/I");
	m11fitpara -> Branch("peak", &energy, "energy/D");
	m11fitpara -> Branch("area1", &m11fitpa1, "m11fitpa1/D");
	m11fitpara -> Branch("mean1", &m11fitpa2, "m11fitpa2/D");
	m11fitpara -> Branch("sigma1", &m11fitpa3, "m11fitpa3/D");
	m11fitpara -> Branch("area2", &m11fitpa4, "m11fitpa4/D");
	m11fitpara -> Branch("mean2", &m11fitpa5, "m11fitpa5/D");
	m11fitpara -> Branch("sigma2", &m11fitpa6, "m11fitpa6/D");
	m11fitpara -> Branch("area3", &m11fitpa7, "m11fitpa7/D");
	m11fitpara -> Branch("mean3", &m11fitpa8, "m11fitpa8/D");
	m11fitpara -> Branch("sigma3", &m11fitpa9, "m11fitpa9/D");
	m11fitpara -> Branch("area4", &m11fitpa10, "m11fitpa10/D");
	m11fitpara -> Branch("mean4", &m11fitpa11, "m11fitpa11/D");
	m11fitpara -> Branch("sigma4", &m11fitpa12, "m11fitpa12/D");
	m11fitpara -> Branch("pol0", &m11fitpa13, "m11fitpa13/D");
	m11fitpara -> Branch("pol1", &m11fitpa14, "m11fitpa14/D");
	m11fitpara -> Branch("chi2", &m11fitchi2, "m11fitchi2/D");

	m11fitpara -> Branch("area1_er", &m11fiter1, "m11fiter1/D");
	m11fitpara -> Branch("mean1_er", &m11fiter2, "m11fiter2/D");
	m11fitpara -> Branch("sigma1_er", &m11fiter3, "m11fiter3/D");
	m11fitpara -> Branch("area2_er", &m11fiter4, "m11fiter4/D");
	m11fitpara -> Branch("mean2_er", &m11fiter5, "m11fiter5/D");
	m11fitpara -> Branch("sigma2_er", &m11fiter6, "m11fiter6/D");
	m11fitpara -> Branch("area3_er", &m11fiter7, "m11fiter7/D");
	m11fitpara -> Branch("mean3_er", &m11fiter8, "m11fiter8/D");
	m11fitpara -> Branch("sigma3_er", &m11fiter9, "m11fiter9/D");
	m11fitpara -> Branch("area4_er", &m11fiter10, "m11fiter10/D");
	m11fitpara -> Branch("mean4_er", &m11fiter11, "m11fiter11/D");
	m11fitpara -> Branch("sigma4_er", &m11fiter12, "m11fiter12/D");
	m11fitpara -> Branch("pol0_er", &m11fiter13, "m11fiter13/D");
	m11fitpara -> Branch("pol1_er", &m11fiter14, "m11fiter14/D");
	m11fitpara -> Branch("pcount1", &m11pcnt1, "m11pcnt1/D");
	m11fitpara -> Branch("pcount2", &m11pcnt2, "m11pcnt2/D");
	m11fitpara -> Branch("pcount1_er", &m11pcnter1, "m11pcnter1/D");
	m11fitpara -> Branch("pcount2_er", &m11pcnter2, "m11pcnter2/D");


	TF1 *m1bg[line_max1-1];
	TF1 *m2bg[line_max1-1];
	char m1bgname[256];
	char m2bgname[256];

	int bin = 4000;// 1bin = 1kev
	//	int bin = 40000;// 1bin = 100ev
	int mul = bin/4000;
	int number_check = 3;
	int bg_range;

	sprintf(hisfile,"/home/kkw/DAQ/ANA500/result/RUN%s/his_%s_M1.root",runnumber3,runnumber6);
	sprintf(res2file,"/home/kkw/DAQ/ANA500/result/RUN%s/fit_quad_%s.C",runnumber3,runnumber6);
	sprintf(res3file,"/home/kkw/DAQ/ANA500/result/RUN%s/fit_quad_%s.root",runnumber3,runnumber6);

	TFile *hf=new TFile(hisfile);
	TH1D * his_temp = new TH1D("his_temp","",bin,0,4000);
	his_temp = (TH1D*)hf->Get("his_tot");
	//		his_temp = (TH1D*)hf->Get("fhis_cday");
	his_temp->SetLineColor(1);



	for(int i = 0;i<(line_max1-1);i++){
		if(i==number_check){	
		number_check = number_check + 4;	
		energy[0] = peak[i-3];
		energy[1] = peak[i-2];
		energy[2] = peak[i-1];
		energy[3] = peak[i];

		xbin[i-3]=his_temp->GetXaxis()->FindBin(peak[i-3]);
		xbin[i-3]=his_temp->GetXaxis()->FindBin(peak[i-2]);
		xbin[i-2]=his_temp->GetXaxis()->FindBin(peak[i-1]);
		xbin[i]=his_temp->GetXaxis()->FindBin(peak[i]);

/*		his_temp->GetXaxis()->SetRange((peak[i-2]-1)*mul,(peak[i-2]+1)*mul);
		maxb1 = his_temp->GetMaximumBin();
		cout<< "max peak1 : "<<maxb1<<endl;
		his_temp->GetXaxis()->SetRange((peak[i]-1)*mul,(peak[i]+1)*mul);
		maxb2 = his_temp->GetMaximumBin();
		cout<< "max peak2 : "<<maxb2<<endl;
*/
		printf("\n");
		printf("\n=========================\n");
		printf("%.2f keV bin : %i	\n",peak[i-3],xbin[i-3]);
		printf("%.2f keV bin : %i	\n",peak[i-2],xbin[i-2]);
		printf("%.2f keV bin : %i	\n",peak[i-1],xbin[i-1]);
		printf("%.2f keV bin : %i	\n",peak[i],xbin[i]);

		c1->cd();
/*		his_temp->GetXaxis()->SetRange(xbin[i-2]-5*mul,xbin[i-2]+5*mul);

		maxb = his_temp->GetMaximumBin();
		maxx = maxb/mul;
		cout<<"maxb(bin) : "<<maxb<<endl;
		cout<<"maxx(kev) : "<<maxx<<endl;
		maxh = his_temp->GetMaximum();
*/
		his_temp->GetXaxis()->SetRange(0,4000*mul);

	
		//efit->SetRange((peak[i-3]-30),(peak[i]+30));
		//p1fit->SetRange((peak[i-3]-15),(peak[i]+15));
		efit->SetRange((peak[i-3]-25),(peak[i]+25));
		p1fit->SetRange((peak[i-3]-25),(peak[i]+25));
		bg_range = 25;

		his_temp->Fit("efit","R+");
		fixpa[1] = efit->GetParameter(0);
		fixpa[2] = efit->GetParameter(1);
		//m2fit->SetRange((maxx-10),(maxx+10));

		//pol1 fit
		his_temp->Fit("p1fit","R+");
		p1fitpa1 = p1fit->GetParameter(0);
		p1fitpa2 = p1fit->GetParameter(1);
		p1fiter1 = p1fit->GetParError(0);
		p1fiter2 = p1fit->GetParError(1);
		p1fitchi2 = (p1fit->GetChisquare())/(p1fit->GetNDF());


		//gaus fit
		cout<<""<<endl;
		cout<<"[gaus 1st peak fitting]"<<endl;
		cout<<""<<endl;


			//1st peak
		gfit->ReleaseParameter(1);
		gfit->SetRange((peak[i-3]-2),(peak[i-3]+2));
		gfit->FixParameter(1,peak[i-3]);
		his_temp->Fit("gfit","R0+");
		gfitpa1[0] = gfit->GetParameter(0);
		gfitpa2[0] = gfit->GetParameter(1);
		gfitpa3[0] = gfit->GetParameter(2);
		gfiter1[0] = gfit->GetParError(0);
		gfiter2[0] = gfit->GetParError(1);
		gfiter3[0] = gfit->GetParError(2);
		gfitchi2[0] = (gfit->GetChisquare())/(gfit->GetNDF());

			//2nd peak
		cout<<""<<endl;
		cout<<"[gaus 2nd peak fitting]"<<endl;

		cout<<""<<endl;
		gfit->ReleaseParameter(1);
		gfit->SetRange((peak[i-2]-2),(peak[i-2]+2));
		gfit->FixParameter(1,peak[i-2]);
		his_temp->Fit("gfit","R0+");
		gfitpa1[1] = gfit->GetParameter(0);
		gfitpa2[1] = gfit->GetParameter(1);
		gfitpa3[1] = gfit->GetParameter(2);
		gfiter1[1] = gfit->GetParError(0);
		gfiter2[1] = gfit->GetParError(1);
		gfiter3[1] = gfit->GetParError(2);
		gfitchi2[1] = (gfit->GetChisquare())/(gfit->GetNDF());

			//3rd peak

		cout<<""<<endl;
		cout<<"[gaus 3rd peak fitting]"<<endl;

		cout<<""<<endl;
		gfit->ReleaseParameter(1);
		gfit->SetRange((peak[i-1]-2),(peak[i-1]+2));
		gfit->FixParameter(1,peak[i-1]);
		his_temp->Fit("gfit","R0+");
		gfitpa1[2] = gfit->GetParameter(0);
		gfitpa2[2] = gfit->GetParameter(1);
		gfitpa3[2] = gfit->GetParameter(2);
		gfiter1[2] = gfit->GetParError(0);
		gfiter2[2] = gfit->GetParError(1);
		gfiter3[2] = gfit->GetParError(2);
		gfitchi2[2] = (gfit->GetChisquare())/(gfit->GetNDF());


		//4th peak


		cout<<""<<endl;
		cout<<"[gaus 4th peak fitting]"<<endl;

		cout<<""<<endl;
		gfit->ReleaseParameter(1);
		gfit->SetRange((peak[i]-2),(peak[i]+2));
		gfit->FixParameter(1,peak[i]);
		his_temp->Fit("gfit","R0+");
		gfitpa1[3] = gfit->GetParameter(0);
		gfitpa2[3] = gfit->GetParameter(1);
		gfitpa3[3] = gfit->GetParameter(2);
		gfiter1[3] = gfit->GetParError(0);
		gfiter2[3] = gfit->GetParError(1);
		gfiter3[3] = gfit->GetParError(2);
		gfitchi2[3] = (gfit->GetChisquare())/(gfit->GetNDF());


		bg1 = 0; bg2 = 0; 
/*		for(int j1=0;j1<(5*mul);j1++){
			bg1 = bg1 + (his_temp->GetBinContent(peak[i-3]-(gfitpa3[0]*3+5)*mul+j1));
			bg2 = bg2 + (his_temp->GetBinContent(peak[i]+(gfitpa3[2]*3+5)*mul-j1));
		}
*/
		for(int j1=0;j1<(bg_range*mul);j1++){
			bg1 = bg1 + (his_temp->GetBinContent(peak[i-3]-(gfitpa3[0]*3+bg_range)*mul+j1));
			bg2 = bg2 + (his_temp->GetBinContent(peak[i]+(gfitpa3[2]*3+bg_range)*mul-j1));
		}

		bg=(bg1+bg2)/(bg_range*2*mul);
		fixpa[0]=bg;


		//area gaus fit

			//1st peak
		cout<<""<<endl;
		cout<<"[area gaus 1st peak fitting]"<<endl;
		cout<<""<<endl;

		agfit->ReleaseParameter(1);
		agfit->SetRange((peak[i-3]-2),(peak[i-3]+2));
		agfit->SetParameters(2.5*gfitpa2[0]*gfitpa3[0],gfitpa2[0],gfitpa3[0]);
		gfit->FixParameter(1,peak[i-3]);

		//agfit->SetParameters(2.5*gfitpa2*gfitpa3,maxx,2);
		his_temp->Fit("agfit","R0+");
		agfitpa1[0] = agfit->GetParameter(0);
		agfitpa2[0] = agfit->GetParameter(1);
		agfitpa3[0] = agfit->GetParameter(2);
		agfiter1[0] = agfit->GetParError(0);
		agfiter2[0] = agfit->GetParError(1);
		agfiter3[0] = agfit->GetParError(2);
		agfitchi2[0] = (agfit->GetChisquare())/(agfit->GetNDF());


			//2nd peak
		cout<<""<<endl;
		cout<<"[area gaus 2nd peak fitting]"<<endl;
		cout<<""<<endl;

		agfit->ReleaseParameter(1);
		agfit->SetRange((peak[i-2]-2),(peak[i-2]+2));
		agfit->SetParameters(2.5*gfitpa2[1]*gfitpa3[1],gfitpa2[1],gfitpa3[1]);
		agfit->FixParameter(1,peak[i-2]);
		his_temp->Fit("agfit","RQ0+");
		agfitpa1[1] = agfit->GetParameter(0);
		agfitpa2[1] = agfit->GetParameter(1);
		agfitpa3[1] = agfit->GetParameter(2);
		agfiter1[1] = agfit->GetParError(0);
		agfiter2[1] = agfit->GetParError(1);
		agfiter3[1] = agfit->GetParError(2);
		agfitchi2[1] = (agfit->GetChisquare())/(agfit->GetNDF());


			//3rd peak
		cout<<""<<endl;
		cout<<"[area gaus 3rd peak fitting]"<<endl;
		cout<<""<<endl;

		agfit->ReleaseParameter(1);
		agfit->SetRange((peak[i-1]-2),(peak[i-1]+2));
		agfit->SetParameters(2.5*gfitpa2[2]*gfitpa3[2],gfitpa2[2],gfitpa3[2]);
		agfit->FixParameter(1,peak[i]);
		his_temp->Fit("agfit","RQ0+");
		agfitpa1[2] = agfit->GetParameter(0);
		agfitpa2[2] = agfit->GetParameter(1);
		agfitpa3[2] = agfit->GetParameter(2);
		agfiter1[2] = agfit->GetParError(0);
		agfiter2[2] = agfit->GetParError(1);
		agfiter3[2] = agfit->GetParError(2);
		agfitchi2[2] = (agfit->GetChisquare())/(agfit->GetNDF());

			//4th peak
		cout<<""<<endl;
		cout<<"[area gaus 4th peak fitting]"<<endl;
		cout<<""<<endl;

		agfit->ReleaseParameter(1);
		agfit->SetRange((peak[i]-2),(peak[i]+2));
		agfit->SetParameters(2.5*gfitpa2[3]*gfitpa3[3],gfitpa2[3],gfitpa3[3]);
		agfit->FixParameter(1,peak[i]);
		his_temp->Fit("agfit","RQ0+");
		agfitpa1[3] = agfit->GetParameter(0);
		agfitpa2[3] = agfit->GetParameter(1);
		agfitpa3[3] = agfit->GetParameter(2);
		agfiter1[3] = agfit->GetParError(0);
		agfiter2[3] = agfit->GetParError(1);
		agfiter3[3] = agfit->GetParError(2);
		agfitchi2[3] = (agfit->GetChisquare())/(agfit->GetNDF());



cout<<"///////m9 fit///////"<<endl;
cout<<endl;
		m9fit_peak->cd((i+1)/4);

		m9fit->ReleaseParameter(0);
		m9fit->ReleaseParameter(1);
		m9fit->ReleaseParameter(2);
		m9fit->ReleaseParameter(3);
		m9fit->ReleaseParameter(4);
		m9fit->ReleaseParameter(5);
		m9fit->ReleaseParameter(6);
		m9fit->ReleaseParameter(7);
		m9fit->ReleaseParameter(8);
		m9fit->ReleaseParameter(9);
		m9fit->ReleaseParameter(10);
		m9fit->ReleaseParameter(11);
		m9fit->ReleaseParameter(12);


		m9fit->SetRange((peak[i-3]-5),(peak[i]+5));

		//m9fit -> SetParameters(agfitpa1[0], peak[i-3], agfitpa3[0], agfitpa1[1], peak[i-2], agfitpa3[1], agfitpa1[2], peak[i-1], agfitpa3[2], agfitpa1[3], peak[i], agfitpa3[3], bg);

		m9fit -> SetParameter(0,agfitpa1[0]);
		m9fit -> SetParameter(1, peak[i-3]);
		m9fit -> SetParameter(2,agfitpa3[0]);
		m9fit -> SetParameter(3, agfitpa1[1]);
		m9fit -> SetParameter(4, peak[i-2]);
		m9fit -> SetParameter(5, agfitpa3[1]);
		m9fit -> SetParameter(6, agfitpa1[2]);
		m9fit -> SetParameter(7, peak[i-1]);
		m9fit -> SetParameter(8, agfitpa3[2]);
		m9fit -> SetParameter(9, agfitpa1[3]);
		m9fit -> SetParameter(10, peak[i]);
		m9fit -> SetParameter(11, agfitpa3[3]);
		m9fit -> SetParameter(12, bg);

		m9fit->FixParameter(1,peak[i-3]);
		m9fit->FixParameter(4,peak[i-2]);
		m9fit->FixParameter(7,peak[i-1]);
		m9fit->FixParameter(10,peak[i]);
/*
		m9fit -> SetParameters(agfitpa1[0], maxb1, agfitpa3[0], agfitpa1[1], maxb2, agfitpa3[0], bg);
		m9fit->FixParameter(1,maxb1);
		m9fit->FixParameter(4,maxb2);
*/
		his_temp -> Fit("m9fit", "R0+"); 
		m9fitpa3 = m9fit->GetParameter(2);		
		m9fitpa6 = m9fit->GetParameter(5);
		m9fitpa9 = m9fit->GetParameter(8);
		m9fitpa12 = m9fit->GetParameter(11);

		m9fit->FixParameter(2,m9fitpa3);
		m9fit->FixParameter(5,m9fitpa6);
		m9fit->FixParameter(8,m9fitpa9);
		m9fit->FixParameter(11,m9fitpa12);
		
		//sigma check
		if(m9fitpa3<0.4){
			cout <<"Too small sigma(1) value"<<endl;
			cout<<">> sigma(1) fix"<<endl;
			m9fit->ReleaseParameter(2);
			if(peak[i-3]<1000){m9fit->FixParameter(2,0.8);}
			else if(peak[i-3]<2000){m9fit->FixParameter(2,1.1);}
			else{m9fit->FixParameter(2,1.6);}
			}

		if(m9fitpa6<0.4){
			cout <<"Too small sigma(2) value"<<endl;
			cout<<">> sigma(2) fix"<<endl;
			m9fit->ReleaseParameter(5);
			if(peak[i-2]<1000){m9fit->FixParameter(5,0.8);}
			else if(peak[i-2]<2000){m9fit->FixParameter(5,1.1);}
			else{m9fit->FixParameter(5,1.6);}
			}

		if(m9fitpa9<0.4){
			cout <<"Too small sigma(3) value"<<endl;
			cout<<">> sigma(3) fix"<<endl;
			m9fit->ReleaseParameter(8);
			if(peak[i-1]<1000){m9fit->FixParameter(8,0.8);}
			else if(peak[i-1]<2000){m9fit->FixParameter(8,1.1);}
			else{m9fit->FixParameter(8,1.6);}
			}

		if(m9fitpa12<0.4){
			cout <<"Too small sigma(4) value"<<endl;
			cout<<">> sigma(4) fix"<<endl;
			m9fit->ReleaseParameter(11);
			if(peak[i]<1000){m9fit->FixParameter(11,0.8);}
			else if(peak[i]<2000){m9fit->FixParameter(11,1.1);}
			else{m9fit->FixParameter(11,1.6);}
			}

		if(m9fitpa3>1.0&&peak[i-3]<1000){
			cout <<"Too big sigma(1) value @ <1000 keV"<<endl;
			cout<<">> sigma fix"<<endl;
			m9fit->ReleaseParameter(2);
			m9fit->FixParameter(2,0.8);
		}	

		if(m9fitpa6>1.0&&peak[i-2]<1000){
			cout <<"Too big sigma(2) value @ <1000 keV"<<endl;
			cout<<">> sigma fix"<<endl;
			m9fit->ReleaseParameter(5);
			m9fit->FixParameter(5,0.8);
		}

		if(m9fitpa9>1.0&&peak[i-1]<1000){
			cout <<"Too big sigma(3) value @ <1000 keV"<<endl;
			cout<<">> sigma fix"<<endl;
			m9fit->ReleaseParameter(8);
			m9fit->FixParameter(8,0.8);
		}

		if(m9fitpa12>1.0&&peak[i]<1000){
			cout <<"Too big sigma(4) value @ <1000 keV"<<endl;
			cout<<">> sigma fix"<<endl;
			m9fit->ReleaseParameter(11);
			m9fit->FixParameter(11,0.8);
		}	
	

		if(m9fitpa3>2.0&&peak[i-3]>=1000&&peak[i-3]<2000){
			cout <<"Too big sigma(1) value @ 1000~2000 keV"<<endl;
			cout<<">> sigma(1) fix"<<endl;
			m9fit->ReleaseParameter(2);
			m9fit->FixParameter(2,1.1);
		}

		if(m9fitpa6>2.0&&peak[i-2]>=1000&&peak[i-2]<2000){
			cout <<"Too big sigma value @ 1000~2000 keV"<<endl;
			cout<<">> sigma(2) fix"<<endl;
			m9fit->ReleaseParameter(5);
			m9fit->FixParameter(5,1.1);
		}
		if(m9fitpa9>2.0&&peak[i-1]>=1000&&peak[i-1]<2000){
			cout <<"Too big sigma value @ 1000~2000 keV"<<endl;
			cout<<">> sigma(3) fix"<<endl;
			m9fit->ReleaseParameter(8);
			m9fit->FixParameter(8,1.1);
		}
		if(m9fitpa9>2.0&&peak[i]>=1000&&peak[i]<2000){
			cout <<"Too big sigma value @ 1000~2000 keV"<<endl;
			cout<<">> sigma(4) fix"<<endl;
			m9fit->ReleaseParameter(11);
			m9fit->FixParameter(11,1.1);
		}

	m9fit->SetParName(0, "Area1");
	m9fit->SetParName(1, "Mean1");
	m9fit->SetParName(2, "Sigma1");
	m9fit->SetParName(3, "Area2");
	m9fit->SetParName(4, "Mean2");
	m9fit->SetParName(5, "Sigma2");
	m9fit->SetParName(6, "Area3");
	m9fit->SetParName(7, "Mean3");
	m9fit->SetParName(8, "Sigma3");
	m9fit->SetParName(9,"Area4");
	m9fit->SetParName(10, "Mean4");
	m9fit->SetParName(11, "Sigma4");
	m9fit->SetParName(12, "Flat");

		his_temp -> Fit("m9fit", "R+"); 
		m9fitpa1 = m9fit->GetParameter(0);
		m9fitpa2 = m9fit->GetParameter(1);
		m9fitpa3 = m9fit->GetParameter(2);
		m9fitpa4 = m9fit->GetParameter(3);
		m9fitpa5 = m9fit->GetParameter(4);
		m9fitpa6 = m9fit->GetParameter(5);
		m9fitpa7 = m9fit->GetParameter(6);
		m9fitpa8 = m9fit->GetParameter(7);
		m9fitpa9 = m9fit->GetParameter(8);
		m9fitpa10 = m9fit->GetParameter(9);
		m9fitpa11 = m9fit->GetParameter(10);
		m9fitpa12 = m9fit->GetParameter(11);
		m9fitpa13 = m9fit->GetParameter(12);

		m9fiter1 = m9fit->GetParError(0);	
		m9fiter2 = m9fit->GetParError(1);	
		m9fiter3 = m9fit->GetParError(2);	
		m9fiter4 = m9fit->GetParError(3);	
		m9fiter5 = m9fit->GetParError(4);	
		m9fiter6 = m9fit->GetParError(5);	
		m9fiter7 = m9fit->GetParError(6);	
		m9fiter8 = m9fit->GetParError(7);
		m9fiter9 = m9fit->GetParError(8);
		m9fiter10 = m9fit->GetParError(9);
		m9fiter11 = m9fit->GetParError(10);
		m9fiter12 = m9fit->GetParError(11);
		m9fiter13 = m9fit->GetParError(12);
		m9fitchi2 = (m9fit->GetChisquare())/(m9fit->GetNDF());
/*
		cout << "[m9fit]" << endl;
		cout << "peak1 : " << m9fitpa2 << endl;
		cout << "m9fit tot. counts (peak1) :" << m9pcnt1 << " +/- " << m9pcnter1 << endl;
//		cout << "m9fit cpd (peak1) :" << m9cpd1 << " +/- " << m9dcpd1 << endl;
		cout << "6 sigma : " << m9fitpa3*6 << endl;
		cout << "peak2 : " << m9fitpa5 << endl;
		cout << "m9fit tot. counts (peak2) :" << m9pcnt2 << " +/- " << m9pcnter2 << endl;
//		cout << "m9fit cpd (peak2) :" << m9cpd2 << " +/- " << m9dcpd2 << endl;
		cout << "6 sigma : " << m9fitpa6*6 << endl;
		cout << "chi2 : " << m9fitchi2 << endl;
*/
		m9fitpara->Fill();
		his_temp->GetXaxis()->SetRange(peak[i-4]-15*mul,peak[i]+15*mul);
		//his_temp->Draw();

		his_temp->SetTitle("double peaks Fitting [Gaus + Exp]");
		his_temp->SetXTitle("Energy[keV]");
		his_temp->SetYTitle("Counts/kev");
		//his_temp->SetYTitle("Counts/kev.day");
		his_temp->SetName("fhis_tot");
		//his_temp->SetName("fit_cday");

		m9fit_peak->Modified();
		m9fit_peak->Update();




//////////m10 fit


cout<<"///////m10 fit///////"<<endl;
cout<<endl;
		m10fit_peak->cd((i+1)/4);

		m10fit->ReleaseParameter(0);
		m10fit->ReleaseParameter(1);
		m10fit->ReleaseParameter(2);
		m10fit->ReleaseParameter(3);
		m10fit->ReleaseParameter(4);
		m10fit->ReleaseParameter(5);
		m10fit->ReleaseParameter(6);
		m10fit->ReleaseParameter(7);
		m10fit->ReleaseParameter(8);
		m10fit->ReleaseParameter(9);
		m10fit->ReleaseParameter(10);
		m10fit->ReleaseParameter(11);
		m10fit->ReleaseParameter(12);
		m10fit->ReleaseParameter(13);


		m10fit->SetRange((peak[i-3]-5),(peak[i]+5));
		//m10fit -> SetParameters(gfitpa1[0], peak[i-3], gfitpa3[0], gfitpa1[1], peak[i-2], gfitpa3[1], gfitpa1[2], peak[i-1], gfitpa3[2], gfitpa1[3], peak[i], gfitpa3[3], fixpa[1], fixpa[2]);

		m10fit -> SetParameter(0,agfitpa1[0]);
		m10fit -> SetParameter(1, peak[i-3]);
		m10fit -> SetParameter(2,agfitpa3[0]);
		m10fit -> SetParameter(3, agfitpa1[1]);
		m10fit -> SetParameter(4, peak[i-2]);
		m10fit -> SetParameter(5, agfitpa3[1]);
		m10fit -> SetParameter(6, agfitpa1[2]);
		m10fit -> SetParameter(7, peak[i-1]);
		m10fit -> SetParameter(8, agfitpa3[2]);
		m10fit -> SetParameter(9, agfitpa1[3]);
		m10fit -> SetParameter(10, peak[i]);
		m10fit -> SetParameter(11, agfitpa3[3]);
		m10fit -> SetParameter(12, fixpa[1]);
		m10fit -> SetParameter(13, fixpa[2]);
		
		m10fit->FixParameter(1,peak[i-3]);
		m10fit->FixParameter(4,peak[i-2]);
		m10fit->FixParameter(7,peak[i-1]);
		m10fit->FixParameter(10,peak[i]);

		his_temp -> Fit("m10fit", "R0+"); 
		m10fitpa3 = m10fit->GetParameter(2);		
		m10fitpa6 = m10fit->GetParameter(5);
		m10fitpa9 = m10fit->GetParameter(8);
		m10fitpa12 = m10fit->GetParameter(11);

		m10fit->FixParameter(2,m10fitpa3);
		m10fit->FixParameter(5,m10fitpa6);
		m10fit->FixParameter(8,m10fitpa9);
		m10fit->FixParameter(11,m10fitpa12);

		//sigma check
		if(m10fitpa3<0.4){
			cout <<"Too small sigma(1) value"<<endl;
			cout<<">> sigma(1) fix"<<endl;
			m10fit->ReleaseParameter(2);
			if(peak[i-3]<1000){m10fit->FixParameter(2,0.8);}
			else if(peak[i-3]<2000){m10fit->FixParameter(2,1.1);}
			else{m10fit->FixParameter(2,1.6);}
			}

		if(m10fitpa6<0.4){
			cout <<"Too small sigma(2) value"<<endl;
			cout<<">> sigma(2) fix"<<endl;
			m10fit->ReleaseParameter(5);
			if(peak[i-2]<1000){m10fit->FixParameter(5,0.8);}
			else if(peak[i-2]<2000){m10fit->FixParameter(5,1.1);}
			else{m10fit->FixParameter(5,1.6);}
			}

		if(m10fitpa9<0.4){
			cout <<"Too small sigma(3) value"<<endl;
			cout<<">> sigma(3) fix"<<endl;
			m10fit->ReleaseParameter(8);
			if(peak[i-1]<1000){m10fit->FixParameter(8,0.8);}
			else if(peak[i-1]<2000){m10fit->FixParameter(8,1.1);}
			else{m10fit->FixParameter(8,1.6);}
			}

		if(m10fitpa12<0.4){
			cout <<"Too small sigma(4) value"<<endl;
			cout<<">> sigma(4) fix"<<endl;
			m10fit->ReleaseParameter(11);
			if(peak[i]<1000){m10fit->FixParameter(11,0.8);}
			else if(peak[i]<2000){m10fit->FixParameter(11,1.1);}
			else{m10fit->FixParameter(11,1.6);}
			}


		if(m10fitpa3>1.0&&peak[i-3]<1000){
			cout <<"Too big sigma(1) value @ <1000 keV"<<endl;
			cout<<">> sigma fix"<<endl;
			m10fit->ReleaseParameter(2);
			m10fit->FixParameter(2,0.8);
		}	

		if(m10fitpa6>1.0&&peak[i-2]<1000){
			cout <<"Too big sigma(2) value @ <1000 keV"<<endl;
			cout<<">> sigma fix"<<endl;
			m10fit->ReleaseParameter(5);
			m10fit->FixParameter(5,0.8);
		}	

		if(m10fitpa9>1.0&&peak[i-1]<1000){
			cout <<"Too big sigma(3) value @ <1000 keV"<<endl;
			cout<<">> sigma fix"<<endl;
			m10fit->ReleaseParameter(8);
			m10fit->FixParameter(8,0.8);
		}

		if(m10fitpa9>1.0&&peak[i-1]<1000){
			cout <<"Too big sigma(4) value @ <1000 keV"<<endl;
			cout<<">> sigma fix"<<endl;
			m10fit->ReleaseParameter(11);
			m10fit->FixParameter(11,0.8);
		}	


		if(m10fitpa3>2.0&&peak[i-3]>=1000&&peak[i-3]<2000){
			cout <<"Too big sigma(1) value @ 1000~2000 keV"<<endl;
			cout<<">> sigma(1) fix"<<endl;
			m10fit->ReleaseParameter(2);
			m10fit->FixParameter(2,1.1);
		}

		if(m10fitpa6>2.0&&peak[i-2]>=1000&&peak[i-2]<2000){
			cout <<"Too big sigma value @ 1000~2000 keV"<<endl;
			cout<<">> sigma(2) fix"<<endl;
			m10fit->ReleaseParameter(5);
			m10fit->FixParameter(5,1.1);
		}

		if(m10fitpa9>2.0&&peak[i-1]>=1000&&peak[i-1]<2000){
			cout <<"Too big sigma value @ 1000~2000 keV"<<endl;
			cout<<">> sigma(3) fix"<<endl;
			m10fit->ReleaseParameter(8);
			m10fit->FixParameter(8,1.1);
		}
		if(m10fitpa9>2.0&&peak[i]>=1000&&peak[i]<2000){
			cout <<"Too big sigma value @ 1000~2000 keV"<<endl;
			cout<<">> sigma(4) fix"<<endl;
			m10fit->ReleaseParameter(11);
			m10fit->FixParameter(11,1.1);
		}

	m10fit->SetParName(0, "Area1");
	m10fit->SetParName(1, "Mean1");
	m10fit->SetParName(2, "Sigma1");
	m10fit->SetParName(3, "Area2");
	m10fit->SetParName(4, "Mean2");
	m10fit->SetParName(5, "Sigma2");
	m10fit->SetParName(6, "Area3");
	m10fit->SetParName(7, "Mean3");
	m10fit->SetParName(8, "Sigma3");
	m10fit->SetParName(9,"Area4");
	m10fit->SetParName(10, "Mean4");
	m10fit->SetParName(11, "Sigma4");
	m10fit->SetParName(12, "exp1");
	m10fit->SetParName(13, "exp2");

		his_temp -> Fit("m10fit", "R+"); 
		m10fitpa1 = m10fit->GetParameter(0);
		m10fitpa2 = m10fit->GetParameter(1);
		m10fitpa3 = m10fit->GetParameter(2);
		m10fitpa4 = m10fit->GetParameter(3);
		m10fitpa5 = m10fit->GetParameter(4);
		m10fitpa6 = m10fit->GetParameter(5);
		m10fitpa7 = m10fit->GetParameter(6);
		m10fitpa8 = m10fit->GetParameter(7);
		m10fitpa9 = m10fit->GetParameter(8);
		m10fitpa10 = m10fit->GetParameter(9);
		m10fitpa11 = m10fit->GetParameter(10);
		m10fitpa12 = m10fit->GetParameter(11);	
		m10fitpa13 = m10fit->GetParameter(12);
		m10fitpa14 = m10fit->GetParameter(13);


		m10fiter1 = m10fit->GetParError(0);	
		m10fiter2 = m10fit->GetParError(1);	
		m10fiter3 = m10fit->GetParError(2);	
		m10fiter4 = m10fit->GetParError(3);	
		m10fiter5 = m10fit->GetParError(4);	
		m10fiter6 = m10fit->GetParError(5);	
		m10fiter7 = m10fit->GetParError(6);	
		m10fiter8 = m10fit->GetParError(7);
		m10fiter9 = m10fit->GetParError(8);
		m10fiter10 = m10fit->GetParError(9);
		m10fiter11 = m10fit->GetParError(10);	
		m10fiter12 = m10fit->GetParError(11);
		m10fiter13 = m10fit->GetParError(12);
		m10fiter14 = m10fit->GetParError(13);
		m10fitchi2 = (m10fit->GetChisquare())/(m10fit->GetNDF());
/*
		cout << "[m10fit]" << endl;
		cout << "peak1 : " << m10fitpa2 << endl;
		cout << "m10fit tot. counts (peak1) :" << m10pcnt1 << " +/- " << m10pcnter1 << endl;
//		cout << "m9fit cpd (peak1) :" << m9cpd1 << " +/- " << m9dcpd1 << endl;
		cout << "6 sigma : " << m10fitpa3*6 << endl;
		cout << "peak2 : " << m10fitpa5 << endl;
		cout << "m9fit tot. counts (peak2) :" << m10pcnt2 << " +/- " << m10pcnter2 << endl;
//		cout << "m9fit cpd (peak2) :" << m9cpd2 << " +/- " << m9dcpd2 << endl;
		cout << "6 sigma : " << m10fitpa6*6 << endl;
		cout << "chi2 : " << m10fitchi2 << endl;
*/
		m10fitpara->Fill();

		//m10fit_peak->cd((i+1)/2);
		his_temp->GetXaxis()->SetRange(peak[i-3]-15*mul,peak[i]+15*mul);
		//his_temp->Draw();

		his_temp->SetTitle("quad peaks Fitting [Gaus + Exp]");
		his_temp->SetXTitle("Energy[keV]");
		his_temp->SetYTitle("Counts/kev");
		//his_temp->SetYTitle("Counts/kev.day");
		his_temp->SetName("fhis_tot");
		//his_temp->SetName("fit_cday");

		m10fit_peak->Modified();
		m10fit_peak->Update();






cout<<"///////m11 fit///////"<<endl;
cout<<endl;
//		m11fit_peak->cd((i+1)/3);
		m11fit_peak->cd();
		m11fit->ReleaseParameter(0);
		m11fit->ReleaseParameter(1);
		m11fit->ReleaseParameter(2);
		m11fit->ReleaseParameter(3);
		m11fit->ReleaseParameter(4);
		m11fit->ReleaseParameter(5);
		m11fit->ReleaseParameter(6);
		m11fit->ReleaseParameter(7);
		m11fit->ReleaseParameter(8);
		m11fit->ReleaseParameter(9);
		m11fit->ReleaseParameter(10);
		m11fit->ReleaseParameter(11);
		m11fit->ReleaseParameter(12);
		m11fit->ReleaseParameter(13);


		m11fit->SetRange((peak[i-3]-5),(peak[i]+5));
//		m11fit -> SetParameters(gfitpa1[0], peak[i-3], gfitpa3[0], gfitpa1[1], peak[i-2], gfitpa3[1], gfitpa1[2], peak[i-1], gfitpa3[2], gfitpa1[3], peak[i], gfitpa3[3], p1fitpa1, p1fitpa2);
	
		m11fit -> SetParameter(0,agfitpa1[0]);
		m11fit -> SetParameter(1, peak[i-3]);
		m11fit -> SetParameter(2,agfitpa3[0]);
		m11fit -> SetParameter(3, agfitpa1[1]);
		m11fit -> SetParameter(4, peak[i-2]);
		m11fit -> SetParameter(5, agfitpa3[1]);
		m11fit -> SetParameter(6, agfitpa1[2]);
		m11fit -> SetParameter(7, peak[i-1]);
		m11fit -> SetParameter(8, agfitpa3[2]);
		m11fit -> SetParameter(9, agfitpa1[3]);
		m11fit -> SetParameter(10, peak[i]);
		m11fit -> SetParameter(11, agfitpa3[3]);
		m11fit -> SetParameter(12, p1fitpa1);
		m11fit -> SetParameter(13, p1fitpa2);
	
		m11fit->FixParameter(1,peak[i-3]);
		m11fit->FixParameter(4,peak[i-2]);
		m11fit->FixParameter(7,peak[i-1]);
		m11fit->FixParameter(10,peak[i]);
/*
		m9fit -> SetParameters(agfitpa1[0], maxb1, agfitpa3[0], agfitpa1[1], maxb2, agfitpa3[0], bg);
		m9fit->FixParameter(1,maxb1);
		m9fit->FixParameter(4,maxb2);
*/

		his_temp -> Fit("m11fit", "R0+"); 
		m11fitpa3 = m11fit->GetParameter(2);		
		m11fitpa6 = m11fit->GetParameter(5);
		m11fitpa9 = m11fit->GetParameter(8);
		m11fitpa12 = m11fit->GetParameter(11);

		m11fit->FixParameter(2,m11fitpa3);
		m11fit->FixParameter(5,m11fitpa6);
		m11fit->FixParameter(8,m11fitpa9);
		m11fit->FixParameter(11,m11fitpa12);

		//sigma check
		if(m11fitpa3<0.4){
			cout <<"Too small sigma(1) value"<<endl;
			cout<<">> sigma(1) fix"<<endl;
			m11fit->ReleaseParameter(2);
			if(peak[i-3]<1000){m11fit->FixParameter(2,0.8);}
			else if(peak[i-3]<2000){m11fit->FixParameter(2,1.1);}
			else{m11fit->FixParameter(2,1.6);}
			}

		if(m11fitpa6<0.4){
			cout <<"Too small sigma(2) value"<<endl;
			cout<<">> sigma(2) fix"<<endl;
			m11fit->ReleaseParameter(5);
			if(peak[i-2]<1000){m11fit->FixParameter(5,0.8);}
			else if(peak[i-2]<2000){m11fit->FixParameter(5,1.1);}
			else{m11fit->FixParameter(5,1.6);}
			}

		if(m11fitpa9<0.4){
			cout <<"Too small sigma(3) value"<<endl;
			cout<<">> sigma(3) fix"<<endl;
			m11fit->ReleaseParameter(8);
			if(peak[i-1]<1000){m11fit->FixParameter(8,0.8);}
			else if(peak[i-1]<2000){m11fit->FixParameter(8,1.1);}
			else{m11fit->FixParameter(8,1.6);}
			}

		if(m11fitpa12<0.4){
			cout <<"Too small sigma(4) value"<<endl;
			cout<<">> sigma(4) fix"<<endl;
			m11fit->ReleaseParameter(11);
			if(peak[i]<1000){m11fit->FixParameter(11,0.8);}
			else if(peak[i]<2000){m11fit->FixParameter(11,1.1);}
			else{m11fit->FixParameter(11,1.6);}
			}


		if(m11fitpa3>1.0&&peak[i-3]<1000){
			cout <<"Too big sigma(1) value @ <1000 keV"<<endl;
			cout<<">> sigma fix"<<endl;
			m11fit->ReleaseParameter(2);
			m11fit->FixParameter(2,0.8);
		}	

		if(m11fitpa6>1.0&&peak[i-2]<1000){
			cout <<"Too big sigma(2) value @ <1000 keV"<<endl;
			cout<<">> sigma fix"<<endl;
			m11fit->ReleaseParameter(5);
			m11fit->FixParameter(5,0.8);
		}	

		if(m11fitpa9>1.0&&peak[i-1]<1000){
			cout <<"Too big sigma(3) value @ <1000 keV"<<endl;
			cout<<">> sigma fix"<<endl;
			m11fit->ReleaseParameter(8);
			m11fit->FixParameter(8,0.8);
		}	
		if(m11fitpa12>1.0&&peak[i]<1000){
			cout <<"Too big sigma(4) value @ <1000 keV"<<endl;
			cout<<">> sigma fix"<<endl;
			m11fit->ReleaseParameter(11);
			m11fit->FixParameter(11,0.8);
		}	


		if(m11fitpa3>2.0&&peak[i-3]>=1000&&peak[i-3]<2000){
			cout <<"Too big sigma(1) value @ 1000~2000 keV"<<endl;
			cout<<">> sigma(1) fix"<<endl;
			m11fit->ReleaseParameter(2);
			m11fit->FixParameter(2,1.1);
		}

		if(m11fitpa6>2.0&&peak[i-2]>=1000&&peak[i-2]<2000){
			cout <<"Too big sigma value @ 1000~2000 keV"<<endl;
			cout<<">> sigma(2) fix"<<endl;
			m11fit->ReleaseParameter(5);
			m11fit->FixParameter(5,1.1);
		}

		if(m11fitpa9>2.0&&peak[i-1]>=1000&&peak[i-1]<2000){
			cout <<"Too big sigma value @ 1000~2000 keV"<<endl;
			cout<<">> sigma(3) fix"<<endl;
			m11fit->ReleaseParameter(8);
			m11fit->FixParameter(8,1.1);
		}

		if(m11fitpa12>2.0&&peak[i]>=1000&&peak[i]<2000){
			cout <<"Too big sigma value @ 1000~2000 keV"<<endl;
			cout<<">> sigma(4) fix"<<endl;
			m11fit->ReleaseParameter(11);
			m11fit->FixParameter(11,1.1);
		}
		if(peak[i-3]>1593&&peak[i-3]<1595){
		m11fit->ReleaseParameter(2); m11fit->FixParameter(2,2.0);
		}
		if(peak[i-2]>1593&&peak[i-2]<1595){
		m11fit->ReleaseParameter(5); m11fit->FixParameter(5,2.0);
		}
		if(peak[i-1]>1593&&peak[i-1]<1595){
		m11fit->ReleaseParameter(8); m11fit->FixParameter(8,2.0);
		}
		if(peak[i]>1593&&peak[i]<1595){
		m11fit->ReleaseParameter(11); m11fit->FixParameter(11,2.0);
		}



	m11fit->SetParName(0, "Area1");
	m11fit->SetParName(1, "Mean1");
	m11fit->SetParName(2, "Sigma1");
	m11fit->SetParName(3, "Area2");
	m11fit->SetParName(4, "Mean2");
	m11fit->SetParName(5, "Sigma2");
	m11fit->SetParName(6, "Area3");
	m11fit->SetParName(7, "Mean3");
	m11fit->SetParName(8, "Sigma3");
	m11fit->SetParName(9,"Area4");
	m11fit->SetParName(10, "Mean4");
	m11fit->SetParName(11, "Sigma4");
	m11fit->SetParName(12, "Intercept");
	m11fit->SetParName(13, "Slope");

	his_temp -> Fit("m11fit", "R+"); 
		m11fitpa1 = m11fit->GetParameter(0);
		m11fitpa2 = m11fit->GetParameter(1);
		m11fitpa3 = m11fit->GetParameter(2);
		m11fitpa4 = m11fit->GetParameter(3);
		m11fitpa5 = m11fit->GetParameter(4);
		m11fitpa6 = m11fit->GetParameter(5);
		m11fitpa7 = m11fit->GetParameter(6);
		m11fitpa8 = m11fit->GetParameter(7);
		m11fitpa9 = m11fit->GetParameter(8);
		m11fitpa10 = m11fit->GetParameter(9);
		m11fitpa11 = m11fit->GetParameter(10);
		m11fitpa12 = m11fit->GetParameter(11);
		m11fitpa13 = m11fit->GetParameter(12);
		m11fitpa14 = m11fit->GetParameter(13);


		m11fiter1 = m11fit->GetParError(0);	
		m11fiter2 = m11fit->GetParError(1);	
		m11fiter3 = m11fit->GetParError(2);	
		m11fiter4 = m11fit->GetParError(3);	
		m11fiter5 = m11fit->GetParError(4);	
		m11fiter6 = m11fit->GetParError(5);	
		m11fiter7 = m11fit->GetParError(6);	
		m11fiter8 = m11fit->GetParError(7);
		m11fiter9 = m11fit->GetParError(8);	
		m11fiter10 = m11fit->GetParError(9);	
		m11fiter11 = m11fit->GetParError(10);	
		m11fiter12 = m11fit->GetParError(11);	
		m11fiter13 = m11fit->GetParError(12);	
		m11fiter14 = m11fit->GetParError(13);	
	
		m11fitchi2 = (m11fit->GetChisquare())/(m11fit->GetNDF());
/*
		cout << "[m11fit]" << endl;
		cout << "peak1 : " << m11fitpa2 << endl;
		cout << "m11fit tot. counts (peak1) :" << m11pcnt1 << " +/- " << m11pcnter1 << endl;
//		cout << "m9fit cpd (peak1) :" << m9cpd1 << " +/- " << m9dcpd1 << endl;
		cout << "6 sigma : " << m11fitpa3*6 << endl;
		cout << "peak2 : " << m11fitpa5 << endl;
		cout << "m9fit tot. counts (peak2) :" << m11pcnt2 << " +/- " << m11pcnter2 << endl;
//		cout << "m9fit cpd (peak2) :" << m9cpd2 << " +/- " << m9dcpd2 << endl;
		cout << "6 sigma : " << m11fitpa6*6 << endl;
		cout << "chi2 : " << m11fitchi2 << endl;
*/
		m11fitpara->Fill();

//		m11fit_peak->cd((i+1)/2);
		his_temp->GetXaxis()->SetRange(peak[i-3]-15*mul,peak[i]+15*mul);
//		his_temp->Draw();

		his_temp->SetTitle("quad peaks Fitting [Gaus + P1]");
		his_temp->SetXTitle("Energy[keV]");
		his_temp->SetYTitle("Counts/kev");
		//his_temp->SetYTitle("Counts/kev.day");
		his_temp->SetName("fhis_tot");
		//his_temp->SetName("fit_cday");

		m11fit_peak->Modified();
		m11fit_peak->Update();




                //
}

		TFile r3f (res3file,"RECREATE");
		r3f.cd();
		gfitpara->Write();
		m9fitpara->Write();
		m10fitpara->Write();
		m11fitpara->Write();
		r3f.Close();
		gfitpara->Reset();
		m9fitpara->Reset();
		m10fitpara->Reset();
		m11fitpara->Reset();
	}
	c1->cd();
	his_temp->GetXaxis()->SetRange(0,4000);
	c1->Modified();
	c1->Update();

	m9fit_peak->SaveAs(res2file);
	m10fit_peak->SaveAs(res2file);
	m11fit_peak->SaveAs(res2file);
	c1->SaveAs(res2file);


}
