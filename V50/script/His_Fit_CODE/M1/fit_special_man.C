#include <iostream>
#include "TMath.h"
#include "TMinuit.h"


/*
Double_t fit1(Double_t * v, Double_t * p){
int peak = 6;
        Double_t x = v[0];
	Double_t A[peak];
	Double_t m[peak];
	Double_t s[peak];
	Double_t exp[peak];

	for(int i =0; i<peak;i++){
	A[i] = p[0+(3*i)];
	m[i] = p[1+(3*i)];
	s[i] = p[2+(3*i)];
	exp[i] = -(x-m[i])*(x-m[i])/2./s[i]/s[i];
	}

        Double_t C = p[2+(3*(peak-1))+1];

	
//        double exp1 = -(x-m1)*(x-m1)/2./s1/s1;
  //      double exp2 = -(x-m2)*(x-m2)/2./s2/s2;
    //    double exp3 = -(x-m3)*(x-m3)/2./s3/s3;

        return A[1]*TMath::Exp(exp[1]) + A[2]*TMath::Exp(exp[2]) + A[3]*TMath::Exp(exp[3]) + A[4]*TMath::Exp(exp[4]) + A[5]*TMath::Exp(exp[5]) + A[6]*TMath::Exp(exp[6]) + C;
}



//Double_t finpart = (par[0]/(sqrt(2*TMath::Pi())*par[2])) + (exp(-((x-par[1])*(x-par[1]))/2/par[2]/par[2])) + par[3] +(par[4]/(sqrt(2*TMath::Pi())*par[6])) + (exp(-((x-par[5])*(x-par[5]))/2/par[6]/par[6])) + par[7] +(par[8]/(sqrt(2*TMath::Pi())*par[10])) + (exp(-((x-par[9])*(x-par[9]))/2/par[10]/par[10])) + par[11] +(par[12]/(sqrt(2*TMath::Pi())*par[14])) + (exp(-((x-par[13])*(x-par[13]))/2/par[14]/par[14])) + par[15] +(par[16]/(sqrt(2*TMath::Pi())*par[18])) + (exp(-((x-par[17])*(x-par[17]))/2/par[18]/par[18])) + par[19] +(par[20]/(sqrt(2*TMath::Pi())*par[22])) + (exp(-((x-par[21])*(x-par[21]))/2/par[22]/par[22])) + par[23];

return finpart;

}

*/


void fit_special_man()
{

	char filename[258];
	cout<< "filename(.root) : ";
	scanf("%s",filename);

	char hisfile[256];
	char rawname[256];
	char res2file[256];//fit his
	char res3file[256];//fit result

	int runnum;
	char runnumber3[256];
	char runnumber6[256];
	char bintag[256];
	char peakfile[258];
	cout<< "run number : ";
	scanf("%i",&runnum);
	sprintf(runnumber6,"%06d",runnum);
	sprintf(runnumber3,"%03d",runnum);

	int binnum;
	cout<<"bin number : ";
	scanf("%i",&binnum);

	int peaknum;
	cout<<"peak number : ";
	scanf("%i",&peaknum);

	double peak[peaknum];

	for (int i = 0; i<peaknum;i++){

	cout<<i<<"/"<<peaknum<<" peak(keV) : ";
	scanf("%lf",&peak[i]);
	}	

/*
	sprintf(peakfile,"/home/kkw/DAQ/ANA500/script/His_Fit_CODE/peak_R%03d/peak_special.dat",runnum);

	//ifstream peakdata1("/home/kkw/DAQ/ANA500/script/His_Fit_CODE/peak_Ta.dat");
	ifstream peakdata1(peakfile);
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
	//ifstream peakdata2("/home/kkw/DAQ/ANA500/script/His_Fit_CODE/peak_Ta.dat");
	ifstream peakdata2(peakfile);
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
*/
	double respa1 = 0.000158;
	double respa2 = 0.547238;
	double respa3 = 0.006808;

	
	//	sprintf(rawhis,"/home/kkw/DAQ/ANA500/result/RUN%s/his_%s.root",runnumber3, runnumber6);
	//	sprintf(rawhis,"his_%s.root",runnumber);
	TF1 * gfit = new TF1("gfit","gaus");
	TF1 * agfit = new TF1("agfit","[0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2])");
	agfit->SetParNames("Area","Mean","Sigma"); //area = sqrt(2pi)*constant*sigma
	//	TF1 * efit = new TF1("efit","expo");
	//	TF1 * e2fit = new TF1("e2fit","expo");
//	TF1 * m1fit = new TF1("m1fit","([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+[3]");
//	m1fit->SetParNames("Area","Mean","Sigma","Flat");
	TF1 * m1fit = new TF1("m1fit","([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+[3]");
	m1fit->SetParNames("Area","Mean","Sigma","Flat");
	TF1 * m2fit = new TF1("m2fit","([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+exp([3]+[4]*x)");
	m2fit->SetParNames("Area","Mean","Sigma","expo1","expo2");
//	TF1 * m3fit = new TF1("m3fit","([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+[3]+exp([4]+[5]*x)");
//	m3fit->SetParNames("Area","Mean","Sigma","pol0","expo1","expo2");
	TF1 * efit = new TF1("efit","expo");
	TF1 * pfit = new TF1("pfit","pol0");

//	TF1 * m1sum = new TF1("m1sum","fit1");


TF1 * m1sum;
TF1 * m2sum;

	if(peaknum==2){
m1sum = new TF1("m1sum","([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+([3]/(sqrt(2*TMath::Pi())*[5])*exp(-((x-[4])*(x-[4]))/2/[5]/[5]))+[6]");
m2sum = new TF1("m2sum","([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+([3]/(sqrt(2*TMath::Pi())*[5])*exp(-((x-[4])*(x-[4]))/2/[5]/[5]))+exp([6]+[7]*x)");
}

	if(peaknum==3){
m1sum = new TF1("m1sum","([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+([3]/(sqrt(2*TMath::Pi())*[5])*exp(-((x-[4])*(x-[4]))/2/[5]/[5]))+([6]/(sqrt(2*TMath::Pi())*[8])*exp(-((x-[7])*(x-[7]))/2/[8]/[8]))+[9]");
m2sum = new TF1("m2sum","([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+([3]/(sqrt(2*TMath::Pi())*[5])*exp(-((x-[4])*(x-[4]))/2/[5]/[5]))+([6]/(sqrt(2*TMath::Pi())*[8])*exp(-((x-[7])*(x-[7]))/2/[8]/[8]))+exp([9]+[10]*x)");
}

	if(peaknum==4){
m1sum = new TF1("m1sum","([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+([3]/(sqrt(2*TMath::Pi())*[5])*exp(-((x-[4])*(x-[4]))/2/[5]/[5]))+([6]/(sqrt(2*TMath::Pi())*[8])*exp(-((x-[7])*(x-[7]))/2/[8]/[8]))+([9]/(sqrt(2*TMath::Pi())*[11])*exp(-((x-[10])*(x-[10]))/2/[11]/[11]))+[12]");
m2sum = new TF1("m2sum","([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+([3]/(sqrt(2*TMath::Pi())*[5])*exp(-((x-[4])*(x-[4]))/2/[5]/[5]))+([6]/(sqrt(2*TMath::Pi())*[8])*exp(-((x-[7])*(x-[7]))/2/[8]/[8]))+([9]/(sqrt(2*TMath::Pi())*[11])*exp(-((x-[10])*(x-[10]))/2/[11]/[11]))+exp([12]+[13])*x");
}

	if(peaknum==6){
TF1 * m1sum = new TF1("m1sum","([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+([3]/(sqrt(2*TMath::Pi())*[5])*exp(-((x-[4])*(x-[4]))/2/[5]/[5]))+([6]/(sqrt(2*TMath::Pi())*[8])*exp(-((x-[7])*(x-[7]))/2/[8]/[8]))+([9]/(sqrt(2*TMath::Pi())*[11])*exp(-((x-[10])*(x-[10]))/2/[11]/[11]))+([12]/(sqrt(2*TMath::Pi())*[14])*exp(-((x-[13])*(x-[13]))/2/[14]/[14]))+([15]/(sqrt(2*TMath::Pi())*[17])*exp(-((x-[16])*(x-[16]))/2/[17]/[17]))+[18]");
}




//TF1 * m1sum = new TF1("m1sum","([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+([3]/(sqrt(2*TMath::Pi())*[5])*exp(-((x-[4])*(x-[4]))/2/[5]/[5]))+([6]/(sqrt(2*TMath::Pi())*[8])*exp(-((x-[7])*(x-[7]))/2/[8]/[8]))+([9]/(sqrt(2*TMath::Pi())*[11])*exp(-((x-[10])*(x-[10]))/2/[11]/[11]))+([12]/(sqrt(2*TMath::Pi())*[14])*exp(-((x-[13])*(x-[13]))/2/[15]/[15]))+([15]/(sqrt(2*TMath::Pi())*[17])*exp(-((x-[16])*(x-[16]))/2/[17]/[17]))+exp([18]+[19]*x)");






	TCanvas * c1 = new TCanvas("c1","fitting",1200,800);
	TCanvas * c2 = new TCanvas("c2","best fitting",1200,800);
//	TCanvas * div_peak1 = new TCanvas("fitting detail 1","fitting detail 1",1200,800);
//	div_peak1->Divide(3,3);

//	TCanvas * div_peak2;
//	TCanvas * div_peak3;
/*
		if((line_max1-1)>=9){
		div_peak2 = new TCanvas("fitting detail 2","fitting detail 2",1200,800);
		div_peak2->Divide(3,3);
		}

		if((line_max1-1)>=18){
		div_peak3 = new TCanvas("fitting detail 3","fitting detail 3",1200,800);
		div_peak3->Divide(3,3);
		}
*/
	//	TCanvas * div_peak4 = new TCanvas("fitting detail 4","fitting detail 4",1200,800);
	//	TCanvas * div_peak5 = new TCanvas("fitting detail 5","fitting detail 5",1200,800);

	//div_peak4->Divide(3,3);

	//	if(line_max1>10){div_peak -> Divide(4,4);}
	/*	if(line_max1<=17&&line_max1>10){div_peak -> Divide(4,4);}
		if(line_max1<=26&&line_max1>17){div_peak -> Divide(5,5);}
		if(line_max1<=31&&line_max1>26){div_peak -> Divide(5,6);}
		if(line_max1<=37&&line_max1>31){div_peak -> Divide(6,6);}
		if(line_max1<=42&&line_max1>37){div_peak -> Divide(6,7);}
		*/
	//	agfit->SetLineColor(9);
	//	gfit->SetLineColor(8);
	m1fit->SetLineColor(4);
	m2fit->SetLineColor(2);
	m1sum->SetLineColor(4);
	m2sum->SetLineColor(2);
	pfit->SetLineColor(5);
	efit->SetLineColor(7);
	pfit->SetLineStyle(7);
	efit->SetLineStyle(7);
	int maxb[peaknum], maxx[peaknum], maxh[peaknum], run;
	double bg1, bg2, bg[peaknum];
	double cnt1, cnt2, cnt3;
	double pa1[peaknum][3];//gfit 
	double cal[peaknum][2]; //integral fit func - m1fit,m2fit
	double epa[6]; //expo fit 
	int xbin[peaknum];
	double energy[peaknum];

	double agfitpa[peaknum][3], agfiter[peaknum][3];
	double gfitpa[peaknum][3], gfiter[peaknum][3];


//	double agfitpa1, agfitpa2, agfitpa3, agfiter1, agfiter2, agfiter3;
//	double gfitpa1, gfitpa2, gfitpa3, gfiter1, gfiter2, gfiter3;
//	double m1fitpa1, m1fitpa2, m1fitpa3, m1fitpa4, m1fiter1, m1fiter2, m1fiter3, m1fiter4;
//	double m2fitpa1, m2fitpa2, m2fitpa3, m2fitpa4, m2fitpa5, m2fiter1, m2fiter2, m2fiter3, m2fiter4, m2fiter5;
	double m2fitpa[peaknum][5], m2fiter[peaknum][5];
	//double fitpa1, fitpa2, fitpa3, fiter1;
	double s1fitpa[peaknum][3], s1fiter[peaknum];
	double s2fitpa[peaknum][3], s2fiter[peaknum];
	double fixpa[3];//pol, efit1, efit2
	double intpa1, intpa2, intpa3, intpa4, intpa5, intpa6, intpa7;
	double m1pcnt, m2pcnt, m3pcnt, intpcnt;
	double m1pcnter, m2pcnter, m3pcnter , intpcnter;
	double m1pcnter1, m1pcnter2, m2pcnter1, m2pcnter2, m3pcnter1, m3pcnter2, intE1, intE2, intE3;
	double gfitchi2,agfitchi2,m1fitchi2,m2fitchi2;
	double m1ymax, m2ymax;




	run=runnum;
	TString m1stat;

	TTree *gfitpara = new TTree("gfitpara","gaus fitting result");
	gfitpara->Branch("run",&run,"run/I");
	gfitpara ->Branch("peak",&energy,"energy/D");
	//gfitpara ->Branch("const",&gfitpa[0],"gfitpa1/D");
	//gfitpara ->Branch("mean",&gfitpa[1],"gfitpa2/D");
	gfitpara ->Branch("fit_par",&gfitpa,"gfitpa/D");
	gfitpara ->Branch("chi2",&gfitchi2,"gfitchi2/D");
	gfitpara ->Branch("fit_er",&gfiter,"gfiter/D");
//	gfitpara ->Branch("const_er",&gfiter1,"gfiter1/D");
//	gfitpara ->Branch("mean_er",&gfiter2,"gfiter2/D");
//	gfitpara ->Branch("sigma_er",&gfiter3,"gfiter3/D");

	TTree *agfitpara = new TTree("agfitpara","area-gaus fitting result");
	agfitpara->Branch("run",&run,"run/I");
	agfitpara ->Branch("peak",&energy,"energy/D");
//	agfitpara ->Branch("area",&agfitpa1,"gfitpa1/D");
//	agfitpara ->Branch("mean",&agfitpa2,"gfitpa2/D");
	agfitpara ->Branch("chi2",&agfitchi2,"gfitchi2/D");
//	agfitpara ->Branch("area_er",&agfiter1,"gfiter1/D");
//	agfitpara ->Branch("mean_er",&agfiter2,"gfiter2/D");
//	agfitpara ->Branch("sigma_er",&agfiter3,"gfiter3/D");
	agfitpara ->Branch("fit_par",&agfitpa,"agfitpa/D");
	agfitpara ->Branch("fit_er",&agfiter,"agfiter/D");
/*
	TTree *m1fitpara = new TTree("m1fitpara","gaus+pol fitting result");
	m1fitpara->Branch("run",&run,"run/I");
	m1fitpara ->Branch("peak",&energy,"energy/D");
	m1fitpara ->Branch("area",&m1fitpa1,"m1fitpa1/D");
	m1fitpara ->Branch("mean",&m1fitpa2,"m1fitpa2/D");
	m1fitpara ->Branch("sigma",&m1fitpa3,"m1fitpa3/D");
	m1fitpara ->Branch("chi2",&m1fitchi2,"m1fitchi2/D");
	m1fitpara ->Branch("pol0",&m1fitpa4,"m1fitpa4/D");
	m1fitpara ->Branch("area_er",&m1fiter1,"m1fiter1/D");
	m1fitpara ->Branch("mean_er",&m1fiter2,"m1fiter2/D");
	m1fitpara ->Branch("sigma_er",&m1fiter3,"m1fiter3/D");
	m1fitpara ->Branch("pol0_er",&m1fiter4,"m1fitper4/D");
	m1fitpara ->Branch("pcount",&m1pcnt,"m1pcnt/D");
	m1fitpara ->Branch("pcount_er",&m1pcnter,"m1pcnter/D");

	TTree *m2fitpara = new TTree("m2fitpara","gaus+expo fitting result");
	m2fitpara->Branch("run",&run,"run/I");
	m2fitpara ->Branch("peak",&energy,"energy/D");
	m2fitpara ->Branch("area",&m2fitpa1,"m2fitpa1/D");
	m2fitpara ->Branch("mean",&m2fitpa2,"m2fitpa2/D");
	m2fitpara ->Branch("sigma",&m2fitpa3,"m2fitpa3/D");
	m2fitpara ->Branch("chi2",&m2fitchi2,"m2fitchi2/D");
	m2fitpara ->Branch("exp1",&m2fitpa4,"m2fitpa4/D");
	m2fitpara ->Branch("exp2",&m2fitpa5,"m2fitpa5/D");
	m2fitpara ->Branch("area_er",&m2fiter1,"m2fiter1/D");
	m2fitpara ->Branch("mean_er",&m2fiter2,"m2fiter2/D");
	m2fitpara ->Branch("sigma_er",&m2fiter3,"m2fiter3/D");
	m2fitpara ->Branch("exp1_er",&m2fiter4,"m2fitper4/D");
	m2fitpara ->Branch("exp2_er",&m2fiter5,"m2fitper5/D");
	m2fitpara ->Branch("pcount",&m2pcnt,"m2pcnt/D");
	m2fitpara ->Branch("pcount_er",&m2pcnter,"m2pcnter/D");
*/
	int fit_type;
	double s1fitpa1, s1fitpa2, s1fitpa3, s1fitchi2, s1fiter1;
	TTree *s1para = new TTree("s1para","Multi Gaus+pol0 fitting result");
	s1para->Branch("run",&run,"run/I");
	s1para ->Branch("mean",&s1fitpa2,"s1fitpa2/D");
	s1para ->Branch("area",&s1fitpa1,"s1fitpa1/D");
	s1para ->Branch("area_er",&s1fiter1,"s1fiter1/D");
	s1para ->Branch("sigma",&s1fitpa3,"s1fitpa3/D");
//	m3fitpara ->Branch("pcount",&m3pcnt,"m3pcnt/D");
//	m3fitpara ->Branch("pcount_er",&m3pcnter,"m3pcnter/D");
	double s2fitpa1, s2fitpa2, s2fitpa3, s2fitchi2, s2fiter1;
	TTree *s2para = new TTree("s2para","Multi Gaus+exp fitting result");
	s2para->Branch("run",&run,"run/I");
	s2para ->Branch("mean",&s2fitpa2,"s2fitpa2/D");
	s2para ->Branch("area",&s2fitpa1,"s2fitpa1/D");
	s2para ->Branch("area_er",&s2fiter1,"s2fiter1/D");
	s2para ->Branch("sigma",&s2fitpa3,"s2fitpa3/D");

	TTree *intpara = new TTree("intpara","integral method result");
	intpara->Branch("run",&run,"run/I");
	intpara ->Branch("peak",&energy,"energy/D");
	intpara ->Branch("e1",&intpa1,"intpa1/D");
	intpara ->Branch("e2",&intpa2,"intpa2/D");
	intpara ->Branch("e3",&intpa3,"intpa3/D");
	intpara ->Branch("e4",&intpa4,"intpa4/D");
	intpara ->Branch("cnt1",&intpa5,"intpa5/D");
	intpara ->Branch("cnt2",&intpa6,"intpa6/D");
	intpara ->Branch("cnt3",&intpa7,"intpa7/D");
	intpara ->Branch("pcount",&intpcnt,"intpcnt/D");
	intpara ->Branch("pcount_er",&intpcnter,"intpcnter/D");

	TF1 *m1bg[peaknum];
	TF1 *m2bg[peaknum];
	char m1bgname[256];
	char m2bgname[256];

	int bin = binnum;// 1bin = 1kev
	//	int bin = 40000;// 1bin = 100ev
	int mul = bin/4000;
	//	for(int hn = 0; hn <hmax; hn++){

	//		run = hn;
	//

	sprintf(hisfile,"/home/kkw/DAQ/ANA500/result/RUN%s/Bin%i/his_%s_M1.root",runnumber3,binnum,runnumber6);
	sprintf(res2file,"/home/kkw/DAQ/ANA500/result/RUN%s/Bin%i/%s.C",runnumber3,binnum,filename);
	sprintf(res3file,"/home/kkw/DAQ/ANA500/result/RUN%s/Bin%i/%s.root",runnumber3,binnum,filename);
	//		sprintf(res2file,"/home/kkw/DAQ/ANA25/histogram/SP02/FIT310_1kev/fit_%s.C",rawhis[hn]);
	//		sprintf(res3file,"/home/kkw/DAQ/ANA25/histogram/SP02/FIT310_1kev/fit_%s.root",rawhis[hn]);
	/*

	//100ev
	sprintf(hisfile,"/home/kkw/DAQ/ANA25/histogram/SP02/FWHM310_100ev/phis_%s.root",rawhis[hn]);
	sprintf(res2file,"/home/kkw/DAQ/ANA25/histogram/SP02/FIT310_100ev/fit_%s.C",rawhis[hn]);
	sprintf(res3file,"/home/kkw/DAQ/ANA25/histogram/SP02/FIT310_100ev/fit_%s.root",rawhis[hn]);
	*/
	TFile *hf=new TFile(hisfile);
	TH1D * his_temp = new TH1D("his_temp","",bin,0,4000);
	TH1D * his_full = new TH1D("his_full","",bin,0,4000);
	his_temp = (TH1D*)hf->Get("his_tot");
	his_full = (TH1D*)hf->Get("his_tot");
	his_temp->SetName("his_temp");
	his_full->SetName("fhis_full");
	//		his_temp = (TH1D*)hf->Get("fhis_cday");
	his_temp->SetLineColor(1);
	//	c1->cd();
	/*
	   efit->SetRange(500,2200);
	   his_temp->Fit("efit","RQ+");
	   efit->GetParameters(&epa[0]);
	   efit->SetRange(2500,2800);
	   his_temp->Fit("efit","RQ+");
	   efit->GetParameters(&epa[2]);
	   efit->SetRange(3000,4000);
	   his_temp->Fit("efit","RQ+");
	   efit->GetParameters(&epa[4]);
	   */


	double sigma_cal[peaknum];
	//for(int i = 0;i<(line_max1-1);i++){
//	for(int i = 0;i<(line_max1-1);i+peaknum){


	his_temp->Fit("efit","RQ0+");
	fixpa[1] = efit->GetParameter(0);
	fixpa[2] = efit->GetParameter(1);

	for(int i = 0; i<peaknum; i++){



		energy[i] = peak[i];
		xbin[i]=his_temp->GetXaxis()->FindBin(peak[i]);
		printf("\n");
		printf("\n=========================\n");
		printf("%.2f keV bin : %i	\n",peak[i],xbin[i]);

		his_temp->GetXaxis()->SetRange(xbin[i]-1*mul,xbin[i]+1*mul);

		sigma_cal[i] = peak[i]*(respa1 + (respa2/peak[i]) + (respa3/(sqrt(peak[i]))));
		maxb[i] = his_temp->GetMaximumBin();
		maxx[i] = maxb[i]/mul;
		cout<<"maxb(bin) : "<<maxb[i]<<endl;
		cout<<"maxx(kev) : "<<maxx[i]<<endl;
		maxh[i] = his_temp->GetMaximum();
	


	his_temp->GetXaxis()->SetRange(0,4000*mul);
//	his_temp->GetXaxis()->SetRange(peak[0]-30,peak[2]+30);
	gfit->SetRange((maxx[i]-5),(maxx[i]+5));
	agfit->SetRange((maxx[i]-2),(maxx[i]+2));
//	m1fit->SetRange((maxx-10),(maxx+10));
//	m2fit->SetRange((maxx-10),(maxx+10));
	efit->SetRange((maxx[i]-30),(maxx[i]+30));
	pfit->SetRange((maxx[i]-15),(maxx[i]+15));

		//gaus fit
		his_temp->Fit("gfit","RQ0+");
		gfitpa[i][0] = gfit->GetParameter(0);
		gfitpa[i][1] = gfit->GetParameter(1);
		gfitpa[i][2] = gfit->GetParameter(2);
		gfiter[i][0] = gfit->GetParError(0);
		gfiter[i][1] = gfit->GetParError(1);
		gfiter[i][2] = gfit->GetParError(2);

		gfitchi2 = (gfit->GetChisquare())/(gfit->GetNDF());

		bg1 = 0; bg2 = 0; 
		for(int j1=0;j1<(5*mul);j1++){
			bg1 = bg1 + (his_temp->GetBinContent((maxb[i])-(gfitpa[i][2]*3+5)*mul+j1));
			bg2 = bg2 + (his_temp->GetBinContent((maxb[i])+(gfitpa[i][2]*3+5)*mul-j1));
			//		cout<<maxx-10+j<<"kev: "<<bg1<<" / "<<maxx+10-j<<"kev : "<<bg2<<endl;
		}
		bg[i]=(bg1+bg2)/(10*mul);
		fixpa[0]=bg[i];


		//area gaus fit
		agfit->ReleaseParameter(0);
		agfit->ReleaseParameter(1);
		agfit->ReleaseParameter(2);

		//	c1->cd();
		//	his_temp->GetXaxis()->SetRange(0,4000*mul);

		//agfit->SetParameters(2.5*gfitpa2*gfitpa3,gfitpa2,gfitpa3);
//		agfit->SetParameters(2.5*gfitpa[i][1]*gfitpa[i][2],gfitpa[i][1],gfitpa[i][2]);
		agfit->SetParameters(2.5*gfitpa[i][1]*gfitpa[i][2],maxx[i],sigma_cal[i]);
		//agfit->SetParameters(2.5*gfitpa2*gfitpa3,maxx,2);

	//	agfit->FixParameter(1,maxx[i]);
	//	agfit->FixParameter(2,sigma_cal[i]);
		his_temp->Fit("agfit","RQ0+");
		agfitpa[i][0] = agfit->GetParameter(0);
		agfitpa[i][1] = agfit->GetParameter(1);
		agfitpa[i][2] = agfit->GetParameter(2);
		agfiter[i][0] = agfit->GetParError(0);
		agfiter[i][1] = agfit->GetParError(1);
		agfiter[i][2] = agfit->GetParError(2);
		agfitchi2 = (agfit->GetChisquare())/(agfit->GetNDF());




		cout<<""<<endl;
		cout<<"[m2 fitting]"<<endl;
		cout<<""<<endl;

		m2fit->ReleaseParameter(0);
		m2fit->ReleaseParameter(1);
		m2fit->ReleaseParameter(2);
		m2fit->ReleaseParameter(3);

	efit->SetRange(peak[0]-30,peak[peaknum-1]+30);
	his_temp->Fit("efit","R0Q+");
	fixpa[1] = efit->GetParameter(0);
	fixpa[2] = efit->GetParameter(1);

		m2fit->SetRange(peak[i]-5,peak[i]+5);
		//m1fit->SetParameters(agfitpa1,agfitpa2,agfitpa3,bg);
		m2fit->SetParameters(agfitpa[i][0],maxx[i],sigma_cal[i],fixpa[1],fixpa[2]);
		m2fit->FixParameter(1,maxx[i]);
		m2fit->FixParameter(2,sigma_cal[i]);
		m2fit->FixParameter(3,fixpa[1]);
		m2fit->FixParameter(4,fixpa[2]);
		his_temp->Fit("m2fit","RQ0+");

		m2fitpa[i][0] = m2fit->GetParameter(0);
		m2fitpa[i][1] = m2fit->GetParameter(1);
		m2fitpa[i][2] = m2fit->GetParameter(2);
		m2fitpa[i][3] = m2fit->GetParameter(3);
		m2fitpa[i][4] = m2fit->GetParameter(4);




	}


		cout<<""<<endl;
		cout<<"[gaus+pol0 fitting] for "<<peaknum<<" peaks"<<endl;
		cout<<""<<endl;


		bg1 = 0; bg2 = 0; 
		for(int j1=0;j1<(5*mul);j1++){
		bg1 = bg1 + (his_temp->GetBinContent((maxb[0])-(gfitpa[0][2]*3+5)*mul+j1));
		bg2 = bg2 + (his_temp->GetBinContent((maxb[peaknum-1])+(gfitpa[peaknum-1][2]*3+5)*mul-j1));
				}
		fixpa[0]=(bg1+bg2)/(10*mul);

		efit->SetRange(peak[0]-30,peak[peaknum-1]+30);
		his_temp->Fit("efit","R0Q+");
		fixpa[1] = efit->GetParameter(0);
		fixpa[2] = efit->GetParameter(1);

		m1sum->SetRange(peak[0]-5,peak[peaknum-1]+5);
		m2sum->SetRange(peak[0]-5,peak[peaknum-1]+5);

		char parname1[256];
		char parname2[256];
		char parname3[256];


			for(int i =0; i<peaknum; i++){
			sprintf(parname1,"area%i",i+1);
			sprintf(parname2,"mean%i",i+1);
			sprintf(parname3,"sigma%i",i+1);

			m1sum->SetParameter(0+(3*i), agfitpa[i][0]);
			m1sum->SetParameter(1+(3*i), maxx[i]);
			m1sum->SetParameter(1+(3*i), peak[i]);
			m1sum->SetParameter(2+(3*i), agfitpa[i][2]);

			m2sum->SetParameter(0+(3*i), agfitpa[i][0]);
			m2sum->SetParameter(1+(3*i), maxx[i]);
			m2sum->SetParameter(1+(3*i), peak[i]);
			m2sum->SetParameter(2+(3*i), agfitpa[i][2]);

			m1sum->SetParName(0+(3*i), parname1);
			m1sum->SetParName(1+(3*i), parname2);
			m1sum->SetParName(2+(3*i), parname3);
			//m1sum->FixParameter(1+(3*i),maxx[i]);
			m1sum->FixParameter(1+(3*i),peak[i]);
			m1sum->FixParameter(2+(3*i),sigma_cal[i]);

			m2sum->SetParName(0+(3*i), parname1);
			m2sum->SetParName(1+(3*i), parname2);
			m2sum->SetParName(2+(3*i), parname3);
			//m2sum->FixParameter(1+(3*i),maxx[i]);
			m2sum->FixParameter(1+(3*i),peak[i]);
			m2sum->FixParameter(2+(3*i),sigma_cal[i]);

			}

//		m1sum->FixParameter(4, 104.5);
//		m1sum->SetParameter(18, fixpa[0]);

		m1sum->SetParameter(peaknum*3,fixpa[0]);
		
		m2sum->SetParameter(peaknum*3,fixpa[1]);
		m2sum->SetParameter((peaknum*3)+1,fixpa[2]);

		m1sum->SetParName(peaknum*3, "bg");
		m2sum->SetParName(peaknum*3,"exp1");
		m2sum->SetParName((peaknum*3)+1,"exp2");

		his_temp->Fit("m1sum","R0+");
		his_temp->Fit("m2sum","R+");

	pfit->SetRange(peak[0]-10,peak[peaknum-1]+10);
	pfit->FixParameter(0,fixpa[0]);
	his_temp->Fit("pfit","R0Q+");

			for(int i =0; i<peaknum; i++){
			s1fitpa[i][0] = m1sum->GetParameter(0+(3*i));
			s1fitpa[i][1] = m1sum->GetParameter(1+(3*i));
			s1fitpa[i][2] = m1sum->GetParameter(2+(3*i));
			s1fiter[i] = m1sum->GetParError(0+(3*i));

			s2fitpa[i][0] = m2sum->GetParameter(0+(3*i));
			s2fitpa[i][1] = m2sum->GetParameter(1+(3*i));
			s2fitpa[i][2] = m2sum->GetParameter(2+(3*i));
			s2fiter[i] = m2sum->GetParError(0+(3*i));

			s1fitpa1= s1fitpa[i][0];
			s1fitpa2= s1fitpa[i][1];
			s1fitpa3= s1fitpa[i][2];
			s1fiter1 = s1fiter[i];
			
			s2fitpa1= s2fitpa[i][0];
			s2fitpa2= s2fitpa[i][1];
			s2fitpa3= s2fitpa[i][2];
			s2fiter1 = s2fiter[i];

			s1para->Fill();
			s2para->Fill();
			}
		

	//m2sum
	



/*

		cout<<""<<endl;
		cout<<"[gaus+exp fitting] for "<<peaknum<<" peaks"<<endl;
		cout<<""<<endl;
		m1sum->SetRange(peak[0]-10,peak[peaknum-1]+10);

		char parname1[256];
		char parname2[256];
		char parname3[256];

	efit->SetRange(peak[0]-30,peak[peaknum-1]+30);
	his_temp->Fit("efit","R+");
	fixpa[1] = efit->GetParameter(0);
	fixpa[2] = efit->GetParameter(1);

		for(int i =0; i<peaknum; i++){
		sprintf(parname1,"area%i",i+1);
		sprintf(parname2,"mean%i",i+1);
		sprintf(parname3,"sigma%i",i+1);

		m1sum->SetParameter(0+(3*i), agfitpa[i][0]);
		cout<<"area"<<i<<" : "<<agfitpa[i][0];
		m1sum->SetParameter(1+(3*i), maxx[i]);
		m1sum->SetParameter(2+(3*i), agfitpa[i][2]);


		m1sum->SetParName(0+(3*i), parname1);
		m1sum->SetParName(1+(3*i), parname2);
		m1sum->SetParName(2+(3*i), parname3);
		m1sum->FixParameter(1+(3*i),maxx[i]);
		m1sum->FixParameter(2+(3*i),sigma_cal[i]);

		}
		cout<<m2fitpa[4][0]<<endl;
		cout<<m2fitpa[5][0]<<endl;

		m1sum->SetParameter(18, fixpa[1]);
		m1sum->SetParameter(19, fixpa[2]);
	//	m1sum->FixParameter(18, fixpa[1]);
	//	m1sum->FixParameter(19, fixpa[2]);
		m1sum->SetParName(18, "exp1");
		m1sum->SetParName(19, "exp2");

//TF1 * total = new TF1("total","m1fit(0)+m1fit(4)+m1fit(8)+m1fit(12)+m1fit(16)+m1fit(20)");

//m1sum->SetParameters(agfitpa[0][0],maxx[0],agfitpa[0][2],bg[0]);

		his_temp->Fit("m1sum","R+");


*/
/*

		cout<<""<<endl;
		cout<<"[m1 fitting]"<<endl;
		cout<<""<<endl;



		//m1 fit
		//m1fit->SetParameters(gfitpa1,peak[i],gfitpa3,bg);
		//m1fit->SetParameters(agfitpa1,peak[i],agfitpa3,bg);
		m1fit->ReleaseParameter(0);
		m1fit->ReleaseParameter(1);
		m1fit->ReleaseParameter(2);
		m1fit->ReleaseParameter(3);
		//m1fit->SetParameters(agfitpa1,agfitpa2,agfitpa3,bg);
		m1fit->SetParameters(agfitpa1,maxx,agfitpa3,bg);

		his_temp->Fit("m1fit","R0+");
		m1fitpa2 = m1fit->GetParameter(1);
		m1fitpa3 = m1fit->GetParameter(2);

		cout<< "m1 fit peak center : "<<m1fitpa2<<endl;
		cout<< "m1 fit peak sigma : "<<m1fitpa3<<endl;

		//peak check
		if(abs(m1fitpa1-peak[i])>3){
			cout <<"fitting mean "<<m2fitpa2<<" : peak "<<peak[i]<<endl;
			cout<<">> peak fix"<<endl;
			m1fit->FixParameter(1,maxx);
		}
		his_temp->Fit("m1fit","R0+");
		m1fitpa3 = m1fit->GetParameter(2);
		//sigma check
		cout<<"cal_sigma : "<<sigma_cal<<endl;
		if(m1fitpa3<(0.9*sigma_cal)){
			cout <<"Too small sigma value"<<endl;
			cout<<">> sigma fix"<<endl;
			m1fit->FixParameter(2,sigma_cal);
		}
		if(m1fitpa3>(6*sigma_cal)){
			cout <<"Too big sigma value"<<endl;
			cout<<">> sigma fix"<<endl;
			m1fit->FixParameter(2,sigma_cal);
		}



		his_temp->Fit("m1fit","R+");
		m1fitpa1 = m1fit->GetParameter(0);
		m1fitpa2 = m1fit->GetParameter(1);
		m1fitpa3 = m1fit->GetParameter(2);
		m1fitpa4 = m1fit->GetParameter(3);
		m1fiter1 = m1fit->GetParError(0);
		m1fiter2 = m1fit->GetParError(1);
		m1fiter3 = m1fit->GetParError(2);
		m1fiter4 = m1fit->GetParError(3);
		m1fitchi2 = (m1fit->GetChisquare())/(m1fit->GetNDF());

		m1ymax = m1fit->Eval(peak[i]);


	//	pfit->ReleaseParameter(0);
	//	pfit->FixParameter(0,m1fitpa4);
	//	his_temp->Fit("pfit","RQ+");
	//	his_full->Fit("pfit","R+");

		sprintf(m1bgname,"m1bg_%i",i);
		m1bg[i] = new TF1(m1bgname,"pol0",peak[i]-15,peak[i]+15);
		m1bg[i]->SetParameter(0,m1fitpa4);
		m1bg[i]->SetLineColor(4);
		m1bg[i]->SetLineStyle(7);
		//m1bg[i]->Draw("SAME");


		printf("BG : %.2f	\n",bg);
		if(m1fitpa3<0){
			cout<<"!!!!m1fit sigma is minus -> will be get ABS!!!"<<endl;
		}
		m1pcnt = m1fitpa1;
		m1pcnter = m1fiter1;
		cout<<"m1fit int.error : "<<m1pcnter<<endl;

		cout<<energy<<endl;
		m1fitpara->Fill();


		cout<<""<<endl;
		cout<<"[m2 fitting]"<<endl;
		cout<<""<<endl;
		//m2fit
		//		m2fit->SetRange((maxx-5),(maxx+5));
		//m2fit->SetParameters(gfitpa1,peak[i],gfitpa3,fixpa[1],fixpa[2]);
		//m2fit->SetParameters(agfitpa1,peak[i],agfitpa3,fixpa[1],fixpa[2]);

		m2fit->ReleaseParameter(0);
		m2fit->ReleaseParameter(1);
		m2fit->ReleaseParameter(2);
		m2fit->ReleaseParameter(3);
		m2fit->ReleaseParameter(4);

		//m2fit->SetParameters(agfitpa1,agfitpa2,agfitpa3,fixpa[1],fixpa[2]);
		m2fit->SetParameters(agfitpa1,maxx,agfitpa3,fixpa[1],fixpa[2]);
		m2fit->FixParameter(3,fixpa[1]);
		m2fit->FixParameter(4,fixpa[2]);

		his_temp->Fit("m2fit","R0+");
		m2fitpa2 = m2fit->GetParameter(1);
		m2fitpa3 = m2fit->GetParameter(2);

		cout<< "m2 fit peak center : "<<m2fitpa2<<endl;
		cout<< "m2 fit peak sigma : "<<m2fitpa3<<endl;

		//peak check
		if(abs(m2fitpa2-peak[i])>3){
			cout <<"fitting mean "<<m2fitpa2<<" : peak "<<peak[i]<<endl;
			cout<<">> peak fix"<<endl;
			m2fit->FixParameter(1,maxx);
		}
		his_temp->Fit("m2fit","RQ0+");
		m2fitpa3 = m2fit->GetParameter(2);
		//sigma check
		if(m2fitpa3<(0.9*sigma_cal)){
			cout <<"Too small sigma value"<<endl;
			cout<<">> sigma fix"<<endl;
			m2fit->FixParameter(2,sigma_cal);
		}
		if(m2fitpa3>(6*sigma_cal)){
			cout <<"Too big sigma value"<<endl;
			cout<<">> sigma fix"<<endl;
			m2fit->FixParameter(2,sigma_cal);
		}	

		his_temp->Fit("m2fit","R+");

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

		m2ymax = m2fit->Eval(peak[i]);

		sprintf(m2bgname,"m2bg_%i",i);
		m2bg[i] = new TF1(m2bgname,"expo",peak[i]-15,peak[i]+15);
		m2bg[i]->SetParameter(0,m2fitpa4);
		m2bg[i]->SetParameter(1,m2fitpa5);
		m2bg[i]->SetLineColor(2);
		m2bg[i]->SetLineStyle(7);
		m2bg[i]->Draw("SAME");
		if(m2fitpa3<0){
			cout<<"!!!!m2fit sigma is minus -> will be get ABS!!!"<<endl;
		}
		m2pcnt = m2fitpa1;
		m2pcnter = m2fiter1;

		m2fitpara->Fill();

		cout<<""<<endl;
		cout<<""<<endl;
		cout<<">>chi squre check<<"<<endl;
		cout<<""<<endl;
		cout<<"m1fit chi^2/NDF = "<<m1fitchi2<<endl;
		cout<<"m2fit chi^2/NDF = "<<m2fitchi2<<endl;
		cout<<""<<endl;

		if(abs(1-m1fitchi2) > abs(1-m2fitchi2)){
		cout<<"bestfit : m2 fit "<<endl;
		his_temp->Fit("m1fit","RQ+");
		his_temp->Fit("m2fit","RQ+");
//		his_full->Fit("m2fit","RQ+");
		fit_type=2;
		fitpa1 = m2fitpa1;
		fitpa2 = m2fitpa2;
		fitpa3 = m2fitpa3;
		fiter1 = m2fiter1;
		fitchi2 = m2fitchi2;

		fitpara->Fill();
		}	
		else if(abs(1-m1fitchi2) < abs(1-m2fitchi2)){
		cout<<"bestfit : m1 fit "<<endl;
		his_temp->Fit("m1fit","RQ+");
		his_temp->Fit("m2fit","RQ+");
//		his_full->Fit("m1fit","RQ+");
		fit_type=1;
		fitpa1 = m1fitpa1;
		fitpa2 = m1fitpa2;
		fitpa3 = m1fitpa3;
		fiter1 = m1fiter1;
		fitchi2 = m1fitchi2;
		fitpara->Fill();
		}	

		cout<<""<<endl;
		cout<<""<<endl;

		//integral.
		cout<<"[auto int.]"<<endl;
		cnt1=0; cnt2=0; cnt3=0;
		int intrange;
		int intrangebin;
		intrange = (int)abs(m1fitpa3*3); //kev
		if(intrange>5){
			if(peak[i]<1000){
	//			cout<<"!!!!int.range(3 sigma) sigma is too big -> will be changed to 2keV(0~1000keV)!!!"<<endl;
				intrange = 2;
			}
			else{
	//			cout<<"!!!!int.range(3 sigma) sigma is too big -> will be changed to 3keV(1000~4000keV)!!!"<<endl;
				intrange = 3;
			}
		}
		intrangebin = intrange*mul;
		//		cout<<"center(bin) : "<<xbin[i]<<endl;
		//		cout<<"center(kev) : "<<his_temp->GetBinCenter(xbin[i])<<endl;
		//		cout<<"intrangebin : "<<intrangebin<<endl;
		//		cout<<"maxx-int.(kev) : "<<his_temp->GetBinCenter(maxx-intrangebin)<<endl;
		for(int j3=0;j3<(intrangebin*2+1);j3++){
			cnt2 = cnt2+(his_temp->GetBinContent(xbin[i]-intrangebin+j3));
		}
		for(int j4=0;j4<(5*mul);j4++){
			cnt1 = cnt1+(his_temp->GetBinContent(xbin[i]-intrangebin-(5*mul)+j4));//5kev
			cnt3 = cnt3+(his_temp->GetBinContent(xbin[i]+intrangebin+(5*mul)-j4));//5kev
		}
		intpa1=(his_temp->GetBinCenter(xbin[i]-intrangebin-(5*mul)))-0.5/mul;
		intpa2=(his_temp->GetBinCenter(xbin[i]-intrangebin))-0.5/mul;
		intpa3=(his_temp->GetBinCenter(xbin[i]+intrangebin))+0.5/mul;
		intpa4=(his_temp->GetBinCenter(xbin[i]+intrangebin+(5*mul)))+0.5/mul;
		intpa5=cnt1/mul;
		intpa6=cnt2/mul;
		intpa7=cnt3/mul;
		intE1 = intpa2-intpa1;
		intE2 = intpa3-intpa2;
		intE3 = intpa4-intpa3;
		intpcnt = intpa6 -((intpa5+intpa7)/10)*intE2;
//		cout<<"auto int. peak range : "<<intrangebin*2+1<<"(bin) / "<<intpa2<<" to "<<intpa3<<" (kev)"<<endl;
//		cout<< "auto int. tot. counts:"<<intpcnt<<" +/- "<<intpcnter<<endl;
//		cout<<""<<endl;

		//intpcnter = sqrt(intpa6 + ((intpa5/(intE1*intE1)+intpa7/(intE3*intE3))*((intE2*intE2)/4)));
		intpcnter = sqrt(intpa6 +(intpa5*((intE2*intE2)/(4*intE1*intE1)))+(intpa7*((intE2*intE2)/(4*intE3*intE3))));

		//		m1fitpara->Fill();
		//		m2fitpara->Fill();
		//		intpara->Fill();



		//c1->cd();
	//	his_temp_f->Draw();

		his_full->SetTitle("Ta Fitting");
		his_full->SetXTitle("Energy[keV]");
		his_full->SetYTitle("Counts/kev");
		//his_temp->SetYTitle("Counts/kev.day");
		//his_temp->SetName("fit_cday");
//	}

		c1->cd();
		his_temp->Draw();
		c2->cd();
		his_full->Draw();
		TFile r3f (res3file,"RECREATE");
		r3f.cd();
	//	his_full->Write("fhis_tot");
		gfitpara->Write();
		m1fitpara->Write();
		m2fitpara->Write();
		fitpara->Write();
		intpara->Write();
	//	c1->Write();
	//	div_peak1->Write();
	//	if((line_max1-1)>=9){div_peak2->Write();}
	//	if((line_max1-1)>=18){div_peak3->Write();}

		r3f.Close();
//	his_full->GetXaxis()->SetRange(0,4000);
//	c2->Modified();
//	c2->Update();
//	c1->Modified();
//	c1->Update();
//	div_peak1->SaveAs(res2file);
//	if((line_max1-1)>=9){div_peak2->SaveAs(res2file);}
//	if((line_max1-1)>=18){div_peak3->SaveAs(res2file);}
	c2->SaveAs(res2file);


*/
		TFile r3f (res3file,"RECREATE");
		r3f.cd();
		s1para->Write();
		s2para->Write();
		r3f.Close();
}
