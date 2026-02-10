#include <iostream>
#include "TMath.h"
#include "TMinuit.h"

Double_t fit_agaus (Double_t *x, Double_t *par){

Double_t gpart1 = par[0]/(sqrt(2*TMath::Pi())*par[2]);
Double_t gpart2 = exp(-((x[0]-par[1])*(x[0]-par[1]))/2/par[2]/par[2]);

Double_t gpart = gpart1*gpart2;

return gpart;
}

void fit_triple_Ta()
{
	char hisfile[256];
	char rawname[256];
	char res2file[256];//fit his
	char res3file[256];//fit result

	int runnum;
	char runnumber3[256];
	char runnumber6[256];
	char bintag[256];
	cout<< "run number : ";
	scanf("%i",&runnum);
	sprintf(runnumber6,"%06d",runnum);
	sprintf(runnumber3,"%03d",runnum);

	char peakfile[256];
	sprintf(peakfile,"/home/kkw/DAQ/ANA500/script/His_Fit_CODE/peak_R%03d/peak_triple_Ta.dat",runnum);

	
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




	TF1 * gfit = new TF1("gfit", "gaus");
	//TF1 * agfit = new TF1("agfit", "[0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2])");
	TF1 * agfit = new TF1("agfit", fit_agaus, 0, 4000, 3);
	agfit->SetParNames("Area", "Mean", "Sigma"); //area = sqrt(2pi)*constant*sigma

	TF1 * efit = new TF1("efit", "expo");
	TF1 * e2fit = new TF1("e2fit", "expo");
	TF1 * p0fit = new TF1("p0fit", "pol0");
	TF1 * p1fit = new TF1("p1fit", "pol1");

	TF1 * m6fit = new TF1("m6fit", "[0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]) + [3]/(sqrt(2*TMath::Pi())*[5])*exp(-((x-[4])*(x-[4]))/2/[5]/[5]) +  [6]/(sqrt(2*TMath::Pi())*[8])*exp(-((x-[7])*(x-[7]))/2/[8]/[8])+[9]");
	m6fit->SetParNames("Area1", "Mean1", "Sigma1", "Area2", "Mean2", "Sigma2", "Area3", "Mean3", "Sigma3","Flat");

	TF1 * m7fit = new TF1("m7fit", "[0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]) + [3]/(sqrt(2*TMath::Pi())*[5])*exp(-((x-[4])*(x-[4]))/2/[5]/[5])+ [6]/(sqrt(2*TMath::Pi())*[8])*exp(-((x-[7])*(x-[7]))/2/[8]/[8]) + exp([9]+[10]*x)");
	m7fit->SetParNames("Area1", "Mean1", "Sigma1", "Area2", "Mean2", "Sigma2", "Area3", "Mean3", "Sigma3", "expo1", "expo2");

	TF1 * m8fit = new TF1("m8fit", "[0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]) + [3]/(sqrt(2*TMath::Pi())*[5])*exp(-((x-[4])*(x-[4]))/2/[5]/[5]) + [6]/(sqrt(2*TMath::Pi())*[8])*exp(-((x-[7])*(x-[7]))/2/[8]/[8]) + [9]+[10]*x");
	m8fit->SetParNames("Area1", "Mean1", "Sigma1", "Area2", "Mean2", "Sigma2", "Area3", "Mean3", "Sigma3", "Intercept", "Slope");

	agfit->SetLineColor(9);
	gfit->SetLineColor(3);
	efit->SetLineColor(5);
	e2fit->SetLineColor(6);
	p1fit->SetLineColor(2);
	m6fit->SetLineColor(kBlue);
	m7fit->SetLineColor(kRed);
	m8fit->SetLineColor(kGreen);

	p0fit->SetLineStyle(7);
	efit->SetLineStyle(7);

	TCanvas * c1 = new TCanvas("c1","fitting",1200,800);
	TCanvas * m6fit_peak = new TCanvas("Gaus+P0 fitting","Gaus+P0 fitting",1200,800);
	TCanvas * m7fit_peak = new TCanvas("Gaus+Exp fitting","Gaus+Exp fitting",1200,800);
//	TCanvas * m8fit_peak = new TCanvas("Gaus+P1 fitting","Gaus+P1 fitting",1200,800);
//	TCanvas * div_peak1 = new TCanvas("fitting detail 1","fitting detail 1",1200,800);

	m6fit_peak->Divide(2,2);	
	m7fit_peak->Divide(2,2);
//	m8fit_peak->Divide(2,2);



	int maxb, maxx, maxh;
	int maxb1, maxb2;
	double bg1, bg2, bg;
	double cnt1, cnt2, cnt3, cnt1_2, cnt3_2;
	double pa1[line_max1-1][3];//gfit 


        double cal[line_max1-1][2]; //integral fit func - m1fit,m2fit
        double epa[6]; //expo fit
        int xbin[line_max1-1];
        double energy[3];

        double agfitpa1[3], agfitpa2[3], agfitpa3[3], agfiter1[3], agfiter2[3], agfiter3[3];

        double gfitpa1[3], gfitpa2[3], gfitpa3[3], gfiter1[3], gfiter2[3], gfiter3[3];

        double p1fitpa1, p1fitpa2, p1fiter1, p1fiter2;

        double m6fitpa1, m6fitpa2, m6fitpa3, m6fitpa4, m6fitpa5, m6fitpa6, m6fitpa7, m6fitpa8, m6fitpa9, m6fitpa10;
        double m6fiter1, m6fiter2, m6fiter3, m6fiter4, m6fiter5, m6fiter6, m6fiter7, m6fiter8, m6fiter9, m6fiter10;
        double m6pcnt1, m6pcnter1, m6cpd1, m6dcpd1;
        double m6pcnt2, m6pcnter2, m6cpd2, m6dcpd2;

        double m7fitpa1, m7fitpa2, m7fitpa3, m7fitpa4, m7fitpa5, m7fitpa6, m7fitpa7, m7fitpa8, m7fitpa9, m7fitpa10, m7fitpa11;
        double m7fiter1, m7fiter2, m7fiter3, m7fiter4, m7fiter5, m7fiter6, m7fiter7, m7fiter8, m7fiter9, m7fiter10, m7fiter11;
        double m7pcnt1, m7pcnter1, m7cpd1, m7dcpd1;
        double m7pcnt2, m7pcnter2, m7cpd2, m7dcpd2;

        double m8fitpa1, m8fitpa2, m8fitpa3, m8fitpa4, m8fitpa5, m8fitpa6, m8fitpa7, m8fitpa8, m8fitpa9, m8fitpa10, m8fitpa11;
        double m8fiter1, m8fiter2, m8fiter3, m8fiter4, m8fiter5, m8fiter6, m8fiter7, m8fiter8, m8fiter9, m8fiter10, m8fiter11;
        double m8pcnt1, m8pcnter1, m8cpd1, m8dcpd1;
        double m8pcnt2, m8pcnter2, m8cpd2, m8dcpd2;
	double fitpa1, fitpa2, fitpa3, fitpa4, fitpa5, fitpa6, fitpa7, fitpa8, fitpa9;
	double fiter1, fiter4, fiter7, fitchi2;
        double gfitchi2[3], agfitchi2[3], p1fitchi2, m6fitchi2, m7fitchi2, m8fitchi2;
        double fixpa[3];
        double m1ymax, m2ymax;

        TString m1stat;

	int fit_type;

        TTree *fitpara = new TTree("fit_type", "best fitting result");
        fitpara -> Branch("area1", &fitpa1, "fitpa1/D");
        fitpara -> Branch("area1_er", &fiter1, "fiter1/D");
        fitpara -> Branch("mean1", &fitpa2, "fitpa2/D");
 

        fitpara -> Branch("area2", &fitpa4, "fitpa4/D");
        fitpara -> Branch("area2_er", &fiter4, "fiter4/D");
        fitpara -> Branch("mean2", &fitpa5, "fitpa5/D");

        fitpara -> Branch("area3", &fitpa7, "fitpa7/D");
        fitpara -> Branch("area3_er", &fiter7, "fiter7/D");
        fitpara -> Branch("mean3", &fitpa8, "fitpa8/D");

        fitpara -> Branch("chi2", &fitchi2, "fitchi2/D");


        TTree *gfitpara = new TTree("gfitpara", "gaus fitting result");
//      gfitpara -> Branch("run", &runnum, "run/I");
        gfitpara -> Branch("peak", &energy, "energy/D");
        gfitpara -> Branch("const", &gfitpa1, "gfitpa1/D");
        gfitpara -> Branch("mean", &gfitpa2, "gfitpa2/D");
        gfitpara -> Branch("chi2", &gfitchi2, "gfitchi2/D");
        gfitpara -> Branch("const_er", &gfiter1, "gfiter1/D");
        gfitpara -> Branch("mean_er", &gfiter2, "gfiter2/D");
        gfitpara -> Branch("sigma_er", &gfiter3, "gfiter3/D");

        TTree *p1fitpara = new TTree("p1fitpara", "pol1 fitting result");
//      p1fitpara -> Branch("run", &runnum, "run/I");
        p1fitpara -> Branch("peak", &energy, "energy/D");
        p1fitpara -> Branch("intercept", &p1fitpa1, "p1fitpa1/D");
        p1fitpara -> Branch("slope", &p1fitpa2, "p1fitpa2/D");
        p1fitpara -> Branch("chi2", &p1fitchi2, "p1fitchi2/D");
        p1fitpara -> Branch("intercept_er", &p1fiter1, "p1fiter1/D");
        p1fitpara -> Branch("slope_er", &p1fiter2, "p1fiter2/D");

        TTree *agfitpara = new TTree("agfitpara", "area-gaus fitting result");
//      agfitpara -> Branch("run", &runnum, "run/I");
        agfitpara -> Branch("peak", &energy, "energy/D");
        agfitpara -> Branch("area", &agfitpa1, "gfitpa1/D");
        agfitpara -> Branch("mean", &agfitpa2, "gfitpa2/D");
        agfitpara -> Branch("chi2", &agfitchi2, "gfitchi2/D");
        agfitpara -> Branch("area_er", &agfiter1, "gfiter1/D");
        agfitpara -> Branch("mean_er", &agfiter2, "gfiter2/D");
        agfitpara -> Branch("sigma_er", &agfiter3, "gfiter3/D");



        TTree *m6fitpara = new TTree("m6fitpara", "gaus+pol0 fitting result");
//      m6fitpara -> Branch("run", &runnum, "run/I");
        m6fitpara -> Branch("peak", &energy, "energy/D");
        m6fitpara -> Branch("area1", &m6fitpa1, "m6fitpa1/D");
        m6fitpara -> Branch("mean1", &m6fitpa2, "m6fitpa2/D");
        m6fitpara -> Branch("sigma1", &m6fitpa3, "m6fitpa3/D");
        m6fitpara -> Branch("area2", &m6fitpa4, "m6fitpa4/D");
        m6fitpara -> Branch("mean2", &m6fitpa5, "m6fitpa5/D");
        m6fitpara -> Branch("sigma2", &m6fitpa6, "m6fitpa6/D");
        m6fitpara -> Branch("area3", &m6fitpa7, "m6fitpa7/D");
        m6fitpara -> Branch("mean3", &m6fitpa8, "m6fitpa8/D");
        m6fitpara -> Branch("sigma3", &m6fitpa9, "m6fitpa9/D");
        m6fitpara -> Branch("pol0", &m6fitpa10, "m6fitpa10/D");
        m6fitpara -> Branch("chi2", &m6fitchi2, "m6fitchi2/D");

        m6fitpara -> Branch("area1_er", &m6fiter1, "m6fiter1/D");
        m6fitpara -> Branch("mean1_er", &m6fiter2, "m6fiter2/D");
        m6fitpara -> Branch("sigma1_er", &m6fiter3, "m6fiter3/D");
        m6fitpara -> Branch("area2_er", &m6fiter4, "m6fiter4/D");
        m6fitpara -> Branch("mean2_er", &m6fiter5, "m6fiter5/D");
        m6fitpara -> Branch("sigma2_er", &m6fiter6, "m6fiter6/D");
        m6fitpara -> Branch("area3_er", &m6fiter7, "m6fiter7/D");
        m6fitpara -> Branch("mean3_er", &m6fiter8, "m6fiter8/D");
        m6fitpara -> Branch("sigma3_er", &m6fiter9, "m6fiter9/D");
        m6fitpara -> Branch("pol0_er", &m6fiter10, "m6fiter10/D");
        m6fitpara -> Branch("pcount1", &m6pcnt1, "m6pcnt1/D");
        m6fitpara -> Branch("pcount2", &m6pcnt2, "m6pcnt2/D");
        m6fitpara -> Branch("pcount1_er", &m6pcnter1, "m6pcnter1/D");
        m6fitpara -> Branch("pcount2_er", &m6pcnter2, "m6pcnter2/D");



        TTree *m7fitpara = new TTree("m7fitpara", "gaus+exp fitting result");
        m7fitpara -> Branch("peak", &energy, "energy/D");
        m7fitpara -> Branch("area1", &m7fitpa1, "m7fitpa1/D");
        m7fitpara -> Branch("mean1", &m7fitpa2, "m7fitpa2/D");
        m7fitpara -> Branch("sigma1", &m7fitpa3, "m7fitpa3/D");
        m7fitpara -> Branch("area2", &m7fitpa4, "m7fitpa4/D");
        m7fitpara -> Branch("mean2", &m7fitpa5, "m7fitpa5/D");
        m7fitpara -> Branch("sigma2", &m7fitpa6, "m7fitpa6/D");
        m7fitpara -> Branch("area3", &m7fitpa7, "m7fitpa7/D");
        m7fitpara -> Branch("mean3", &m7fitpa8, "m7fitpa8/D");
        m7fitpara -> Branch("sigma3", &m7fitpa9, "m7fitpa9/D");
        m7fitpara -> Branch("exp1", &m7fitpa10, "m7fitpa10/D");
        m7fitpara -> Branch("exp2", &m7fitpa11, "m7fitpa11/D");
        m7fitpara -> Branch("chi2", &m7fitchi2, "m7fitchi2/D");

        m7fitpara -> Branch("area1_er", &m7fiter1, "m7fiter1/D");
        m7fitpara -> Branch("mean1_er", &m7fiter2, "m7fiter2/D");
        m7fitpara -> Branch("sigma1_er", &m7fiter3, "m7fiter3/D");
        m7fitpara -> Branch("area2_er", &m7fiter4, "m7fiter4/D");
        m7fitpara -> Branch("mean2_er", &m7fiter5, "m7fiter5/D");
        m7fitpara -> Branch("sigma2_er", &m7fiter6, "m7fiter6/D");
        m7fitpara -> Branch("area3_er", &m7fiter7, "m7fiter7/D");
        m7fitpara -> Branch("mean3_er", &m7fiter8, "m7fiter8/D");
        m7fitpara -> Branch("sigma3_er", &m7fiter9, "m7fiter9/D");
        m7fitpara -> Branch("exp1_er", &m7fiter10, "m7fiter10/D");
        m7fitpara -> Branch("exp2_er", &m7fiter11, "m7fiter11/D");
        m7fitpara -> Branch("pcount1", &m7pcnt1, "m7pcnt1/D");
        m7fitpara -> Branch("pcount2", &m7pcnt2, "m7pcnt2/D");
        m7fitpara -> Branch("pcount1_er", &m7pcnter1, "m7pcnter1/D");
        m7fitpara -> Branch("pcount2_er", &m7pcnter2, "m7pcnter2/D");


        TF1 *m1bg[line_max1-1];
        TF1 *m2bg[line_max1-1];
        char m1bgname[256];
        char m2bgname[256];

        int bin = 4000;// 1bin = 1kev
        //      int bin = 40000;// 1bin = 100ev
        int mul = bin/4000;
        int number_check = 2;

        sprintf(hisfile,"/home/kkw/DAQ/ANA500/result/RUN%s/Bin8000/his_%s_M1.root",runnumber3,runnumber6);
        sprintf(res2file,"/home/kkw/DAQ/ANA500/result/RUN%s/Bin8000/fit_triple_%s.C",runnumber3,runnumber6);
        sprintf(res3file,"/home/kkw/DAQ/ANA500/result/RUN%s/Bin8000/fit_triple_%s.root",runnumber3,runnumber6);

        TFile *hf=new TFile(hisfile);
        TH1D * his_temp = new TH1D("his_temp","",bin,0,4000);
        TH1D * his_temp1 = new TH1D("his_temp1","",bin,0,4000);
        TH1D * his_temp2 = new TH1D("his_temp2","",bin,0,4000);
        TH1D * his_temp3 = new TH1D("his_temp3","",bin,0,4000);
        his_temp = (TH1D*)hf->Get("his_tot");
        his_temp1 = (TH1D*)hf->Get("his_tot");
        his_temp2 = (TH1D*)hf->Get("his_tot");
        his_temp3 = (TH1D*)hf->Get("his_tot");
        //              his_temp = (TH1D*)hf->Get("fhis_cday");
        his_temp->SetLineColor(1);

        double respa1 = 0.000158;
        double respa2 = 0.547238;
        double respa3 = 0.006808;



double sigma_cal[3];

        for(int i = 0;i<(line_max1-1);i++){
                if(i==number_check){
                number_check = number_check + 3;

                energy[0] = peak[i-2];
                energy[1] = peak[i-1];
                energy[2] = peak[i];






                sigma_cal[0] = peak[i-2]*(respa1 + (respa2/peak[i-2]) + (respa3/(sqrt(peak[i-2]))));
                sigma_cal[1] = peak[i-1]*(respa1 + (respa2/peak[i-1]) + (respa3/(sqrt(peak[i-1]))));
                sigma_cal[2] = peak[i]*(respa1 + (respa2/peak[i]) + (respa3/(sqrt(peak[i]))));






                xbin[i-2]=his_temp->GetXaxis()->FindBin(peak[i-2]);
                xbin[i-1]=his_temp->GetXaxis()->FindBin(peak[i-1]);
                xbin[i]=his_temp->GetXaxis()->FindBin(peak[i]);


                printf("\n");
                printf("\n=========================\n");
                printf("%.2f keV bin : %i       \n",peak[i-2],xbin[i-2]);
                printf("%.2f keV bin : %i       \n",peak[i-1],xbin[i-1]);
                printf("%.2f keV bin : %i       \n",peak[i],xbin[i]);

                c1->cd();




		his_temp->GetXaxis()->SetRange(0,4000*mul);

	
		efit->SetRange((peak[i-2]-30),(peak[i]+30));
		p1fit->SetRange((peak[i-2]-15),(peak[i]+15));

		his_temp->Fit("efit","R+");
		fixpa[1] = efit->GetParameter(0);
		fixpa[2] = efit->GetParameter(1);
		//m2fit->SetRange((maxx-10),(maxx+10));

		//pol1 fit
		his_temp->Fit("p1fit","RQ0+");
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
		gfit->SetRange((peak[i-2]-2),(peak[i-2]+2));
		gfit->FixParameter(1,peak[i-2]);
		his_temp->Fit("gfit","RQ0+");
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
		gfit->SetRange((peak[i-1]-2),(peak[i-1]+2));
		gfit->FixParameter(1,peak[i-1]);
		his_temp->Fit("gfit","RQ0+");
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
		gfit->SetRange((peak[i]-2),(peak[i]+2));
		gfit->FixParameter(1,peak[i]);
		his_temp->Fit("gfit","RQ0+");
		gfitpa1[2] = gfit->GetParameter(0);
		gfitpa2[2] = gfit->GetParameter(1);
		gfitpa3[2] = gfit->GetParameter(2);
		gfiter1[2] = gfit->GetParError(0);
		gfiter2[2] = gfit->GetParError(1);
		gfiter3[2] = gfit->GetParError(2);
		gfitchi2[2] = (gfit->GetChisquare())/(gfit->GetNDF());

		bg1 = 0; bg2 = 0; 
		for(int j1=0;j1<(5*mul);j1++){
			bg1 = bg1 + (his_temp->GetBinContent(peak[i-2]-(gfitpa3[0]*3+5)*mul+j1));
			bg2 = bg2 + (his_temp->GetBinContent(peak[i]+(gfitpa3[2]*3+5)*mul-j1));
		}

		bg=(bg1+bg2)/(10*mul);
		fixpa[0]=bg;


		//area gaus fit

			//1st peak
		cout<<""<<endl;
		cout<<"[area gaus 1st peak fitting]"<<endl;
		cout<<""<<endl;

		agfit->ReleaseParameter(1);
		agfit->SetRange((peak[i-2]-2),(peak[i-2]+2));
		agfit->SetParameters(2.5*gfitpa2[0]*gfitpa3[0],gfitpa2[0],gfitpa3[0]);
		gfit->FixParameter(1,peak[i-2]);

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
		agfit->SetRange((peak[i-1]-2),(peak[i-1]+2));
		agfit->SetParameters(2.5*gfitpa2[1]*gfitpa3[1],gfitpa2[1],gfitpa3[1]);
		agfit->FixParameter(1,peak[i-1]);
		his_temp->Fit("agfit","R0+");
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
		agfit->SetRange((peak[i]-2),(peak[i]+2));
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


/*
		cout<<""<<endl;
		cout<<"[m1 fitting]"<<endl;
		cout<<""<<endl;
*/

cout<<"///////m6 fit///////"<<endl;
cout<<endl;
		c1->cd();

		m6fit->ReleaseParameter(0);
		m6fit->ReleaseParameter(1);
		m6fit->ReleaseParameter(2);
		m6fit->ReleaseParameter(3);
		m6fit->ReleaseParameter(4);
		m6fit->ReleaseParameter(5);
		m6fit->ReleaseParameter(6);
		m6fit->ReleaseParameter(7);
		m6fit->ReleaseParameter(8);
		m6fit->ReleaseParameter(9);

		m6fit->SetRange((peak[i-2]-5),(peak[i]+5));

		m6fit -> SetParameters(agfitpa1[0], peak[i-2], agfitpa3[0], agfitpa1[1], peak[i-1], agfitpa3[1], agfitpa1[2], peak[i], agfitpa3[2], bg);

		if(agfit)
		m6fit->FixParameter(1,peak[i-2]);
		m6fit->FixParameter(4,peak[i-1]);
		m6fit->FixParameter(7,peak[i]);
/*
		m6fit -> SetParameters(agfitpa1[0], maxb1, agfitpa3[0], agfitpa1[1], maxb2, agfitpa3[0], bg);
		m6fit->FixParameter(1,maxb1);
		m6fit->FixParameter(4,maxb2);
*/
		his_temp -> Fit("m6fit", "R0+"); 
		m6fitpa3 = m6fit->GetParameter(2);		
		m6fitpa6 = m6fit->GetParameter(5);
		m6fitpa9 = m6fit->GetParameter(8);


		//sigma check

                cout<<"1st cal_sigma : "<<sigma_cal[0]<<endl;
                if(m6fitpa3<(0.9*sigma_cal[0])){
                        cout <<"Too small sigma value"<<endl;
                        cout<<">> sigma fix"<<endl;
                        m6fit->FixParameter(2,sigma_cal[0]);
                }
                if(m6fitpa3>(6*sigma_cal[0])){
                        cout <<"Too big sigma value"<<endl;
                        cout<<">> sigma fix"<<endl;
                        m6fit->FixParameter(2,sigma_cal[0]);
                }

                cout<<"2nd cal_sigma : "<<sigma_cal[1]<<endl;
                if(m6fitpa6<(0.9*sigma_cal[1])){
                        cout <<"Too small sigma value"<<endl;
                        cout<<">> sigma fix"<<endl;
                        m6fit->FixParameter(5,sigma_cal[1]);
                }
                if(m6fitpa6>(6*sigma_cal[1])){
                        cout <<"Too big sigma value"<<endl;
                        cout<<">> sigma fix"<<endl;
                        m6fit->FixParameter(5,sigma_cal[1]);
                }

                cout<<"3rd cal_sigma : "<<sigma_cal[2]<<endl;
                if(m6fitpa9<(0.9*sigma_cal[2])){
                        cout <<"Too small sigma value"<<endl;
                        cout<<">> sigma fix"<<endl;
                        m6fit->FixParameter(8,sigma_cal[2]);
                }
                if(m6fitpa9>(6*sigma_cal[2])){
                        cout <<"Too big sigma value"<<endl;
                        cout<<">> sigma fix"<<endl;
                        m6fit->FixParameter(8,sigma_cal[2]);
                }


		c1->cd();

		his_temp -> Fit("m6fit", "R+"); 
		m6fitpa1 = m6fit->GetParameter(0);
		m6fitpa2 = m6fit->GetParameter(1);
		m6fitpa3 = m6fit->GetParameter(2);
		m6fitpa4 = m6fit->GetParameter(3);
		m6fitpa5 = m6fit->GetParameter(4);
		m6fitpa6 = m6fit->GetParameter(5);
		m6fitpa7 = m6fit->GetParameter(6);
		m6fitpa8 = m6fit->GetParameter(7);
		m6fitpa9 = m6fit->GetParameter(8);
		m6fitpa10 = m6fit->GetParameter(9);

		m6fiter1 = m6fit->GetParError(0);	
		m6fiter2 = m6fit->GetParError(1);	
		m6fiter3 = m6fit->GetParError(2);	
		m6fiter4 = m6fit->GetParError(3);	
		m6fiter5 = m6fit->GetParError(4);	
		m6fiter6 = m6fit->GetParError(5);	
		m6fiter7 = m6fit->GetParError(6);	
		m6fiter8 = m6fit->GetParError(7);
		m6fiter9 = m6fit->GetParError(8);
		m6fiter10 = m6fit->GetParError(9);
		m6fitchi2 = (m6fit->GetChisquare())/(m6fit->GetNDF());

		cout << "[m6fit]" << endl;
		cout << "peak1 : " << m6fitpa2 << endl;
		cout << "m6fit tot. counts (peak1) :" << m6pcnt1 << " +/- " << m6pcnter1 << endl;
//		cout << "m6fit cpd (peak1) :" << m6cpd1 << " +/- " << m6dcpd1 << endl;
		cout << "6 sigma : " << m6fitpa3*6 << endl;
		cout << "peak2 : " << m6fitpa5 << endl;
		cout << "m6fit tot. counts (peak2) :" << m6pcnt2 << " +/- " << m6pcnter2 << endl;
//		cout << "m6fit cpd (peak2) :" << m6cpd2 << " +/- " << m6dcpd2 << endl;
		cout << "6 sigma : " << m6fitpa6*6 << endl;
		cout << "chi2 : " << m6fitchi2 << endl;

		m6fitpara->Fill();

		m6fit_peak->cd((i+1)/3);
		his_temp1 -> Fit("m6fit", "RQ"); 
		his_temp1->GetXaxis()->SetRange(peak[i-2]-15*mul,peak[i]+15*mul);
		//his_temp->Draw();

		his_temp1->SetTitle("double peaks Fitting [Gaus + Exp]");
		his_temp1->SetXTitle("Energy[keV]");
		his_temp1->SetYTitle("Counts/kev");
		//his_temp->SetYTitle("Counts/kev.day");
		his_temp1->SetName("fhis_tot");
		//his_temp->SetName("fit_cday");

		m6fit_peak->Modified();
		m6fit_peak->Update();




//////////m7 fit


cout<<"///////m7 fit///////"<<endl;
cout<<endl;
		c1->cd();

		m7fit->ReleaseParameter(0);
		m7fit->ReleaseParameter(1);
		m7fit->ReleaseParameter(2);
		m7fit->ReleaseParameter(3);
		m7fit->ReleaseParameter(4);
		m7fit->ReleaseParameter(5);
		m7fit->ReleaseParameter(6);
		m7fit->ReleaseParameter(7);
		m7fit->ReleaseParameter(8);
		m7fit->ReleaseParameter(9);
		m7fit->ReleaseParameter(10);



		m7fit->SetRange((peak[i-2]-5),(peak[i]+5));
		m7fit -> SetParameters(gfitpa1[0], peak[i-2], gfitpa3[0], gfitpa1[1], peak[i-1], gfitpa3[1], gfitpa1[2], peak[i], gfitpa3[2], fixpa[1], fixpa[2]);
		
		m7fit->FixParameter(1,peak[i-2]);
		m7fit->FixParameter(4,peak[i-1]);
		m7fit->FixParameter(7,peak[i]);

		his_temp -> Fit("m7fit", "R0+"); 
		m7fitpa3 = m7fit->GetParameter(2);		
		m7fitpa6 = m7fit->GetParameter(5);
		m7fitpa9 = m7fit->GetParameter(8);

		//sigma check


                cout<<"1st cal_sigma : "<<sigma_cal[0]<<endl;
                if(m7fitpa3<(0.9*sigma_cal[0])){
                        cout <<"Too small sigma value"<<endl;
                        cout<<">> sigma fix"<<endl;
                        m7fit->FixParameter(2,sigma_cal[0]);
                }
                if(m7fitpa3>(6*sigma_cal[0])){
                        cout <<"Too big sigma value"<<endl;
                        cout<<">> sigma fix"<<endl;
                        m7fit->FixParameter(2,sigma_cal[0]);
                }

                cout<<"2nd cal_sigma : "<<sigma_cal[1]<<endl;
                if(m7fitpa6<(0.9*sigma_cal[1])){
                        cout <<"Too small sigma value"<<endl;
                        cout<<">> sigma fix"<<endl;
                        m7fit->FixParameter(5,sigma_cal[1]);
                }
                if(m7fitpa6>(6*sigma_cal[1])){
                        cout <<"Too big sigma value"<<endl;
                        cout<<">> sigma fix"<<endl;
                        m7fit->FixParameter(5,sigma_cal[1]);
                }

                cout<<"3rd cal_sigma : "<<sigma_cal[2]<<endl;
                if(m7fitpa9<(0.9*sigma_cal[2])){
                        cout <<"Too small sigma value"<<endl;
                        cout<<">> sigma fix"<<endl;
                        m7fit->FixParameter(8,sigma_cal[2]);
                }
                if(m7fitpa9>(6*sigma_cal[2])){
                        cout <<"Too big sigma value"<<endl;
                        cout<<">> sigma fix"<<endl;
                        m7fit->FixParameter(8,sigma_cal[2]);
                }

		c1->cd();
		his_temp -> Fit("m7fit", "R+"); 
		m7fitpa1 = m7fit->GetParameter(0);
		m7fitpa2 = m7fit->GetParameter(1);
		m7fitpa3 = m7fit->GetParameter(2);
		m7fitpa4 = m7fit->GetParameter(3);
		m7fitpa5 = m7fit->GetParameter(4);
		m7fitpa6 = m7fit->GetParameter(5);
		m7fitpa7 = m7fit->GetParameter(6);
		m7fitpa8 = m7fit->GetParameter(7);
		m7fitpa9 = m7fit->GetParameter(8);
		m7fitpa10 = m7fit->GetParameter(9);
		m7fitpa11 = m7fit->GetParameter(10);

		m7fiter1 = m7fit->GetParError(0);	
		m7fiter2 = m7fit->GetParError(1);	
		m7fiter3 = m7fit->GetParError(2);	
		m7fiter4 = m7fit->GetParError(3);	
		m7fiter5 = m7fit->GetParError(4);	
		m7fiter6 = m7fit->GetParError(5);	
		m7fiter7 = m7fit->GetParError(6);	
		m7fiter8 = m7fit->GetParError(7);
		m7fiter9 = m7fit->GetParError(8);
		m7fiter10 = m7fit->GetParError(9);
		m7fiter11 = m7fit->GetParError(10);	
		m7fitchi2 = (m7fit->GetChisquare())/(m7fit->GetNDF());

		cout << "[m7fit]" << endl;
		cout << "peak1 : " << m7fitpa2 << endl;
		cout << "m7fit tot. counts (peak1) :" << m7pcnt1 << " +/- " << m7pcnter1 << endl;
//		cout << "m6fit cpd (peak1) :" << m6cpd1 << " +/- " << m6dcpd1 << endl;
		cout << "6 sigma : " << m7fitpa3*6 << endl;
		cout << "peak2 : " << m7fitpa5 << endl;
		cout << "m6fit tot. counts (peak2) :" << m7pcnt2 << " +/- " << m7pcnter2 << endl;
//		cout << "m6fit cpd (peak2) :" << m6cpd2 << " +/- " << m6dcpd2 << endl;
		cout << "6 sigma : " << m7fitpa6*6 << endl;
		cout << "chi2 : " << m7fitchi2 << endl;

		m7fitpara->Fill();

		m7fit_peak->cd((i+1)/3);
		his_temp2 -> Fit("m7fit", "R"); 
		//m7fit_peak->cd((i+1)/2);
		his_temp2->GetXaxis()->SetRange(peak[i-2]-15*mul,peak[i]+15*mul);
		//his_temp->Draw();

		his_temp2->SetTitle("triple Ta peaks Fitting [Gaus + Exp]");
		his_temp2->SetXTitle("Energy[keV]");
		his_temp2->SetYTitle("Counts/kev");
		//his_temp->SetYTitle("Counts/kev.day");
		his_temp2->SetName("fhis_tot");
		//his_temp->SetName("fit_cday");

		m7fit_peak->Modified();
		m7fit_peak->Update();

/*




cout<<"///////m8 fit///////"<<endl;
cout<<endl;
		c1->cd();
		m8fit->ReleaseParameter(0);
		m8fit->ReleaseParameter(1);
		m8fit->ReleaseParameter(2);
		m8fit->ReleaseParameter(3);
		m8fit->ReleaseParameter(4);
		m8fit->ReleaseParameter(5);
		m8fit->ReleaseParameter(6);
		m8fit->ReleaseParameter(7);
		m8fit->ReleaseParameter(8);
		m8fit->ReleaseParameter(9);
		m8fit->ReleaseParameter(10);


		m8fit->SetRange((peak[i-2]-5),(peak[i]+5));
		m8fit -> SetParameters(gfitpa1[0], peak[i-2], gfitpa3[0], gfitpa1[1], peak[i-1], gfitpa3[1], gfitpa1[2], peak[i], gfitpa3[2], p1fitpa1, p1fitpa2);
		
		m8fit->FixParameter(1,peak[i-2]);
		m8fit->FixParameter(4,peak[i-1]);
		m8fit->FixParameter(7,peak[i]);

		his_temp -> Fit("m8fit", "R0+"); 
		m8fitpa3 = m8fit->GetParameter(2);		
		m8fitpa6 = m8fit->GetParameter(5);
		m8fitpa9 = m8fit->GetParameter(8);

		//sigma check
		if(m8fitpa3<0.4){
			cout <<"Too small sigma(1) value"<<endl;
			cout<<">> sigma(1) fix"<<endl;
			if(peak[i-2]<1000){m8fit->FixParameter(2,0.8);}
			else if(peak[i-2]<2000){m8fit->FixParameter(2,1.1);}
			else{m8fit->FixParameter(2,1.6);}
			}

		if(m8fitpa6<0.4){
			cout <<"Too small sigma(2) value"<<endl;
			cout<<">> sigma(2) fix"<<endl;
			if(peak[i-1]<1000){m8fit->FixParameter(5,0.8);}
			else if(peak[i-1]<2000){m8fit->FixParameter(5,1.1);}
			else{m8fit->FixParameter(5,1.6);}
			}

		if(m8fitpa9<0.4){
			cout <<"Too small sigma(3) value"<<endl;
			cout<<">> sigma(3) fix"<<endl;
			if(peak[i]<1000){m8fit->FixParameter(8,0.8);}
			else if(peak[i]<2000){m8fit->FixParameter(8,1.1);}
			else{m8fit->FixParameter(8,1.6);}
			}

		if(m8fitpa3>1.0&&peak[i-2]<1000){
			cout <<"Too big sigma(1) value @ <1000 keV"<<endl;
			cout<<">> sigma fix"<<endl;
			m8fit->FixParameter(2,0.8);
		}	

		if(m8fitpa6>1.0&&peak[i-1]<1000){
			cout <<"Too big sigma(2) value @ <1000 keV"<<endl;
			cout<<">> sigma fix"<<endl;
			m8fit->FixParameter(5,0.8);
		}	

		if(m8fitpa9>1.0&&peak[i]<1000){
			cout <<"Too big sigma(3) value @ <1000 keV"<<endl;
			cout<<">> sigma fix"<<endl;
			m8fit->FixParameter(8,0.8);
		}	


		if(m8fitpa3>2.0&&peak[i-2]>=1000&&peak[i-2]<2000){
			cout <<"Too big sigma(1) value @ 1000~2000 keV"<<endl;
			cout<<">> sigma(1) fix"<<endl;
			m8fit->FixParameter(2,1.1);
		}

		if(m8fitpa6>2.0&&peak[i-1]>=1000&&peak[i-1]<2000){
			cout <<"Too big sigma value @ 1000~2000 keV"<<endl;
			cout<<">> sigma(2) fix"<<endl;
			m8fit->FixParameter(5,1.1);
		}

		if(m8fitpa9>2.0&&peak[i]>=1000&&peak[i]<2000){
			cout <<"Too big sigma value @ 1000~2000 keV"<<endl;
			cout<<">> sigma(3) fix"<<endl;
			m8fit->FixParameter(8,1.1);
		}


	his_temp -> Fit("m8fit", "R+"); 
		m8fitpa1 = m8fit->GetParameter(0);
		m8fitpa2 = m8fit->GetParameter(1);
		m8fitpa3 = m8fit->GetParameter(2);
		m8fitpa4 = m8fit->GetParameter(3);
		m8fitpa5 = m8fit->GetParameter(4);
		m8fitpa6 = m8fit->GetParameter(5);
		m8fitpa7 = m8fit->GetParameter(6);
		m8fitpa8 = m8fit->GetParameter(7);
		m8fitpa9 = m8fit->GetParameter(8);
		m8fitpa10 = m8fit->GetParameter(9);
		m8fitpa11 = m8fit->GetParameter(10);

		m8fiter1 = m8fit->GetParError(0);	
		m8fiter2 = m8fit->GetParError(1);	
		m8fiter3 = m8fit->GetParError(2);	
		m8fiter4 = m8fit->GetParError(3);	
		m8fiter5 = m8fit->GetParError(4);	
		m8fiter6 = m8fit->GetParError(5);	
		m8fiter7 = m8fit->GetParError(6);	
		m8fiter8 = m8fit->GetParError(7);
		m8fiter9 = m8fit->GetParError(8);	
		m8fiter10 = m8fit->GetParError(9);	
		m8fiter11 = m8fit->GetParError(10);	
	
		m8fitchi2 = (m8fit->GetChisquare())/(m8fit->GetNDF());

		cout << "[m8fit]" << endl;
		cout << "peak1 : " << m8fitpa2 << endl;
		cout << "m8fit tot. counts (peak1) :" << m8pcnt1 << " +/- " << m8pcnter1 << endl;
//		cout << "m6fit cpd (peak1) :" << m6cpd1 << " +/- " << m6dcpd1 << endl;
		cout << "6 sigma : " << m8fitpa3*6 << endl;
		cout << "peak2 : " << m8fitpa5 << endl;
		cout << "m6fit tot. counts (peak2) :" << m8pcnt2 << " +/- " << m8pcnter2 << endl;
//		cout << "m6fit cpd (peak2) :" << m6cpd2 << " +/- " << m6dcpd2 << endl;
		cout << "6 sigma : " << m8fitpa6*6 << endl;
		cout << "chi2 : " << m8fitchi2 << endl;

		m8fitpara->Fill();
		m8fit_peak->cd((i+1)/3);
//		m8fit_peak->cd((i+1)/2);
		his_temp3 -> Fit("m8fit", "R"); 
		his_temp3->GetXaxis()->SetRange(peak[i-2]-15*mul,peak[i]+15*mul);
//		his_temp->Draw();

		his_temp3->SetTitle("double peaks Fitting [Gaus + P1]");
		his_temp3->SetXTitle("Energy[keV]");
		his_temp3->SetYTitle("Counts/kev");
		//his_temp->SetYTitle("Counts/kev.day");
		his_temp3->SetName("fhis_tot");
		//his_temp->SetName("fit_cday");

		m8fit_peak->Modified();
		m8fit_peak->Update();
*/


                cout<<""<<endl;
                cout<<""<<endl;
                cout<<">>chi squre check<<"<<endl;
                cout<<""<<endl;
                cout<<"m6fit chi^2/NDF = "<<m6fitchi2<<endl;
                cout<<"m7fit chi^2/NDF = "<<m7fitchi2<<endl;
                cout<<""<<endl;

                if(abs(1-m6fitchi2) > abs(1-m7fitchi2)){
                cout<<"bestfit : m7 fit "<<endl;
//                his_temp->Fit("m7fit","RQ+");
//              his_full->Fit("m2fit","RQ+");
                fit_type=7;
                fitpa1 = m7fitpa1;
                fitpa2 = m7fitpa2;
                fitpa3 = m7fitpa3;
                fiter1 = m7fiter1;
                fitchi2 = m7fitchi2;

                }
                else if(abs(1-m6fitchi2) < abs(1-m6fitchi2)){
                cout<<"bestfit : m6 fit "<<endl;

                fit_type=6;
                fitpa1 = m6fitpa1;
                fitpa2 = m6fitpa2;
                fitpa3 = m6fitpa3;
                fiter1 = m6fiter1;
                fitchi2 = m6fitchi2;
                fitpara->Fill();
                }

                cout<<""<<endl;



		c1->Modified();
}
	}
		TFile r3f (res3file,"RECREATE");
		r3f.cd();
		gfitpara->Write();
		m6fitpara->Write();
		m7fitpara->Write();
		//m8fitpara->Write();
		fitpara->Write();
		r3f.Close();
	c1->cd();
	his_temp->GetXaxis()->SetRange(0,4000);
	c1->Modified();
	c1->Update();

	m6fit_peak->SaveAs(res2file);
	m7fit_peak->SaveAs(res2file);
//	m8fit_peak->SaveAs(res2file);
	c1->SaveAs(res2file);


}
