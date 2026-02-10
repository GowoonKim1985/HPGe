#include <iostream>
#include "TMath.h"
#include "TMinuit.h"

void fit_double()
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


        int binnum;
        cout<<"bin number : ";
        scanf("%i",&binnum);

	char peakfile[256];
	sprintf(peakfile,"/home/kkw/DAQ/ANA500/script/His_Fit_CODE/peak_R%03d/peak_double.dat",runnum);

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

	agfit->SetLineColor(9);
	gfit->SetLineColor(3);
	efit->SetLineColor(5);
	e2fit->SetLineColor(6);
	p1fit->SetLineColor(2);
	//m3fit->SetLineColor(kBlue);
	m3fit->SetLineColor(kRed);
	m4fit->SetLineColor(kRed);
	m5fit->SetLineColor(kGreen);

	p0fit->SetLineStyle(7);
	efit->SetLineStyle(7);

	TCanvas * c1 = new TCanvas("c1","fitting",1200,800);
//	TCanvas * m3fit_peak = new TCanvas("Gaus+P0 fitting","Gaus+P0 fitting",1200,800);
//	TCanvas * m4fit_peak = new TCanvas("Gaus+Exp fitting","Gaus+Exp fitting",1200,800);
//	TCanvas * m5fit_peak = new TCanvas("Gaus+P1 fitting","Gaus+P1 fitting",1200,800);
//	TCanvas * div_peak1 = new TCanvas("fitting detail 1","fitting detail 1",1200,800);

//	m3fit_peak->Divide(4,3);	
//	m4fit_peak->Divide(4,3);
//	m5fit_peak->Divide(4,3);



	int maxb, maxh;
	double maxx[2];
	int maxb1, maxb2;
	double bg1, bg2, bg;
	double cnt1, cnt2, cnt3, cnt1_2, cnt3_2;
	double pa1[line_max1-1][3];//gfit 

	double cal[line_max1-1][2]; //integral fit func - m1fit,m2fit
	double epa[6]; //expo fit 
	int xbin[line_max1-1];
	double energy[2];

	double agfitpa1[2], agfitpa2[2], agfitpa3[2], agfiter1[2], agfiter2[2], agfiter3[2];

	double gfitpa1[2], gfitpa2[2], gfitpa3[2], gfiter1[2], gfiter2[2], gfiter3[2];

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

	double gfitchi2[2], agfitchi2[2], p1fitchi2, m3fitchi2, m4fitchi2, m5fitchi2;
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


	TF1 *m1bg[line_max1-1];
	TF1 *m2bg[line_max1-1];
	char m1bgname[256];
	char m2bgname[256];

	int bin = binnum;// 1bin = 1kev
	//	int bin = 40000;// 1bin = 100ev
	int mul = bin/4000;
	int odd_check = 1;

	sprintf(hisfile,"/home/kkw/DAQ/ANA500/result/RUN%s/Bin%i/his_%s_M1.root",runnumber3,binnum,runnumber6);
	sprintf(res2file,"/home/kkw/DAQ/ANA500/result/RUN%s/Bin%i/fit_double_%s.C",runnumber3,binnum,runnumber6);
	sprintf(res3file,"/home/kkw/DAQ/ANA500/result/RUN%s/Bin%i/fit_double_%s.root",runnumber3,binnum,runnumber6);

	TFile *hf=new TFile(hisfile);
	TH1D * his_temp = new TH1D("his_temp","",bin,0,4000);
//	TH1D * his_temp1 = new TH1D("his_temp1","",bin,0,4000);
//	TH1D * his_temp2 = new TH1D("his_temp2","",bin,0,4000);
//	TH1D * his_temp3 = new TH1D("his_temp3","",bin,0,4000);
	his_temp = (TH1D*)hf->Get("his_tot");
//	his_temp1 = (TH1D*)hf->Get("his_tot");
//	his_temp2 = (TH1D*)hf->Get("his_tot");
//	his_temp3 = (TH1D*)hf->Get("his_tot");
	//		his_temp = (TH1D*)hf->Get("fhis_cday");
	his_temp->SetLineColor(1);

        double respa1 = 0.000158;
        double respa2 = 0.547238;
        double respa3 = 0.006808;

	double sigma_cal[2];

	for(int i = 0;i<(line_max1-1);i++){
		if(i==odd_check){	
		odd_check = odd_check + 2;	

		energy[0] = peak[i-1];
		energy[1] = peak[i];

		sigma_cal[0] = peak[i-1]*(respa1 + (respa2/peak[i-1]) + (respa3/(sqrt(peak[i-1]))));
		sigma_cal[1] = peak[i]*(respa1 + (respa2/peak[i]) + (respa3/(sqrt(peak[i]))));


		xbin[i-1]=his_temp->GetXaxis()->FindBin(peak[i-1]);
		xbin[i]=his_temp->GetXaxis()->FindBin(peak[i]);

		his_temp->GetXaxis()->SetRange((peak[i-1]-2),(peak[i-1]+2));
		maxb1 = his_temp->GetMaximumBin();
		cout<< "max peak1 : "<<maxb1<<endl;
		his_temp->GetXaxis()->SetRange((peak[i]-2),(peak[i]+2));
		maxb2 = his_temp->GetMaximumBin();
		cout<< "max peak2 : "<<maxb2<<endl;

		maxx[i-1]=maxb1;
		maxx[i]=maxb2;

		printf("\n");
		printf("\n=========================\n");
		printf("%.2f keV bin : %i	\n",peak[i-1],xbin[i-1]);
		printf("%.2f keV bin : %i	\n",peak[i],xbin[i]);

		c1->cd();
/*		his_temp->GetXaxis()->SetRange(xbin[i-1]-5*mul,xbin[i-1]+5*mul);

		maxb = his_temp->GetMaximumBin();
		maxx = maxb/mul;
		cout<<"maxb(bin) : "<<maxb<<endl;
		cout<<"maxx(kev) : "<<maxx<<endl;
		maxh = his_temp->GetMaximum();
*/
		his_temp->GetXaxis()->SetRange(0,4000*mul);

	
		efit->SetRange((peak[i-1]-30),(peak[i]+30));
		p1fit->SetRange((peak[i-1]-15),(peak[i]+15));

		his_temp->Fit("efit","R0+");
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
		gfit->SetRange((peak[i-1]-2),(peak[i-1]+2));
		gfit->FixParameter(1,peak[i-1]);
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
		gfit->SetRange((peak[i]-2),(peak[i]+2));
		gfit->FixParameter(1,peak[i]);
		his_temp->Fit("gfit","R0+");
		gfitpa1[1] = gfit->GetParameter(0);
		gfitpa2[1] = gfit->GetParameter(1);
		gfitpa3[1] = gfit->GetParameter(2);
		gfiter1[1] = gfit->GetParError(0);
		gfiter2[1] = gfit->GetParError(1);
		gfiter3[1] = gfit->GetParError(2);
		gfitchi2[1] = (gfit->GetChisquare())/(gfit->GetNDF());

		bg1 = 0; bg2 = 0; 
		for(int j1=0;j1<(5*mul);j1++){
			bg1 = bg1 + (his_temp->GetBinContent(peak[i-1]-(gfitpa3[0]*3+5)*mul+j1));
			bg2 = bg2 + (his_temp->GetBinContent(peak[i]+(gfitpa3[1]*3+5)*mul-j1));
		}

		bg=(bg1+bg2)/(10*mul);
		fixpa[0]=bg;


		//area gaus fit
			//1st peak
		cout<<""<<endl;
		cout<<"[area gaus 1st peak fitting]"<<endl;
		cout<<""<<endl;

		agfit->ReleaseParameter(1);
		agfit->SetRange((peak[i-1]-2),(peak[i-1]+2));
		agfit->SetParameters(2.5*gfitpa2[0]*gfitpa3[0],gfitpa2[0],gfitpa3[0]);
		gfit->FixParameter(1,peak[i-1]);

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
		agfit->SetRange((peak[i]-2),(peak[i]+2));
		agfit->SetParameters(2.5*gfitpa2[1]*gfitpa3[1],gfitpa2[1],gfitpa3[1]);
		agfit->FixParameter(1,peak[i]);
		his_temp->Fit("agfit","RQ0+");
		agfitpa1[1] = agfit->GetParameter(0);
		agfitpa2[1] = agfit->GetParameter(1);
		agfitpa3[1] = agfit->GetParameter(2);
		agfiter1[1] = agfit->GetParError(0);
		agfiter2[1] = agfit->GetParError(1);
		agfiter3[1] = agfit->GetParError(2);
		agfitchi2[1] = (agfit->GetChisquare())/(agfit->GetNDF());


/*
		cout<<""<<endl;
		cout<<"[m1 fitting]"<<endl;
		cout<<""<<endl;
*/

cout<<"///////m3 fit///////"<<endl;
cout<<endl;

		//m3fit_peak->cd((i+1)/2);
//	c1->cd();

		m3fit->ReleaseParameter(0);
		m3fit->ReleaseParameter(1);
		m3fit->ReleaseParameter(2);
		m3fit->ReleaseParameter(3);
		m3fit->ReleaseParameter(4);
		m3fit->ReleaseParameter(5);
		m3fit->ReleaseParameter(6);



		m3fit->SetRange((peak[i-1]-5),(peak[i]+5));
//		m1fit->SetParameters(agfitpa1,agfitpa2,agfitpa3,bg);


		//m3fit
		//m3fit -> SetParameters(gfitpa1, peak[i], gfitpa3, gfitpa1, peak[i]+10, gfitpa3, bg);
		//m3fit -> SetParameters(gpa1[i-1], peak[i-1], gpa2[i-1], gpa1[i], peak[i], gpa2[i], bg);
		//m3fit -> SetParameters(gpa1[i-1], peak[i-1], gpa2[i-1], gpa1[i], peak[i], gpa2[i], bg);

		m3fit -> SetParameters(agfitpa1[0], peak[i-1], agfitpa3[0], agfitpa1[1], peak[i], agfitpa3[1], bg);

//		m3fit->FixParameter(1,maxx[i-1]);
//		m3fit->FixParameter(4,maxx[i]);
		m3fit->FixParameter(1,peak[i-1]);
		m3fit->FixParameter(4,peak[i]);
		if(peak[i]==338.53){m3fit->FixParameter(4,338.2);}
/*
		m3fit -> SetParameters(agfitpa1[0], maxb1, agfitpa3[0], agfitpa1[1], maxb2, agfitpa3[0], bg);
		m3fit->FixParameter(1,maxb1);
		m3fit->FixParameter(4,maxb2);
*/
		his_temp -> Fit("m3fit", "R0+");
		m3fitpa3 = m3fit->GetParameter(2);		
		m3fitpa6 = m3fit->GetParameter(5);

		//sigma check
		if(m3fitpa3<0.4){
			cout <<"Too small sigma(1) value"<<endl;
			cout<<">> sigma(1) fix"<<endl;
			m3fit->FixParameter(2, sigma_cal[i-1]);
//			if(peak[i-1]<1000){m3fit->FixParameter(2,0.8);}
//			else if(peak[i-1]<2000){m3fit->FixParameter(2,1.1);}
//			else{m3fit->FixParameter(2,1.6);}
			}

		if(m3fitpa6<0.4){
			cout <<"Too small sigma(2) value"<<endl;
			cout<<">> sigma(2) fix"<<endl;
			//if(peak[i]<1000){m3fit->FixParameter(5,0.8);}
			//else if(peak[i]<2000){m3fit->FixParameter(5,1.1);}
			//else{m3fit->FixParameter(5,1.6);}
			m3fit->FixParameter(5, sigma_cal[i]);
			}

		if(m3fitpa3>(sigma_cal[0]*1.64)){
			cout <<"Too big sigma(1) value @ <1000 keV"<<endl;
			cout<<">> sigma fix"<<endl;
//			m3fit->FixParameter(2,0.8);
			m3fit->FixParameter(2, sigma_cal[0]);
		}	

		if(m3fitpa6>(sigma_cal[1]*1.64)){
			cout <<"Too big sigma(2) value @ <1000 keV"<<endl;
			cout<<">> sigma fix"<<endl;
//			m3fit->FixParameter(5,0.8);
			m3fit->FixParameter(5, sigma_cal[1]);
		}	



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
/*
		m3fit_peak->cd((i+1)/2);
		his_temp1->GetXaxis()->SetRange(peak[i-1]-5,peak[i]+5);
		his_temp1 -> Fit("m3fit", "RQ+"); 
//		his_temp1->GetXaxis()->SetRange(peak[i-1]-15*mul,peak[i]+15*mul);
		//his_temp->Draw();

		his_temp1->SetTitle("double peaks Fitting [Gaus + Exp]");
		his_temp1->SetXTitle("Energy[keV]");
		his_temp1->SetYTitle("Counts/kev");
		//his_temp->SetYTitle("Counts/kev.day");
		his_temp1->SetName("fhis_tot");
		//his_temp->SetName("fit_cday");

		m3fit_peak->Modified();
		m3fit_peak->Update();
*/



//////////m4 fit


cout<<"///////m4 fit///////"<<endl;
cout<<endl;

//	c1->cd();
		m4fit->ReleaseParameter(0);
		m4fit->ReleaseParameter(1);
		m4fit->ReleaseParameter(2);
		m4fit->ReleaseParameter(3);
		m4fit->ReleaseParameter(4);
		m4fit->ReleaseParameter(5);
		m4fit->ReleaseParameter(6);
		m4fit->ReleaseParameter(7);



		m4fit->SetRange((peak[i-1]-5),(peak[i]+5));
		m4fit -> SetParameters(gfitpa1[0], peak[i-1], gfitpa3[0], gfitpa1[1], peak[i], gfitpa3[1], fixpa[1], fixpa[2]);
		
		//m4fit->FixParameter(1,peak[i-1]);
		//m4fit->FixParameter(4,peak[i]);
		m4fit->FixParameter(1,maxx[i-1]);
		m4fit->FixParameter(4,maxx[i]);

		his_temp -> Fit("m4fit", "R0+"); 
		m4fitpa3 = m4fit->GetParameter(2);		
		m4fitpa6 = m4fit->GetParameter(5);

		//sigma check

		//sigma check
		if(m4fitpa3<0.4){
			cout <<"Too small sigma(1) value"<<endl;
			cout<<">> sigma(1) fix"<<endl;
			m4fit->FixParameter(2, sigma_cal[i-1]);
//			if(peak[i-1]<1000){m3fit->FixParameter(2,0.8);}
//			else if(peak[i-1]<2000){m3fit->FixParameter(2,1.1);}
//			else{m3fit->FixParameter(2,1.6);}
			}

		if(m4fitpa6<0.4){
			cout <<"Too small sigma(2) value"<<endl;
			cout<<">> sigma(2) fix"<<endl;
			//if(peak[i]<1000){m3fit->FixParameter(5,0.8);}
			//else if(peak[i]<2000){m3fit->FixParameter(5,1.1);}
			//else{m3fit->FixParameter(5,1.6);}
			m4fit->FixParameter(5, sigma_cal[i]);
			}

		if(m4fitpa3>(sigma_cal[0]*1.64)){
			cout <<"Too big sigma(1) value"<<endl;
			cout<<">> sigma fix"<<endl;
//			m3fit->FixParameter(2,0.8);
			m4fit->FixParameter(2, sigma_cal[0]);
		}	

		if(m4fitpa6>(sigma_cal[1]*1.64)){
			cout <<"Too big sigma(2) value"<<endl;
			cout<<">> sigma fix"<<endl;
//			m3fit->FixParameter(5,0.8);
			m4fit->FixParameter(5, sigma_cal[1]);
		}	




		his_temp -> Fit("m4fit", "R0+"); 

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

		cout << "[m4fit]" << endl;
		cout << "peak1 : " << m4fitpa2 << endl;
		cout << "m4fit tot. counts (peak1) :" << m4pcnt1 << " +/- " << m4pcnter1 << endl;
//		cout << "m3fit cpd (peak1) :" << m3cpd1 << " +/- " << m3dcpd1 << endl;
		cout << "6 sigma : " << m4fitpa3*6 << endl;
		cout << "peak2 : " << m4fitpa5 << endl;
		cout << "m3fit tot. counts (peak2) :" << m4pcnt2 << " +/- " << m4pcnter2 << endl;
//		cout << "m3fit cpd (peak2) :" << m3cpd2 << " +/- " << m3dcpd2 << endl;
		cout << "6 sigma : " << m4fitpa6*6 << endl;
		cout << "chi2 : " << m4fitchi2 << endl;

		m4fitpara->Fill();
/*
		//m4fit_peak->cd((i+1)/2);
		m4fit_peak->cd((i+1)/2);
		his_temp2->GetXaxis()->SetRange(peak[i-1]-5,peak[i]+5);
		his_temp2 -> Fit("m4fit", "RQ");  
//		his_temp2->GetXaxis()->SetRange(peak[i-1]-15*mul,peak[i]+15*mul);
		//his_temp->Draw();

		his_temp2->SetTitle("double peaks Fitting [Gaus + Exp]");
		his_temp2->SetXTitle("Energy[keV]");
		his_temp2->SetYTitle("Counts/kev");
		//his_temp->SetYTitle("Counts/kev.day");
		his_temp2->SetName("fhis_tot");
		//his_temp->SetName("fit_cday");

		m4fit_peak->Modified();
		m4fit_peak->Update();



*/
/*

cout<<"///////m5 fit///////"<<endl;
cout<<endl;
		//m5fit_peak->cd((i+1)/2);
//	c1->cd();
		m5fit->ReleaseParameter(0);
		m5fit->ReleaseParameter(1);
		m5fit->ReleaseParameter(2);
		m5fit->ReleaseParameter(3);
		m5fit->ReleaseParameter(4);
		m5fit->ReleaseParameter(5);
		m5fit->ReleaseParameter(6);
		m5fit->ReleaseParameter(7);



		m5fit->SetRange((peak[i-1]-5),(peak[i]+5));
		m5fit -> SetParameters(gfitpa1[0], peak[i-1], gfitpa3[0], gfitpa1[1], peak[i], gfitpa3[1], p1fitpa1, p1fitpa2);
		
		m5fit->FixParameter(1,peak[i-1]);
		m5fit->FixParameter(4,peak[i]);

		his_temp -> Fit("m5fit", "R+"); 
		m5fitpa3 = m5fit->GetParameter(2);		
		m5fitpa6 = m5fit->GetParameter(5);


		//sigma check
		if(m5fitpa3<0.4){
			cout <<"Too small sigma(1) value"<<endl;
			cout<<">> sigma(1) fix"<<endl;
			m5fit->FixParameter(2, sigma_cal[i-1]);
//			if(peak[i-1]<1000){m3fit->FixParameter(2,0.8);}
//			else if(peak[i-1]<2000){m3fit->FixParameter(2,1.1);}
//			else{m3fit->FixParameter(2,1.6);}
			}

		if(m5fitpa6<0.4){
			cout <<"Too small sigma(2) value"<<endl;
			cout<<">> sigma(2) fix"<<endl;
			//if(peak[i]<1000){m3fit->FixParameter(5,0.8);}
			//else if(peak[i]<2000){m3fit->FixParameter(5,1.1);}
			//else{m3fit->FixParameter(5,1.6);}
			m5fit->FixParameter(5, sigma_cal[i]);
			}

		if(m5fitpa3>(sigma_cal[i-1]*1.64)){
			cout <<"Too big sigma(1) value"<<endl;
			cout<<">> sigma fix"<<endl;
//			m3fit->FixParameter(2,0.8);
			m5fit->FixParameter(2, sigma_cal[i-1]);
		}	

		if(m5fitpa6>(sigma_cal[i]*1.64)){
			cout <<"Too big sigma(2) value"<<endl;
			cout<<">> sigma fix"<<endl;
//			m3fit->FixParameter(5,0.8);
			m5fit->FixParameter(5, sigma_cal[i]);
		}	



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

		cout << "[m5fit]" << endl;
		cout << "peak1 : " << m5fitpa2 << endl;
		cout << "m5fit tot. counts (peak1) :" << m5pcnt1 << " +/- " << m5pcnter1 << endl;
//		cout << "m3fit cpd (peak1) :" << m3cpd1 << " +/- " << m3dcpd1 << endl;
		cout << "6 sigma : " << m5fitpa3*6 << endl;
		cout << "peak2 : " << m5fitpa5 << endl;
		cout << "m3fit tot. counts (peak2) :" << m5pcnt2 << " +/- " << m5pcnter2 << endl;
//		cout << "m3fit cpd (peak2) :" << m3cpd2 << " +/- " << m3dcpd2 << endl;
		cout << "6 sigma : " << m5fitpa6*6 << endl;
		cout << "chi2 : " << m5fitchi2 << endl;

		m5fitpara->Fill();

*/
/*
		m5fit_peak->cd((i+1)/2);
		his_temp3->GetXaxis()->SetRange(peak[i-1]-5,peak[i]+5);
		his_temp3 -> Fit("m5fit", "RQ"); 
//		his_temp3->GetXaxis()->SetRange(peak[i-1]-15*mul,peak[i]+15*mul);
//		m5fit_peak->cd((i+1)/2);

//		his_temp->Draw();

		his_temp3->SetTitle("double peaks Fitting [Gaus + P1]");
		his_temp3->SetXTitle("Energy[keV]");
		his_temp3->SetYTitle("Counts/kev");
		//his_temp->SetYTitle("Counts/kev.day");
		his_temp3->SetName("fhis_tot");
		//his_temp->SetName("fit_cday");

		m5fit_peak->Modified();
		m5fit_peak->Update();
*/


                //

	}}

		TFile r3f (res3file,"RECREATE");
		r3f.cd();
		gfitpara->Write();
		m3fitpara->Write();
		m4fitpara->Write();
		m5fitpara->Write();
	//	intpara->Write();
		r3f.Close();




}
