#include <iostream>
#include "TMath.h"
#include "TMinuit.h"
#include <string>

void fit_1peak_ch_251023()
{
	string workdir("/data/HPGe/USERS/kkw/DAQ/ANA125");

	char hisfile[256];
	char rawname[256];
	char resfile[256];//fit result

	char runnumber3[256];
	char runnumber6[256];
	char bintag[256];
	char peakfile[258];
	
	int runnum = 650;
//	int runnum = 624;

	int expnum = 2;
	int ich = 1;
	sprintf(runnumber6,"%06d",runnum);
	sprintf(runnumber3,"%03d",runnum);

	int binnum = 16000;// 1bin = 1kev
	int bw = binnum/binnum; //bin width, keV

	int const N = 1; //peaks

	//	double peak[N] = {583.2, 609.3, 768.36, 794.95, 1001.04, 1553.8};
	//	double peak[N] = {609.3, 911, 1120};
	double peak[N] = {1553.77};
	//General resolution


	//exp 2
	double respa[3];

	double respa_ch[14][3]={
	  {0.363128,	0.000667908,	-8.40E-08},
	  {0.370658,	0.000637953,	-7.62E-08},
	  {0.445986,	0.000680498,	-9.62E-08},
	  {1.99356,	0.00228603,	-9.97E-07},
	  {0.436251,	0.000626869,	-7.70E-08},
	  {0.402777,	0.000545669,	-4.38E-08},
	  {0.582696,	0.000821223,	-1.54E-07},
	  {0.495129,	0.000449309,	8.03E-09},
	  {0.511883,	0.000459608,	-2.06E-08},
	  {0 ,0 ,0},
	  {0.482192,	0.000504369,	-2.71E-08},
	  {0.432383,	0.000486415,	-3.44E-08},
	  {0.429384,	0.000619621,	-8.22E-08},
	  {0.448204, 	0.000548424,	-4.75E-08} 	
	};

	respa[0]=respa_ch[ich-1][0];
	respa[1]=respa_ch[ich-1][1];
	respa[2]=respa_ch[ich-1][2];
	cout<<"respa "<<respa[0]<<" / "<<respa[1]<<" / "<<respa[2]<<endl;

	TF1 * gfit = new TF1("gfit","gaus");
	TF1 * agfit = new TF1("agfit","[0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2])");
	agfit->SetParNames("Area","Mean","Sigma"); //area = sqrt(2pi)*constant*sigma

	TF1 * p0fit = new TF1("p0fit","pol0");
	TF1 * p1fit = new TF1("p1fit","pol1");

	TF1 * fit = new TF1("fit", "([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+[3]+[4]*x");
	fit->SetParNames("Area", "Mean", "Sigma", "Intercept", "Slope");



	TCanvas * c1 = new TCanvas("c1","fitting",1200,800);

	fit->SetLineColor(kGreen);
	//	p0fit->SetLineColor(kBlue);
	p1fit->SetLineColor(kGreen);
	//	p0fit->SetLineStyle(7);
	p1fit->SetLineStyle(7);

	int maxb, maxh, run;
	double maxx;
	double bg1, bg2, bg;
	double pa1[N][3];//gfit 
	double cal[N][2]; //integral fit func - m1fit,m2fit

	int xbin[N];
	double energy;
	double gfitpa[3], gfiter[3], gfitchi2;
	double agfitpa[3], agfiter[3], agfitchi2;
	double fitpa[5], fiter[5], fitmean, fitsig, fitchi2;

	//	double fitpa[3], fiter[3], fitchi2;
	double fixbg[2];//bg line pa. for m1, m2, m3

	TTree *gfitpara = new TTree("gfitpara","gaus fitting result");
	gfitpara->Branch("run",&runnum,"run/I");
	gfitpara ->Branch("peak",&energy,"energy/D");
	gfitpara ->Branch("const",&gfitpa[0],"gfitpa0/D");
	gfitpara ->Branch("mean",&gfitpa[1],"gfitpa1/D");
	gfitpara ->Branch("sigma",&gfitpa[2],"gfitpa2/D");
	gfitpara ->Branch("chi2",&gfitchi2,"gfitchi2/D");
	gfitpara ->Branch("const_er",&gfiter[0],"gfiter0/D");
	gfitpara ->Branch("mean_er",&gfiter[1],"gfiter1/D");
	gfitpara ->Branch("sigma_er",&gfiter[2],"gfiter2/D");

	TTree *agfitpara = new TTree("agfitpara","area-gaus fitting result");
	agfitpara->Branch("run",&runnum,"run/I");
	agfitpara ->Branch("peak",&energy,"energy/D");
	agfitpara ->Branch("area",&agfitpa[0],"gfitpa0/D");
	agfitpara ->Branch("mean",&agfitpa[1],"gfitpa1/D");
	agfitpara ->Branch("sigma",&agfitpa[2],"gfitpa2/D");
	agfitpara ->Branch("chi2",&agfitchi2,"gfitchi2/D");
	agfitpara ->Branch("area_er",&agfiter[0],"gfiter0/D");
	agfitpara ->Branch("mean_er",&agfiter[1],"gfiter1/D");
	agfitpara ->Branch("sigma_er",&agfiter[2],"gfiter2/D");


	TTree *fitpara = new TTree("fitpara", "gaus+pol1 fitting result");
	fitpara -> Branch("run", &runnum, "run/I");
	fitpara -> Branch("peak", &energy, "energy/D");
	fitpara -> Branch("area", &fitpa[0], "fitpa0/D");
	fitpara -> Branch("mean", &fitpa[1], "fitpa1/D");
	fitpara -> Branch("sigma", &fitpa[2], "fitpa2/D");
	fitpara -> Branch("pol0", &fitpa[3], "fitpa3/D");
	fitpara -> Branch("pol1", &fitpa[4], "fitpa4/D");
	fitpara -> Branch("chi2", &fitchi2, "fitchi2/D");
	fitpara -> Branch("area_er", &fiter[0], "fiter0/D");
	fitpara -> Branch("mean_er", &fiter[1], "fiter1/D");
	fitpara -> Branch("sigma_er", &fiter[2], "fiter2/D");
	fitpara -> Branch("pol0_er", &fiter[3], "fitper3/D");
	fitpara -> Branch("pol1_er", &fiter[4], "fitper4/D");


	TF1 *m1bg[N];
	TF1 *m2bg[N];
	char m1bgname[256];
	char m2bgname[256];


	if(expnum==1){
	//run624(run1) rrs on
	  //	  sprintf(hisfile,"%s/result/V50/Bin%i/his_%s.v0.rrs_on.root", workdir.c_str(), binnum, runnumber6);
	  //	  sprintf(resfile,"%s/result/V50/Bin%i/res_%s.v0.rrs_on.root", workdir.c_str(), binnum, runnumber6);
	  sprintf(hisfile,"%s/result/V50/Bin%i/his_run1.root", workdir.c_str(), binnum);
	  sprintf(resfile,"%s/result/V50/Bin%i/res_run1.root", workdir.c_str(), binnum);



	}
	if(expnum==2){
	//run650 run2
	  //	  sprintf(hisfile,"%s/result/V50/Bin%i/his_000650.v1cut.root", workdir.c_str(), binnum);
	  //	  sprintf(resfile,"%s/result/V50/Bin%i/res_000650.v1cut.root", workdir.c_str(), binnum);

	  sprintf(hisfile,"%s/result/V50/Bin%i/his_000650.v1.root", workdir.c_str(), binnum);
	  sprintf(resfile,"%s/result/V50/Bin%i/res_000650.v1.ch%i.root", workdir.c_str(), binnum, ich);

	}


	TFile *hf=new TFile(hisfile);
	TH1D * his_temp = new TH1D("his_temp","",binnum,0,4000);

//	his_temp = (TH1D*)hf->Get("his_tot1");//13det
//	his_temp = (TH1D*)hf->Get("his_tot3");//12det(4,10)
	char hisname[245];
	sprintf(hisname,"his%i",ich);
      	his_temp = (TH1D*)hf->Get(hisname);//11det
	his_temp->Draw();
	his_temp->SetName("his_temp");
	his_temp->SetLineColor(1);

	double sigma_cal;

	for(int i = 0;i<N;i++){

		c1->cd();
		energy = peak[i];
		//		xbin[i]=his_temp->GetXaxis()->FindBin(peak[i]);
		xbin[i]=his_temp->FindBin(peak[i]);

		printf("\n");
		printf("\n=========================\n");
		printf("%.2f keV bin : %i	\n",peak[i],xbin[i]);

		//		his_temp->GetXaxis()->SetRange(xbin[i]-1*mul,xbin[i]+1*mul);
		//		sigma_cal = peak[i]*(respa[0] + (respa[1]/peak[i]) + (respa[2]/(sqrt(peak[i]))));		
		sigma_cal = sqrt(respa[0]*respa[0]+respa[1]*peak[i]+respa[2]*peak[i]*peak[i]);

		cout<<"Cal. Sigma : "<<sigma_cal<<" (3 sigma : "<<3*sigma_cal<<endl;
		cout<<"3 sigma range : "<<his_temp->FindBin(peak[i]-(3*sigma_cal))<<" ~ "<<his_temp->FindBin(peak[i]+(3*sigma_cal))<<endl;
		//		his_temp->GetXaxis()->SetRange(his_temp->Getxbin[i]-1*mul,xbin[i]+1*mul);
		his_temp->GetXaxis()->SetRange(his_temp->FindBin(peak[i]-(3*sigma_cal)),his_temp->FindBin(peak[i]+(3*sigma_cal)));


		maxb = his_temp->GetMaximumBin();
		maxx = (his_temp->GetBinCenter(maxb));
		cout<<"maxb(bin) : "<<maxb<<endl;
		cout<<"maxx(kev) : "<<maxx<<endl;
		maxh = his_temp->GetMaximum();

		his_temp->GetXaxis()->SetRange(0,binnum);
		gfit->SetRange((maxx-3*sigma_cal),(maxx+3*sigma_cal));
		agfit->SetRange((maxx-3*sigma_cal),(maxx+3*sigma_cal));
		fit->SetRange((maxx-6*sigma_cal),(maxx+6*sigma_cal));
//		fit->SetRange((maxx-24*sigma_cal),(maxx+24*sigma_cal));

//		p0fit->SetRange((maxx-30),(maxx+30));
		p1fit->SetRange((maxx-30),(maxx+30));

		//		his_temp->Fit("p0fit","RQ0+");
		//		fixbg1 = p0fit->GetParameter(0);

		his_temp->Fit("p1fit","RQ+");
		fixbg[0] = p1fit->GetParameter(0);
		fixbg[1] = p1fit->GetParameter(1);

		//gaus fit
		his_temp->Fit("gfit","RQ0+");
		gfitpa[0] = gfit->GetParameter(0);
		gfitpa[1] = gfit->GetParameter(1);
		gfitpa[2] = gfit->GetParameter(2);
		gfiter[0] = gfit->GetParError(0);
		gfiter[1] = gfit->GetParError(1);
		gfiter[2] = gfit->GetParError(2);

		gfitchi2 = (gfit->GetChisquare())/(gfit->GetNDF());


		//area gaus fit

		//agfit->SetParameters(2.5*gfitpa[1]*gfitpa[2],gfitpa[1],gfitpa[2]);
		agfit->SetParameters(2.5*gfitpa[1]*gfitpa[2],peak[i],sigma_cal);
		//agfit->SetParameters(2.5*gfitpa2*gfitpa3,maxx,2);
		his_temp->Fit("agfit","R0Q+");
		//		his_temp->Fit("agfit","R+");
		agfitpa[0] = agfit->GetParameter(0);
		agfitpa[1] = agfit->GetParameter(1);
		agfitpa[2] = agfit->GetParameter(2);
		agfiter[0] = agfit->GetParError(0);
		agfiter[1] = agfit->GetParError(1);
		agfiter[2] = agfit->GetParError(2);
		agfitchi2 = (agfit->GetChisquare())/(agfit->GetNDF());

		cout<<""<<endl;
		cout<<"[fin fitting]"<<endl;
		cout<<""<<endl;
		//fit

		fit->ReleaseParameter(0);
		fit->ReleaseParameter(1);
		fit->ReleaseParameter(2);
		fit->ReleaseParameter(3);
		fit->ReleaseParameter(4);
		cout<<"maxx"<<maxx<<endl;
		//		fit->SetParameters(agfitpa[0],maxx,sigma_cal,fixbg[0],fixbg[1]);

		fit->SetParameter(0,agfitpa[0]);
		fit->SetParameter(1,maxx);
		fit->SetParLimits(1,maxx-3*sigma_cal, maxx+3*sigma_cal);
		fit->SetParameter(2,sigma_cal);
		fit->SetParLimits(2,0.5*sigma_cal, 1.5*sigma_cal);
		fit->SetParameter(3,fixbg[0]);
		fit->SetParameter(4,fixbg[1]);

//		fit->FixParameter(3,fixbg[0]);
//		fit->FixParameter(4,fixbg[1]);

		/*
		//peak check
		his_temp->Fit("fit","R+");
		fitpa[1] = fit->GetParameter(1);
		cout<<"fitpa[1]-peak[i] "<<fitpa[1]-peak[i]<<endl;
		cout<<"3*sigma_cal "<<3*sigma_cal<<endl;
		
		if(abs(fitpa[1]-peak[i])>3*(sigma_cal)){
			cout <<"test fitting mean "<<fitpa[1]<<" for peak "<<peak[i]<<endl;
			cout<<">> peak fixed at"<<peak[i]<<endl;
			fit->FixParameter(1,maxx);
		}
		*/
		/*
		//sigma check
		his_temp->Fit("fit","R0+");
		
		fitpa[2] = fit->GetParameter(2);
		if((fitpa[2]<(0.8*sigma_cal))||(fitpa[2]>(3*sigma_cal))){
			cout <<"m3 test fitting sigma "<<fitpa[2]<<" for peak "<<peak[i]<<endl;
			cout<<">> sigma fixed at "<<sigma_cal<<endl;
			fit->FixParameter(2,sigma_cal);

		}
		*/
		/*
		//peak re-check
		his_temp->Fit("fit","R0+");
		fitpa[1] = fit->GetParameter(1);
		if(abs(fitpa[1]-peak[i])>0.3*(sigma_cal)){
			cout <<"m3 test fitting mean "<<fitpa[1]<<" for peak "<<peak[i]<<endl;
			cout<<">> peak fixed at"<<peak[i]<<endl;
			fit->FixParameter(1,maxx);
			
		}
		*/
			his_temp->Fit("fit","R+");

		fitpa[0] = fit->GetParameter(0);
		fitpa[1] = fit->GetParameter(1);
		fitpa[2] = fit->GetParameter(2);
		fitpa[3] = fit->GetParameter(3);
		fitpa[4] = fit->GetParameter(4);
		fiter[0] = fit->GetParError(0);
		fiter[1] = fit->GetParError(1);
		fiter[2] = fit->GetParError(2);
		fiter[3] = fit->GetParError(3);
		fiter[4] = fit->GetParError(4);
		fitchi2 = (fit->GetChisquare())/(fit->GetNDF());

		fitpara->Fill();
		//find best fit

		cout<<""<<endl;
		cout<<"fit chi^2/NDF = "<<fitchi2<<endl;
		cout<<""<<endl;
		/*
		if((abs(1.-m1fitchi2)<abs(1.-m2fitchi2)) && (abs(1.-m1fitchi2)<abs(1.-fitchi2))){
			cout<<"bestfit : m1 fit "<<endl;
			fit_type = 1;
			fitpa[0] = m1fitpa[0];
			fitpa[1] = m1fitpa[1];
			fitpa[2] = m1fitpa[2];
			fiter[0] = m1fiter[0];
			fiter[1] = m1fiter[1];
			fiter[2] = m1fiter[2];
			fitchi2 = m1fitchi2;
//			c2->cd();
		}
		else if((abs(1.-m2fitchi2)<abs(1.-m1fitchi2)) && (abs(1.-m2fitchi2)<abs(1.-fitchi2))){
			cout<<"bestfit : m2 fit "<<endl;
			fit_type = 2;
			fitpa[0] = m2fitpa[0];
			fitpa[1] = m2fitpa[1];
			fitpa[2] = m2fitpa[2];
			fiter[0] = m2fiter[0];
			fiter[1] = m2fiter[1];
			fiter[2] = m2fiter[2];
			fitchi2 = m2fitchi2;
//			c2->cd();
		}
		else if((abs(1.-fitchi2)<abs(1.-m1fitchi2)) && (abs(1.-fitchi2)<abs(1.-m2fitchi2))){
			cout<<"bestfit : m3 fit "<<endl;
			fit_type = 3;
			fitpa[0] = fitpa[0];
			fitpa[1] = fitpa[1];
			fitpa[2] = fitpa[2];
			fiter[0] = fiter[0];
			fiter[1] = fiter[1];
			fiter[2] = fiter[2];
			fitchi2 = fitchi2;
//			c2->cd();
		}
		*/

		cout<<fitpa[0]<<endl;
		cout<<fitpa[1]<<endl;
		cout<<fitpa[2]<<endl;
		//		fitpara->Fill();
		cout<<""<<endl;
		cout<<""<<endl;

	}

	char histitle[256];
	sprintf(histitle,"source peaks Fitting");
	his_temp->SetTitle(histitle);
	his_temp->SetXTitle("Energy[keV]");
	his_temp->SetYTitle("Counts/kev");

	c1->cd();

	his_temp->GetXaxis()->SetRangeUser(100, 3000);

	his_temp->Draw();
	TFile rf (resfile,"RECREATE");
	rf.cd();
	gfitpara->Write();
	fitpara->Write();
	c1->Write();

	rf.Close();
}
