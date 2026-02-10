#include <iostream>
#include "TMath.h"
#include "TMinuit.h"
#include <string>

void fit_mpeakBG_func()
{
	string workdir("/data/HPGe/USERS/kkw/DAQ/ANA125");

	char hisfile[256];
	char rawname[256];

	char runnumber3[256];
	char runnumber6[256];
	
	int runnum = 624;
	int binnum = 8000;

	sprintf(runnumber6,"%06d",runnum);
	sprintf(runnumber3,"%03d",runnum);

	int bin = binnum;// 1bin = 1kev
	//	int bin = 40000;// 1bin = 100ev
	int mul = bin/8000;
	//	for(int hn = 0; hn <hmax; hn++){

	int const N=2;
	double peak[N]={583, 609};

	//General resolution
	double respa[3];
	respa[0] = 0.000158;
	respa[1] = 0.547238;
	respa[2] = 0.006808;

/*
	// run511 resolution
		
        double respa[3];
        respa[0] = 0.000548398;
        respa[1] = 1.843895463;
        respa[2] = -0.024315872;

*/


	TF1 * gfit = new TF1("gfit","gaus");
	TF1 * agfit = new TF1("agfit","[0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2])");
	agfit->SetParNames("Area","Mean","Sigma"); //area = sqrt(2pi)*constant*sigma
	TF1 * p1fit = new TF1("p1fit","pol1");


	TString terms[N+1];

	for(int i=0; i<N; i++){
        terms[i] = Form("([%i]/(sqrt(2*TMath::Pi())*[%i])*exp(-((x-[%i])*(x-[%i]))/(2*[%i]*[%i])))", 3*i, 3*i+2, 3*i+1, 3*i+1, 3*i+2, 3*i+2);
	//	cout<<"term "<<i<<" "<<terms[i]<<endl;
	}
	
	terms[N]=Form("[%i]+[%i]*x",3*N, 3*N+1);
	//	cout<<"term "<<N<<" "<<terms[N]<<endl;

	TString model=terms[0];
	for(int i=1; i<=N; i++){
	  model=model+"+"+terms[i];
	  //	  cout<<model<<endl;
	  }

	
	TF1 * mfit = new TF1("mfit", model);

	for(int i=0; i<N; i++){
	  mfit->SetParName(3*i, Form("Area%i",i+1));
	  mfit->SetParName(3*i+1, Form("Mean%i",i+1));
	  mfit->SetParName(3*i+2, Form("Sigma%i",i+1));	  
	}
	mfit->SetParName(3*N, "Intercept");
	mfit->SetParName(3*N+1, "Slope");


	TCanvas * cbg = new TCanvas("cbg","BG fitting",1000,800);
	TCanvas * cfit = new TCanvas("cfit","peak fitting",1000,800);

	mfit->SetLineColor(kRed);
	p1fit->SetLineColor(kBlue);
	p1fit->SetLineStyle(7);

	int maxb[N], maxh[N];
	double maxx[N];

	double bg1, bg2, bg;
	double cnt1, cnt2, cnt3;
	double pa1[N][3];//gfit 
	double cal[N][2]; //integral fit func - m1fit,m2fit
	int xbin[N];
	double energy[N];
	double gfitpa[N][3], gfiter[N][3], gfitchi2[N];
	double agfitpa[N][3], agfiter[N][3], agfitchi2[N];
	double mfitpa[3*N+2], mfiter[3*N+2], mfitchi2;

	double fixbg[2];//bg line pa. for m1, m2, m3

	TFile *hf = new TFile(Form("%s/result/RUN%s/Bin%i/his_%s.rrs_on.root", workdir.c_str(), runnumber3, binnum, runnumber6));
	TH1D * his_temp = new TH1D("his_temp","",bin,0,4000);
      	his_temp = (TH1D*)hf->Get("his_tot3");//11det
	his_temp->Draw();
	his_temp->SetName("his_temp");
	his_temp->SetLineColor(1);
	//	his_temp->SetTitle(histitle);
	his_temp->SetXTitle("Energy[keV]");
	his_temp->SetYTitle("Counts/kev");


	double sigma_cal[N];

	for(int i = 0;i<N;i++){
	  energy[i] = peak[i];
	  xbin[i]=his_temp->GetXaxis()->FindBin(peak[i]);
	  printf("%.2f keV bin : %i	\n",peak[i],xbin[i]);

	  sigma_cal[i] = peak[i]*(respa[0] + (respa[1]/peak[i]) + (respa[2]/(sqrt(peak[i]))));
	  cout<<"Cal. Sigma : "<<sigma_cal[i]<<endl;

	  his_temp->GetXaxis()->SetRange(xbin[i]-1*mul,xbin[i]+1*mul);
	  
	  maxb[i] = his_temp->GetMaximumBin();
	  maxx[i] = (his_temp->GetBinCenter(maxb[i]));
	  cout<<"maxb(bin) : "<<maxb[i]<<endl;
	  cout<<"maxx(kev) : "<<maxx[i]<<endl;
	  maxh[i] = his_temp->GetMaximum();

	}

	cbg->cd();
	p1fit->SetRange((energy[0]-50),(energy[N-1]+50));
	his_temp->GetXaxis()->SetRangeUser(energy[0]-50, energy[N-1]+50);
	his_temp->Draw();
	his_temp->Fit("p1fit","R");
	fixbg[0] = p1fit->GetParameter(0);
	fixbg[1] = p1fit->GetParameter(1);
	cbg->Update();

	
	cfit->cd();
	his_temp->GetXaxis()->SetRange(0,4000*mul);	
	his_temp->Fit("p1fit","RQ+");

	for(int i=0; i<N; i++){
	gfit->SetRange((maxx[i]-3*sigma_cal[i]),(maxx[i]+3*sigma_cal[i]));
	agfit->SetRange((maxx[i]-3*sigma_cal[i]),(maxx[i]+3*sigma_cal[i]));

		//gaus fit
	his_temp->Fit("gfit","RQ0+");
	gfitpa[i][0] = gfit->GetParameter(0);
	gfitpa[i][1] = gfit->GetParameter(1);
	gfitpa[i][2] = gfit->GetParameter(2);
	//	gfiter[0] = gfit->GetParError(0);
	//	gfiter[1] = gfit->GetParError(1);
	//	gfiter[2] = gfit->GetParError(2);

	//area gaus fit

	agfit->SetParameters(2.5*gfitpa[i][1]*gfitpa[i][2],peak[i],sigma_cal[i]);

	his_temp->Fit("agfit","R0Q+");
	agfitpa[i][0] = agfit->GetParameter(0);
	agfitpa[i][1] = agfit->GetParameter(1);
	agfitpa[i][2] = agfit->GetParameter(2);
	//	agfiter[i][0] = agfit->GetParError(0);
	//	agfiter[i][1] = agfit->GetParError(1);
	//	agfiter[2] = agfit->GetParError(2);
	//	agfitchi2 = (agfit->GetChisquare())/(agfit->GetNDF());

	}

	mfit->SetRange(maxx[0]-5,maxx[N-1]+5);


	//		gfitchi2 = (gfit->GetChisquare())/(gfit->GetNDF());
		
		//gaus+ bg(pol1) fitting

	cout<<""<<endl;
	cout<<"[m-gauss fitting]"<<endl;
	cout<<""<<endl;
		//m2fit

	for(int i=0; i<3*N+2; i++){mfit->ReleaseParameter(i);}

	for(int i=0; i<N; i++){
	 mfit->SetParameter(3*i,agfitpa[i][0]);
	 mfit->SetParameter(3*i+1,peak[i]);
	 mfit->SetParameter(3*i+2,sigma_cal[i]);

	 mfit->SetParLimits(3*i+1, peak[i]-3*sigma_cal[i], peak[i]+3*sigma_cal[i]);
	 mfit->SetParLimits(3*i+2, 0.1*sigma_cal[i], 2*sigma_cal[i]);
	}

	mfit->FixParameter(3*N,fixbg[0]);
	mfit->FixParameter(3*N+1,fixbg[1]);


		//fitting
	his_temp->Fit("mfit","R+");

	for(int i=0; i<3*N+2; i++){
	  mfitpa[i] = mfit->GetParameter(i);
	  mfiter[i] = mfit->GetParError(i);
	}
	mfitchi2 = (mfit->GetChisquare())/(mfit->GetNDF());

		//		m2fitpara->Fill();


	cout<<">>chi squre check<<"<<endl;
	cout<<"mfit chi^2/NDF = "<<mfitchi2<<endl;
	cout<<""<<endl;
	
	his_temp->GetXaxis()->SetRangeUser(energy[0]-5, energy[N-1]+5);
	his_temp->Draw();
	cfit->Update();
	
	


}
