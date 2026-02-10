#include <iostream>
#include "TMath.h"
#include "TMinuit.h"
#include <string>

void fit_1peakBG_func()
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

	int const N=1;
	double peak[N]={1553.8};

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
	//	TF1 * efit = new TF1("efit","expo");
	//	TF1 * p0fit = new TF1("p0fit","pol0");
	TF1 * p1fit = new TF1("p1fit","pol1");

	//	TF1 * m1fit = new TF1("m1fit","([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+[3]");
	//	m1fit->SetParNames("Area","Mean","Sigma","Flat");

	//	TF1 * m2fit = new TF1("m2fit","([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+exp([3]+[4]*x)");
	//	m2fit->SetParNames("Area","Mean","Sigma","expo1","expo2");

	TF1 * mfit = new TF1("mfit", "([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+[3]+[4]*x");
	mfit->SetParNames("Area", "Mean", "Sigma", "Intercept", "Slope");


	TCanvas * cbg = new TCanvas("cbg","BG fitting",1000,800);
	TCanvas * cfit = new TCanvas("cfit","peak fitting",1000,800);

	mfit->SetLineColor(kRed);
	p1fit->SetLineColor(kBlue);
	p1fit->SetLineStyle(7);

	int maxb, maxh, run;
	double maxx;

	double bg1, bg2, bg;
	double cnt1, cnt2, cnt3;
	double pa1[N][3];//gfit 
	double cal[N][2]; //integral fit func - m1fit,m2fit
	int xbin[N];
	double energy;
	double gfitpa[3], gfiter[3], gfitchi2;
	double agfitpa[3], agfiter[3], agfitchi2;
	double mfitpa[5], mfiter[5], mfitmean, mfitsig, mfitchi2;

	double fitpa[3], fiter[3], fitchi2;
	double fixbg[2];//bg line pa. for m1, m2, m3

	/*
	TTree *m3fitpara = new TTree("m3fitpara", "gaus+pol1 fitting result");
	m3fitpara -> Branch("run", &runnum, "run/I");
	m3fitpara -> Branch("peak", &energy, "energy/D");
	m3fitpara -> Branch("area", &mfitpa[0], "mfitpa0/D");
	m3fitpara -> Branch("mean", &mfitpa[1], "mfitpa1/D");
	m3fitpara -> Branch("sigma", &mfitpa[2], "mfitpa2/D");
	m3fitpara -> Branch("pol0", &mfitpa[3], "mfitpa3/D");
	m3fitpara -> Branch("pol1", &mfitpa[4], "mfitpa4/D");
	m3fitpara -> Branch("chi2", &mfitchi2, "mfitchi2/D");
	m3fitpara -> Branch("area_er", &mfiter[0], "mfiter0/D");
	m3fitpara -> Branch("mean_er", &mfiter[1], "mfiter1/D");
	m3fitpara -> Branch("sigma_er", &mfiter[2], "mfiter2/D");
	m3fitpara -> Branch("pol0_er", &mfiter[3], "mfitper3/D");
	m3fitpara -> Branch("pol1_er", &mfiter[4], "mfitper4/D");
	*/

	TFile *hf = new TFile(Form("%s/result/RUN%s/Bin%i/his_%s.rrs_on.root", workdir.c_str(), runnumber3, binnum, runnumber6));
	TH1D * his_temp = new TH1D("his_temp","",bin,0,4000);
      	his_temp = (TH1D*)hf->Get("his_tot3");//11det
	his_temp->Draw();
	his_temp->SetName("his_temp");
	his_temp->SetLineColor(1);
	//	his_temp->SetTitle(histitle);
	his_temp->SetXTitle("Energy[keV]");
	his_temp->SetYTitle("Counts/kev");


	double sigma_cal;

	for(int i = 0;i<N;i++){

		energy = peak[i];
		xbin[i]=his_temp->GetXaxis()->FindBin(peak[i]);
		printf("\n");
		printf("\n=========================\n");
		printf("%.2f keV bin : %i	\n",peak[i],xbin[i]);

		his_temp->GetXaxis()->SetRange(xbin[i]-1*mul,xbin[i]+1*mul);

		sigma_cal = peak[i]*(respa[0] + (respa[1]/peak[i]) + (respa[2]/(sqrt(peak[i]))));
		cout<<"Cal. Sigma : "<<sigma_cal<<endl;

		maxb = his_temp->GetMaximumBin();
		//maxx = maxb/mul;
		maxx = (his_temp->GetBinCenter(maxb));
		cout<<"maxb(bin) : "<<maxb<<endl;
		cout<<"maxx(kev) : "<<maxx<<endl;
		maxh = his_temp->GetMaximum();

		his_temp->GetXaxis()->SetRange(0,4000*mul);
//		gfit->SetRange((maxx-5),(maxx+5));
//		agfit->SetRange((maxx-5),(maxx+5));
//		m1fit->SetRange((maxx-10),(maxx+10));
//		m2fit->SetRange((maxx-10),(maxx+10));
//		m3fit->SetRange((maxx-10),(maxx+10));
		gfit->SetRange((maxx-3*sigma_cal),(maxx+3*sigma_cal));
		agfit->SetRange((maxx-3*sigma_cal),(maxx+3*sigma_cal));
		mfit->SetRange((maxx-3*sigma_cal),(maxx+3*sigma_cal));
//		m3fit->SetRange((maxx-24*sigma_cal),(maxx+24*sigma_cal));


		p1fit->SetRange((maxx-50),(maxx+50));

		cbg->cd();
		his_temp->Fit("p1fit","R");
		fixbg[0] = p1fit->GetParameter(0);
		fixbg[1] = p1fit->GetParameter(1);
		his_temp->GetXaxis()->SetRangeUser(energy-50, energy+50);
		his_temp->Draw();
		cbg->Update();

		cfit->cd();
		his_temp->Fit("p1fit","RQ+");
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

		agfit->SetParameters(2.5*gfitpa[1]*gfitpa[2],peak[i],sigma_cal);

		his_temp->Fit("agfit","R0Q+");
		agfitpa[0] = agfit->GetParameter(0);
		agfitpa[1] = agfit->GetParameter(1);
		agfitpa[2] = agfit->GetParameter(2);
		agfiter[0] = agfit->GetParError(0);
		agfiter[1] = agfit->GetParError(1);
		agfiter[2] = agfit->GetParError(2);
		agfitchi2 = (agfit->GetChisquare())/(agfit->GetNDF());


		//gaus+ bg(pol1) fitting

		cout<<""<<endl;
		cout<<"[m2 fitting]"<<endl;
		cout<<""<<endl;
		//m2fit

		mfit->ReleaseParameter(0);
		mfit->ReleaseParameter(1);
		mfit->ReleaseParameter(2);
		mfit->ReleaseParameter(3);
		mfit->ReleaseParameter(4);

		mfit->SetParameters(agfitpa[0],maxx,sigma_cal,fixbg[0],fixbg[1]);
		mfit->FixParameter(3,fixbg[0]);
		mfit->FixParameter(4,fixbg[1]);
		mfit->SetParLimits(1, energy-5, energy+5);
		mfit->SetParLimits(2, 0.1*sigma_cal, 2*sigma_cal);
		/*
		//peak check
		his_temp->Fit("m2fit","RQ0+");
		m2fitpa[1] = m2fit->GetParameter(1);
		if(abs(m2fitpa[1]-peak[i])>3*(sigma_cal)){
			cout <<"m2 test fitting mean "<<m2fitpa[1]<<" for peak "<<peak[i]<<endl;
			cout<<">> peak fixed at"<<peak[i]<<endl;
			m2fit->FixParameter(1,maxx);
		}

		//sigma check
		his_temp->Fit("m2fit","RQ0+");
		m2fitpa[2] = m2fit->GetParameter(2);
		if((m2fitpa[2]<(0.8*sigma_cal))||(m2fitpa[2]>(3*sigma_cal))){
			cout <<"m2 test fitting sigma "<<m2fitpa[2]<<" for peak "<<peak[i]<<endl;
			cout<<">> sigma fixed at "<<sigma_cal<<endl;
			m2fit->FixParameter(2,sigma_cal);
		}

		//peak re-check
		his_temp->Fit("m2fit","RQ0+");
		m2fitpa[1] = m2fit->GetParameter(1);
		if(abs(m2fitpa[1]-peak[i])>3*(sigma_cal)){
			cout <<"m2 test fitting mean "<<m2fitpa[1]<<" for peak "<<peak[i]<<endl;
			cout<<">> peak fixed at"<<peak[i]<<endl;
			m2fit->FixParameter(1,maxx);
		}
		*/

		//fitting
		his_temp->Fit("mfit","R+");

		mfitpa[0] = mfit->GetParameter(0);
		mfitpa[1] = mfit->GetParameter(1);
		mfitpa[2] = mfit->GetParameter(2);
		mfitpa[3] = mfit->GetParameter(3);
		mfitpa[4] = mfit->GetParameter(4);
		mfiter[0] = mfit->GetParError(0);
		mfiter[1] = mfit->GetParError(1);
		mfiter[2] = mfit->GetParError(2);
		mfiter[3] = mfit->GetParError(3);
		mfiter[4] = mfit->GetParError(4);
		mfitchi2 = (mfit->GetChisquare())/(mfit->GetNDF());

		//		m2fitpara->Fill();


		cout<<">>chi squre check<<"<<endl;
		cout<<"mfit chi^2/NDF = "<<mfitchi2<<endl;
		cout<<""<<endl;

		his_temp->GetXaxis()->SetRangeUser(energy-10, energy+10);
		his_temp->Draw();
		cfit->Update();

		/*
	char histitle[256];
	sprintf(histitle,"source peaks Fitting");

	c1->cd();
	his_temp->Draw();
		*/
}
}
