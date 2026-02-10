#include <iostream>
#include "TMath.h"
#include "TMinuit.h"
#include <string>

void fit_5peak_multi_250514()
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
	int expnum=2;

	sprintf(runnumber6,"%06d",runnum);
	sprintf(runnumber3,"%03d",runnum);

	int binnum = 8000;// 1bin = 1kev
	int bw = 4000/binnum; //bin width, keV

	int const N = 5; //peaks

	double peak[N] = {768.36, 772.26, 783.3, 785.4, 794.95};
	//	double peak[N] = {768.36, 772.26, 783.3, 785.96, 794.95};

	double pr1, pr2;
	//4000bin rrs on 
	//	double peakfix[N] = {767.506, 771.525, 782.902, 785.066, 785.643, 794.908};//rrs on 4000
	//	double peakfix[N] = {767.42, 771.419, 782.739, 784.892, 785.466, 794.684};//tot 4000
	//	double peakfix[N] = {767.326, 771.373, 782.828, 785.007, 785.589, 794.917};//rrs on 8000
	//		double peakfix[N] = {767.403, 771.402, 782.723, 784.877, 785.451, 794.67};//tot 8000
		double peakfix[N] = {768.745, 772.26, 783.3, 786.235, 795.723};//peak 6, v1 new cal, run1

		//run2 v1. 609, 911 16000bin
		//	pr1=0.3;		pr2=0.3; // v1 new cal. run1

	//General resolution


	double respa[3];

	if(expnum==1){
	  respa[0] = 0.498182;
	  respa[1] = 0.000373447;
          respa[2] = 4.40e-8;
	}

	if(expnum==2){
	  respa[0] = 0.446;
	  respa[1] = 6.06E-4;
	  respa[2] = -6.31E-8;
	}

	TF1 * gfit = new TF1("gfit","gaus");
	TF1 * agfit = new TF1("agfit","[0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2])");
	agfit->SetParNames("Area","Mean","Sigma"); //area = sqrt(2pi)*constant*sigma

	TF1 * p0fit = new TF1("p0fit","pol0");
	TF1 * p1fit = new TF1("p1fit","pol1");

	TString term1 = "([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/(2*[2]*[2])))";
	TString term2 = "([3]/(sqrt(2*TMath::Pi())*[5])*exp(-((x-[4])*(x-[4]))/(2*[5]*[5])))";
	TString term3 = "([6]/(sqrt(2*TMath::Pi())*[8])*exp(-((x-[7])*(x-[7]))/(2*[8]*[8])))";
	TString term4 = "([9]/(sqrt(2*TMath::Pi())*[11])*exp(-((x-[10])*(x-[10]))/(2*[11]*[11])))";
	TString term5 = "([12]/(sqrt(2*TMath::Pi())*[14])*exp(-((x-[13])*(x-[13]))/(2*[14]*[14])))";	
	TString term6 = "([15]/(sqrt(2*TMath::Pi())*[17])*exp(-((x-[16])*(x-[16]))/(2*[17]*[17])))";
	TString linear = "([15] + [16]*x)";
	//	TString linear = "([18] + [19]*x)";
	TString formula = term1 + " + " + term2 + " + " + term3 + " + " + term4 + " + " + term5 + " + " + linear;
	//TString formula = term1 + " + " + term2 + " + " + term3 + " + " + term4 + " + " + term5 + " + " + linear;

	TF1 *fit = new TF1("fit", formula);

	fit->SetParName(0,"Area1");
	fit->SetParName(1,"Mean1");
	fit->SetParName(2,"Sigma1");

	fit->SetParName(3,"Area2");
	fit->SetParName(4,"Mean2");
	fit->SetParName(5,"Sigma2");

	fit->SetParName(6,"Area3");
	fit->SetParName(7,"Mean3");
	fit->SetParName(8,"Sigma3");

	fit->SetParName(9,"Area4");
	fit->SetParName(10,"Mean4");
	fit->SetParName(11,"Sigma4");

	fit->SetParName(12,"Area5");
	fit->SetParName(13,"Mean5");
	fit->SetParName(14,"Sigma5");
	/*
	fit->SetParName(15,"Area6");
	fit->SetParName(16,"Mean6");
	fit->SetParName(17,"Sigma6");
	*/
	fit->SetParName(15,"Intercept");
	fit->SetParName(16,"Slope");


	TCanvas * c1 = new TCanvas("c1","fitting",1200,800);

	fit->SetLineColor(kRed);
	gfit->SetLineColor(kBlue);
	p1fit->SetLineColor(kGreen);
	p1fit->SetLineStyle(7);

	int maxb, maxh, run;
	double maxx;
	double pa1[N][3];//gfit 

	int xbin[N];
	double energy;
	double gfitpa[N][3], gfiter[N][3], gfitchi2;
	double agfitpa[N][3], agfiter[N][3], agfitchi2;
	double fitpa[3*N+2], fiter[3*N+2], fitchi2;

	//	double fitpa[3], fiter[3], fitchi2;
	double fixbg[2];//bg line pa. for m1, m2, m3

	TTree *fitpara = new TTree("fitpara", "gaus+pol1 fitting result");
	fitpara -> Branch("run", &runnum, "run/I");
	fitpara -> Branch("peak", &energy, "energy/D");
	
	for(int i=0; i<N; i++){
	  fitpara->Branch(Form("area_p%i",i),&fitpa[3*i], Form("fitpa[%i]/D",3*i));
	  fitpara->Branch(Form("mean_p%i",i),&fitpa[3*i+1], Form("fitpa[%i]/D",3*i+1));
	  fitpara->Branch(Form("sigma_p%i",i),&fitpa[3*i+2], Form("fitpa[%i]/D",3*i+2));

	  fitpara->Branch(Form("area_err%i",i),&fiter[3*i], Form("fiter[%i]/D",3*i));
	  fitpara->Branch(Form("mean_err%i",i),&fiter[3*i+1], Form("fiter[%i]/D",3*i+1));
	  fitpara->Branch(Form("sigma_err%i",i),&fiter[3*i+2], Form("fiter[%i]/D",3*i+2));
	}
	
	//6 peaks
	fitpara -> Branch("pol0", &fitpa[15], "fitpa[15]/D");
	fitpara -> Branch("pol1", &fitpa[16], "fitpa[16]/D");

	fitpara -> Branch("pol0_err", &fiter[15], "fitper[15]/D");
	fitpara -> Branch("pol1_err", &fiter[16], "fitper[16]/D");


	TF1 *m1bg[N];
	TF1 *m2bg[N];
	char m1bgname[256];
	char m2bgname[256];


	if(expnum==1){
	  //	  sprintf(hisfile,"%s/result/V50/Bin%i/his_000624.v0.rrs_on.root", workdir.c_str(), binnum);
	  //	  sprintf(resfile,"%s/result/V50/Bin%i/mfit_p5_000624.v0.rrs_on.root", workdir.c_str(), binnum);
	  	  sprintf(hisfile,"%s/result/V50/Bin%i/his_run1.root", workdir.c_str(), binnum);
	  	  sprintf(resfile,"%s/result/V50/Bin%i/mfit_p5_run1.root", workdir.c_str(), binnum);

	}

	if(expnum==2){
	  sprintf(hisfile,"%s/result/V50/Bin%i/his_000650.v1cut.root", workdir.c_str(), binnum);
	  sprintf(resfile,"%s/result/V50/Bin%i/mfit_p5_000650_v1cut.root", workdir.c_str(), binnum);
	  //	  sprintf(hisfile,"%s/result/V50/Bin%i/his_000650.v1.root", workdir.c_str(), binnum);
	  //	  sprintf(resfile,"%s/result/V50/Bin%i/mfit_000650_v1.root", workdir.c_str(), binnum);
	}


	TFile *hf=new TFile(hisfile);
	TH1D * his_temp = new TH1D("his_temp","",binnum,0,4000);

//	his_temp = (TH1D*)hf->Get("his_tot1");//13det
//	his_temp = (TH1D*)hf->Get("his_tot3");//12det(4,10)
//      	his_temp = (TH1D*)hf->Get("his_tot4");//11det rrs_on.
//      	his_temp = (TH1D*)hf->Get("his_run1");//11det 2peak
      	his_temp = (TH1D*)hf->Get("his_tot3");//11det 2peak
	his_temp->Draw();
	his_temp->SetName("his_temp");
	his_temp->SetLineColor(1);

	double sigma_cal[N];

	for(int i = 0; i<N; i++){

		c1->cd();
		energy = peak[i];
		//		xbin[i]=his_temp->GetXaxis()->FindBin(peak[i]);
		xbin[i]=his_temp->FindBin(peak[i]);

		printf("\n");
		printf("\n=========================\n");
		printf("%.2f keV bin : %i	\n",peak[i],xbin[i]);

		//		his_temp->GetXaxis()->SetRange(xbin[i]-1*mul,xbin[i]+1*mul);
		//		sigma_cal[i] = peak[i]*(respa[0] + (respa[1]/peak[i]) + (respa[2]/(sqrt(peak[i]))));
		sigma_cal[i] = sqrt(respa[0]*respa[0]+respa[1]*peak[i]+respa[2]*peak[i]*peak[i]);

		cout<<"Cal. Sigma : "<<sigma_cal[i]<<" (3 sigma : "<<3*sigma_cal[i]<<endl;
		cout<<"3 sigma range : "<<his_temp->FindBin(peak[i]-(3*sigma_cal[i]))<<" ~ "<<his_temp->FindBin(peak[i]+(3*sigma_cal[i]))<<endl;
		//		his_temp->GetXaxis()->SetRange(his_temp->Getxbin[i]-1*mul,xbin[i]+1*mul);
		his_temp->GetXaxis()->SetRange(his_temp->FindBin(peak[i]-(3*sigma_cal[i])),his_temp->FindBin(peak[i]+(3*sigma_cal[i])));


		maxb = his_temp->GetMaximumBin();
		maxx = (his_temp->GetBinCenter(maxb));
		cout<<"maxb(bin) : "<<maxb<<endl;
		cout<<"maxx(kev) : "<<maxx<<endl;
		maxh = his_temp->GetMaximum();

		his_temp->GetXaxis()->SetRange(0,binnum);
		gfit->SetRange((maxx-3*sigma_cal[i]),(maxx+3*sigma_cal[i]));
		agfit->SetRange((maxx-3*sigma_cal[i]),(maxx+3*sigma_cal[i]));


		//gaus fit
		//		his_temp->Fit("gfit","RQ0+");
		his_temp->Fit("gfit","R0+");
		gfitpa[i][0] = gfit->GetParameter(0);
		gfitpa[i][1] = gfit->GetParameter(1);
		gfitpa[i][2] = gfit->GetParameter(2);
		gfiter[i][0] = gfit->GetParError(0);
		gfiter[i][1] = gfit->GetParError(1);
		gfiter[i][2] = gfit->GetParError(2);

		gfitchi2 = (gfit->GetChisquare())/(gfit->GetNDF());


		//area gaus fit

		//agfit->SetParameters(2.5*gfitpa[1]*gfitpa[2],gfitpa[1],gfitpa[2]);
		agfit->SetParameters(2.5*gfitpa[i][1]*gfitpa[i][2],peak[i],sigma_cal[i]);
		//agfit->SetParameters(2.5*gfitpa2*gfitpa3,maxx,2);
		his_temp->Fit("agfit","R0Q+");
		//		his_temp->Fit("agfit","R+");
		agfitpa[i][0] = agfit->GetParameter(0);
		agfitpa[i][1] = agfit->GetParameter(1);
		agfitpa[i][2] = agfit->GetParameter(2);
		agfiter[i][0] = agfit->GetParError(0);
		agfiter[i][1] = agfit->GetParError(1);
		agfiter[i][2] = agfit->GetParError(2);
		agfitchi2 = (agfit->GetChisquare())/(agfit->GetNDF());
	}


	//bg line fit
	//	p1fit->SetRange((his_temp->FindBin(peak[0]-50)),(his_temp->FindBin(peak[4]+50)));
	p1fit->SetRange(peak[0]-50,peak[N-1]+50);
	//		p1fit->SetRange((maxx-30),(maxx+30));

		//		his_temp->Fit("p0fit","RQ0+");
		//		fixbg1 = p0fit->GetParameter(0);

	his_temp->Fit("p1fit","R+");
	fixbg[0] = p1fit->GetParameter(0);
	fixbg[1] = p1fit->GetParameter(1);


	//main fit
	//	fit->SetRange((his_temp->FindBin(peak[0]-50)),(his_temp->FindBin(peak[4]+50)));
	fit->SetRange(peak[0]-50,peak[4]+50);

	for(int i=0; i<(3*N+2); i++){fit->ReleaseParameter(i);}

	//	double pr; //par limit range
	//	for(int i=0; i<5; i++){
	for(int i=0; i<N; i++){
	  //	  fit->SetParameter(i*3,agfitpa[i][0]);
	  //	  fit->SetParameter(i*3,0);

	  //peak cnt preset
	  //before
	  /*
	  if(i==0){fit->SetParameter(i*3,agfitpa[i][0]);}
	  if(i==5){fit->SetParameter(i*3,agfitpa[i][0]);}
	  else{	  fit->SetParameter(i*3,0);}
	  */
	  
	  //2nd run

	  if(expnum==1){
	    fit->SetParameter(0,63.*bw);
	    fit->SetParameter(3,10.*bw);
	    fit->SetParameter(6,0.*bw);
	    fit->SetParameter(9,24.*bw);
	    //	    fit->SetParameter(12,6.*bw);
	    fit->SetParameter(15,31.*bw);
	  }

	  if(expnum==2){
	    fit->SetParameter(0,387.*bw);
	    fit->SetParameter(3,12.*bw);
	    fit->SetParameter(6,0.*bw);
	    fit->SetParameter(9,(26.+13.)*bw);
	    //	    fit->SetParameter(12,13.*bw);
	    fit->SetParameter(15,36.*bw);
	  }
	  //	  fit->SetParLimits(9,-10,100);
	  //	  fit->SetParLimits(12,0,13.*bw*3);
	  //	  fit->SetParameter(i*3+1,his_temp->FindBin(peak[i]));
	  
	  fit->SetParameter(i*3+1,peak[i]);
	  //	  fit->SetParLimits(i*3+1,his_temp->FindBin(peak[i])-3*sigma_cal[i], his_temp->FindBin(peak[i])+3*sigma_cal[i]);
	  //	  pr=2*sigma_cal[i];
	  //	  	  pr=1.22;//rrs on hist
	  //	  pr=0.91;//new cal 8000 
	  //		  pr=0.82;//new cal 4000
	  
	  //his1
	  
	  if(expnum==1){pr1=0.3; pr2=0.3;}
	  if(expnum==2){pr1=0.3; pr2=0.3;}
	  // 	  pr1=sigma_cal[i]*3;	  	  	  pr2=sigma_cal[i]*3;


		  //	  fit->SetParLimits(i*3+1,(peak[i]-pr), (peak[i]+pr));
	  //	  fit->SetParLimits(i*3+1,(peak[i]-pr), (peak[i]));
	  fit->SetParLimits(i*3+1,(peak[i]-pr1), (peak[i]+pr2));
	  //	  cout<<"pr "<<pr<<endl;
	  //	  cout<<"peak "<<i<<" par limit "<<peak[i]-pr<<" - "<<peak[i]<<endl;
	  cout<<"peak "<<i<<" par limit "<<peak[i]-pr1<<" - "<<peak[i]+pr2<<endl;
  	  fit->SetParameter(i*3+2,sigma_cal[i]);
	  fit->SetParLimits(i*3+2,0.8*sigma_cal[i],1.2*sigma_cal[i]);
	}


	//6peak
	
	//5peaks
	fit->SetParameter(15,fixbg[0]);
	fit->SetParameter(16,fixbg[1]);

	fit->FixParameter(15,fixbg[0]);
	fit->FixParameter(16,fixbg[1]);
	

		
		//sigma check

	his_temp->Fit("fit","RL0+");
	for(int i=0; i<N; i++){
	  fitpa[i*3+2]=fit->GetParameter(i*3+2);
 	  if((fitpa[i*3+2]<(0.8*sigma_cal[i]))||(fitpa[i*3+2]>(3*sigma_cal[i]))){
	    cout <<"peak"<<peak[i]<<" sigma"<<i<<" : "<<fitpa[i*3+2]<<endl;
	    cout<<">> sigma fixed at "<<sigma_cal[i]<<endl;
	    fit->FixParameter(i*3+2, sigma_cal[i]);
	    
	  }
	}



	his_temp->Fit("fit","RL+");
	
	

	

	//	for(int ip=0; ip<=19; ip=ip+3){ 
	for(int ip=0; ip<=16; ip=ip+3){ 
	  if(ip==15){cout<<fit->GetParameter(ip)<<" * "<<fit->GetParError(ip)<<" * "<<fit->GetParameter(ip+1)<<" * "<<fit->GetParError(ip+1)<<endl;}
	  else{cout<<fit->GetParameter(ip)<<" * "<<fit->GetParError(ip)<<" * "<<fit->GetParameter(ip+1)<<" * "<<fit->GetParError(ip+1)<< " * " <<fit->GetParameter(ip+2)<<" * "<<fit->GetParError(ip+2)<<endl;}
	  
	}
	
			/*
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

		cout<<fitpa[0]<<endl;
		cout<<fitpa[1]<<endl;
		cout<<fitpa[2]<<endl;
		fitpara->Fill();
		cout<<""<<endl;
		cout<<""<<endl;



	


		fit->SetRange((maxx-6*sigma_cal),(maxx+6*sigma_cal));
//		fit->SetRange((maxx-24*sigma_cal),(maxx+24*sigma_cal));

//		p0fit->SetRange((maxx-30),(maxx+30));

		cout<<""<<endl;
		cout<<"[fin fitting]"<<endl;
		cout<<""<<endl;
		//fit

		*/
	

	char histitle[256];
	sprintf(histitle,"peaks Fitting");
	his_temp->SetTitle(histitle);
	his_temp->SetXTitle("Energy[keV]");
	his_temp->SetYTitle("Counts/kev");

	c1->cd();
	his_temp->GetXaxis()->SetRangeUser(750,815);
	//	his_temp->GetYaxis()->SetRangeUser(20,150);//4000
		his_temp->GetYaxis()->SetRangeUser(10,80);//8000
	his_temp->Draw();
	TFile rf (resfile,"RECREATE");
	rf.cd();
	//	gfitpara->Write();
	fitpara->Write();
	c1->Write();

	rf.Close();
	
}
