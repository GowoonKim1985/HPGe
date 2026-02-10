#include <iostream>
#include "TMath.h"
#include "TMinuit.h"
#include <string>

void fit_7peak_multi_251024()
{
//  string workdir("/data/HPGe/USERS/kkw/DAQ/ANA125");
  string workdir("/home/kkw/study/V50");



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

  int binnum = 4000;// 1bin = 1kev
  double bw = 4000./((double)binnum); //bin width, keV
  
  int const N = 7; //peaks
  double peak[N] = {768.36, 772.29, 783.29, 785.4, 785.96, 794.95, 806.179};
  
  double pr1, pr2;
  //4000bin rrs on 
  //	double peakfix[N] = {767.506, 771.525, 782.902, 785.066, 785.643, 794.908};//rrs on 4000
  //	double peakfix[N] = {767.42, 771.419, 782.739, 784.892, 785.466, 794.684};//tot 4000
  //	double peakfix[N] = {767.326, 771.373, 782.828, 785.007, 785.589, 794.917};//rrs on 8000
  //		double peakfix[N] = {767.403, 771.402, 782.723, 784.877, 785.451, 794.67};//tot 8000
  //  double peakfix[N] = {768.745, 772.26, 783.3, 785.4, 786.235, 795.723};//peak 6, v1 new cal, run1
  
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
  TString term7 = "([18]/(sqrt(2*TMath::Pi())*[20])*exp(-((x-[19])*(x-[19]))/(2*[20]*[20])))";
  //	TString linear = "([15] + [16]*x)";
  TString linear = "([21] + [22]*x)";
  TString formula = term1 + " + " + term2 + " + " + term3 + " + " + term4 + " + " + term5 + " + " + term6 + " + " + term7 " + " + linear;
  
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
  
  fit->SetParName(15,"Area6");
  fit->SetParName(16,"Mean6");
  fit->SetParName(17,"Sigma6");

  fit->SetParName(18,"Area7");
  fit->SetParName(19,"Mean7");
  fit->SetParName(20,"Sigma7");

  fit->SetParName(21,"Intercept");
  fit->SetParName(22,"Slope");
  
  
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
  fitpara -> Branch("pol0", &fitpa[21], "fitpa[21]/D");
  fitpara -> Branch("pol1", &fitpa[22], "fitpa[22]/D");
  
  fitpara -> Branch("pol0_err", &fiter[21], "fitper[21]/D");
  fitpara -> Branch("pol1_err", &fiter[22], "fitper[22]/D");
  

  TF1 *m1bg[N];
  TF1 *m2bg[N];
  char m1bgname[256];
  char m2bgname[256];
  
  if(expnum==1){
	  //	  sprintf(hisfile,"%s/result/V50/Bin%i/his_000624.v0.rrs_on.root", workdir.c_str(), binnum);
	  //	  sprintf(resfile,"%s/result/V50/Bin%i/mfit_000624.v0.rrs_on.root", workdir.c_str(), binnum);
    sprintf(hisfile,"%s/result/V50/Bin%i/his_run1.root", workdir.c_str(), binnum);
    sprintf(resfile,"%s/result/V50/Bin%i/mfit_run1.v0.root", workdir.c_str(), binnum);
  }

  if(expnum==2){
    //    sprintf(hisfile,"%s/result/V50/Bin%i/his_000650.v1cut.root", workdir.c_str(), binnum);
    //    sprintf(resfile,"%s/result/V50/Bin%i/mfit_p6_000650_v1cut.root", workdir.c_str(), binnum);
    sprintf(hisfile,"%s/result/Bin%i/his_000650.v1.root", workdir.c_str(), binnum);
    sprintf(resfile,"%s/result/Bin%i/mfit_000650_v1.root", workdir.c_str(), binnum);
  }


  TFile *hf=new TFile(hisfile);
  TH1D * his_temp = new TH1D("his_temp","",binnum,0,4000);

  his_temp = (TH1D*)hf->Get("his_tot3");//11det
  his_temp->Draw();
  his_temp->SetName("his_temp");
  his_temp->SetLineColor(1);

  double sigma_cal[N];

  //pre fitting
  for(int i = 0; i<N; i++){

    c1->cd();
    energy = peak[i];
    xbin[i]=his_temp->FindBin(peak[i]);
    
    printf("\n");
    printf("\n=========================\n");
    printf("%.2f keV bin : %i	\n",peak[i],xbin[i]);
    
    //sigma calculation
    sigma_cal[i] = sqrt(respa[0]*respa[0]+respa[1]*peak[i]+respa[2]*peak[i]*peak[i]);
    
    cout<<"Cal. Sigma : "<<sigma_cal[i]<<" (3 sigma : "<<3*sigma_cal[i]<<endl;
    cout<<"3 sigma range : "<<his_temp->FindBin(peak[i]-(3*sigma_cal[i]))<<" ~ "<<his_temp->FindBin(peak[i]+(3*sigma_cal[i]))<<endl;
    his_temp->GetXaxis()->SetRange(his_temp->FindBin(peak[i]-(3*sigma_cal[i])),his_temp->FindBin(peak[i]+(3*sigma_cal[i])));
    
    // fitting range set : btw x@maxH +/- 3sigma
    maxb = his_temp->GetMaximumBin();
    maxx = (his_temp->GetBinCenter(maxb));
    cout<<"maxb(bin) : "<<maxb<<endl;
    cout<<"maxx(kev) : "<<maxx<<endl;
    maxh = his_temp->GetMaximum();
    
    his_temp->GetXaxis()->SetRange(0,binnum);
    gfit->SetRange((maxx-3*sigma_cal[i]),(maxx+3*sigma_cal[i]));
    agfit->SetRange((maxx-3*sigma_cal[i]),(maxx+3*sigma_cal[i]));
    
    
    //gaus fit : to get pre set values
    //		his_temp->Fit("gfit","RQ0+");
    his_temp->Fit("gfit","R0+");
    gfitpa[i][0] = gfit->GetParameter(0);
    gfitpa[i][1] = gfit->GetParameter(1);
    gfitpa[i][2] = gfit->GetParameter(2);
    gfiter[i][0] = gfit->GetParError(0);
    gfiter[i][1] = gfit->GetParError(1);
    gfiter[i][2] = gfit->GetParError(2);
    
    gfitchi2 = (gfit->GetChisquare())/(gfit->GetNDF());
    
    
    //area gaus fit : to get pre set values
    
    agfit->SetParameters(2.5*gfitpa[i][1]*gfitpa[i][2],peak[i],sigma_cal[i]);
    his_temp->Fit("agfit","R0Q+");
    agfitpa[i][0] = agfit->GetParameter(0);
    agfitpa[i][1] = agfit->GetParameter(1);
    agfitpa[i][2] = agfit->GetParameter(2);
    agfiter[i][0] = agfit->GetParError(0);
    agfiter[i][1] = agfit->GetParError(1);
    agfiter[i][2] = agfit->GetParError(2);
    agfitchi2 = (agfit->GetChisquare())/(agfit->GetNDF());
  }
  
  //bg line fit : pre set values, pol1
  p1fit->SetRange(peak[0]-50,peak[N-1]+50);
  his_temp->Fit("p1fit","R+");
  fixbg[0] = p1fit->GetParameter(0);
  fixbg[1] = p1fit->GetParameter(1);
  
  
  //main fit

  for(int i=0; i<(3*N+2); i++){fit->ReleaseParameter(i);} //parameter reset
  //  fit->SetRange(peak[0]-50,peak[5]+50);
    fit->SetRange(peak[0]-50,peak[5]+50);

  //preset parameter, limit
  if(expnum==1){
    fit->SetParameter(0,63.*bw);
    fit->SetParameter(3,10.*bw);
    fit->SetParameter(6,0.*bw);
    fit->SetParameter(9,24.*bw);
    fit->SetParameter(12,6.*bw);
    fit->SetParameter(15,31.*bw);
  }

  if(expnum==2){
    //area pre set
    fit->SetParameter(0, 365.0*bw);//768.36kev
    fit->SetParameter(3, 6.32*bw);//772.29kev
    fit->SetParameter(6, 24.2*bw);//783.29 kev,v50
 //   fit->SetParameter(6, 0);//783.29 kev,v50
    fit->SetParameter(9, 26.2*bw);//785.4 kev
    fit->SetParameter(12, 102.1*bw);//785.96kev
    fit->SetParameter(15, 18.0*bw);//794.95kev
    
    
    //area limit
    fit->SetParLimits(0, 234.*bw, 495.9*bw);
    fit->SetParLimits(3, 0., 13.6*bw);
    fit->SetParLimits(6, 0., 24.2*bw);
//    fit->SetParLimits(6, 0.*bw, 1.*bw);
//    fit->SetParLimits(9, 1.*bw, 51.1*bw);
    fit->SetParLimits(9, 1.*bw, 51.1*bw);
    fit->SetParLimits(12, 59.*bw, 145.5*bw);
    fit->SetParLimits(15, 0., 38.6*bw);
     
    }

  //bg par set
  fit->SetParameter(18,fixbg[0]);
  fit->SetParameter(19,fixbg[1]);

  for(int i=0; i<N; i++){
    fit->SetParameter(i*3+1,peak[i]);//peak preset
	  
    //peak range limit
    if(expnum==1){pr1=0.3; pr2=0.3;}
    if(expnum==2){pr1=0.3; pr2=0.3;}
    
    fit->SetParLimits(i*3+1,(peak[i]-pr1), (peak[i]+pr2));
    cout<<"peak "<<i<<" par limit "<<peak[i]-pr1<<" - "<<peak[i]+pr2<<endl;

    //sigma set + limit
    fit->SetParameter(i*3+2,sigma_cal[i]);
    fit->SetParLimits(i*3+2,0.8*sigma_cal[i],1.2*sigma_cal[i]);

//    fit->SetParLimits(i*3+2,0.7*sigma_cal[i],1.3*sigma_cal[i]);

/*
    //sigma fix for low peak 
    fit->FixParameter(5,sigma_cal[1]);
    fit->FixParameter(8,sigma_cal[2]);
    fit->FixParameter(11,sigma_cal[3]);
    fit->FixParameter(17,sigma_cal[5]);
*/

    cout<<""<<endl;
    cout<<"peak "<<i<<" "<<peak[i]<<" keV"<<endl;
    cout<<"peak range : "<<peak[i]-pr1<<" ~ "<<peak[i]+pr2<<endl;
    cout<<"sigma set : "<<sigma_cal[i]<< " ( "<<0.8*sigma_cal[i]<<" ~ "<<1.2*sigma_cal[i]<<" )"<<endl;

  }
  
  /*		
  //fix peak for 0 count
  fit->FixParameter(5, sigma_cal[1]);// peak2 sigma 772.29 kev
  fit->FixParameter(8, sigma_cal[2]);// peak3 sigma 783.29 keV  
  fit->FixParameter(11, sigma_cal[3]);// peak4 sigma 785.4 keV
  fit->FixParameter(17, sigma_cal[5]); //peam6 sigma 795.95 keV

  fit->FixParameter(4, peak[1]);// peak2 sigma 772.29 kev
  fit->FixParameter(7, peak[2]);// peak3 sigma 783.29 keV  
  fit->FixParameter(10, peak[3]);// peak4 sigma 785.4 keV
  fit->FixParameter(16, peak[5]); //peam6 sigma 795.95 keV
  */
	
  his_temp->Fit("fit","R+L");

	

  for(int ip=0; ip<=19; ip=ip+3){ 
    //	for(int ip=0; ip<=16; ip=ip+3){ 
    if(ip==18){cout<<fit->GetParameter(ip)<<" * "<<fit->GetParError(ip)<<" * "<<fit->GetParameter(ip+1)<<" * "<<fit->GetParError(ip+1)<<endl;}
    else{cout<<fit->GetParameter(ip)<<" * "<<fit->GetParError(ip)<<" * "<<fit->GetParameter(ip+1)<<" * "<<fit->GetParError(ip+1)<< " * " <<fit->GetParameter(ip+2)<<" * "<<fit->GetParError(ip+2)<<endl;}
    
  }
	
  
  char histitle[256];
  char ytitle[256];
  sprintf(histitle,"peaks Fitting");
  sprintf(ytitle,"Counts/%.1f keV", bw);
  cout<<"bw "<<bw<<endl;
  his_temp->SetTitle(histitle);
  his_temp->SetXTitle("Energy[keV]");
  his_temp->SetYTitle(ytitle);
  //	his_temp->SetYTitle("Counts/kev");
  
  c1->cd();
//  his_temp->GetXaxis()->SetRangeUser(750,815);
  his_temp->GetXaxis()->SetRangeUser(peak[0]-50,peak[5]+50);
  //	his_temp->GetYaxis()->SetRangeUser(20,150);//4000
  //	his_temp->GetYaxis()->SetRangeUser(10,80);//8000
  his_temp->Draw();
  TFile rf (resfile,"RECREATE");
  rf.cd();
  //	gfitpara->Write();
  fitpara->Write();
  c1->Write();
  
  rf.Close();
	
}
