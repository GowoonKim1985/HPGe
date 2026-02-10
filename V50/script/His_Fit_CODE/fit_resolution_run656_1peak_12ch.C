#include <iostream>
#include "TMath.h"
#include "TMinuit.h"

void fit_resolution_run656_ch() {

  int const N=17;
  //    double peak[N] = {661.66, 1173.23, 1332.49, 1460};
  //  double peak[N] = {583.19, 727.33, 794.95, 911.20, 1588.19};
  //  double peak[N] = {583.19, 911.20, 1588.19};
  //  double peak[N] = {238, 338, 583, 727, 860, 911,1588, 2614};
  double peak = 727;

    //  double peak[N] = {391.7, 661.66, 1173.23, 1332.49, 1460, 1836.05};
  //  double ref_adc[N]= {};
  
  int runnum=656;
  int binnum=16000;

  char hisfile[256];
  sprintf(hisfile, "/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN%i/Bin%i/his_%06d.v1.root", runnum, binnum, runnum);


  TCanvas *cpeak[16];
  for(int i=1; i<N; i++){
    cpeak[i-1] = new TCanvas(Form("ch%i",i), Form("ch%i",i), 600, 450);
  }


  TFile *hf=new TFile(hisfile);
  TH1D * his_temp[17];
  char hisname[256];
  for(int i=1; i<17; i++){
    sprintf(hisname,"his%i",i);
    his_temp[i] = new TH1D("his_temp","",binnum,0,4000);
    his_temp[i] = (TH1D*)hf->Get(hisname);
    his_temp[i]->SetLineColor(1);
    cpeak[i-1]->cd();
    his_temp[i]->Draw();
  }


  //  his_temp->Draw();
  //  his_temp->SetName("his_temp");




  //  TCanvas *cres = new TCanvas("cres", "cres", 800, 600);

  TGraphErrors *gr1[16], *gr2[16], *gr3[16], *final_gr[16];

  TF1 * gfit = new TF1("gfit", "gaus");
  //  TF1 * rfit = new TF1("rfit","(sqrt([0]+[1]*x+[2]*x*x))/2.355"); //fwhm(E)=sqrt(a+b*E+c*E*E). sigma = fwhm/2.355 / sigma-E fucntion
  TF1 * rfit = new TF1("rfit","sqrt([0]*[0]+[1]*x+[2]*x*x)"); //fwhm(E)=sqrt(a+b*E+c*E*E). sigma = fwhm/2.355 / sigma-E fucntion


  double sigma[N], sigma_err[N], mean[N], mean_err[N];
  char a;

  for(int i=1; i<N; i++){
    //    peak_err[i]=0;
    cpeak[i-1]->cd();
    his_temp[i]->GetXaxis()->SetRange(his_temp[i]->FindBin(peak-15), his_temp[i]->FindBin(peak+15));
    gfit->SetRange(peak-5, peak+5);
    his_temp[i]->Fit("gfit", "R0");
    sigma[i]=gfit->GetParameter(2);
    gfit->SetRange(peak-sigma[i]*3, peak+sigma[i]*3);
    his_temp[i]->Fit("gfit", "R+");
    cpeak[i-1]->cd();
    his_temp[i]->Draw();
    sigma[i]=gfit->GetParameter(2);
    sigma_err[i]=gfit->GetParError(2);
    mean[i]=gfit->GetParameter(1);
    mean_err[i]=gfit->GetParError(1);

    cout<<"sigma "<<i<<" "<<sigma[i]<<endl;
    cout<<"sigma err"<<i<<" "<<sigma_err[i]<<endl;
		cpeak[i-1]->Modified();
		cpeak[i-1]->Update();
		//                cin>>a;

		
  }
    cout<<"ch * mean * mean_er * sigma * sigma_er"<<endl;
  for(int i=1; i<17; i++){

    cout<<i<<" * "<<mean[i]<<" * "<<mean_err[i]<<" * "<<sigma[i]<<" * "<<sigma_err[i]<<endl;
    
  }

  /*
  cres->cd();
      TGraph *hres = new TGraph(N, peak, sigma);
      //  TGraphErrors *hres = new TGraphErrors(N, peak, sigma, peak_err, sigma_err);


    hres->SetTitle("Resolution Fit;Energy (keV);Sigma (keV)");
    hres->SetMarkerStyle(20);
    //    hres->GetXaxis()->SetLimits(0,4000);
    //    rfit->SetParameters(1,0.001, 0,00001);
    rfit->SetParameters(0.5,0.005, 0,00001);
    //    rfit->SetRange(600,1500);
    rfit->SetRange(100,3000);
    hres->Fit("rfit","R+");                               


    double chi2 = (rfit->GetChisquare())/(rfit->GetNDF());
    cout<<"chi2 "<<chi2<<endl;

    hres->Draw("AP");
    rfit->Draw("SAME");
		cres->Modified();
		cres->Update();
  */
}
