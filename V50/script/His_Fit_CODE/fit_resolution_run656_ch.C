#include <iostream>
#include "TMath.h"
#include "TMinuit.h"

void fit_resolution_run656_ch() {

  int const N=8;
  //    double peak[N] = {661.66, 1173.23, 1332.49, 1460};
  //    double peak[N] = {583.19, 727.33, 794.95, 911.20, 1588.19};
  //  double peak[N] = {583.19, 911.20, 1588.19};
  double peak[N] = {238, 338, 583, 727, 860, 911,1588, 2614};
    //  double peak[N] = {391.7, 661.66, 1173.23, 1332.49, 1460, 1836.05};
  //  double ref_adc[N]= {};
  
  int runnum=656;
  int binnum=16000;

  char hisfile[256];
  sprintf(hisfile, "/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN%i/Bin%i/his_%06d.v1.root", runnum, binnum, runnum);

  int ich=3;
  char hisname[256];
  sprintf(hisname,"his%i",ich);
  TFile *hf=new TFile(hisfile);
  TH1D * his_temp = new TH1D("his_temp","",binnum,0,4000);

  //  his_temp = (TH1D*)hf->Get("his_tot3");
  his_temp = (TH1D*)hf->Get(hisname);
  his_temp->Draw();
  his_temp->SetName("his_temp");
  his_temp->SetLineColor(1);


  TCanvas *cpeak[N];
  for(int i=0; i<N; i++){
    cpeak[i] = new TCanvas(Form("cpeak%i",i), Form("cpeak%i",i), 600, 450);
  }
  TCanvas *cres = new TCanvas("cres", "cres", 800, 600);

  TGraphErrors *gr1[16], *gr2[16], *gr3[16], *final_gr[16];

  TF1 * gfit = new TF1("gfit", "gaus");
  //  TF1 * rfit = new TF1("rfit","(sqrt([0]+[1]*x+[2]*x*x))/2.355"); //fwhm(E)=sqrt(a+b*E+c*E*E). sigma = fwhm/2.355 / sigma-E fucntion
    TF1 * rfit = new TF1("rfit","sqrt([0]*[0]+[1]*x+[2]*x*x)"); //fwhm(E)=sqrt(a+b*E+c*E*E). sigma = fwhm/2.355 / sigma-E fucntion


  double sigma[N], sigma_err[N], peak_err[N];
  char a;

  for(int i=0; i<N; i++){
    peak_err[i]=0;
    cpeak[i]->cd();
    //    his_temp->GetXaxis()->SetRange(peak[i]-15, peak[i]+15);
    his_temp->GetXaxis()->SetRange(his_temp->FindBin(peak[i]-15), his_temp->FindBin(peak[i]+15));
    gfit->SetRange(peak[i]-5, peak[i]+5);
    his_temp->Fit("gfit", "R0");
    sigma[i]=gfit->GetParameter(2);
    gfit->SetRange(peak[i]-sigma[i]*3, peak[i]+sigma[i]*3);
    his_temp->Fit("gfit", "R+");
    cpeak[i]->cd();
    his_temp->Draw();
    sigma[i]=gfit->GetParameter(2);
    sigma_err[i]=gfit->GetParError(2);

    cout<<"sigma "<<i<<" "<<sigma[i]<<endl;
    cout<<"sigma err"<<i<<" "<<sigma_err[i]<<endl;
		cpeak[i]->Modified();
		cpeak[i]->Update();
		//                cin>>a;

		
  }
  cres->cd();
  
  TGraph *hres;
  if(ich==11||ich==8){hres = new TGraph(N-1, peak, sigma);    }
  else{hres = new TGraph(N, peak, sigma);}
      //  TGraphErrors *hres = new TGraphErrors(N, peak, sigma, peak_err, sigma_err);

hres->SetTitle(Form("Resolution Fit his%i;Energy (keV);Sigma (keV)", ich));

    hres->SetMarkerStyle(20);
    //    hres->GetXaxis()->SetLimits(0,4000);
    //    rfit->SetParameters(1,0.001, 0,00001);
    rfit->SetParameters(0.5,0.005, 0,00001);
    //        rfit->SetRange(600,1500);
	    rfit->SetRange(100,3000);
    hres->Fit("rfit","R+");                               


    double chi2 = (rfit->GetChisquare())/(rfit->GetNDF());
    cout<<"chi2 "<<chi2<<endl;

    hres->Draw("AP");
    rfit->Draw("SAME");
    cres->Modified();
    cres->Update();

    cout<<hisname<<" * "<<rfit->GetParameter(0)<<" * "<<rfit->GetParError(0)
<<" * "<<rfit->GetParameter(1)<<" * "<<rfit->GetParError(1)
	<<" * "<<rfit->GetParameter(2)<<" * "<<rfit->GetParError(2)<<endl;
}
