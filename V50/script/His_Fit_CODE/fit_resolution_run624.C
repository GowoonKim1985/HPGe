#include <iostream>
#include "TMath.h"
#include "TMinuit.h"

void fit_resolution_run624() {

  int const N=4;
    double peak[N] = {661.66, 1173.23, 1332.49, 1460};
    //  double peak[N] = {391.7, 661.66, 1173.23, 1332.49, 1460, 1836.05};
  //  double ref_adc[N]= {};
  
    double sigfit[N];
    double sigcal[N];

  int runnum=620;
  int binnum=16000;
  double bw = 4000/binnum;

  char hisfile[256];
  sprintf(hisfile, "/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN%i/Bin%i/his_%06d.root", runnum, binnum, runnum);



  TFile *hf=new TFile(hisfile);
  TH1D * his_temp = new TH1D("his_temp","",binnum,0,4000);

  his_temp = (TH1D*)hf->Get("his_tot3");
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
  TF1 * rfit = new TF1("rfit","(sqrt([0]+[1]*x+[2]*x*x))/2.355"); //fwhm(E)=sqrt(a+b*E+c*E*E). sigma = fwhm/2.355 / sigma-E fucntion
  double sigma[N];
  char a;
  int peakbin[N];

  for(int i=0; i<N; i++){
    cpeak[i]->cd();

    //    his_temp->GetXaxis()->SetRange(peak[i]-15, peak[i]+15);
    //    gfit->SetRange(peak[i]-5, peak[i]+5);
    cout << "peakbin "<< his_temp->FindBin(peak[i]-15)<<endl;
    his_temp->GetXaxis()->SetRange(his_temp->FindBin(peak[i]-15), his_temp->FindBin(peak[i]+15));
    gfit->SetRange(peak[i]-5, peak[i]+5);

    his_temp->Fit("gfit", "R");
    sigma[i]=gfit->GetParameter(2);
    cout << "sigma "<<sigma[i]<<endl;
    gfit->SetRange(peak[i]-sigma[i]*3, peak[i]+sigma[i]*3);
    his_temp->Fit("gfit", "R+");
    cpeak[i]->cd();
    his_temp->Draw();
    sigma[i]=gfit->GetParameter(2);
    sigfit[i]=sigma[i];

    cout<<"sigma "<<i<<" "<<sigma[i]<<endl;
		cpeak[i]->Modified();
		cpeak[i]->Update();
		//                cin>>a;

		
  }
  cres->cd();
  double respa[3];
  TGraph *hres = new TGraph(N, peak, sigma);
  hres->SetTitle("Resolution Fit;Energy (keV);Sigma (keV)");
  hres->SetMarkerStyle(20);
    //    hres->GetXaxis()->SetLimits(0,4000);
    //    rfit->SetParameters(1,0.001, 0,00001);
  rfit->SetParameters(0.5,0.005, 0,00001);
  rfit->SetRange(600,1500);
  hres->Fit("rfit","R+");                               
  respa[0]=rfit->GetParameter(0);
  respa[1]=rfit->GetParameter(1);
  respa[2]=rfit->GetParameter(2);

  double chi2 = (rfit->GetChisquare())/(rfit->GetNDF());
  cout<<"chi2 "<<chi2<<endl;

  hres->Draw("AP");
  rfit->Draw("SAME");
  cres->Modified();
  cres->Update();

  cout<<"peak * sigfit * sigcal"<<endl;
  for(int i=0; i<N; i++){
    cout<<peak[i]<<" * "<<sigfit[i]<<" * "<<(sqrt(respa[0]+respa[1]*peak[i]+respa[2]*peak[i]*peak[i]))/2.355<<endl;
  }
		


}
