#include <iostream>
#include "TMath.h"
#include "TMinuit.h"

void fit_resolution_run624_v1() {

  int const N=7;
  double peak[N] = {609,911,1120, 1174, 1332, 1460,1553};
  //  double peak[N] = {609};
    //  double peak[N] = {391.7, 661.66, 1173.23, 1332.49, 1460, 1836.05};
  //  double ref_adc[N]= {};
  
  int runnum=624;
  int binnum=16000;

  char hisfile[256];
  sprintf(hisfile, "/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN%i/Bin%i/his_%06d.v1.root", runnum, binnum, runnum);



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
        

  TString term1 = "([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/(2*[2]*[2])))";

//  TF1 * gfit = new TF1("gfit", "gaus");
 TF1 * gfit = new TF1("gfit", term1);
        gfit->SetParName(0,"Area");
        gfit->SetParName(1,"Mean");
        gfit->SetParName(2,"Sigma");




  TF1 * rfit = new TF1("rfit","(sqrt([0]+[1]*x+[2]*x*x))/2.355"); //fwhm(E)=sqrt(a+b*E+c*E*E). sigma = fwhm/2.355 / sigma-E fucntion
  double sigma[N];
  char a;

  for(int i=0; i<N; i++){

    cpeak[i]->cd();
    //    his_temp->GetXaxis()->SetRange(peak[i]-15, peak[i]+15);
    his_temp->GetXaxis()->SetRange(his_temp->FindBin(peak[i]-15), his_temp->FindBin(peak[i]+15));
    gfit->SetRange(peak[i]-5, peak[i]+5);
gfit->ReleaseParameter(2);           
  gfit->SetParameter(0,1);
  gfit->SetParameter(1,peak[i]);
  gfit->SetParameter(2,0.9);


    if(peak[i]==351){
      gfit->SetParLimits(2,0.5,1);
    }

    if(peak[i]==911){
      gfit->SetParLimits(2,0.7,1);
    }
      cout<<"test"<<endl;
      //          }    
his_temp->Fit("gfit", "R0");
    sigma[i]=gfit->GetParameter(2);
    gfit->SetRange(peak[i]-sigma[i]*3, peak[i]+sigma[i]*3);
    his_temp->Fit("gfit", "R+");
    cpeak[i]->cd();
    his_temp->Draw();
    sigma[i]=gfit->GetParameter(2);


    cout<<"sigma "<<i<<" "<<sigma[i]<<endl;
		cpeak[i]->Modified();
		cpeak[i]->Update();
		//                cin>>a;

		
  }
  cres->cd();
    TGraph *hres = new TGraph(N, peak, sigma);
    hres->SetTitle("Resolution Fit;Energy (keV);Sigma (keV)");
    hres->SetMarkerStyle(20);
    //    hres->GetXaxis()->SetLimits(0,4000);
    //    rfit->SetParameters(1,0.001, 0,00001);
    rfit->SetParameters(0.5,0.005, 0,00001);
    rfit->SetRange(100,3000);
    hres->Fit("rfit","R+");                               


    double chi2 = (rfit->GetChisquare())/(rfit->GetNDF());
    cout<<"chi2 "<<chi2<<endl;

    hres->Draw("AP");
    rfit->Draw("SAME");
		cres->Modified();
		cres->Update();
}
