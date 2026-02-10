#include <iostream>
#include "TMath.h"
#include "TMinuit.h"

void fit_resolution_241126() {


  char a;
  int const N=4;
    double peak[N] = {661.66, 1173.23, 1332.49, 1460};
    //  double peak[N] = {391.7, 661.66, 1173.23, 1332.49, 1460, 1836.05};
  //  double ref_adc[N]= {};
  
  int runnum=620;
  int binnum=4000;

  char hisfile[256];
  sprintf(hisfile, "/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN%i/Bin%i/his_%06d.root", runnum, binnum, runnum);
  char resfile[256];
  sprintf(resfile,"/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN%i/Bin%i/res_%06d.root", runnum, binnum, runnum);
  TFile rf (resfile, "RECREATE");

  //  char a; 

  TCanvas *cpeak[N];
  for(int i=0; i<N; i++){
    cpeak[i] = new TCanvas(Form("cpeak%i",i), Form("cpeak%i",i), 1200, 1000);
    cpeak[i]->Divide(4,4);
  }

  TCanvas *cres = new TCanvas("cres", "cres", 1200, 1000);
  cres->Divide(4,4);
  
  TFile *hf=new TFile(hisfile);
  TH1D * his_temp = new TH1D("his_temp","",binnum,0,4000);
  
  his_temp = (TH1D*)hf->Get("his_tot3");
  //  test->cd();
  //  his_temp->Draw();
  his_temp->SetName("his_tot");
  his_temp->SetLineColor(1);

  
  TH1D * his_det[14];
  char hisname[256];
  for(int i=0; i<14; i++){
    sprintf(hisname,"his%i",i+1);
    //    his_det[i] = new TH1D(hisname,"",binnum,0,4000);
    his_det[i] = (TH1D*)hf->Get(hisname);
    //  test->cd();
    //  his_det[i]->Draw();



  }
  


  //TGraphErrors *gr1[16], *gr2[16], *gr3[16], *final_gr[16];

  TF1 * gfit = new TF1("gfit", "gaus");
  TF1 * rfit = new TF1("rfit","(sqrt([0]+[1]*x+[2]*x*x))/2.355"); //fwhm(E)=sqrt(a+b*E+c*E*E). sigma = fwhm/2.355 / sigma-E fucntion
  double sigma[15][N], reschi2[15], respa[15][3];


  for(int i=0; i<N; i++){

    //tot
    cpeak[i]->cd(16);
    his_temp->GetXaxis()->SetRange(his_temp->FindBin(peak[i]-15), his_temp->FindBin(peak[i]+15));
    gfit->SetRange(peak[i]-5, peak[i]+5);
    his_temp->Fit("gfit", "R0");
    sigma[14][i]=gfit->GetParameter(2);
    gfit->SetRange(peak[i]-sigma[14][i]*3, peak[i]+sigma[14][i]*3);
    //    his_temp->Fit("gfit", "R+");
    his_temp->Fit("gfit", "R+");
    cpeak[i]->cd(16);
    his_temp->Draw();
    sigma[14][i]=gfit->GetParameter(2);


    //each det.
    for(int j=0; j<14; j++){

      cpeak[i]->cd(j+1);
      his_det[j]->GetXaxis()->SetRange(his_det[j]->FindBin(peak[i]-15), his_det[j]->FindBin(peak[i]+15));
      gfit->SetRange(peak[i]-5, peak[i]+5);
      his_det[j]->Fit("gfit", "R0");
      sigma[j][i]=gfit->GetParameter(2);
      gfit->SetRange(peak[i]-sigma[j][i]*3, peak[i]+sigma[j][i]*3);
      his_det[j]->Fit("gfit", "R+");
      //      his_det[j]->Fit("gfit", "R+");
      //      cpeak[i]->cd(j+1);
      his_det[j]->Draw();
      sigma[j][i]=gfit->GetParameter(2);
    }

    //    cout<<"sigma "<<i<<" "<<sigma[0][i]<<endl;
		cpeak[i]->Modified();
		cpeak[i]->Update();

  rf.cd();
  cpeak[i]->Write();	                
  //cin>>a;
  
  }

  //tot
  cres->cd(16);
  TGraph *hres = new TGraph(N, peak, sigma[14]);
  hres->SetTitle("11det res Fit;Energy (keV);Sigma (keV)");
  hres->SetMarkerStyle(20);
  rfit->SetParameters(0.5,0.005, 0,00001);
  rfit->SetRange(600,1500);
  hres->Fit("rfit","R+");                               

  respa[14][0]=rfit->GetParameter(0);
  respa[14][1]=rfit->GetParameter(1);
  respa[14][2]=rfit->GetParameter(2);
  reschi2[14] = (rfit->GetChisquare())/(rfit->GetNDF());

  hres->Draw("AP");
  //  rfit->Draw("SAME");

  //tot
  TGraph *hres_det[14];
  for(int i=0; i<14; i++){
    cres->cd(i+1);
    hres_det[i] = new TGraph(N, peak, sigma[i]);
    hres_det[i]->SetTitle(Form("res%i Fit;Energy (keV);Sigma (keV)",i+1));
    hres_det[i]->SetMarkerStyle(20);
    rfit->SetParameters(0.5,0.005, 0,0001);
    rfit->SetRange(600,1500);
    hres_det[i]->Fit("rfit","R+");                               
  respa[i][0]=rfit->GetParameter(0);
  respa[i][1]=rfit->GetParameter(1);
  respa[i][2]=rfit->GetParameter(2);

    reschi2[i] = (rfit->GetChisquare())/(rfit->GetNDF());

  hres_det[i]->Draw("AP");
  //  rfit->Draw("SAME");

  }

  cres->Modified();
  cres->Update();


  double tpeak, fitsig, calsig, trespa[3];
  int hisnum;
  TTree * fittree = new TTree("fittree","peak fitting");
  fittree->Branch("hisnum",&hisnum,"hisnum/I");
  fittree->Branch("peak",&tpeak,"peak/D");
  fittree->Branch("fitsig",&fitsig,"fitsig/D");
  fittree->Branch("calsig",&calsig,"calsig/D");

  TTree * restree = new TTree("restree","resolution fn parameter");
  restree->Branch("hisnum",&hisnum,"hisnum/I");
  restree->Branch("respa",trespa,"respa[3]/D");


  for(int i=0; i<15; i++){
      hisnum=i+1;
    for(int j=0; j<N; j++){
      tpeak=peak[j];
      //      cout<<"tpeak"<<tpeak<<endl;
      if(i==14){hisnum=999;}
      fitsig=sigma[i][j];
      calsig=sqrt(respa[i][0]+respa[i][1]*tpeak+respa[i][2]*tpeak*tpeak)/2.355;
      fittree->Fill();
    }
    
    trespa[0]=respa[i][0];
    trespa[1]=respa[i][1];
    trespa[2]=respa[i][2];
    restree->Fill();
  }


  rf.cd();
  fittree->Write();
  restree->Write();
  //  for(int i=0; i<N; i++){  cpeak[i]->Write();}
  cres->Write();
  rf.Close();

}
