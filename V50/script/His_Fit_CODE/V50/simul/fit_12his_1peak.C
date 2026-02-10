#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooPlot.h"
#include "RooFitResult.h"

using namespace RooFit;

void fit_12his_1peak()

{

  // 1)histogram load    

  string workdir("/data/HPGe/USERS/kkw/DAQ/ANA125");

  char hisfile[256];
  char rawname[256];
  char resfile[256];//fit result
  char runnumber3[256];
  char runnumber6[256];

  int runnum = 650;
  int expnum=2;
  int binnum = 16000;// 1bin = 1kev
  double bw = 4000./((double)binnum); //bin width, keV

  sprintf(runnumber6,"%06d",runnum);
  sprintf(runnumber3,"%03d",runnum);

  sprintf(hisfile,"%s/result/V50/Bin%i/his_000650.v1.root", workdir.c_str(), binnum);
  sprintf(resfile,"%s/result/V50/Bin%i/sim_000650.v1.root", workdir.c_str(), binnum);
  TFile *hf = new TFile(hisfile);

  TH1* h[14]; char hisname[14];
  for(int i=0; i<14; i++){
    sprintf(hisname,"his%i",i+1);
    h[i] = (TH1D*)hf->Get(hisname);
  }

  //create roodatahist
  
  //peak, sigma pre set
  double peak = 1553.77;//1 peak
  double respa[14][3]={
    {0.363128,      0.000667908,    -8.40E-08},
    {0, 0, 0},
    {0.445986,      0.000680498,    -9.62E-08},
    {0, 0 ,0},
    {0.436251,      0.000626869,    -7.70E-08},
    {0.402777,      0.000545669,    -4.38E-08},
    {0.582696,      0.000821223,    -1.54E-07},
    {0.495129,      0.000449309,    8.03E-09},
    {0.511883,      0.000459608,    -2.06E-08},
    {0, 0, 0},
    {0.482192,      0.000504369,    -2.71E-08},
    {0.432383,      0.000486415,    -3.44E-08},
    {0.429384,      0.000619621,    -8.22E-08},
    {0.448204,      0.000548424,    -4.75E-08}
  };

  /*
  double respa[3];
  respa[0] = 0.446;
  respa[1] = 6.06E-4;
  respa[2] = -6.31E-8;
  */

  //  double sigma_cal = std::sqrt(respa[0]*respa[0] + respa[1]*peak + respa[2]*peak*peak);//sigma calculation
  double sigma_cal[14];
  for(int i=0; i<14; i++){
    sigma_cal[i] = std::sqrt(respa[i][0]*respa[i][0] + respa[i][1]*peak + respa[i][2]*peak*peak);
  }

  RooRealVar E("E", "Energy (keV)", peak-50, peak+50); // E range 

  std::vector<RooDataHist*> data(14);
  for(int i=0; i<14; i++){
    TString name  = Form("data%d", i+1);
    TString title = Form("dataset %d", i+1);
    data[i] = new RooDataHist(name.Data(), title.Data(), E, Import(*h[i]));
  }

  //BG model (p1)

  std::vector<RooPolynomial*> bkg(14);
  std::vector<RooRealVar*> p0(14);
  std::vector<RooRealVar*> p1(14);

  for(int i=0; i<14; i++){
  //  RooRealVar p0_1("p0_1", "const1", 0.0);  // no range → free
  //  RooRealVar p1_1("p1_1", "slope1", 0.0);  // no range → free
    TString p0name = Form("p0_%d", i+1);
    TString p1name = Form("p1_%d", i+1);
    TString bkgname = Form("bkg%d", i+1);
    TString bkgtitle = Form("background%d", i+1);

    p0[i] = new RooRealVar(p0name, Form("const%d", i+1), 0.0);
    p1[i] = new RooRealVar(p1name, Form("slope%d", i+1), 0.0);

    //    bkg[i] = new RooPolynomial(bkgname, bkgtitle, E, RooArgList(*p0[i], *p1[i]));
    bkg[i] = new RooPolynomial(bkgname, bkgtitle, E, RooArgList(*p0[i]));
  }
  // gaus model

  RooRealVar mean("mean", "Peak mean", peak, peak-10, peak+10);

  std::vector<RooRealVar*> sigma(14);
  std::vector<RooGaussian*> gaus(14);

  for (int i = 0; i < 14; i++) {
    TString sname = Form("sigma%d", i+1);
    TString gname = Form("gaus%d", i+1);
    TString gtitle = Form("gaussian%d", i+1);
    
    sigma[i] = new RooRealVar(sname, Form("Sigma (hist%d)", i+1), sigma_cal[i], 0.7*sigma_cal[i], 1.3*sigma_cal[i]);
    gaus[i] = new RooGaussian(gname, gtitle, E, mean, *sigma[i]);
  }

  //extended yields (area parameters)

  double n1[14];
  for(int i=0; i<14; i++){n1[i]=h[i]->Integral(h[i]->FindBin(1500), h[i]->FindBin(1600));}

  // fitting model

std::vector<RooRealVar*> Nsig(14);
std::vector<RooRealVar*> Nbkg(14);
std::vector<RooAddPdf*>  model(14);

for (int i = 0; i < 14; i++) {

    TString nSname = Form("Nsig%d", i+1);
    TString nBname = Form("Nbkg%d", i+1);
    TString mname  = Form("model%d", i+1);
    TString mtitle = Form("sig+bkg%d", i+1);

    Nsig[i] = new RooRealVar(nSname, Form("Signal yield (hist%d)", i+1),
                             0.2 * n1[i], 0, 10.0 * n1[i]);
    Nbkg[i] = new RooRealVar(nBname, Form("Bkg yield (hist%d)", i+1),
                             0.8 * n1[i], 0, 10.0 * n1[i]);

    model[i] = new RooAddPdf(mname, mtitle,
                             RooArgList(*gaus[i], *bkg[i]),
                             RooArgList(*Nsig[i], *Nbkg[i]));
 }
  // simultaneous setting
 RooCategory sample("sample", "sample");

  sample.defineType("hist1");
  sample.defineType("hist3");
  sample.defineType("hist5");
  sample.defineType("hist6");
  sample.defineType("hist7");
  sample.defineType("hist8");
  sample.defineType("hist9");
  sample.defineType("hist11");
  sample.defineType("hist12");
  sample.defineType("hist13");
  sample.defineType("hist14");


 // do not use 2,4,11 det. using 11dets.
 RooDataHist combData("combData", "combined data", E,
		      Index(sample),
		      Import("hist1", *data[0]),
		      Import("hist3", *data[2]),
		      Import("hist5", *data[4]),
		      Import("hist6", *data[5]),
		      Import("hist7", *data[6]),
		      Import("hist8", *data[7]),
		      Import("hist9", *data[8]),
		      Import("hist11", *data[10]),
		      Import("hist12", *data[11]),
		      Import("hist13", *data[12]),
		      Import("hist14", *data[13])
		      );
 
 RooSimultaneous simPdf("simPdf", "simultaneous pdf", sample);
 for (int i = 0; i < 14; i++) {
   if(i!=1&&i!=3&&i!=9){
     simPdf.addPdf(*model[i], Form("hist%d", i+1));
   }
 }


// fitting
 RooFitResult *fitres = simPdf.fitTo(combData, Save(), Extended(true), PrintLevel(-1));

  // plot drawing

 TCanvas *tc; TCanvas *bc;

 tc = new TCanvas("tc", "Top Simultaneous Fit (1553.8 keV)", 1500, 1000);
 tc->Divide(4,2);
 bc = new TCanvas("bc", "Bot Simultaneous Fit (1553.8 keV)", 1500, 1000);
 bc->Divide(4,2);

 for(int i=0; i<7; i++){
   if(i!=1&&i!=3){
   tc->cd(1+i);
   TString hname = Form("hist%d", i+1);
   RooPlot *frame = E.frame(Title(Form("Dataset %d", i+1)));
   combData.plotOn(frame, Cut(Form("sample==sample::%s", hname.Data())));
   simPdf.plotOn(frame, Slice(sample, hname), ProjWData(sample, combData));
   frame->Draw();
   }
 }
 for(int i=8; i<14; i++){
   if(i!=9){
   bc->cd(i-7);
   TString hname = Form("hist%d", i+1);
   RooPlot *frame = E.frame(Title(Form("Dataset %d", i+1)));
   combData.plotOn(frame, Cut(Form("sample==sample::%s", hname.Data())));
   simPdf.plotOn(frame, Slice(sample, hname), ProjWData(sample, combData));
   frame->Draw();
   }
 }
  
  //fitting result 
  fitres->Print("v");

  /*
  // ----- totals (after fit) -----
  auto &pars = fitres->floatParsFinal();   // order of floating params
  int i1 = pars.index("Nsig1");
  int i2 = pars.index("Nsig2");
  int j1 = pars.index("Nbkg1");
  int j2 = pars.index("Nbkg2");

  const TMatrixDSym& C = fitres->covarianceMatrix();

  // signal total and uncertainty (with covariance)
  double Nsig_tot = Nsig1.getVal() + Nsig2.getVal();
  double var_sig  = std::pow(Nsig1.getError(),2) + std::pow(Nsig2.getError(),2)
    + 2.0 * C(i1,i2); // sqrt((e1)^2 + (e2)^2 + 2Cov)
  double err_sig  = (var_sig>0) ? std::sqrt(var_sig) : 0.0; // true->sqrt(var_sig), false->0

  // background total and uncertainty
  double Nbkg_tot = Nbkg1.getVal() + Nbkg2.getVal();
  double var_bkg  = std::pow(Nbkg1.getError(),2) + std::pow(Nbkg2.getError(),2)
    + 2.0 * C(j1,j2);
  double err_bkg  = (var_bkg>0) ? std::sqrt(var_bkg) : 0.0;
  
  // grand total (signal+background)
  double N_all = Nsig_tot + Nbkg_tot;
  
  printf("\n[Totals]\n");
  printf("  Nsig_total = %.3f ± %.3f\n", Nsig_tot, err_sig);
  printf("  Nbkg_total = %.3f ± %.3f\n", Nbkg_tot, err_bkg);
  printf("  N_all      = %.3f\n\n",      N_all);
  */
}
