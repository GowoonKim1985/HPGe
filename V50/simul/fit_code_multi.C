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

void fit_code_multi()

{

  // 1)histogram load    

  string workdir("/data/HPGe/USERS/kkw/DAQ/ANA125");

  char hisfile[256];
  char rawname[256];
  char resfile[256];//fit result
  char runnumber3[256];
  char runnumber6[256];
  const int N=9;
  int runnum = 650;
  int expnum=2;
  int binnum = 16000;// 1bin = 1kev
  double bw = 4000./((double)binnum); //bin width, keV

  sprintf(runnumber6,"%06d",runnum);
  sprintf(runnumber3,"%03d",runnum);

  sprintf(hisfile,"%s/result/V50/Bin%i/his_000650.v1.root", workdir.c_str(), binnum);
  sprintf(resfile,"%s/result/V50/Bin%i/sim_000650.v1.root", workdir.c_str(), binnum);
  TFile *hf = new TFile(hisfile);
  TH1* h[N]; char hisname[235];

  for(int i=0; i<N; i++){
    sprintf(hisname,"his%i",i+1);
  h[i] = (TH1D*)hf->Get(hisname);
  }

  //  h[0] = (TH1D*)hf->Get("his1");

  //create roodatahist
  double peak = 1553.77;//1 peak
  RooRealVar E("E", "Energy (keV)", peak-10, peak+10); // E range
  E.setRange("fitRange", peak-10, peak+10);
  //  RooDataHist data1("data1", "dataset 1", E, Import(*h1));


  //  RooRealVar E("E", "Energy (keV)", 1500, 1600); // E range 
  /*
  RooDataHist data1("data1", "dataset 1", E, Import(*h[0]));
  RooDataHist data2("data2", "dataset 2", E, Import(*h[1]));  
  RooDataHist data3("data3", "dataset 3", E, Import(*h[2]));
  */
RooDataHist* data[N];
for (int i = 0; i < N; i++) {
    data[i] = new RooDataHist(Form("data%d", i+1),
                              Form("dataset %d", i+1),
                              E,
                              Import(*h[i]));
}



  //peak, sigma pre set


  double respa[14][3]={
    {0.363128,      0.000667908,    -8.40E-08},
    {0.4,      0.00066,    -8E-08},
    {0.445986,      0.000680498,    -9.62E-08},
    {0.4,      0.00066,    -8E-08},
    {0.436251,      0.000626869,    -7.70E-08},
    {0.402777,      0.000545669,    -4.38E-08},
    {0.582696,      0.000821223,    -1.54E-07},
    {0.495129,      0.000449309,    8.03E-09},
    {0.511883,      0.000459608,    -2.06E-08},
    {0.4,      0.00066,    -8E-08},
    {0.482192,      0.000504369,    -2.71E-08},
    {0.432383,      0.000486415,    -3.44E-08},
    {0.429384,      0.000619621,    -8.22E-08},
    {0.448204,      0.000548424,    -4.75E-08}
  };

  double sigma_cal[14];
  for(int i=0; i<14; i++){
sigma_cal[i] = std::sqrt(respa[i][0]*respa[i][0] + respa[i][1]*peak + respa[i][2]*peak*peak);//sigma calculation

}
 
  //BG model (p1)
  
  RooRealVar* p0[N];
  RooRealVar* p1[N];
  //  RooPolynomial* bkg[N];
  RooChebychev* bkg[N];

  /*
  for(int i=0; i<N; i++){
    p0[i] = new RooRealVar(Form("p0_%d", i+1), Form("const%d", i+1), 1, 0, 2);
    p1[i] = new RooRealVar(Form("p1_%d", i+1), Form("slope%d", i+1), 0.0, -1e-4, 1e-4);
    //    bkg[i] = new RooPolynomial(Form("bkg%d", i+1), Form("background%d", i+1), E, RooArgList(*p0[i], *p1[i]));
    bkg[i] = new RooPolynomial(Form("bkg%d", i+1), Form("background%d", i+1), E, RooArgList(*p0[i]));
  }
  */

 for(int i=0; i<N; i++){
    p0[i] = new RooRealVar(Form("p0_%d", i+1), Form("const%d", i+1), 0, -1, 1);
    p1[i] = new RooRealVar(Form("p1_%d", i+1), Form("slope%d", i+1), 0, -1, 1);
    //    bkg[i] = new RooChebychev(Form("bkg%d", i+1), Form("background%d", i+1), E, RooArgList(*p0[i]));
    bkg[i] = new RooChebychev(Form("bkg%d", i+1), Form("background%d", i+1), E, RooArgList(*p0[i], *p1[i]));
  }
    // f(E) = 1 + p0*T1(E) + p1*T2(E)
    // f(E) = 1 + p0*T1(E)

  // gaus model
  
  RooRealVar mean("mean", "Peak mean", peak, peak-1, peak+1);

  RooRealVar *sigma[N];
  RooGaussian *gaus[N];




  for(int i=0; i<N; i++){

    sigma[i] = new RooRealVar(Form("sigma%i",i+1), Form("Sigma (hist%i)",i+1), sigma_cal[i], 0.7*sigma_cal[i], 1.3*sigma_cal[i]);
    gaus[i] = new RooGaussian(Form("gaus%i",i+1), Form("gaussian%i",i+1), E, mean, *sigma[i]);
  }

  //extended yields (area parameters)

  double totN[N];
  RooRealVar* Nsig[N];
  RooRealVar* Nbkg[N];
  RooAddPdf* model[N];

  for(int i=0; i<N; i++){
    totN[i] = h[i]->Integral(h[i]->FindBin(peak-10), h[i]->FindBin(peak+10));
    Nsig[i] = new RooRealVar(Form("Nsig%i",i+1), Form("Signal yield (hist%i)",i+1), 0.2*totN[i], 0, 10.0*totN[i]);
    Nbkg[i] = new RooRealVar(Form("Nbkg%i",i+1), Form("Bkg yield (hist%i)",i+1),   0.8*totN[i], 0, 10.0*totN[i]);
    model[i] = new RooAddPdf(Form("model%i",i+1), Form("sig+bkg%i",i+1), RooArgList(*gaus[i], *bkg[i]), RooArgList(*Nsig[i], *Nbkg[i]));

}

  // fitting model

  // simultaneous setting

  RooCategory sample("sample", "sample");
  char samname[246];
  for(int i=0; i<N; i++){  
    if(i!=1&&i!=3&&i!=9){
      sprintf(samname,"hist%i",i+1);
      sample.defineType(samname);
    }
  }
  //  RooDataHist combData("combData", "combined data", E, Index(sample), Import("hist1", *data[0]), Import("hist2", *data[1]), Import("hist3", *data[2]));
  




  //  RooDataHist combData("combData", "combined data", E, Index(sample), Import("hist1", *data[0]), Import("hist3", *data[2]), Import("hist5", *data[4]), Import("hist6", *data[5]), Import("hist7", *data[6]), Import("hist8", *data[7]), Import("hist9", *data[8]), Import("hist11", *data[10]));
  //  RooDataHist combData("combData", "combined data", E, Index(sample), Import("hist1", *data[0]), Import("hist3", *data[2]), Import("hist5", *data[4]), Import("hist6", *data[5]), Import("hist7", *data[6]), Import("hist8", *data[7]), Import("hist9", *data[8]), Import("hist11", *data[10]), Import("hist12", *data[11]), Import("hist13", *data[12]), Import("hist14", *data[13]));
 

RooDataHist combData("combData", "combined data", E, Index(sample));
sample.setLabel("hist1"); combData.add(*data[0]);
sample.setLabel("hist3"); combData.add(*data[2]);
sample.setLabel("hist5"); combData.add(*data[4]);
sample.setLabel("hist6"); combData.add(*data[5]);
sample.setLabel("hist7"); combData.add(*data[6]);
sample.setLabel("hist8"); combData.add(*data[7]);
sample.setLabel("hist9"); combData.add(*data[8]);


  /*
RooDataHist combData("combData", "combined data", E, Index(sample));
combData.add(*data[0],  "hist1");
combData.add(*data[2],  "hist3");
combData.add(*data[4],  "hist5");
combData.add(*data[5],  "hist6");
combData.add(*data[6],  "hist7");
combData.add(*data[7],  "hist8");
combData.add(*data[8],  "hist9");
  */ 
  
  RooSimultaneous simPdf("simPdf", "simultaneous pdf", sample);
  //  simPdf.addPdf(model1, "hist1");
  //  simPdf.addPdf(model2, "hist2");
  //  simPdf.addPdf(model2, "hist3");
  for(int i=0; i<N; i++){
    if(i!=1&&i!=3&&i!=9){
    simPdf.addPdf(*model[i], Form("hist%i",i+1));
    }
  }
  /*
  simPdf.addPdf(*model[0], "hist1");
  simPdf.addPdf(*model[1], "hist2");
  simPdf.addPdf(*model[2], "hist3");
  */
  // fitting
  RooFitResult *fitres = simPdf.fitTo(combData, Save(), Range("fitRange"), Extended(true), PrintLevel(-1));

  // plot drawing

  TCanvas *tc = new TCanvas("tc", "Top array", 1800, 1000);
  tc->Divide(4,2);
  for(int i=0; i<7; i++){
    tc->cd(i+1);
    RooPlot *frame = E.frame(Title(Form("Histogram %i",i+1)));
    combData.plotOn(frame, Cut(Form("sample==sample::hist%i",i+1)));
    
    simPdf.plotOn(frame, Slice(sample, Form("hist%i",i+1)), ProjWData(sample, combData), Range("fitRange"), NormRange("fitRange"));
    frame->Draw();
  }

  TCanvas *bc = new TCanvas("bc", "Bot array", 1800, 1000);
  bc->Divide(4,2);
  for(int i=7; i<N; i++){
    bc->cd(i-6);
    RooPlot *frame = E.frame(Title(Form("Histogram %i",i+1)));
    combData.plotOn(frame, Cut(Form("sample==sample::hist%i",i+1)));    
    simPdf.plotOn(frame, Slice(sample, Form("hist%i",i+1)), ProjWData(sample, combData), Range("fitRange"), NormRange("fitRange"));
    frame->Draw();
  }

  /*
  TCanvas *c = new TCanvas("c", "Simultaneous Fit (1553.8 keV)", 1200, 500);
  c->Divide(3,1);

  c->cd(1);
  RooPlot *frame1 = E.frame(Title("Histogram 1"));
  combData.plotOn(frame1, Cut("sample==sample::hist1"));
  simPdf.plotOn(frame1, Slice(sample, "hist1"), ProjWData(sample, combData));
  frame1->Draw();
  
  c->cd(2);
  RooPlot *frame2 = E.frame(Title("Histogram 2"));
  combData.plotOn(frame2, Cut("sample==sample::hist2"));
  simPdf.plotOn(frame2, Slice(sample, "hist2"), ProjWData(sample, combData));
  frame2->Draw();

  c->cd(3);
  RooPlot *frame3 = E.frame(Title("Histogram 3"));
  combData.plotOn(frame3, Cut("sample==sample::hist3"));
  simPdf.plotOn(frame3, Slice(sample, "hist3"), ProjWData(sample, combData));
  frame3->Draw();
  */
  //fitting result 
  fitres->Print("v");

  
  // ----- totals (after fit) -----
  auto &pars = fitres->floatParsFinal();   // order of floating params
  int cnt[N]; int cnt_bg[N];
  char signame[256], bkgname[256];
  for(int i=0; i<N; i++){
    if(i!=1&&i!=3&&i!=9){
      sprintf(signame,"Nsig%i",i+1);
      sprintf(bkgname,"Nbkg%i",i+1);
      cnt[i] = pars.index(signame);
      cnt_bg[i] = pars.index(bkgname);
    }
  }
  /*
  cnt[0] = pars.index("Nsig1");
  cnt[1] = pars.index("Nsig2");
  cnt_bg[0] = pars.index("Nbkg1");
  cnt_bg[1] = pars.index("Nbkg2");
  */
  
  const TMatrixDSym& C = fitres->covarianceMatrix();
  double Nsig_tot =0;
  double var_sig =0;
  double cov_sig = 0;

  for(int i=0; i<N; i++){
    if(i!=1&&i!=3&&i!=9){
      Nsig_tot += Nsig[i]->getVal();
      var_sig  += std::pow(Nsig[i]->getError(),2);
      for(int j=i+1; j<N; j++){
	if(j!=1&&j!=3&&j!=9){
	  cov_sig += 2.0 * C(cnt[i], cnt[j]);
	}
      }
    }
  }
    var_sig = var_sig+cov_sig;
    double err_sig  = (var_sig>0) ? std::sqrt(var_sig) : 0.0; // true->sqrt(var_sig), false->
  
  // background total and uncertainty

    double Nbkg_tot=0;
    double var_bkg=0;
    double cov_bkg = 0;

  for(int i=0; i<N; i++){
    if(i!=1&&i!=3&&i!=9){
      Nbkg_tot += Nbkg[i]->getVal();
      var_bkg  += std::pow(Nbkg[i]->getError(),2);
      for(int j=i+1; j<N; j++){
	if(j!=1&&j!=3&&j!=9){
	  cov_bkg += 2.0 * C(cnt_bg[i], cnt_bg[j]);
	}
      }
    }
  }
    var_bkg = var_bkg+cov_bkg;
  double err_bkg  = (var_bkg>0) ? std::sqrt(var_bkg) : 0.0;
  
  // grand total (signal+background)
  double N_all = Nsig_tot + Nbkg_tot;
  
  printf("\n[Totals]\n");
  printf("  Nsig_total = %.3f ± %.3f\n", Nsig_tot, err_sig);
  printf("  Nbkg_total = %.3f ± %.3f\n", Nbkg_tot, err_bkg);
  printf("  N_all      = %.3f\n\n",      N_all);
  
}
