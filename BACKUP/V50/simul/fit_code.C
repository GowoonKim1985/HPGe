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

void fit_code()

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

  TH1* h1 = (TH1D*)hf->Get("his1");
  TH1* h2 = (TH1D*)hf->Get("his2");

  //create roodatahist

  RooRealVar E("E", "Energy (keV)", 1500, 1600); // E range 
  RooDataHist data1("data1", "dataset 1", E, Import(*h1));
  RooDataHist data2("data2", "dataset 2", E, Import(*h2));


  //peak, sigma pre set
  double peak = 1553.77;//1 peak

  double respa[3];
  respa[0] = 0.446;
  respa[1] = 6.06E-4;
  respa[2] = -6.31E-8;
  double sigma_cal = std::sqrt(respa[0]*respa[0] + respa[1]*peak + respa[2]*peak*peak);//sigma calculation

  //BG model (p1)

  //  RooRealVar p0_1("p0_1", "const1", 0.0);  // no range → free
  //  RooRealVar p1_1("p1_1", "slope1", 0.0);  // no range → free
  RooRealVar p0_1("p0_1", "const1", -4.0);  // no range → free
  RooRealVar p1_1("p1_1", "slope1", 0.03);  // no range → free
  RooPolynomial bkg1("bkg1", "background1", E, RooArgList(p0_1, p1_1));
  
  RooRealVar p0_2("p0_2", "const2", 0.0);
  RooRealVar p1_2("p1_2", "slope2", 0.0);
  RooPolynomial bkg2("bkg2", "background2", E, RooArgList(p0_2, p1_2));

  // gaus model
  
  RooRealVar mean("mean", "Peak mean", peak, peak-10, peak+10);
  RooRealVar sigma1("sigma1", "Sigma (hist1)", sigma_cal, 0.7*sigma_cal, 1.3*sigma_cal);
  RooRealVar sigma2("sigma2", "Sigma (hist2)", sigma_cal, 0.7*sigma_cal, 1.3*sigma_cal);
  RooGaussian gaus1("gaus1", "gaussian1", E, mean, sigma1);
  RooGaussian gaus2("gaus2", "gaussian2", E, mean, sigma2);

  //extended yields (area parameters)

  double n1 = h1->Integral(h1->FindBin(1500), h1->FindBin(1600));
  double n2 = h2->Integral(h2->FindBin(1500), h2->FindBin(1600));

  // fitting model

  RooRealVar Nsig1("Nsig1", "Signal yield (hist1)", 0.2*n1, 0, 10.0*n1);
  RooRealVar Nbkg1("Nbkg1", "Bkg yield (hist1)",   0.8*n1, 0, 10.0*n1);
  RooAddPdf  model1("model1", "sig+bkg1", RooArgList(gaus1, bkg1), RooArgList(Nsig1, Nbkg1));

  RooRealVar Nsig2("Nsig2", "Signal yield (hist2)", 0.2*n2, 0, 10.0*n2);
  RooRealVar Nbkg2("Nbkg2", "Bkg yield (hist2)",   0.8*n2, 0, 10.0*n2);
  RooAddPdf  model2("model2", "sig+bkg2", RooArgList(gaus2, bkg2), RooArgList(Nsig2, Nbkg2));

  // simultaneous setting

  RooCategory sample("sample", "sample");
  sample.defineType("hist1");
  sample.defineType("hist2");
  
  RooDataHist combData("combData", "combined data", E, Index(sample), Import("hist1", data1), Import("hist2", data2));
  
  RooSimultaneous simPdf("simPdf", "simultaneous pdf", sample);
  simPdf.addPdf(model1, "hist1");
  simPdf.addPdf(model2, "hist2");
  
  // fitting
  RooFitResult *fitres = simPdf.fitTo(combData, Save(), Extended(true), PrintLevel(-1));

  // plot drawing
  TCanvas *c = new TCanvas("c", "Simultaneous Fit (1553.8 keV)", 1200, 500);
  c->Divide(2,1);

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
  
  //fitting result 
  fitres->Print("v");


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

}
