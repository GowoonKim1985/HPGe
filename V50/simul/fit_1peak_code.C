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

void fit_1peak_code()

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

  TH1* h1 = (TH1D*)hf->Get("his_tot3");
  //  TH1* h2 = (TH1D*)hf->Get("his2");


  //peak, sigma pre set
  double peak = 1553.77;//1 peak

  double respa[3];
  respa[0] = 0.446;
  respa[1] = 6.06E-4;
  respa[2] = -6.31E-8;
  double sigma_cal = std::sqrt(respa[0]*respa[0] + respa[1]*peak + respa[2]*peak*peak);//sigma calculation

  //create roodatahist

  RooRealVar E("E", "Energy (keV)", peak-50, peak+50); // E range 
  E.setRange("fitRange", peak-(sigma_cal*6), peak+(sigma_cal*6));
  RooDataHist data1("data1", "dataset 1", E, Import(*h1));


  //BG model (p1)
  /*
  RooRealVar p0_1("p0_1", "const1", 0.0);  // no range → free
  RooRealVar p1_1("p1_1", "slope1", 0.0);  // no range → free
  //  RooPolynomial bkg1("bkg1", "background1", E, RooArgList(p0_1, p1_1)); //poly 1+p0*E+p1*E*E
  RooPolynomial bkg1("bkg1", "background1", E, RooArgList(p0_1));
  */

  RooRealVar c1("c1", "Chebychev slope", 0.0, -1.0, 1.0);
  RooChebychev bkg1("bkg1", "Chebyshev background", E, RooArgList(c1));



  // gaus model
  
  RooRealVar mean("mean", "Peak mean", peak, peak-10, peak+10);
  RooRealVar sigma1("sigma1", "Sigma (hist1)", sigma_cal, 0.7*sigma_cal, 1.3*sigma_cal);
  RooGaussian gaus1("gaus1", "gaussian1", E, mean, sigma1);

  //extended yields (area parameters)

  double n1 = h1->Integral(h1->FindBin(peak-30), h1->FindBin(peak+30));
  
  // fitting model

    RooRealVar Nsig1("Nsig1", "Signal yield (hist1)", 0.9*n1, 0, 10.0*n1);
    RooRealVar Nbkg1("Nbkg1", "Bkg yield (hist1)",   0.1*n1, 0, 10.0*n1);
  //  RooRealVar Nsig1("Nsig1", "Signal yield (hist1)", ns1, 0, 2.0*n1);
  //  RooRealVar Nbkg1("Nbkg1", "Bkg yield (hist1)",   nb1, 0, 50);
  RooAddPdf  model1("model1", "sig+bkg1", RooArgList(gaus1, bkg1), RooArgList(Nsig1, Nbkg1));

  // simultaneous setting

  RooCategory sample("sample", "sample");
  sample.defineType("hist1");
  RooDataHist combData("combData", "combined data", E, Index(sample), Import("hist1", data1));
  
  RooSimultaneous simPdf("simPdf", "simultaneous pdf", sample);
  simPdf.addPdf(model1, "hist1");
  RooFitResult *fitres = simPdf.fitTo(combData, Save(), Range("fitRange"), Extended(true), PrintLevel(-1));

  // plot drawing
  TCanvas *c = new TCanvas("c", "Simultaneous Fit (1553.8 keV)", 800, 500);
  //  c->Divide(2,1);
  
  c->cd();
  RooPlot *frame1 = E.frame(Title("Histogram 1"), Range(1500,1600));
  combData.plotOn(frame1, Cut("sample==sample::hist1"));
  simPdf.plotOn(frame1, Slice(sample, "hist1"), ProjWData(sample, combData), Range("fitRange"), NormRange("fitRange"));
  frame1->Draw();
  fitres->Print("v");


double c1_val = c1.getVal();
double Emin = E.getMin();
double Emax = E.getMax();

double A = 1 - c1_val - (2*c1_val * Emin) / (Emax - Emin);
double B = (2*c1_val) / (Emax - Emin);

printf("\n[Converted background parameters]\n");
printf(" A (const term) = %.6f\n", A);
printf(" B (slope term) = %.6f\n", B);
printf(" -> Equivalent linear form: f(E) = %.6f + %.6f * E\n", A, B);

  /*
  double slope = p0_1.getVal();   // slope term in RooPolynomial
  double N_bkg = Nbkg1.getVal();  // total background yield
  double E_min = peak - (sigma_cal*6);
  double E_max = peak + (sigma_cal*6);
  double L = E_max - E_min;       // range width

  double norm = (E_max - E_min) + slope * 0.5 * (E_max*E_max - E_min*E_min);

  double A = (N_bkg / norm);          // normalization
  double B = A * slope;               // linear term in counts/keV per keV


  printf("\n[Background reconstruction]\n");
  printf("  Range         : %.2f - %.2f keV\n", E_min, E_max);
  printf("  Nbkg (fit)    : %.3f counts\n", N_bkg);
  printf("  slope (fit)   : %.6f\n", slope);
  printf("  A (const)     : %.6e counts/keV\n", A);
  printf("  B (slope term): %.6e counts/keV^2\n", B);
  printf("  -> B(E) = A + B*E\n");

  TF1 *f_bkg = new TF1("f_bkg", [=](double *x, double *p){ return A + B * x[0]; }, E_min, E_max, 0);
  f_bkg->SetLineColor(kRed);
  f_bkg->Draw("same");
  */

}
