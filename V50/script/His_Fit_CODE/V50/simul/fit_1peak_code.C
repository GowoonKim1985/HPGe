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

  string workdir("/home/kkw/study/V50");

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

  sprintf(hisfile,"%s/result/Bin%i/his_000650.v1.root", workdir.c_str(), binnum);
  sprintf(resfile,"%s/result/V50/Bin%i/sim_000650.v1.root", workdir.c_str(), binnum);
  TFile *hf = new TFile(hisfile);

  TH1* h1 = (TH1D*)hf->Get("his_tot3");
  //  TH1* h2 = (TH1D*)hf->Get("his2");


  //peak, sigma pre set
//  double peak = 1553.77;//1 peak
//  double peak = 1001.4;//1 peak
  double peak = 726.86;//1 peak
  //  double peak = 768.36;//1 peak

  double respa[3];
  respa[0] = 0.446;
  respa[1] = 6.06E-4;
  respa[2] = -6.31E-8;
  double sigma_cal = std::sqrt(respa[0]*respa[0] + respa[1]*peak + respa[2]*peak*peak);//sigma calculation

  //create roodatahist
  /*
  RooRealVar E("E", "Energy (keV)", peak-50, peak+50); // E range 
  E.setRange("fitRange", peak-(sigma_cal*6), peak+(sigma_cal*6));
  RooDataHist data1("data1", "dataset 1", E, Import(*h1));
  */
  double w_range = 10;
  //  double f_range = sigma_cal*6;
  double f_range = 10;  
  int ibin = h1->FindBin(peak-w_range);
  int fbin = h1->FindBin(peak+w_range);
  int nb = fbin-ibin+1;
  double e1 = h1->GetXaxis()->GetBinLowEdge(ibin);
  double e2 = h1->GetXaxis()->GetBinUpEdge (fbin);
 RooRealVar E("E","Energy (keV)", e1, e2);
 E.setBins(nb);
   E.setRange("fitRange", peak-f_range, peak+f_range);  
   //  E.setRange("fitRange", e1, e2);
  RooDataHist data1("data1","data set 1",E,h1);



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
  RooRealVar sigma1("sigma1", "Sigma (hist1)", sigma_cal, 0.8*sigma_cal, 1.2*sigma_cal);
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
//  RooPlot *frame1 = E.frame(Title("Histogram 1"), Range(1500,1600));
//  RooPlot *frame1 = E.frame(Title("Histogram 1"), Range(peak-20, peak+20));
  RooPlot *frame1 = E.frame(Title("Histogram 1"), Range(peak-w_range, peak+w_range));
  
  //  combData.plotOn(frame1, Cut("sample==sample::hist1"));



combData.plotOn(frame1, Cut("sample==sample::hist1"), Range("fitRange"), CutRange("fitRange"),DataError(RooAbsData::Poisson)); // 권장: 에러 일치
  
  simPdf.plotOn(frame1, Slice(sample, "hist1"), ProjWData(sample, combData), Range("fitRange"), NormRange("fitRange"));
  frame1->Draw();
  fitres->Print("v");


  cout<<"goodness of fit"<<endl;
  std::cout << "[Goodness of Fit] Chi2/NDF = " << frame1->chiSquare() << std::endl;


  TCanvas *c_res = new TCanvas("c_res", "Pull and residual (side-by-side)", 1000, 800);
  c_res->Divide(2,2);

  c_res->cd(1);
  RooHist* hpull = frame1->pullHist();
  // RooHist* hpull  = frame->pullHist("data1","model1");   // (data-fit)/σ
  RooPlot* f_pull = E.frame(Title("Pull"));
  f_pull->addPlotable(hpull, "P");
  f_pull->Draw();

c_res->cd(2);
RooHist* hresid = frame1->residHist();
//RooHist* hresid = frame->residHist("data1","model1");  // (data-fit)

RooPlot* f_resid = E.frame(Title("Residual"));
f_resid->addPlotable(hresid, "P");
f_resid->Draw();

// (2,1) Pull y-projection (가우시안 확인)
c_res->cd(3);
TH1D* hPullDist = new TH1D("hPullDist","Pull Distribution;(Data-Fit)/#sigma;Entries",30,-5,5);
for (int i=0; i<hpull->GetN(); ++i) { double x,y; hpull->GetPoint(i,x,y); hPullDist->Fill(y); }
hPullDist->Fit("gaus","Q");
hPullDist->Draw();

// (2,2) Residual y-projection
c_res->cd(4);
TH1D* hResidDist = new TH1D("hResidDist","Residual Distribution;Data - Fit;Entries",40,-50,50);
for (int i=0; i<hresid->GetN(); ++i) { double x,y; hresid->GetPoint(i,x,y); hResidDist->Fill(y); }
hResidDist->Fit("gaus","Q");
hResidDist->Draw();



  
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
  */

}
