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

void fit_mpeak()

{

  // 1)histogram load    

  //  string workdir("/data/HPGe/USERS/kkw/DAQ/ANA125");
  string workdir("/home/kkw/study/V50");

  const int pn = 6;
  double peak[pn] = {768.36, 772.29, 783.29, 785.4, 785.96, 794.95};
 
  
  char hisfile[256];
  char rawname[256];
  char resfile[256];//fit result
  char runnumber3[256];
  char runnumber6[256];
  const int N=1;
  int runnum = 650;
  int expnum=2;
  int version=2;
  int binnum = 16000;// 1bin = 1kev
  double bw = 4000./((double)binnum); //bin width, keV

  sprintf(runnumber6,"%06d",runnum);
  sprintf(runnumber3,"%03d",runnum);

  //  sprintf(hisfile,"%s/result/V50/Bin%i/his_000650.v1.root", workdir.c_str(), binnum);
  //  sprintf(resfile,"%s/result/V50/Bin%i/sim_000650.v1.root", workdir.c_str(), binnum);
  sprintf(hisfile,"%s/result/Bin%i/his_000650.v1.root", workdir.c_str(), binnum);
  sprintf(resfile,"%s/result/Bin%i/sim_000650.v1.root", workdir.c_str(), binnum);
  TFile *hf = new TFile(hisfile);
  TH1* h[N]; char hisname[235];
  /*
  for(int i=0; i<N; i++){
    sprintf(hisname,"his%i",i+1);
  h[i] = (TH1D*)hf->Get(hisname);
  }
  */
    h[0] = (TH1D*)hf->Get("his_tot3");

  //create roodatahist
  //  double peak = 1553.77;//1 peak
  //  RooRealVar E("E", "Energy (keV)", peak-20, peak+20); // E range
  //  E.setRange("fitRange", peak-10, peak+10);

    //  RooRealVar E("E", "Energy (keV)", 750, 810); // E range
    //  E.setRange("fitRange", 750, 810);
    double range = 50;
  RooRealVar E("E", "Energy (keV)", peak[0]-range, peak[5]+range); // E range
  E.setRange("fitRange", peak[0]-range, peak[5]+range);

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

/*
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
*/

 double respa[14][3];
 respa[0][0] = 0.446;
 respa[0][1] = 6.06E-4;
 respa[0][2] = -6.31E-8;
 
  double sigma_cal[N][pn];
  for(int i=0; i<N; i++){
    for(int j=0; j<pn; j++){
      sigma_cal[i][j] = std::sqrt(respa[i][0]*respa[i][0] + respa[i][1]*peak[j] + respa[i][2]*peak[j]*peak[j]);//sigma calculation
    }
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
        bkg[i] = new RooChebychev(Form("bkg%d", i+1), Form("background%d", i+1), E, RooArgList(*p0[i]));
    //      bkg[i] = new RooChebychev(Form("bkg%d", i+1), Form("background%d", i+1), E, RooArgList(*p0[i], *p1[i]));
  }
    // f(E) = 1 + p0*T1(E) + p1*T2(E)
    // f(E) = 1 + p0*T1(E)

  // gaus model

 //peak5 mean get to distinguish mean4


 //  RooRealVar mean("mean", "Peak mean", peak, peak-0.5, peak+0.5);
  RooRealVar* mean[pn];
  RooRealVar *sigma[N][pn];
  RooGaussian *gaus[N][pn];

  for(int j=0; j<pn; j++){
    //        if(j==3){mean[j] = new RooRealVar(Form("mean_%i",j+1), Form("Peak mean %i", j+1), peak[j], peak[j]-0.25, peak[j]+0.25);}
	//        if(j==4){mean[j] = new RooRealVar(Form("mean_%i",j+1), Form("Peak mean %i", j+1), peak[j], peak[j]-0.25, peak[j]+0.25);}
    //        if(j==4){mean[j] = new RooRealVar(Form("mean_%i",j+1), Form("Peak mean %i", j+1), 786.05, 786.05-0.27, 786.05+0.27);}	
    //        else{mean[j] = new RooRealVar(Form("mean_%i",j+1), Form("Peak mean %i", j+1), peak[j], peak[j]-0.3, peak[j]+0.3);}
    mean[j] = new RooRealVar(Form("mean_%i",j+1), Form("Peak mean %i", j+1), peak[j], peak[j]-0.3, peak[j]+0.3);
    //   mean[j] = new RooRealVar(Form("mean_%i",j+1), Form("Peak mean %i", j+1), peak[j], peak[j]-0.5, peak[j]+0.5);
  }

 
  for(int i=0; i<N; i++){
    for(int j=0; j<pn; j++){
      sigma[i][j] = new RooRealVar(Form("sigma%i_%i",i+1,j+1), Form("Sigma (hist%i peak%i)",i+1,j+1), sigma_cal[i][j], 0.8*sigma_cal[i][j], 1.2*sigma_cal[i][j]);
      gaus[i][j] = new RooGaussian(Form("gaus%i_%i",i+1, j+1), Form("gaussian%i_%i",i+1, j+1), E, *mean[j], *sigma[i][j]);
    }
  }

  //extended yields (area parameters)

  double totN[N];
  RooRealVar* Nsig[N][pn];
  RooRealVar* Nbkg[N];
  RooAddPdf* model[N];
  /*

  if(version==1){
    double Nsig_init[pn]  = {365.0, 6.32, 1e-4, 26.2, 102.1, 18.0};
    double Nsig_min[pn]   = {234.0, 0.0, 0.0, 1.0, 59.0, 0.0};
    double Nsig_max[pn]   = {495.9, 13.6, 24.2, 51.1, 145.5, 38.6};
  }
  else if(version==2){
    double Nsig_init[pn]  = {381.8, 6.61, 1e-4, 25.3, 106.9, 18.9};
    double Nsig_min[pn]   = {244.0, 0.0, 0.0, 1.0, 62.0, 0.0};
    double Nsig_max[pn]   = {518.75, 14.1, 25.3, 53.4, 152.22, 40.3};
  }
  //  double Nsig_min[pn]   = {234.0, 0.0, 0.0, 17, 87, 0.0};
  //  double Nsig_max[pn]   = {495.9, 14.0, 24.2, 35, 117, 38.5};
  */
  double Nsig_init[pn], Nsig_min[pn], Nsig_max[pn];
if(version==1){
  double temp_init[6] = {365.0, 6.32, 1e-4, 26.2, 102.1, 18.0};
  double temp_min[6]  = {234.0, 0.0, 0.0, 1.0, 59.0, 0.0};
  double temp_max[6]  = {495.9, 13.6, 24.2, 51.1, 145.5, 38.6};
  for(int i=0;i<pn;i++){ Nsig_init[i]=temp_init[i]; Nsig_min[i]=temp_min[i]; Nsig_max[i]=temp_max[i]; }
}
else{
  double temp_init[6] = {381.8, 6.61, 1e-4, 25.3, 106.9, 18.9};
  double temp_min[6]  = {244.0, 0.0, 0.0, 1.0, 62.0, 0.0};
  double temp_max[6]  = {518.75, 14.1, 25.3, 53.4, 152.22, 40.3};
  for(int i=0;i<pn;i++){ Nsig_init[i]=temp_init[i]; Nsig_min[i]=temp_min[i]; Nsig_max[i]=temp_max[i]; }
}

  
  for(int i=0; i<N; i++){
    totN[i] = h[i]->Integral(h[i]->FindBin(peak[0]-range), h[i]->FindBin(peak[pn-1]+range));
    for(int j=0; j<pn; j++){

      //      Nsig[i][j] = new RooRealVar(Form("Nsig%i_%i",i+1, j+1), Form("Signal yield (hist%i peak%i)",i+1, j+1), 0.2*totN[i], 0, 1.0*totN[i]);
      Nsig[i][j] = new RooRealVar(Form("Nsig%i_%i",i+1, j+1), Form("Signal yield (hist%i peak%i)",i+1, j+1), Nsig_init[j], Nsig_min[j], Nsig_max[j]);
    }
      Nbkg[i] = new RooRealVar(Form("Nbkg%i",i+1), Form("Bkg yield (hist%i)",i+1),   0.2*totN[i], 0, 1.0*totN[i]);

    //    model[i] = new RooAddPdf(Form("model%i",i+1), Form("sig+bkg%i",i+1), RooArgList(*gaus[i], *bkg[i]), RooArgList(*Nsig[i], *Nbkg[i]));
}

 for(int i=0; i<N; i++){
   RooArgList gausList;
   RooArgList yieldList;

   for (int k=0; k<pn; k++) {
     gausList.add(*gaus[i][k]);
     yieldList.add(*Nsig[i][k]);
   }
   RooAddPdf* signal = new RooAddPdf(Form("signal_%d", i+1), Form("signal hist%d", i+1), gausList, yieldList);

   model[i] = new RooAddPdf(Form("model%i",i+1), Form("sig+bkg%i",i+1), RooArgList(*signal, *bkg[i]), RooArgList(*Nbkg[i]));
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
  

std::map<std::string, RooDataHist*> dhmap;
dhmap["hist1"]  = data[0];
/*
 dhmap["hist3"]  = data[2];
dhmap["hist5"]  = data[4];
dhmap["hist6"]  = data[5];
dhmap["hist7"]  = data[6];
dhmap["hist8"]  = data[7];
dhmap["hist9"]  = data[8];
dhmap["hist11"] = data[10];
dhmap["hist12"] = data[11];
dhmap["hist13"] = data[12];
dhmap["hist14"] = data[13];
*/
RooDataHist combData("combData", "combined data", RooArgList(E), sample, dhmap);

  
  RooSimultaneous simPdf("simPdf", "simultaneous pdf", sample);
  for(int i=0; i<N; i++){
    if(i!=1&&i!=3&&i!=9){
    simPdf.addPdf(*model[i], Form("hist%i",i+1));
    }
  }
  // fitting
  RooFitResult *fitres = simPdf.fitTo(combData, Save(), Range("fitRange"), Extended(true), PrintLevel(-1));


  
    TCanvas *c = new TCanvas("tc", "tc", 1200, 800);
    c->cd();
    RooPlot *frame = E.frame(Title("Histogram 1"));
    combData.plotOn(frame, Cut("sample==sample::hist1"));
    
    simPdf.plotOn(frame, Slice(sample, "hist1"), ProjWData(sample, combData), Range("fitRange"), NormRange("fitRange"));
    frame->Draw();
  // plot drawing
  /*
  TCanvas *tc = new TCanvas("tc", "Top array", 1800, 1000);
  tc->Divide(4,2);
  for(int i=0; i<7; i++){
    if(i!=1&&i!=3){
    tc->cd(i+1);
    RooPlot *frame = E.frame(Title(Form("Histogram %i",i+1)));
    combData.plotOn(frame, Cut(Form("sample==sample::hist%i",i+1)));
    
    simPdf.plotOn(frame, Slice(sample, Form("hist%i",i+1)), ProjWData(sample, combData), Range("fitRange"), NormRange("fitRange"));
    frame->Draw();
    }
  }
  */
  /*
  TCanvas *bc = new TCanvas("bc", "Bot array", 1800, 1000);
  bc->Divide(4,2);
  for(int i=7; i<N; i++){
    if(i!=9){
    bc->cd(i-6);
    RooPlot *frame = E.frame(Title(Form("Histogram %i",i+1)));
    combData.plotOn(frame, Cut(Form("sample==sample::hist%i",i+1)));    
    simPdf.plotOn(frame, Slice(sample, Form("hist%i",i+1)), ProjWData(sample, combData), Range("fitRange"), NormRange("fitRange"));
    frame->Draw();
    }
    }
  */

  //fitting result 
  fitres->Print("v");

  //fitting check
  fitres->correlationMatrix().Print();


// === 리밋에 걸린 파라미터만 출력 ===
std::cout << "\n[Limit Hit Parameters]" << std::endl;
std::cout << "--------------------------------------------" << std::endl;

const RooArgList& fitPars = fitres->floatParsFinal();
bool anyLimit = false;

for (int i = 0; i < fitPars.getSize(); ++i) {
    RooRealVar* p = (RooRealVar*)&fitPars[i];
    double val = p->getVal();
    double minv = p->getMin();
    double maxv = p->getMax();
    double tol = 1e-6; // 허용 오차 (경계 근처 판정용)

    bool lowHit  = std::fabs(val - minv) < tol * std::max(1.0, std::fabs(minv));
    bool highHit = std::fabs(val - maxv) < tol * std::max(1.0, std::fabs(maxv));

    if (lowHit || highHit) {
        anyLimit = true;
        std::cout << "⚠️  " << p->GetName()
                  << " reached " << (lowHit ? "LOWER" : "UPPER")
                  << " limit (" << val << " = "
                  << (lowHit ? minv : maxv) << ")\n";
    }
}

if (!anyLimit)
    std::cout << "✅ No parameters reached their limits." << std::endl;

std::cout << "--------------------------------------------" << std::endl;

  
  
  // ----- totals (after fit) -----
  auto &pars = fitres->floatParsFinal();   // order of floating params
  int cnt[N][pn]; int cnt_bg[N];
  char signame[256], bkgname[256];
  for(int i=0; i<N; i++){
    if(i!=1&&i!=3&&i!=9){
      for(int j=0; j<pn; j++){
	sprintf(signame,"Nsig%i_%i",i+1, j+1);
	cnt[i][j] = pars.index(signame);
      }
      sprintf(bkgname,"Nbkg%i",i+1);
      cnt_bg[i] = pars.index(bkgname);
    }
  }
  
  const TMatrixDSym& C = fitres->covarianceMatrix();
  double Nsig_tot =0;
  double var_sig =0;
  double cov_sig = 0;
int tp=2; //target peak 2 :v50

  for(int i=0; i<N; i++){
    if(i!=1&&i!=3&&i!=9){
      Nsig_tot += Nsig[i][tp]->getVal();
      var_sig  += std::pow(Nsig[i][tp]->getError(),2);
      for(int j=i+1; j<N; j++){
	if(j!=1&&j!=3&&j!=9){
	  cov_sig += 2.0 * C(cnt[i][tp], cnt[j][tp]);
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


// ===== 피팅 결과 요약 출력 =====
std::cout << "\n[Parameter Results]" << std::endl;
std::cout << "--------------------------------------------" << std::endl;

const RooArgList& finalPars = fitres->floatParsFinal();
for (int i = 0; i < finalPars.getSize(); ++i) {
    RooRealVar* par = (RooRealVar*)&finalPars[i];
    std::cout << par->GetName()
              << " * " << par->getVal()
              << " * " << par->getError()
              << std::endl;
}

std::cout << "--------------------------------------------" << std::endl;
  
  printf("\n[Totals]\n");
  printf("  Nsig_total = %.3f ± %.3f\n", Nsig_tot, err_sig);
  printf("  Nbkg_total = %.3f ± %.3f\n", Nbkg_tot, err_bkg);
  printf("  N_all      = %.3f\n\n",      N_all);
  
}
