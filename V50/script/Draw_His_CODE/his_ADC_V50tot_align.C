R__LOAD_LIBRARY(libGui)
R__LOAD_LIBRARY(libFramework)
//R__LOAD_LIBRARY(libRawObjs)
//R__LOAD_LIBRARY(libCupObjs)
#include <fstream>
void his_ADC_V50tot_align()
{

  double tottime = 11618056/(60*60*24); //sec

  int corr[14];//corr0 - corr
  corr[0]=0;
  corr[1]=0;
  corr[2]=41;
  corr[3]=0;
  corr[4]=4;
  corr[5]=-25;
  corr[6]=-40;
  corr[7]=-37;
  corr[8]=-14;
  corr[9]=0;
  corr[10]=-14;
  corr[11]=-66;
  corr[12]=2;
  corr[13]=-70;

  int hisbin = 16000;

	TH1D * rawhis[14];
	TH1D * rawhis_cpd[14];
	TH1D * his[14];
	TH1D * his_cpd[14];


	TFile * hisfile = new TFile("/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN624/hisADC_V50_tot.root");
	char hisname[256];

	for(int i = 0;i<14;i++){
	  sprintf(hisname,"his%i",i+1);
	  rawhis[i]=(TH1D*)hisfile->Get(hisname);
	  his[i]= new TH1D (hisname, "", hisbin, 0, hisbin);
	}


	int binc, rebin;
	for(int ihis=0; ihis<14; ihis++){
	  cout<<"hist"<<ihis<<endl;
	  for(int ibin=1; ibin<16000; ibin++){
	    rebin = ibin + corr[0] - corr[ihis];

	    if(rebin>0&&rebin<=16000){
	      //    cout<<"rebin "<<ibin<<" -> "<<rebin<<endl;
	      binc = rawhis[ihis]->GetBinContent(ibin);
	      //	      binc2 = rawhis_cpd[ihis]->GetBinContent(ibin);
	      //	      if(ihis==0){cout<<"8000bin binc2 : " <<binc2<<endl;}
	      //	      if(ibin==5000){cout<<"8000bin binc2 : " <<binc2<<endl;}
	      his[ihis]->SetBinContent(rebin, binc);
	      //	      his_cpd[ihis]->SetBinContent(rebin, binc2);
	    }
	  }
	}

	//	his_cpd[0]->Draw();
	TString rawhisfile = "hisADC_V50_tot_align.root";//outfile name
	TString outfile = "/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN624/" + rawhisfile;
TFile f(outfile,"RECREATE");

char title[256];
char ytitle1[256];
char ytitle2[256];

for(int c = 0; c<14; c++){

	sprintf(title,"Spectrum ch%i",c+1);
	sprintf(ytitle1,"Counts/ADC ch");

	his[c]->SetTitle(title);
	his[c]->SetYTitle(ytitle1);
	his[c]->SetXTitle("ADC ch");

	his[c]->SetStats(kFALSE);
	his[c]->Write(0);

}

 char hiscpdname[256];

for(int c = 0; c<14; c++){

  his_cpd[c] = (TH1D*)his[c]->Clone();
  sprintf(hiscpdname,"his_cpd%i",c+1);
  his_cpd[c]->SetName(hiscpdname);
  his_cpd[c]->Scale(1./tottime);
  his_cpd[c]->Sumw2(0);

	sprintf(title,"Spectrum ch%i (cpd)",c+1);
	sprintf(ytitle1,"Counts/ADC ch / day");

	his_cpd[c]->SetTitle(title);
	his_cpd[c]->SetYTitle(ytitle1);
	his_cpd[c]->SetXTitle("ADC ch");

	his_cpd[c]->SetStats(kFALSE);
	his_cpd[c]->Write(0);

}

TH1D * tothis = new TH1D("tothis","w/o 2,4,10 E hist",hisbin,0,hisbin);
TH1D * tothis_cpd = new TH1D("tothis_cpd","w/o 2,4,10 E hist cpd",hisbin,0,hisbin);

 for(int c = 0; c<14; c++){
   if(c!=1&&c!=3&&c!=9){
     cout<<"hist "<<c+1<<" added..."<<endl;
     tothis->Add(his[c]);
     //   tothis_cpd->Add(his_cpd[c]);
   }
 }

  tothis -> SetName("tothis");
  tothis -> SetTitle("w/o 2,4,10 E hist");
  tothis -> SetYTitle("Counts/ADC ch / day");

 tothis->Write();

  tothis_cpd = (TH1D*)tothis->Clone();
  tothis_cpd->Scale(1./tottime);
  tothis_cpd->Sumw2(0);

  tothis_cpd -> SetName("tothis_cpd");
  tothis_cpd -> SetTitle("w/o 2,4,10 E hist cpd");
  tothis_cpd -> SetYTitle("Counts/ADC ch / day");



 tothis_cpd->Write();
f.Close();

}
