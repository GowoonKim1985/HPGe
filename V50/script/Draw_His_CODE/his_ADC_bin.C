R__LOAD_LIBRARY(libGui)
R__LOAD_LIBRARY(libFramework)
//R__LOAD_LIBRARY(libRawObjs)
//R__LOAD_LIBRARY(libCupObjs)
#include <fstream>
void his_ADC_bin()
{

	int runnum;
	char bintag[256];
	int tch = 16;
	cout<< "run number : ";
	scanf("%i", &runnum);
	TString runnumber6 = Form("%06d",runnum);	
	TString runnumber3 = Form("%03d",runnum);

	int binnum = 16000; // fadc125


	TString rawhis = "hisADC_"+runnumber6+".root";
	TString workdir = "/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN"+runnumber3;
	TString hisfile = workdir + "/" + rawhis;
	TString anafile = workdir + "/Gana_" + runnumber6 + ".root";
	cout<< hisfile<<endl;

	char hisname[256];

	int hisbin = binnum;
	//	double binsize = 4000. / ((double) binnum);
	double binsize = 1;
	TH1D * his[tch];
	for(int i = 0;i<tch;i++){
		sprintf(hisname,"his%i",i+1);
		his[i] = new TH1D(hisname,"",hisbin,0,hisbin);
	}

	TChain * fitpara = new TChain("fitpara");
	fitpara->Add(anafile); //file name

	int ch, evt, file, run, mul;
	double fit_h, time;

	int tot = fitpara->GetEntries();
	cout<<"total event : "<<tot<<endl;

	TLeaf *runleaf = fitpara->GetLeaf("run");
	TLeaf *fileleaf = fitpara->GetLeaf("file");
	TLeaf *evtleaf = fitpara->GetLeaf("evt");
	TLeaf *timeleaf = fitpara->GetLeaf("time");
	TLeaf *chleaf = fitpara->GetLeaf("ch");
	TLeaf *mulhleaf = fitpara->GetLeaf("mul");
	TLeaf *fithleaf = fitpara->GetLeaf("fit_h");

	fitpara->GetEntry(tot-1);

	time = timeleaf->GetValue();
	int fsubfile = fileleaf->GetValue();
	int fsubevt = evtleaf->GetValue();
	cout<<"final sub file : "<<fsubfile<<endl;

	for(int i = 0; i < tot; i++){ //all event

	  fitpara->GetEntry(i);
	  file = fileleaf->GetValue();
	  evt = evtleaf->GetValue();
	  fitpara->GetEntry(i);
	  fit_h = fithleaf->GetValue(); 
	  ch = (int) chleaf->GetValue();
	  //	  cout<<"ch:"<<ch<<endl;
	  his[ch-1]->Fill(fit_h);

	}


TFile f(hisfile,"RECREATE");
char title[256];
char ytitle1[256];
char ytitle2[256];

for(int c = 0; c<tch; c++){

	sprintf(title,"Spectrum ch%i",c+1);
	sprintf(ytitle1,"Counts/ADC Ch");

	his[c]->SetTitle(title);
	his[c]->SetYTitle(ytitle1);
	his[c]->SetXTitle("ADC Ch");

	his[c]->SetStats(kFALSE);
	his[c]->Write(0);

}

f.Close();

}
