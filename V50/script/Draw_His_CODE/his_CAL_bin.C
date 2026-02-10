R__LOAD_LIBRARY(libGui)
R__LOAD_LIBRARY(libFramework)
//R__LOAD_LIBRARY(libRawObjs)
//R__LOAD_LIBRARY(libCupObjs)
#include <fstream>
void his_CAL_bin()
{
	/*
	///////////////////fitting check canvas//////////////////////
	TCanvas * can = new TCanvas("can", "", 1200, 1000);
	can->Divide(4, 4); //wave form fitting check (all chn)

	TCanvas * can_d = new TCanvas("can_d", "", 700, 500); // (each chn)
	/////////////////////////////////////////////////////////////
	*/
	int runnum;
	char runnumber3[256];
	char runnumber6[256];
	char bintag[256];
	cout<< "run number : ";
	scanf("%i",&runnum);
	sprintf(runnumber6,"%06d",runnum);
	sprintf(runnumber3,"%03d",runnum);

	int binnum;
	cout<< "x-bin number (default:4000(1ADC)) : ";
	scanf("%i",&binnum);


	char hisfile[256];
	char calfile[256];
	char hisname[256];
	//char hisname_day[256];



	//cupcluster, personal dir
	//sprintf(calfile,"/home/hpge/ANA500/RUN%s/Cal_%s.root",runnumber3,runnumber6); 
//	sprintf(calfile,"/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN%s/Cal_%s.root",runnumber3,runnumber6); 
//	sprintf(hisfile,"/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN%s/Bin%i/his_%s.root",runnumber3,binnum,runnumber6); 
	sprintf(calfile,"/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN%s/Cal_%s.root",runnumber3,runnumber6); 
	sprintf(hisfile,"/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN%s/Bin%i/his_%s.root",runnumber3,binnum,runnumber6); 

	//sprintf(pedfile_name,"/home/kkw/DAQ/ANA/data/ped_%s.dat",runnumber); //cupcluster, personal dir

	int hisbin = binnum;
	double binsize = 4000. / ((double) binnum);
	TH1D * his[16];
	for(int i = 0;i<16;i++){
		sprintf(hisname,"his%i",i+1);
		//	sprintf(hisname_perday,"his%i_perday",i+1);
		his[i] = new TH1D(hisname,"",hisbin,0,4000);
		//	his_perday[i] = new TH1D(hisname_perday,hisname_perday,hisbin,0,4000);
	}

	TH1D * his_tot1 = new TH1D("his_tot1","14det",hisbin,0,4000);
	TH1D * his_tot2 = new TH1D("his_tot2","13det",hisbin,0,4000);
	TH1D * his_tot3 = new TH1D("his_tot3","12det",hisbin,0,4000);
	TH1D * his_tot1_perday = new TH1D("his_perday1","14det counts per day",hisbin,0,4000);
	TH1D * his_tot2_perday = new TH1D("his_perday2","13det counts per day",hisbin,0,4000);
	TH1D * his_tot3_perday = new TH1D("his_perday3","12det counts per day",hisbin,0,4000);


	TChain * calpara = new TChain("calpara");
	// chain->Add("./161130/161130_test06.root"); //file name
	calpara->Add(calfile); //file name

	double time_sec, time_day, energy;
	int ch, evt, file, run, mul;
	int tot = calpara->GetEntries();
	cout<<"total event : "<<tot<<endl;


	TLeaf *timesecleaf = calpara->GetLeaf("time_sec");
	TLeaf *timedayleaf = calpara->GetLeaf("time_day");
	TLeaf *energyleaf = calpara->GetLeaf("energy");
	TLeaf *chleaf = calpara->GetLeaf("ch");
	TLeaf *evtleaf = calpara->GetLeaf("evt");
	TLeaf *fileleaf = calpara->GetLeaf("file");
	TLeaf *runleaf = calpara->GetLeaf("run");
	TLeaf *mulhleaf = calpara->GetLeaf("mul");


	calpara->GetEntry(tot-1);
	time_sec = timesecleaf->GetValue();
	time_day = timedayleaf->GetValue();
//	time_sec = 6592986 + 3013930;
//	time_day = time_sec/(60*60*24);
	cout<<"total time : "<<time_sec<<" sec / "<<time_day<<" day"<<endl;


	for(int j = 0; j < tot; j++){ //all event


        calpara->GetEntry(j);
        file = fileleaf->GetValue();

		//	for(int j = 0; j < 10; j++){ //all event
//		if(file>=85&&file<=2753){
		calpara->GetEntry(j);
		energy = energyleaf->GetValue(); 
		ch = (int) chleaf->GetValue();
		his[ch-1]->Fill(energy);
		his_tot1->Fill(energy);
		if(ch!=4){his_tot2->Fill(energy);}
		if(ch!=4&&ch!=10){his_tot3->Fill(energy);}
//	}
}

TFile f(hisfile,"RECREATE");
char title[256];
char ytitle1[256];
char ytitle2[256];
//	can1->cd();
for(int c = 0; c<16; c++){
	//      his[c]->Draw(0);

	sprintf(hisname,"his%i_day",c+1);
	TH1D *his_temp = (TH1D*)his[c]->Clone();
	his_temp->Sumw2();
	his_temp->Scale((1./time_day));
	his_temp->SetName(hisname);

	sprintf(title,"Spectrum ch%i",c+1);
	sprintf(ytitle1,"Counts/(%.1f keV)", binsize);
	sprintf(ytitle2,"Counts/day/(%.1f keV)", binsize);

	his[c]->SetTitle(title);
	his[c]->SetYTitle(ytitle1);
	his[c]->SetXTitle("Energy[keV]");


	his_temp->SetTitle(title);
	his_temp->SetYTitle(ytitle2);
	his_temp->SetXTitle("Energy[keV]");

	his[c]->SetStats(kFALSE);
	his_temp->SetStats(kFALSE);
	his_temp->Write(0);
	his[c]->Write(0);

}
TH1D *his_tottemp1 = (TH1D*)his_tot1->Clone();
his_tottemp1->Sumw2();
his_tottemp1->Scale((1./time_day));
his_tottemp1->SetName("his_totday1");

his_tot1->SetTitle("Spectrum 14det.");
his_tot1->SetXTitle("Energy[keV]");
his_tot1->SetYTitle(ytitle1);
his_tottemp1->SetTitle("Spectrum 14det.");
his_tottemp1->SetXTitle("Energy[keV]");
his_tottemp1->SetYTitle(ytitle2);

his_tot1->SetStats(kFALSE);
his_tottemp1->SetStats(kFALSE);
his_tottemp1->Write(0);
his_tot1->Write(0);

TH1D *his_tottemp2 = (TH1D*)his_tot2->Clone();
his_tottemp2->Sumw2();
his_tottemp2->Scale((1./time_day));
his_tottemp2->SetName("his_totday2");

his_tot2->SetTitle("Spectrum 13det.");
his_tot2->SetXTitle("Energy[keV]");
his_tot2->SetYTitle(ytitle1);
his_tottemp2->SetTitle("Spectrum 13det.");
his_tottemp2->SetXTitle("Energy[keV]");
his_tottemp2->SetYTitle(ytitle2);

his_tot2->SetStats(kFALSE);
his_tottemp2->SetStats(kFALSE);
his_tottemp2->Write(0);
his_tot2->Write(0);


TH1D *his_tottemp3 = (TH1D*)his_tot3->Clone();
his_tottemp3->Sumw2();
his_tottemp3->Scale((1./time_day));
his_tottemp3->SetName("his_totday3");

his_tot3->SetTitle("Spectrum 12det.");
his_tot3->SetXTitle("Energy[keV]");
his_tot3->SetYTitle(ytitle1);
his_tottemp3->SetTitle("Spectrum 12det.");
his_tottemp3->SetXTitle("Energy[keV]");
his_tottemp3->SetYTitle(ytitle2);

his_tot3->SetStats(kFALSE);
his_tottemp3->SetStats(kFALSE);
his_tottemp3->Write(0);
his_tot3->Write(0);


f.Close();

}
