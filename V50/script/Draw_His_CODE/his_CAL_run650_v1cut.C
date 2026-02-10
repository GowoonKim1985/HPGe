R__LOAD_LIBRARY(libGui)
R__LOAD_LIBRARY(libFramework)
//R__LOAD_LIBRARY(libRawObjs)
//R__LOAD_LIBRARY(libCupObjs)
#include <fstream>
void his_CAL_run650_v1cut()
{
	/*
	///////////////////fitting check canvas//////////////////////
	TCanvas * can = new TCanvas("can", "", 1200, 1000);
	can->Divide(4, 4); //wave form fitting check (all chn)

	TCanvas * can_d = new TCanvas("can_d", "", 700, 500); // (each chn)
	/////////////////////////////////////////////////////////////
	*/
  int runnum = 650;
	char runnumber3[256];
	char runnumber6[256];
	char bintag[256];

	sprintf(runnumber6,"%06d",runnum);
	sprintf(runnumber3,"%03d",runnum);

	int binnum = 400;
	
		//rrs on period
	int isub[3], fsub[3];
	isub[0]=0;
	fsub[0]=114;
	isub[1]=205;
	fsub[1]=470;
	isub[2]=527;
	fsub[2]=3308;
	
	char hisfile[256];
	char calfile[256];
	char hisname[256];
	//char hisname_day[256];



	//cupcluster, personal dir
	//sprintf(calfile,"/home/hpge/ANA500/RUN%s/Cal_%s.root",runnumber3,runnumber6); 
//	sprintf(calfile,"/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN%s/Cal_%s.root",runnumber3,runnumber6); 
//	sprintf(hisfile,"/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN%s/Bin%i/his_%s.root",runnumber3,binnum,runnumber6); 
//	sprintf(calfile,"/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN%s/Cal_%s_2peak.root",runnumber3,runnumber6); 
//	sprintf(hisfile,"/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN%s/Bin%i/his_%s_2peak_rrs_on.root",runnumber3,binnum,runnumber6); 

	//tot run, source calibration
	sprintf(calfile,"/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN%s/Cal_v1/Cal_%s.v1.root",runnumber3,runnumber6); 
	sprintf(hisfile,"/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN%s/Bin%i/his_%s.v1cut.root",runnumber3,binnum,runnumber6); 

	//	sprintf(calfile,"/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN%s/Cal_%s.root",runnumber3,runnumber6); 
	//	sprintf(hisfile,"/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN%s/Bin%i/his_%s.rrs_on.root",runnumber3,binnum,runnumber6); 


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

	//	TH1D * his_tot1 = new TH1D("his_tot1","14det",hisbin,0,4000);
	TH1D * his_tot1 = new TH1D("his_tot1","13det w/o 10",hisbin,0,4000);
	TH1D * his_tot2 = new TH1D("his_tot2","12det w/o 2, 10",hisbin,0,4000);
	TH1D * his_tot3 = new TH1D("his_tot3","12det w/o 2, 4, 10",hisbin,0,4000);
	TH1D * his_tot1_perday = new TH1D("his_perday1","13det counts per day",hisbin,0,4000);
	TH1D * his_tot2_perday = new TH1D("his_perday2","12det counts per day",hisbin,0,4000);
	TH1D * his_tot3_perday = new TH1D("his_perday3","11det counts per day",hisbin,0,4000);

	TChain * calpara = new TChain("calpara");
	calpara->Add(calfile);


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


	
	//rrs on

	double itime[3], ftime[3];
	int ifile, icheck;

	for(int i=0; i<tot; i++){
	  calpara->GetEntry(i);	  
	  ifile=fileleaf->GetValue();
	  if(ifile==isub[0]){
	    itime[0]=timesecleaf->GetValue();
	    cout<<itime[0]<<endl;
	    icheck=i;
	    break;
	  }
	}

	for(int i=icheck; i<tot; i++){
	  calpara->GetEntry(i);	  
	  ifile=fileleaf->GetValue();
	  if(ifile==(fsub[0]+1)){
	    calpara->GetEntry(i-1);
	    ftime[0]=timesecleaf->GetValue();
	    cout<<ftime[0]<<endl;
	    icheck=i-1;
	    break;
	  }
	}

	for(int i=icheck; i<tot; i++){
	  calpara->GetEntry(i);	  
	  ifile=fileleaf->GetValue();
	  if(ifile==isub[1]){
	    itime[1]=timesecleaf->GetValue();
	    cout<<itime[1]<<endl;
	    icheck=i;
	    break;
	  }
	}

	for(int i=icheck; i<tot; i++){
	  calpara->GetEntry(i);	  
	  ifile=fileleaf->GetValue();
	  if(ifile==(fsub[1]+1)){
	    calpara->GetEntry(i-1);
	    ftime[1]=timesecleaf->GetValue();
	    cout<<ftime[1]<<endl;
	    icheck=i-1;
	    break;
	  }
	}


	for(int i=icheck; i<tot; i++){
	  calpara->GetEntry(i);	  
	  ifile=fileleaf->GetValue();
	  if(ifile==isub[2]){
	    itime[2]=timesecleaf->GetValue();
	    cout<<itime[2]<<endl;
	    icheck=i;
	    break;
	  }
	}

	//fin file end time
	calpara->GetEntry(tot-1);
	ftime[2]=timesecleaf->GetValue();
	cout<<ftime[2]<<endl;
	

	cout<<"time 1 : "<<itime[0]<<" to "<<ftime[0]<<endl;
	cout<<"time 1 : "<<ftime[0]-itime[0]<<endl;
	cout<<"time 2 : "<<itime[1]<<" to "<<ftime[1]<<endl;
	cout<<"time 2 : "<<ftime[1]-itime[1]<<endl;
	cout<<"time 3 : "<<itime[2]<<" to "<<ftime[2]<<endl;
	cout<<"time 3 : "<<ftime[2]-itime[2]<<endl;


	time_sec = (ftime[0]-itime[0])+(ftime[1]-itime[1])+(ftime[2]-itime[2]);
	time_day = time_sec/(60*60*24);

	cout<<"total time : "<<time_sec<<" sec / "<<time_day<<" day"<<endl;
	

	/*

	//tot time

	calpara->GetEntry(tot-1);
	time_sec = timesecleaf->GetValue();
	time_day = time_sec/(60*60*24);
	
	cout<<"total time : "<<time_sec<<" sec / "<<time_day<<" day"<<endl;
	*/

	
	for(int j = 0; j < tot; j++){ //all event

	  //	  cout<<j<<endl;
	  calpara->GetEntry(j);

	  file = fileleaf->GetValue();
	  
	  if(((file>=isub[0])&&(file<=fsub[0]))||((file>=isub[1])&&(file<=fsub[1]))||((file>=isub[2])&&(file<=fsub[2]))){
	    
	    calpara->GetEntry(j);
	    energy = energyleaf->GetValue(); 
	    ch = (int) chleaf->GetValue();
	    his[ch-1]->Fill(energy);
	    his_tot1->Fill(energy);
	    
	    if(ch!=2){his_tot2->Fill(energy);}
	    if(ch!=2&&ch!=4){his_tot3->Fill(energy);}

	   	  }
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

his_tot1->SetTitle("Spectrum 13det.");
his_tot1->SetXTitle("Energy[keV]");
his_tot1->SetYTitle(ytitle1);
his_tottemp1->SetTitle("Spectrum 13det.");
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

his_tot2->SetTitle("Spectrum 12det.");
his_tot2->SetXTitle("Energy[keV]");
his_tot2->SetYTitle(ytitle1);
his_tottemp2->SetTitle("Spectrum 12det.");
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

his_tot3->SetTitle("Spectrum 11det.");
his_tot3->SetXTitle("Energy[keV]");
his_tot3->SetYTitle(ytitle1);
his_tottemp3->SetTitle("Spectrum 11det.");
his_tottemp3->SetXTitle("Energy[keV]");
his_tottemp3->SetYTitle(ytitle2);

his_tot3->SetStats(kFALSE);
his_tottemp3->SetStats(kFALSE);
his_tottemp3->Write(0);
his_tot3->Write(0);


f.Close();
	
}
