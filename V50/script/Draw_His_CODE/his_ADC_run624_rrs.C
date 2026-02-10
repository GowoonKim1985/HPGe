R__LOAD_LIBRARY(libGui)
R__LOAD_LIBRARY(libFramework)
//R__LOAD_LIBRARY(libRawObjs)
//R__LOAD_LIBRARY(libCupObjs)
#include <fstream>
void his_ADC_run624_rrs()
{
	/*
	///////////////////fitting check canvas//////////////////////
	TCanvas * can = new TCanvas("can", "", 1200, 1000);
	can->Divide(4, 4); //wave form fitting check (all chn)

	TCanvas * can_d = new TCanvas("can_d", "", 700, 500); // (each chn)
	/////////////////////////////////////////////////////////////
	*/
	int runnum = 624;
	char runnumber3[256];
	char runnumber6[256];
	char bintag[256];

	sprintf(runnumber6,"%06d",runnum);
	sprintf(runnumber3,"%03d",runnum);

	int binnum = 16000;

		//rrs on period
	int isub[2], fsub[2];
	isub[0]=4;
	fsub[0]=1067;
	isub[1]=1688;
	fsub[1]=2695;

	char hisfile[256];
	char anafile[256];
	char hisname[256];
	//char hisname_day[256];



	//cupcluster, personal dir
	//sprintf(calfile,"/home/hpge/ANA500/RUN%s/Cal_%s.root",runnumber3,runnumber6); 
//	sprintf(calfile,"/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN%s/Cal_%s.root",runnumber3,runnumber6); 
//	sprintf(hisfile,"/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN%s/Bin%i/his_%s.root",runnumber3,binnum,runnumber6); 
//	sprintf(calfile,"/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN%s/Cal_%s_2peak.root",runnumber3,runnumber6); 
//	sprintf(hisfile,"/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN%s/Bin%i/his_%s_2peak_rrs_on.root",runnumber3,binnum,runnumber6); 

	//tot run, source calibration
	sprintf(anafile,"/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN%s/Gana/Gana_%s.root",runnumber3,runnumber6); 

	sprintf(hisfile,"/data/HPGe/USERS/kkw/DAQ/ANA125/result/RUN%s/hisADC_%s_rrs_on.root",runnumber3,runnumber6); 


	int hisbin = binnum;
	double binsize = 1.;

	TH1D * his[16];
	for(int i = 0;i<16;i++){
		sprintf(hisname,"his%i",i+1);
		//	sprintf(hisname_perday,"his%i_perday",i+1);
		his[i] = new TH1D(hisname,"",hisbin,0,binnum);
		//	his_perday[i] = new TH1D(hisname_perday,hisname_perday,hisbin,0,4000);
	}

	//	TH1D * his_tot1 = new TH1D("his_tot1","14det",hisbin,0,4000);

	TChain * fitpara = new TChain("fitpara");
	fitpara->Add(anafile);


	double fit_h, time, time_sec, time_day;
	int ch, evt, file, run, mul;
	int tot = fitpara->GetEntries();
	cout<<"total event : "<<tot<<endl;
	TLeaf *runleaf = fitpara->GetLeaf("run");
	TLeaf *fileleaf = fitpara->GetLeaf("file");
	TLeaf *evtleaf = fitpara->GetLeaf("evt");
	TLeaf *timeleaf = fitpara->GetLeaf("time");
	TLeaf *chleaf = fitpara->GetLeaf("ch");
	TLeaf *mulhleaf = fitpara->GetLeaf("mul");
	TLeaf *fithleaf = fitpara->GetLeaf("fit_h");

	
	//rrs on

	double itime[2], ftime[2];
	int ifile, icheck;

	for(int i=0; i<tot; i++){
	  fitpara->GetEntry(i);	  
	  ifile=fileleaf->GetValue();
	  if(ifile==isub[0]){
	    itime[0]=(timeleaf->GetValue())/1000000000.;
	    cout<<itime[0]<<endl;
	    icheck=i;
	    break;
	  }
	}

	for(int i=icheck; i<tot; i++){
	  fitpara->GetEntry(i);	  
	  ifile=fileleaf->GetValue();
	  if(ifile==(fsub[0]+1)){
	    fitpara->GetEntry(i-1);
	    ftime[0]=(timeleaf->GetValue())/1000000000.;
	    cout<<ftime[0]<<endl;
	    icheck=i-1;
	    break;
	  }
	}

	for(int i=icheck; i<tot; i++){
	  fitpara->GetEntry(i);	  
	  ifile=fileleaf->GetValue();
	  if(ifile==isub[1]){
	    itime[1]=(timeleaf->GetValue())/1000000000.;
	    cout<<itime[1]<<endl;
	    icheck=i;
	    break;
	  }
	}

	for(int i=icheck; i<tot; i++){
	  fitpara->GetEntry(i);	  
	  ifile=fileleaf->GetValue();
	  if(ifile==(fsub[1]+1)){
	    fitpara->GetEntry(i-1);
	    ftime[1]=(timeleaf->GetValue())/1000000000.;
	    cout<<ftime[1]<<endl;
	    icheck=i-1;
	    break;
	  }
	}

	cout<<"time 1 : "<<itime[0]<<" to "<<ftime[0]<<endl;
	cout<<"time 1 : "<<ftime[0]-itime[0]<<endl;
	cout<<"time 2 : "<<itime[1]<<" to "<<ftime[1]<<endl;
	cout<<"time 2 : "<<ftime[1]-itime[1]<<endl;


	time_sec = (ftime[0]-itime[0])+(ftime[1]-itime[1]);
	time_day = time_sec/(60*60*24);

//	cout<<"total time : "<<time_sec<<" sec / "<<time_day<<" day"<<endl;
	printf("total time : %.1f sec / %.1f day", time_sec, time_day);
//	cout<<"total time : "<<time_sec<<" sec / "<<time_day<<" day"<<endl;
	

	
	/*
	//tot time

	fitpara->GetEntry(tot-1);
	time_sec = (timeleaf->GetValue())/1000000000.;
	time_day = time_sec/(60*60*24);

	cout<<"total time : "<<time_sec<<" sec / "<<time_day<<" day"<<endl;
	*/

	
	for(int j = 0; j < tot; j++){ //all event

	  //	  cout<<j<<endl;
	  fitpara->GetEntry(j);
	  file = fileleaf->GetValue();
	  
	  if((file>=isub[0]&&file<=fsub[0])||((file>=isub[1])&&(file<=fsub[1]))){
	    
	  //	    fitpara->GetEntry(j);
	    fit_h = fithleaf->GetValue(); 
	    ch = (int) chleaf->GetValue();
	    his[ch-1]->Fill(fit_h);

	  }
}

TFile f(hisfile,"RECREATE");
char title[256];
char ytitle[256];

//	can1->cd();
for(int c = 0; c<16; c++){
	//      his[c]->Draw(0);


	sprintf(title,"Spectrum ch%i",c+1);

	sprintf(ytitle,"Counts");

	his[c]->SetTitle(title);
	his[c]->SetYTitle(ytitle);
	his[c]->SetXTitle("ADC Channel");


	his[c]->SetStats(kFALSE);
	his[c]->Write(0);

}

f.Close();
	
}
