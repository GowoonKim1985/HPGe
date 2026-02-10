R__LOAD_LIBRARY(libGui)
R__LOAD_LIBRARY(libFramework)
//R__LOAD_LIBRARY(libRawObjs)
//R__LOAD_LIBRARY(libCupObjs)
#include <fstream>
void his_CAL_M2()
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

	int binN;
	cout<< "x-bin number (default:4000(1ADC)) : ";
	scanf("%i",&binN);

	char hisfile[256];
	char calfile[256];
	char hisname[256];

	//cupcluster, personal dir
	//sprintf(calfile,"/home/hpge/ANA500/RUN%s/Cal_%s.root",runnumber3,runnumber6); 
	sprintf(calfile,"/home/kkw/DAQ/ANA500/result/RUN%s/Cal_%s.root",runnumber3,runnumber6); 
	sprintf(hisfile,"/home/kkw/DAQ/ANA500/result/RUN%s/his_%s_M2.root",runnumber3,binnum,runnumber6); 

	//sprintf(pedfile_name,"/home/kkw/DAQ/ANA/data/ped_%s.dat",runnumber); //cupcluster, personal dir

	double MinE = 0;
	double MaxE = 4000;
	double e1, e2, e3;

	// 180mTa coincidence
	e1=215.3;
	e2=332.3;
	e3= e1+e2;

//	char 1dimdetname[256];


	//mul2 graph
	TH1F *mul2_1dim;
	mul2_1dim = new TH1F("mul2_1dim", "", binN, MinE, MaxE);
	mul2_1dim->GetXaxis()->SetTitle("Energy(keV)");
	mul2_1dim->GetYaxis()->SetTitle("Count");
	
	TH2F *mul2_2dim; // E1 >  E2
	mul2_2dim = new TH2F("mul2_2dim","", binN, MinE, MaxE, binN, MinE, MaxE);
	mul2_2dim->GetXaxis()->SetTitle("Energy1(keV)");
	mul2_2dim->GetYaxis()->SetTitle("Energy2(kev)");

/*
	// mul2 dec graph
	TH1F *mul2_1dimdet[14];
	TH2F *mul2_2dimdet[14];

		for(int j = 0; j<14; j++){
		sprintf(2dimdetname, "mul2_1d_det%i", j);
		mul2_1dimdet[j] =new TH1F(2dimdetname, "", 14, 0, 14);
		mul2_1dimdet[j]->GetXaxis()->SetTitle("Det.");
		mul2_1dimdet[j]->GetYaxis()->SetTitle("Count");
		}

		for(int j = 0; j<14; j++){
		sprintf(mul2detname, "mul2_2dim_det%i", j);
		mul2_2dimdet[j] =new TH1F(mul2detname, "", 14, 0, 14);
		mul2_2dimdet[j]->GetXaxis()->SetTitle("Det.");
		mul2_2dimdet[j]->GetYaxis()->SetTitle("Count");
		}

	//mul2 coin graph
	TH1F *coin_det[14];
	for(int k = 0; k<14; k++){
		sprintf(coindetname, "coin_det%i", k);
		coin_det[k] =new TH1F(coindetname, "", 14, 0, 14);
		coin_det[k]->GetXaxis()->SetTitle("Det.");
		coin_det[k]->GetYaxis()->SetTitle("Count");
	}
*/


	// mul2 dec ntuple
	//double mul2_e1, mul2_e2, mul2_d1, mul2_d2, P_info;
	//TNtuple * mul2_para = new TNtuple("mul2_para", "mul2", "mul2_e1:mul2_e2:mul2_d1:mul2_d2:P_info");

	double mul2_e1, mul2_e2, mul2_ch1, mul2_ch2, totE;
	TNtuple * mul2_para = new TNtuple("mul2_para", "mul2", "evt:file:mul2_e1:mul2_e2:mul2_ch1:mul2_ch2:totE");



	//file load
	TChain * calpara = new TChain("calpara");
	calpara->Add(calfile); //file name

	
	double time_sec, time_day, energy;
	int ch, evt, file, run, mul;
	int tot = calpara->GetEntries();
	cout<<"total event : "<<tot<<endl;

	//for m2 parameter
	double m2evt[2], m2file[2], m2ch[2], m2energy[2];

	m2energy[0]=0; m2energy[1]=0;
	m2evt[0]=0; m2evt[1]=0;
	m2file[0]=0; m2file[1]=0;
	m2ch[0]=0; m2ch[1]=0;

	TLeaf *timesecleaf = calpara->GetLeaf("time_sec");
	TLeaf *timedayleaf = calpara->GetLeaf("time_day");
	TLeaf *energyleaf = calpara->GetLeaf("energy");
	TLeaf *chleaf = calpara->GetLeaf("ch");
	TLeaf *evtleaf = calpara->GetLeaf("evt");
	TLeaf *fileleaf = calpara->GetLeaf("file");
	TLeaf *runleaf = calpara->GetLeaf("run");
	TLeaf *mulleaf = calpara->GetLeaf("mul");


	calpara->GetEntry(tot-1);
	cout<< timesecleaf->GetValue()<<endl;
	time_sec = timesecleaf->GetValue();
	time_day = timedayleaf->GetValue();
	cout<<"total time : "<<time_sec<<" sec / "<<time_day<<" day"<<endl;


	for(int j = 0; j < tot; j++){ // for 1 dim

		calpara->GetEntry(j);
		mul = mulleaf->GetValue();

		if(mul==2){
		energy = energyleaf->GetValue(); 
		ch = (int) chleaf->GetValue();
//		mul2_det[ch-1]->Fill(energy);
		mul2_1dim->Fill(energy);
		}
	}
	
	cout<<""<<endl;
	cout<<"1dim finished"<<endl;
	cout<<""<<endl;

	for(int j = 0; j < tot; j++){ // for 2 dim

		calpara->GetEntry(j);
		mul = mulleaf->GetValue();

		if(mul==2){
			if(m2energy[0]==0){ // first peak	
			m2energy[0] = energyleaf->GetValue(); 
			m2ch[0] = (int) chleaf->GetValue();
			m2evt[0] = evtleaf->GetValue();
			m2file[0] = fileleaf->GetValue();
			}
			else{
			m2energy[1] = energyleaf->GetValue(); 
			m2ch[1] = (int) chleaf->GetValue();
			m2evt[1] = evtleaf->GetValue();
			m2file[1] = fileleaf->GetValue();
	
				if(m2evt[0]==m2evt[1]&&m2file[0]==m2file[1]){

					if(m2energy[0]>m2energy[1]){
					mul2_2dim->Fill(m2energy[0],m2energy[1]);
					mul2_para->Fill(m2evt[0], m2file[0], m2energy[0],m2energy[1],m2ch[0],m2ch[1],(m2energy[0]+m2energy[1]));
					}
					if(m2energy[0]<m2energy[1]){
					mul2_2dim->Fill(m2energy[1],m2energy[0]);
					mul2_para->Fill(m2evt[0], m2file[0], m2energy[1],m2energy[0],m2ch[1],m2ch[0],(m2energy[0]+m2energy[1]));
					}
				}
				else{cout<<"code check"<<endl;}

			//cout<<""<<endl;
			//cout<<m2file[0]<<" file "<<m2evt[0]<<" evt set was recorded"<<endl;

			m2energy[0]=0; m2energy[1]=0;
			m2evt[0]=0; m2evt[1]=0;
			m2file[0]=0; m2file[1]=0;
			m2ch[0]=0; m2ch[1]=0;			
			}	
		}
	}

	cout<<""<<endl;
	cout<<"2dim finished"<<endl;
	cout<<""<<endl;

TFile f(hisfile,"RECREATE");
char title[256];
/*
for(int c = 0; c<16; c++){

	sprintf(hisname,"his%i_day",c+1);
	TH1D *his_temp = (TH1D*)his[c]->Clone();
	his_temp->Sumw2();
	his_temp->Scale((1./time_day));
	his_temp->SetName(hisname);

	sprintf(title,"Spectrum ch%i",c+1);
	his[c]->SetTitle(title);
	his[c]->SetXTitle("Energy[keV]");
	his[c]->SetYTitle("Counts");
	his_temp->SetTitle(title);
	his_temp->SetXTitle("Energy[keV]");
	his_temp->SetYTitle("Counts/day");

	his[c]->SetStats(kFALSE);
	his_temp->SetStats(kFALSE);
	his_temp->Write(0);
	his[c]->Write(0);

}
*/

mul2_1dim-> Write(0);
mul2_2dim->Write(0);
mul2_para->Write(0);
/*
TH1D *his_tottemp = (TH1D*)his_tot->Clone();
//his_tottemp->Scale((1./time_day));
his_tottemp->Sumw2();
his_tottemp->Scale((1./time_day));
//his_tottemp->Scale(dayscale);
his_tottemp->SetName("his_totday");

his_tot->SetTitle("Spectrum 14det.");
his_tot->SetXTitle("Energy[keV]");
his_tot->SetYTitle("Counts");
his_tottemp->SetTitle("Spectrum 14det.");
his_tottemp->SetXTitle("Energy[keV]");
his_tottemp->SetYTitle("Counts/day");

his_tot->SetStats(kFALSE);
his_tottemp->SetStats(kFALSE);
his_tottemp->Write(0);
his_tot->Write(0);
*/

f.Close();

}
