#include <iostream>
#include "TMath.h"
#include "TMinuit.h"

void count_M2_RUN151()
{

	char hisfile[256];
	char calfile[256];
	char res2file[256];//fit his
	char res3file[256];//fit result

	int runnum;
	char runnumber3[256];
	char runnumber6[256];
	char bintag[256];
	cout<< "run number : ";
	scanf("%i",&runnum);
	sprintf(runnumber6,"%06d",runnum);
	sprintf(runnumber3,"%03d",runnum);
	int binnum;
	cout<< "bin number : ";
	scanf("%i",&binnum);

	sprintf(calfile,"/home/kkw/DAQ/ANA500/result/RUN%s/Cal_%s.root",runnumber3,runnumber6);
	sprintf(hisfile,"/home/kkw/DAQ/ANA500/result/RUN%s/Bin%i/his_%s_M2.root",runnumber3,binnum, runnumber6);
	sprintf(res2file,"/home/kkw/DAQ/ANA500/result/RUN%s/Bin%i/count_M2_%s.txt",runnumber3,binnum,runnumber6);
	sprintf(res3file,"/home/kkw/DAQ/ANA500/result/RUN%s/Bin%i/count_M2_%s.root",runnumber3,binnum,runnumber6);


	//energy set 
//	double fpeak[3]; //sigma*3, 
//	fpeak[0] = 0.8 * 3; //0~1000keV, sig. avg. 0.8
//	fpeak[1] = 1.1 * 3; //0~2000keV, sig. avg. 1.1	int set = 1;

	double fpeak[3];
	int set = 2;
	int peak = 2;
	double e1[peak], e2[peak], e3[peak];
	double e1min[peak], e1max[peak], e2min[peak], e2max[peak], e3min[peak], e3max[peak];

	int e1bin[peak], e1bin_min[peak], e1bin_max[peak];
	int e2bin[peak], e2bin_min[peak], e2bin_max[peak];
	int e3bin[peak], e3bin_min[peak], e3bin_max[peak];

	double Eset[peak][3];

	//EC
	Eset[0][0] = 332.3;
	Eset[0][1] = 215.3;
	//Beta
	Eset[1][0] = 350.9;
	Eset[1][1] = 234.0;
	


	//e1[peak] = 332.3;
	//e2[peak] = 215.3;


	//time check
	double time_sec, time_day;

        TChain * calpara = new TChain("calpara");
	calpara->Add(calfile); //file name
        TLeaf *timesecleaf = calpara->GetLeaf("time_sec");
        TLeaf *timedayleaf = calpara->GetLeaf("time_day");

	int tot = calpara->GetEntries();
	calpara->GetEntry(tot-1);
//	time_sec = timesecleaf->GetValue();
//	time_day = timedayleaf->GetValue();


        time_sec = 450044 + 1904769;//151
//        time_sec = 6592986+3013930;//160
        time_day = time_sec/(24*60*60);


	printf("\nMeasured time (sec) : %0.2f sec \n", time_sec);		
	printf("\nMeasured time (day) : %0.2f day \n", time_day);		

        TChain * mul2_para = new TChain("mul2_para");
	mul2_para->Add(hisfile); //file name

	int m2cnt, area;
	double evt, eset, file, E1, E2, ch1, ch2, totE; 
	int p1, p2, p3, p4, p5, p3_double, p4_double;

	double respa1 = 0.000158;
	double respa2 = 0.547238;
	double respa3 = 0.006808;

        TLeaf *evtleaf = mul2_para->GetLeaf("evt");
        TLeaf *fileleaf = mul2_para->GetLeaf("file");
        TLeaf *esetleaf = mul2_para->GetLeaf("eset");
        TLeaf *e1leaf = mul2_para->GetLeaf("mul2_e1");
        TLeaf *e2leaf = mul2_para->GetLeaf("mul2_e2");
        TLeaf *ch1leaf = mul2_para->GetLeaf("mul2_ch1");
        TLeaf *ch2leaf = mul2_para->GetLeaf("mul2_ch2");
	TLeaf *toteleaf = mul2_para->GetLeaf("totE");

	TNtuple * areapara = new TNtuple("areapara", "area", "Eset:area:evt:file:mul2_e1:mul2_e2:mul2_ch1:mul2_ch2:totE");

	for(int pc=0 ; pc<peak ; pc++){

	p1=0; p2=0; p3=0; p4=0; p5=0;
	e1[pc] = Eset[pc][0];
	e2[pc] = Eset[pc][1];
	e3[pc] = e1[pc]+e2[pc];

	fpeak[0] = 3 * e1[pc] * (respa1 + (respa2 / e1[pc]) + (respa3/(sqrt(e1[pc]))));
	fpeak[1] = 3 * e2[pc] * (respa1 + (respa2 / e2[pc]) + (respa3/(sqrt(e2[pc]))));
	//fpeak[2] = 3 * e3[pc] * (respa1 + (respa2 / e3[pc]) + (respa3/(sqrt(e3[pc]))));
	fpeak[2] = fpeak[0] + fpeak[1];

cout<< "sigma for peak 1,2,3 : "<<fpeak[0]/3<<" : "<<fpeak[1]/3<<" : "<<fpeak[2]/3<<" : "<<endl;
	e1min[pc] = e1[pc]-fpeak[0];
	e1max[pc] = e1[pc]+fpeak[0];

	e2min[pc] = e2[pc]-fpeak[1];
	e2max[pc] = e2[pc]+fpeak[1];

//	e3min[pc] = e1min[pc]+e2min[pc];
//	e3max[pc] = e1max[pc]+e2max[pc];
	e3min[pc] = e3[pc]-fpeak[2];
	e3max[pc] = e3[pc]+fpeak[2];

	printf("\n------------------------------------\n");
	cout << endl;
	cout<< "[Energy Set "<<pc+1<<"]"<<endl;
	cout << endl;
	cout << endl;
	cout << "E1 = " << e1[pc] << " keV, " << "E2 = " << e2[pc] << " keV" << endl;
	cout << endl;
	cout << "min x of "<<e1[pc]<<"keV = " << e1[pc] - fpeak[0] << " keV" << endl;
	cout << "max x of "<<e1[pc]<<"keV = " << e1[pc] + fpeak[0] << " keV" << endl;
	//cout << "bin of min x of "<<e1[pc]<<" keV = " << e1min[pc] << endl;
	//cout << "bin of "<<e1[pc]<<"keV = " << e1bin[pc] << endl;
	//cout << "bin of max x of "<<e1[pc]<<"keV = " << e1max[pc] << endl;

	cout << endl;
	cout << "min x of "<<e2[pc]<<"keV = " << e2[pc] - fpeak[1] << " keV" << endl;
	cout << "max x of "<<e2[pc]<<"keV = " << e2[pc] + fpeak[1] << " keV" << endl;
	//cout << "bin of min x of "<<e2[pc]<<" keV = " << e2min[pc] << endl;
	//cout << "bin of "<<e2[pc]<<"keV = " << e2bin[pc] << endl;
	//cout << "bin of max x of "<<e2[pc]<<"keV = " << e2max[pc] << endl;

	cout << endl;
//	cout << "min x of "<<e3[3]<<"keV (E1+E2) = " << e1[1] - fpeak[3] << " keV" << endl;
//	cout << "max x of "<<e3[3]<<"keV (E1+E2) = " << e1[1] + fpeak[3] << " keV" << endl;
//	cout << "bin of minest x of "<<e1[pc]<<" + "<<e2[pc]<<" keV = " << e3min[pc] << endl;
//	cout << "bin of maxext x of "<<e1[pc]<<" + "<<e2[pc]<<" keV = " << e3max[pc] << endl;

	printf("\n------------------------------------\n");

	m2cnt = mul2_para->GetEntries();	
	//P1: all perpect coin., no energy loss


	        for(int j = 0; j < m2cnt; j++){ //all event

		mul2_para->GetEntry(j);

		evt = evtleaf->GetValue();

		file = fileleaf->GetValue();
		E1 = e1leaf->GetValue();
		E2 = e2leaf->GetValue();
		ch1 = ch1leaf->GetValue();
		ch2 = ch2leaf->GetValue();
		totE = toteleaf->GetValue();
		if(totE<e3max[pc]){
			if(totE>e3min[pc]&&totE<e3max[pc]){
				if(E1>e1min[pc] && E1<e1max[pc] && E2>e2min[pc] && E2<e2max[pc]){
				p1 = p1+1;
				areapara->Fill(pc,1,evt,file,E1,E2,ch1,ch2,totE);
				}
				else{
				p2 = p2+1;
				areapara->Fill(pc,2,evt,file,E1,E2,ch1,ch2,totE);
				}
			}

			else{
				if((E1>e1min[pc]&&E1<e1max[pc])||
				(E2>e1min[pc]&&E2<e2max[pc])){
				p4 = p4+1;
				areapara->Fill(pc,4,evt,file,E1,E2,ch1,ch2,totE);
				}
				else if((E1>e2min[pc]&&E1<e2max[pc])||
				(E2>e2min[pc]&&E2<e2max[pc])){
				p3 = p3+1;
				areapara->Fill(pc,3,evt,file,E1,E2,ch1,ch2,totE);
				}
				else{
				p5 = p5+1;
				areapara->Fill(pc,5,evt,file,E1,E2,ch1,ch2,totE);
				}
			}
			}
		}

cout << "P1 counts(Perfect coin.) : "<<p1<<endl;
cout << "P2 counts(coin. w/ compton effect, w/o P1) : "<<p2<<endl;
cout << "P3 counts(w/o P1) : "<<p3<<endl;
cout << "P4 counts(w/p P1) : "<<p4<<endl;
cout << "others : "<<p5<<endl;
cout <<"total analized count check : "<<p1+p2+p3+p4+p5<<endl;;
	 }

TFile f(res3file, "RECREATE");
areapara->Write(0);
f. Close();
			

}
