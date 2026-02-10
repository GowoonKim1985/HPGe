#include <iostream>
#include "TMath.h"
#include "TMinuit.h"
#include <string>

void integral_source_ch()
{
	string workdir("/data/HPGe/USERS/kkw/DAQ/ANA500");

	char hisfile[256];
	char rawname[256];
	char res2file[256];//fit his
	char res3file[256];//fit result

	int runnum;
	char runnumber3[256];
	char runnumber6[256];
	char bintag[256];
	char peakfile[258];
	cout<< "run number : ";
	scanf("%i",&runnum);
	sprintf(runnumber6,"%06d",runnum);
	sprintf(runnumber3,"%03d",runnum);

	int binnum;
	cout<<"bin number : ";
	scanf("%i",&binnum);

	int chnum;
	cout<<"ch number : ";
	scanf("%i",&chnum);
	char chname[258];
	sprintf(chname,"his%i", chnum);
//	char nuname[256];
//	cout<<"Ana Nu ( 238U / 232Th / 40K / 182Ta ): ";
//	scanf("%s", nuname);

	sprintf(peakfile,"%s/script/His_Fit_CODE/peak/peak_source.dat", workdir.c_str());

	ifstream peakdata1(peakfile);
	string line_n1;
	string line_n2;
	int line_max1 = 0;
	int line_max2 = 0;
	double file, evt;
	char rtemp1[256];//nuclear info.
	double rtemp2;//peak energy
	while (peakdata1.good()){
		getline(peakdata1,line_n1);
		peakdata1>>rtemp1>>rtemp2;
		++line_max1;
	}
	cout<<"line max : "<<line_max1<<endl;

	char nu[line_max1-1][256];
	double peak[line_max1-1];

	double readfile, readevt;
	ifstream peakdata2(peakfile);
	while (peakdata2.good()){
		if(line_max2<(line_max1-1)){
			getline(peakdata2, line_n2);
			peakdata2>>rtemp1>>rtemp2;
			sprintf(nu[line_max2],"%s",rtemp1);
			peak[line_max2]=rtemp2;
			cout<<"line "<<line_max2<<" : "<<nu[line_max2]<<" , "<<peak[line_max2]<<endl;
		}
		else{break;}
		++line_max2;
	}

	double respa1 = 0.000158;
	double respa2 = 0.547238;
	double respa3 = 0.006808;


	int maxb, maxh, run;
	double maxx;
	int N = line_max1-1;
	double bg1, bg2, bg;
	double cnt1, cnt2, cnt3;
	double pa1[N][3];//gfit 
	double cal[N][2]; //integral fit func - m1fit,m2fit
	int xbin[N];
	double energy;


	double sigma_cal;
	double rawcnt, peakcnt;
	int peakbin, sidebin;
	TTree *intpara = new TTree("intpara","integral result");
	intpara->Branch("run",&runnum,"run/I");
	intpara ->Branch("peak",&energy,"energy/D");
	intpara ->Branch("raw_count",&rawcnt,"rawcnt/D");
	intpara ->Branch("peak_count",&peakcnt,"peakcnt/D");


	TF1 *m1bg[N];
	TF1 *m2bg[N];
	char m1bgname[256];
	char m2bgname[256];

	int bin = binnum;// 1bin = 1kev
	int mul = bin/4000;
	//	for(int hn = 0; hn <hmax; hn++){

	//		run = hn;
	//
	sprintf(hisfile,"%s/result/RUN%s/Bin%i/his_%s.root", workdir.c_str(), runnumber3, binnum, runnumber6);
	sprintf(res2file,"%s/result/RUN%s/Bin%i/fit_source_%s.C", workdir.c_str(), runnumber3,binnum,runnumber6);
	sprintf(res3file,"%s/result/RUN%s/Bin%i/int_source_%s.root", workdir.c_str(), runnumber3,binnum,runnumber6);
	TFile *hf=new TFile(hisfile);
	TH1D * his_temp = new TH1D("his_temp","",bin,0,4000);
	his_temp = (TH1D*)hf->Get(chname);
	his_temp->SetName("his_temp");
	his_temp->SetLineColor(1);

	sidebin = 10;

	for(int i = 0;i<N;i++){
		energy = peak[i];
		xbin[i]=his_temp->GetXaxis()->FindBin(peak[i]);
		printf("\n");
		printf("\n=========================\n");
		printf("%.2f keV bin : %i	\n",peak[i],xbin[i]);

		his_temp->GetXaxis()->SetRange(xbin[i]-1*mul,xbin[i]+1*mul);

		sigma_cal = peak[i]*(respa1 + (respa2/peak[i]) + (respa3/(sqrt(peak[i]))));
		cout<<"Cal. Sigma : "<<sigma_cal<<endl;
		peakbin = (int) (sigma_cal*3);
			if(peakbin<2){peakbin=2;}
		
		cout <<"peak bin (3 sigma or 2) : "<<peakbin<<endl;

		maxb = his_temp->GetMaximumBin();
		//maxx = maxb/mul;
		maxx = (his_temp->GetBinCenter(maxb))/mul;
		cout<<"maxb(bin) : "<<maxb<<endl;
		cout<<"maxx(kev) : "<<maxx<<endl;
		maxh = his_temp->GetMaximum();

		rawcnt = his_temp->Integral(maxb-peakbin*mul, maxb+peakbin*mul);
		cout<<"raw peak counts("<<maxb-peakbin<<"~"<<maxb+peakbin<<" bin) : "<<rawcnt<<endl;

		bg1 = 0; bg2 = 0; 
		for(int j1=1;j1<=(sidebin*mul);j1++){
			bg1 = bg1 + (his_temp->GetBinContent((maxb)-peakbin*mul-j1));
			bg2 = bg2 + (his_temp->GetBinContent((maxb)+peakbin*mul+j1));
		//	cout<<"test count : "<<j1<<endl;
			//		cout<<maxx-10+j<<"kev: "<<bg1<<" / "<<maxx+10-j<<"kev : "<<bg2<<endl;
		}
		cout <<"left BG counts("<<maxb-peakbin*mul-sidebin*mul<<"~"<<maxb-peakbin*mul-1<<" bin) : "<<bg1<<endl;
		cout <<"right BG counts("<<maxb+peakbin*mul+1<<"~"<<maxb+peakbin*mul+sidebin*mul<<" bin) : "<<bg2<<endl;
		bg = (bg1+bg2)/(sidebin*2);
		cout<<"avg BG height : "<<bg<<endl;
		peakcnt = rawcnt - bg*(peakbin*2+1);
		cout<<"BG sub. peak counts("<<maxb-peakbin<<"~"<<maxb+peakbin<<" bin) : "<<peakcnt<<endl;
		intpara->Fill();
	}


	TFile r3f (res3file,"RECREATE");
	r3f.cd();
	intpara->Write();
	//		intpara->Write();
	//	c1->Write();
	//	div_peak1->Write();
	//	if((line_max1-1)>=9){div_peak2->Write();}
	//	if((line_max1-1)>=18){div_peak3->Write();}

	r3f.Close();
	//	c2->Modified();
	//	c2->Update();
	//	c1->Modified();
	//	c1->Update();
	//	div_peak1->SaveAs(res2file);
	//	if((line_max1-1)>=9){div_peak2->SaveAs(res2file);}
	//	if((line_max1-1)>=18){div_peak3->SaveAs(res2file);}
}
