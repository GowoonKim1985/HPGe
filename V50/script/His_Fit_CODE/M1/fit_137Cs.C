#include <iostream>
#include "TMath.h"
#include "TMinuit.h"

void fit_137Cs()
{

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

	sprintf(peakfile,"/home/kkw/DAQ/ANA500/script/His_Fit_CODE/peak_R%03d/peak_137Cs.dat",runnum);

	//ifstream peakdata1("/home/kkw/DAQ/ANA500/script/His_Fit_CODE/peak_40K.dat");
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
	//ifstream peakdata2("/home/kkw/DAQ/ANA500/script/His_Fit_CODE/peak_40K.dat");
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



	
	//	sprintf(rawhis,"/home/kkw/DAQ/ANA500/result/RUN%s/his_%s.root",runnumber3, runnumber6);
	//	sprintf(rawhis,"his_%s.root",runnumber);
	TF1 * gfit = new TF1("gfit","gaus");
	TF1 * agfit = new TF1("agfit","[0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2])");
	agfit->SetParNames("Area","Mean","Sigma"); //area = sqrt(2pi)*constant*sigma
	//	TF1 * efit = new TF1("efit","expo");
	//	TF1 * e2fit = new TF1("e2fit","expo");
	TF1 * m1fit = new TF1("m1fit","([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+[3]");
	m1fit->SetParNames("Area","Mean","Sigma","Flat");
	TF1 * m2fit = new TF1("m2fit","([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+exp([3]+[4]*x)");
	m2fit->SetParNames("Area","Mean","Sigma","expo1","expo2");
	TF1 * m3fit = new TF1("m3fit","([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+[3]+exp([4]+[5]*x)");
	m3fit->SetParNames("Area","Mean","Sigma","pol0","expo1","expo2");
	TF1 * efit = new TF1("efit","expo");
	TF1 * pfit = new TF1("pfit","pol0");


	TCanvas * c1 = new TCanvas("c1","fitting",1200,800);
	TCanvas * div_peak1 = new TCanvas("fitting detail 1","fitting detail 1",1200,800);
	div_peak1->Divide(3,3);

	TCanvas * div_peak2;
	TCanvas * div_peak3;

		if((line_max1-1)>=9){
		div_peak2 = new TCanvas("fitting detail 2","fitting detail 2",1200,800);
		div_peak2->Divide(3,3);
		}

		if((line_max1-1)>=18){
		div_peak3 = new TCanvas("fitting detail 3","fitting detail 3",1200,800);
		div_peak3->Divide(3,3);
		}
	//	TCanvas * div_peak4 = new TCanvas("fitting detail 4","fitting detail 4",1200,800);
	//	TCanvas * div_peak5 = new TCanvas("fitting detail 5","fitting detail 5",1200,800);

	//div_peak4->Divide(3,3);

	//	if(line_max1>10){div_peak -> Divide(4,4);}
	/*	if(line_max1<=17&&line_max1>10){div_peak -> Divide(4,4);}
		if(line_max1<=26&&line_max1>17){div_peak -> Divide(5,5);}
		if(line_max1<=31&&line_max1>26){div_peak -> Divide(5,6);}
		if(line_max1<=37&&line_max1>31){div_peak -> Divide(6,6);}
		if(line_max1<=42&&line_max1>37){div_peak -> Divide(6,7);}
		*/
	//	agfit->SetLineColor(9);
	//	gfit->SetLineColor(8);
	m1fit->SetLineColor(4);
	m2fit->SetLineColor(2);
	pfit->SetLineColor(5);
	efit->SetLineColor(7);
	pfit->SetLineStyle(7);
	efit->SetLineStyle(7);
	int maxb, maxx, maxh, run;
	double bg1, bg2, bg;
	double cnt1, cnt2, cnt3;
	double pa1[line_max1-1][3];//gfit 
	double cal[line_max1-1][2]; //integral fit func - m1fit,m2fit
	double epa[6]; //expo fit 
	int xbin[line_max1-1];
	double energy;
	double agfitpa1, agfitpa2, agfitpa3, agfiter1, agfiter2, agfiter3;
	double gfitpa1, gfitpa2, gfitpa3, gfiter1, gfiter2, gfiter3;
	double m1fitpa1, m1fitpa2, m1fitpa3, m1fitpa4, m1fiter1, m1fiter2, m1fiter3, m1fiter4;
	double m2fitpa1, m2fitpa2, m2fitpa3, m2fitpa4, m2fitpa5, m2fiter1, m2fiter2, m2fiter3, m2fiter4, m2fiter5;
	double m3fitpa[6], m3fiter[6];
	double fixpa[3];//pol, efit1, efit2
	double intpa1, intpa2, intpa3, intpa4, intpa5, intpa6, intpa7;
	double m1pcnt, m2pcnt, m3pcnt, intpcnt;
	double m1pcnter, m2pcnter, m3pcnter , intpcnter;
	double m1pcnter1, m1pcnter2, m2pcnter1, m2pcnter2, m3pcnter1, m3pcnter2, intE1, intE2, intE3;
	double gfitchi2,agfitchi2,m1fitchi2,m2fitchi2,m3fitchi2;
	double m1ymax, m2ymax;

	run=runnum;
	TString m1stat;

	TTree *gfitpara = new TTree("gfitpara","gaus fitting result");
	gfitpara->Branch("run",&run,"run/I");
	gfitpara ->Branch("peak",&energy,"energy/D");
	gfitpara ->Branch("const",&gfitpa1,"gfitpa1/D");
	gfitpara ->Branch("mean",&gfitpa2,"gfitpa2/D");
	gfitpara ->Branch("chi2",&gfitchi2,"gfitchi2/D");
	gfitpara ->Branch("const_er",&gfiter1,"gfiter1/D");
	gfitpara ->Branch("mean_er",&gfiter2,"gfiter2/D");
	gfitpara ->Branch("sigma_er",&gfiter3,"gfiter3/D");

	TTree *agfitpara = new TTree("agfitpara","area-gaus fitting result");
	agfitpara->Branch("run",&run,"run/I");
	agfitpara ->Branch("peak",&energy,"energy/D");
	agfitpara ->Branch("area",&agfitpa1,"gfitpa1/D");
	agfitpara ->Branch("mean",&agfitpa2,"gfitpa2/D");
	agfitpara ->Branch("chi2",&agfitchi2,"gfitchi2/D");
	agfitpara ->Branch("area_er",&agfiter1,"gfiter1/D");
	agfitpara ->Branch("mean_er",&agfiter2,"gfiter2/D");
	agfitpara ->Branch("sigma_er",&agfiter3,"gfiter3/D");

	TTree *m1fitpara = new TTree("m1fitpara","gaus+pol fitting result");
	m1fitpara->Branch("run",&run,"run/I");
	m1fitpara ->Branch("peak",&energy,"energy/D");
	m1fitpara ->Branch("area",&m1fitpa1,"m1fitpa1/D");
	m1fitpara ->Branch("mean",&m1fitpa2,"m1fitpa2/D");
	m1fitpara ->Branch("sigma",&m1fitpa3,"m1fitpa3/D");
	m1fitpara ->Branch("chi2",&m1fitchi2,"m1fitchi2/D");
	m1fitpara ->Branch("pol0",&m1fitpa4,"m1fitpa4/D");
	m1fitpara ->Branch("area_er",&m1fiter1,"m1fiter1/D");
	m1fitpara ->Branch("mean_er",&m1fiter2,"m1fiter2/D");
	m1fitpara ->Branch("sigma_er",&m1fiter3,"m1fiter3/D");
	m1fitpara ->Branch("pol0_er",&m1fiter4,"m1fitper4/D");
	m1fitpara ->Branch("pcount",&m1pcnt,"m1pcnt/D");
	m1fitpara ->Branch("pcount_er",&m1pcnter,"m1pcnter/D");

	TTree *m2fitpara = new TTree("m2fitpara","gaus+expo fitting result");
	m2fitpara->Branch("run",&run,"run/I");
	m2fitpara ->Branch("peak",&energy,"energy/D");
	m2fitpara ->Branch("area",&m2fitpa1,"m2fitpa1/D");
	m2fitpara ->Branch("mean",&m2fitpa2,"m2fitpa2/D");
	m2fitpara ->Branch("sigma",&m2fitpa3,"m2fitpa3/D");
	m2fitpara ->Branch("chi2",&m2fitchi2,"m2fitchi2/D");
	m2fitpara ->Branch("exp1",&m2fitpa4,"m2fitpa4/D");
	m2fitpara ->Branch("exp2",&m2fitpa5,"m2fitpa5/D");
	m2fitpara ->Branch("area_er",&m2fiter1,"m2fiter1/D");
	m2fitpara ->Branch("mean_er",&m2fiter2,"m2fiter2/D");
	m2fitpara ->Branch("sigma_er",&m2fiter3,"m2fiter3/D");
	m2fitpara ->Branch("exp1_er",&m2fiter4,"m2fitper4/D");
	m2fitpara ->Branch("exp2_er",&m2fiter5,"m2fitper5/D");
	m2fitpara ->Branch("pcount",&m2pcnt,"m2pcnt/D");
	m2fitpara ->Branch("pcount_er",&m2pcnter,"m2pcnter/D");

	TTree *m3fitpara = new TTree("m3fitpara","gaus+pol0+expo fitting result");
	m3fitpara->Branch("run",&run,"run/I");
	m3fitpara ->Branch("peak",&energy,"energy/D");
	m3fitpara ->Branch("area",&m3fitpa[0],"m2fitpa[0]/D");
	m3fitpara ->Branch("mean",&m3fitpa[1],"m2fitpa[1]/D");
	m3fitpara ->Branch("sigma",&m3fitpa[2],"m2fitpa[2]/D");
	m3fitpara ->Branch("chi2",&m3fitchi2,"m2fitchi2/D");
	m3fitpara ->Branch("pol0",&m3fitpa[3],"m3fitpa[3]/D");
	m3fitpara ->Branch("exp1",&m3fitpa[4],"m3fitpa[4]4/D");
	m3fitpara ->Branch("exp2",&m3fitpa[5],"m3fitpa[5]/D");
	m3fitpara ->Branch("area_er",&m3fiter[0],"m3fiter[0]/D");
	m3fitpara ->Branch("mean_er",&m3fiter[1],"m3fiter[1]/D");
	m3fitpara ->Branch("sigma_er",&m3fiter[2],"m3fiter[2]/D");
	m3fitpara ->Branch("pol0_er",&m3fiter[3],"m3fiter[3]/D");
	m3fitpara ->Branch("exp1_er",&m3fiter[4],"m3fitper[4]/D");
	m3fitpara ->Branch("exp2_er",&m3fiter[5],"m3fitper[5]/D");
	m3fitpara ->Branch("pcount",&m3pcnt,"m3pcnt/D");
	m3fitpara ->Branch("pcount_er",&m3pcnter,"m3pcnter/D");

	TTree *intpara = new TTree("intpara","integral method result");
	intpara->Branch("run",&run,"run/I");
	intpara ->Branch("peak",&energy,"energy/D");
	intpara ->Branch("e1",&intpa1,"intpa1/D");
	intpara ->Branch("e2",&intpa2,"intpa2/D");
	intpara ->Branch("e3",&intpa3,"intpa3/D");
	intpara ->Branch("e4",&intpa4,"intpa4/D");
	intpara ->Branch("cnt1",&intpa5,"intpa5/D");
	intpara ->Branch("cnt2",&intpa6,"intpa6/D");
	intpara ->Branch("cnt3",&intpa7,"intpa7/D");
	intpara ->Branch("pcount",&intpcnt,"intpcnt/D");
	intpara ->Branch("pcount_er",&intpcnter,"intpcnter/D");

	TF1 *m1bg[line_max1-1];
	TF1 *m2bg[line_max1-1];
	char m1bgname[256];
	char m2bgname[256];

	int bin = 4000;// 1bin = 1kev
	//	int bin = 40000;// 1bin = 100ev
	int mul = bin/4000;
	//	for(int hn = 0; hn <hmax; hn++){

	//		run = hn;
	//

	sprintf(hisfile,"/home/kkw/DAQ/ANA500/result/RUN%s/his_%s_M1.root",runnumber3,runnumber6);
	sprintf(res2file,"/home/kkw/DAQ/ANA500/result/RUN%s/fit_137Cs_%s.C",runnumber3,runnumber6);
	sprintf(res3file,"/home/kkw/DAQ/ANA500/result/RUN%s/fit_137Cs_%s.root",runnumber3,runnumber6);
	//		sprintf(res2file,"/home/kkw/DAQ/ANA25/histogram/SP02/FIT310_1kev/fit_%s.C",rawhis[hn]);
	//		sprintf(res3file,"/home/kkw/DAQ/ANA25/histogram/SP02/FIT310_1kev/fit_%s.root",rawhis[hn]);
	/*

	//100ev
	sprintf(hisfile,"/home/kkw/DAQ/ANA25/histogram/SP02/FWHM310_100ev/phis_%s.root",rawhis[hn]);
	sprintf(res2file,"/home/kkw/DAQ/ANA25/histogram/SP02/FIT310_100ev/fit_%s.C",rawhis[hn]);
	sprintf(res3file,"/home/kkw/DAQ/ANA25/histogram/SP02/FIT310_100ev/fit_%s.root",rawhis[hn]);
	*/
	TFile *hf=new TFile(hisfile);
	TH1D * his_temp = new TH1D("his_temp","",bin,0,4000);
	his_temp = (TH1D*)hf->Get("his_tot");
	//		his_temp = (TH1D*)hf->Get("fhis_cday");
	his_temp->SetLineColor(1);
	//	c1->cd();
	/*
	   efit->SetRange(500,2200);
	   his_temp->Fit("efit","RQ+");
	   efit->GetParameters(&epa[0]);
	   efit->SetRange(2500,2800);
	   his_temp->Fit("efit","RQ+");
	   efit->GetParameters(&epa[2]);
	   efit->SetRange(3000,4000);
	   his_temp->Fit("efit","RQ+");
	   efit->GetParameters(&epa[4]);
	   */




	for(int i = 0;i<(line_max1-1);i++){
		energy = peak[i];
		xbin[i]=his_temp->GetXaxis()->FindBin(peak[i]);
		printf("\n");
		printf("\n=========================\n");
		printf("%.2f keV bin : %i	\n",peak[i],xbin[i]);

		c1->cd();
		his_temp->GetXaxis()->SetRange(xbin[i]-5*mul,xbin[i]+5*mul);

		maxb = his_temp->GetMaximumBin();
		maxx = maxb/mul;
		cout<<"maxb(bin) : "<<maxb<<endl;
		cout<<"maxx(kev) : "<<maxx<<endl;
		maxh = his_temp->GetMaximum();

		his_temp->GetXaxis()->SetRange(0,4000*mul);
		gfit->SetRange((maxx-5),(maxx+5));
		agfit->SetRange((maxx-5),(maxx+5));
		m1fit->SetRange((maxx-10),(maxx+10));
		m2fit->SetRange((maxx-10),(maxx+10));
		//	m3fit->SetRange((maxx-10),(maxx+10));

		efit->SetRange((maxx-30),(maxx+30));
		pfit->SetRange((maxx-15),(maxx+15));

		his_temp->Fit("efit","R+");
		fixpa[1] = efit->GetParameter(0);
		fixpa[2] = efit->GetParameter(1);
		//m2fit->SetRange((maxx-10),(maxx+10));

		//gaus fit
		his_temp->Fit("gfit","RQ0+");
		//his_temp->Fit("gfit","R+");
		gfitpa1 = gfit->GetParameter(0);
		gfitpa2 = gfit->GetParameter(1);
		gfitpa3 = gfit->GetParameter(2);
		gfiter1 = gfit->GetParError(0);
		gfiter2 = gfit->GetParError(1);
		gfiter3 = gfit->GetParError(2);

		gfitchi2 = (gfit->GetChisquare())/(gfit->GetNDF());

		bg1 = 0; bg2 = 0; 
		for(int j1=0;j1<(5*mul);j1++){
			bg1 = bg1 + (his_temp->GetBinContent((maxb)-(gfitpa3*3+5)*mul+j1));
			bg2 = bg2 + (his_temp->GetBinContent((maxb)+(gfitpa3*3+5)*mul-j1));
			//		cout<<maxx-10+j<<"kev: "<<bg1<<" / "<<maxx+10-j<<"kev : "<<bg2<<endl;
		}
		bg=(bg1+bg2)/(10*mul);
		fixpa[0]=bg;

cout<<gfitpa2<<endl;

		//area gaus fit
		//	c1->cd();
		//	his_temp->GetXaxis()->SetRange(0,4000*mul);

		agfit->SetParameters(2.5*gfitpa2*gfitpa3,gfitpa2,gfitpa3);
		//agfit->SetParameters(2.5*gfitpa2*gfitpa3,maxx,2);
		his_temp->Fit("agfit","R0+");
		agfitpa1 = agfit->GetParameter(0);
		agfitpa2 = agfit->GetParameter(1);
		agfitpa3 = agfit->GetParameter(2);
		agfiter1 = agfit->GetParError(0);
		agfiter2 = agfit->GetParError(1);
		agfiter3 = agfit->GetParError(2);
		agfitchi2 = (agfit->GetChisquare())/(agfit->GetNDF());

		cout<<""<<endl;
		cout<<"[m1 fitting]"<<endl;
		cout<<""<<endl;



		//m1 fit
		//m1fit->SetParameters(gfitpa1,peak[i],gfitpa3,bg);
		//m1fit->SetParameters(agfitpa1,peak[i],agfitpa3,bg);
		m1fit->ReleaseParameter(0);
		m1fit->ReleaseParameter(1);
		m1fit->ReleaseParameter(2);
		m1fit->ReleaseParameter(3);
		m1fit->SetParameters(agfitpa1,agfitpa2,agfitpa3,bg);

		his_temp->Fit("m1fit","R0+");
		m1fitpa2 = m1fit->GetParameter(1);

		cout<< "m1 fit peak center : "<<m1fitpa2<<endl;
		cout<< "m1 fit peak sigma : "<<m1fitpa3<<endl;

		//peak check
		if(abs(m1fitpa1-peak[i])>3){
			cout <<"fitting mean "<<m2fitpa2<<" : peak "<<peak[i]<<endl;
			cout<<">> peak fix"<<endl;
			m1fit->FixParameter(1,peak[i]);
		}
		his_temp->Fit("m1fit","R0+");
		m1fitpa3 = m1fit->GetParameter(2);
		//sigma check
		if(m1fitpa3<0.4){
			cout <<"Too small sigma value"<<endl;
			cout<<">> sigma fix"<<endl;
			if(peak[i]<1000){m1fit->FixParameter(2,0.8);}
			else if(peak[i]<2000){m1fit->FixParameter(2,1.1);}
			else{m1fit->FixParameter(2,1.6);}
		}
		if(m1fitpa3<0){
			cout <<"minus sigma value "<<endl;
			cout<<">> sigma fix"<<endl;
			if(peak[i]<1000){m1fit->FixParameter(2,0.8);}
			if(peak[i]>=1000&&peak[i]<2000){m1fit->FixParameter(2,1.1);}
		}
		if(m1fitpa3>1.0&&peak[i]<1000){
			cout <<"Too big sigma value @ <1000 keV"<<endl;
			cout<<">> sigma fix"<<endl;
			m1fit->FixParameter(2,0.8);
		}	
		if(m1fitpa3>2.0&&peak[i]>=1000&&peak[i]<2000){
			cout <<"Too big sigma value @ 1000~2000 keV"<<endl;
			cout<<">> sigma fix"<<endl;
			m1fit->FixParameter(2,1.1);
		}	

		his_temp->Fit("m1fit","R+");
		m1fitpa1 = m1fit->GetParameter(0);
		m1fitpa2 = m1fit->GetParameter(1);
		m1fitpa3 = m1fit->GetParameter(2);
		m1fitpa4 = m1fit->GetParameter(3);
		m1fiter1 = m1fit->GetParError(0);
		m1fiter2 = m1fit->GetParError(1);
		m1fiter3 = m1fit->GetParError(2);
		m1fiter4 = m1fit->GetParError(3);
		m1fitchi2 = (m1fit->GetChisquare())/(m1fit->GetNDF());

		m1ymax = m1fit->Eval(peak[i]);


		pfit->ReleaseParameter(0);
		pfit->FixParameter(0,m1fitpa4);
		his_temp->Fit("pfit","R+");

		sprintf(m1bgname,"m1bg_%i",i);
		m1bg[i] = new TF1(m1bgname,"pol0",peak[i]-15,peak[i]+15);
		m1bg[i]->SetParameter(0,m1fitpa4);
		m1bg[i]->SetLineColor(4);
		m1bg[i]->SetLineStyle(7);
		//m1bg[i]->Draw("SAME");


		printf("BG : %.2f	\n",bg);
		if(m1fitpa3<0){
			cout<<"!!!!m1fit sigma is minus -> will be get ABS!!!"<<endl;
		}
		m1pcnt = m1fitpa1;
		m1pcnter = m1fiter1;
		cout<<"m1fit int.error : "<<m1pcnter<<endl;

		cout<<energy<<endl;
		m1fitpara->Fill();


		cout<<""<<endl;
		cout<<"[m2 fitting]"<<endl;
		cout<<""<<endl;
		//m2fit
		/*			if(peak[i]<2200){
					m2fit->SetParameters(gfitpa1,peak[i],gfitpa3,epa[0],epa[1]);
					}
					else if(peak[i]>2500&&peak[i]<2700){
					m2fit->SetParameters(gfitpa1,peak[i],gfitpa3,epa[2],epa[3]);
					}
					else if(peak[i]>3000){
					m2fit->SetParameters(gfitpa1,peak[i],gfitpa3,epa[4],epa[5]);
					}
					*/
		//		m2fit->SetRange((maxx-5),(maxx+5));
		//m2fit->SetParameters(gfitpa1,peak[i],gfitpa3,fixpa[1],fixpa[2]);
		//m2fit->SetParameters(agfitpa1,peak[i],agfitpa3,fixpa[1],fixpa[2]);

		m2fit->ReleaseParameter(0);
		m2fit->ReleaseParameter(1);
		m2fit->ReleaseParameter(2);
		m2fit->ReleaseParameter(3);
		m2fit->ReleaseParameter(4);

		m2fit->SetParameters(agfitpa1,agfitpa2,agfitpa3,fixpa[1],fixpa[2]);
		m2fit->FixParameter(3,fixpa[1]);
		m2fit->FixParameter(4,fixpa[2]);

		his_temp->Fit("m2fit","R0+");
		m2fitpa2 = m2fit->GetParameter(1);

		cout<< "m2 fit peak center : "<<m2fitpa2<<endl;
		cout<< "m2 fit peak sigma : "<<m2fitpa3<<endl;

		//peak check
		if(abs(m2fitpa2-peak[i])>3){
			cout <<"fitting mean "<<m2fitpa2<<" : peak "<<peak[i]<<endl;
			cout<<">> peak fix"<<endl;
			m2fit->FixParameter(1,peak[i]);
		}
		his_temp->Fit("m2fit","R0+");
		m2fitpa3 = m2fit->GetParameter(2);
		//sigma check
		if(m2fitpa3<0){
			cout <<"minus sigma value "<<endl;
			cout<<">> sigma fix"<<endl;
			if(peak[i]<1000){m2fit->FixParameter(2,0.8);}
			if(peak[i]>=1000&&peak[i]<2000){m2fit->FixParameter(2,1.1);}
		}
		if(m2fitpa3<0.4){
			cout <<"Too small sigma value"<<endl;
			cout<<">> sigma fix"<<endl;
			if(peak[i]<1000){m2fit->FixParameter(2,0.8);}
			else if(peak[i]<2000){m2fit->FixParameter(2,1.1);}
			else{m2fit->FixParameter(2,1.6);}
		}
		if(m2fitpa3>1.0&&peak[i]<1000){
			cout <<"Too big sigma value @ <1000 keV"<<endl;
			cout<<">> sigma fix"<<endl;
			m2fit->FixParameter(2,0.8);
		}	
		if(m2fitpa3>2.0&&peak[i]>=1000&&peak[i]<2000){
			cout <<"Too big sigma value @ 1000~2000 keV"<<endl;
			cout<<">> sigma fix"<<endl;
			m2fit->FixParameter(2,1.1);
		}	
		his_temp->Fit("m2fit","R+");

		m2fitpa1 = m2fit->GetParameter(0);
		m2fitpa2 = m2fit->GetParameter(1);
		m2fitpa3 = m2fit->GetParameter(2);
		m2fitpa4 = m2fit->GetParameter(3);
		m2fitpa5 = m2fit->GetParameter(4);
		m2fiter1 = m2fit->GetParError(0);
		m2fiter2 = m2fit->GetParError(1);
		m2fiter3 = m2fit->GetParError(2);
		m2fiter4 = m2fit->GetParError(3);
		m2fiter5 = m2fit->GetParError(4);
		m2fitchi2 = (m2fit->GetChisquare())/(m2fit->GetNDF());

		m2ymax = m2fit->Eval(peak[i]);

		sprintf(m2bgname,"m2bg_%i",i);
		m2bg[i] = new TF1(m2bgname,"expo",peak[i]-15,peak[i]+15);
		m2bg[i]->SetParameter(0,m2fitpa4);
		m2bg[i]->SetParameter(1,m2fitpa5);
		m2bg[i]->SetLineColor(2);
		m2bg[i]->SetLineStyle(7);
		m2bg[i]->Draw("SAME");
		if(m2fitpa3<0){
			cout<<"!!!!m2fit sigma is minus -> will be get ABS!!!"<<endl;
		}
		m2pcnt = m2fitpa1;
		m2pcnter = m2fiter1;

		m2fitpara->Fill();
		/*
		//m3 fit

		m3fit->SetParameters(gfitpa1,peak[i],gfitpa3,bg,0,0);
		his_temp->Fit("m3fit","RQ0+");

		cout<<"[m1fit]"<<endl;
		cout<< "m1fit tot. counts:"<<m1pcnt<<" +/- "<<m2pcnter<<endl;
		cout<<"6 sigma : "<<m1fitpa3*6<<endl;
		cout<<""<<endl;
		cout<<"[m2fit]"<<endl;
		cout<< "m2fit tot. counts:"<<m2pcnt<<" +/- "<<m2pcnter<<endl;
		cout<<"6 sigma : "<<m2fitpa3*6<<endl;
		cout<<""<<endl;
		*/
		/*
		//special fit 2614kev

		his_temp->Fit(sp1fit,"R0+");
		his_temp->Fit(sp2fit,"R0+");
		his_temp->Fit(sp3fit,"R0+");
		sp1fit->GetParameters(&pa5[0]);
		sp1fit->GetParameters(&pa5[2]);
		sp1fit->GetParameters(&pa5[5]);
		spfit->SetParameters(pa5);
		his_temp->Fit(spfit,"R+");
		*/
		//integral.
		cout<<"[auto int.]"<<endl;
		cnt1=0; cnt2=0; cnt3=0;
		int intrange;
		int intrangebin;
		intrange = (int)abs(m1fitpa3*3); //kev
		if(intrange>5){
			if(peak[i]<1000){
				cout<<"!!!!int.range(3 sigma) sigma is too big -> will be changed to 2keV(0~1000keV)!!!"<<endl;
				intrange = 2;
			}
			else{
				cout<<"!!!!int.range(3 sigma) sigma is too big -> will be changed to 3keV(1000~4000keV)!!!"<<endl;
				intrange = 3;
			}
		}
		intrangebin = intrange*mul;
		//		cout<<"center(bin) : "<<xbin[i]<<endl;
		//		cout<<"center(kev) : "<<his_temp->GetBinCenter(xbin[i])<<endl;
		//		cout<<"intrangebin : "<<intrangebin<<endl;
		//		cout<<"maxx-int.(kev) : "<<his_temp->GetBinCenter(maxx-intrangebin)<<endl;
		for(int j3=0;j3<(intrangebin*2+1);j3++){
			cnt2 = cnt2+(his_temp->GetBinContent(xbin[i]-intrangebin+j3));
		}
		for(int j4=0;j4<(5*mul);j4++){
			cnt1 = cnt1+(his_temp->GetBinContent(xbin[i]-intrangebin-(5*mul)+j4));//5kev
			cnt3 = cnt3+(his_temp->GetBinContent(xbin[i]+intrangebin+(5*mul)-j4));//5kev
		}
		intpa1=(his_temp->GetBinCenter(xbin[i]-intrangebin-(5*mul)))-0.5/mul;
		intpa2=(his_temp->GetBinCenter(xbin[i]-intrangebin))-0.5/mul;
		intpa3=(his_temp->GetBinCenter(xbin[i]+intrangebin))+0.5/mul;
		intpa4=(his_temp->GetBinCenter(xbin[i]+intrangebin+(5*mul)))+0.5/mul;
		intpa5=cnt1/mul;
		intpa6=cnt2/mul;
		intpa7=cnt3/mul;
		intE1 = intpa2-intpa1;
		intE2 = intpa3-intpa2;
		intE3 = intpa4-intpa3;
		intpcnt = intpa6 -((intpa5+intpa7)/10)*intE2;
		cout<<"auto int. peak range : "<<intrangebin*2+1<<"(bin) / "<<intpa2<<" to "<<intpa3<<" (kev)"<<endl;
		cout<< "auto int. tot. counts:"<<intpcnt<<" +/- "<<intpcnter<<endl;
		cout<<""<<endl;

		//intpcnter = sqrt(intpa6 + ((intpa5/(intE1*intE1)+intpa7/(intE3*intE3))*((intE2*intE2)/4)));
		intpcnter = sqrt(intpa6 +(intpa5*((intE2*intE2)/(4*intE1*intE1)))+(intpa7*((intE2*intE2)/(4*intE3*intE3))));

		//		m1fitpara->Fill();
		//		m2fitpara->Fill();
		//		intpara->Fill();


		if(i<9){div_peak1->cd(i+1);}
		if(i>=9&&i<18){div_peak2->cd(i-9+1);}
		if(i>=18&&i<27){div_peak3->cd(i-18+1);}
		//		if(i>=27&&i<36){div_peak4->cd(i-36+1);}
		//	div_peak->cd(i+1);
		his_temp->GetXaxis()->SetRange(xbin[i]-15*mul,xbin[i]+15*mul);

		cout<<"m1 y : m2 y = "<< m1ymax<<" : "<<m2ymax<<endl;
		if(m1ymax>m2ymax){his_temp->SetMaximum(m1ymax+m1ymax*0.2);}
		if(m1ymax<m2ymax){his_temp->SetMaximum(m2ymax+m2ymax*0.2);}

		his_temp->Draw();

		his_temp->SetTitle("137Cs Fitting");
		his_temp->SetXTitle("Energy[keV]");
		his_temp->SetYTitle("Counts/kev");
		//his_temp->SetYTitle("Counts/kev.day");
		his_temp->SetName("fhis_tot");
		//his_temp->SetName("fit_cday");

		if(i<9){
			div_peak1->Modified();
			div_peak1->Update();
		}
		if(i>=9&&i<18){
			div_peak2->Modified();
			div_peak2->Update();
		}
		if(i>=18&&i<27){
			div_peak3->Modified();
			div_peak3->Update();
		}
		/*	if(i>=27&&i<36){
			div_peak4->Modified();
			div_peak4->Update();
			}
			*/
		//		c1->Modified();
		//		c1->Update();
		//		c1->SaveAs(res2file);

	}
		TFile r3f (res3file,"RECREATE");
		r3f.cd();
		gfitpara->Write();
		m1fitpara->Write();
		m2fitpara->Write();
		intpara->Write();
	//	c1->Write();
	//	div_peak1->Write();
	//	if((line_max1-1)>=9){div_peak2->Write();}
	//	if((line_max1-1)>=18){div_peak3->Write();}

		r3f.Close();
	c1->cd();
	his_temp->GetXaxis()->SetRange(0,4000);
	c1->Modified();
	c1->Update();
//	div_peak1->SaveAs(res2file);
//	if((line_max1-1)>=9){div_peak2->SaveAs(res2file);}
//	if((line_max1-1)>=18){div_peak3->SaveAs(res2file);}
	c1->SaveAs(res2file);
}
