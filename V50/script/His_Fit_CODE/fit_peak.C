#include <iostream>
#include "TMath.h"
#include "TMinuit.h"
#include <string>

void fit_peak()
{
	string workdir("/data/HPGe/USERS/kkw/DAQ/ANA125");

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
/*
	int chnum;
	cout<<"ch number : ";
	scanf("%i",&chnum);
	char chname[258];
	sprintf(chname,"his%i", chnum);
*/
//	char nuname[256];
//	cout<<"Ana Nu ( 238U / 232Th / 40K / 182Ta ): ";
//	scanf("%s", nuname);

	sprintf(peakfile,"%s/script/His_Fit_CODE/peak/peak_BG.dat", workdir.c_str());

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


	//General resolution
	double respa[3];
	respa[0] = 0.000158;
	respa[1] = 0.547238;
	respa[2] = 0.006808;

/*
	// run511 resolution
		
        double respa[3];
        respa[0] = 0.000548398;
        respa[1] = 1.843895463;
        respa[2] = -0.024315872;

*/


	TF1 * gfit = new TF1("gfit","gaus");
	TF1 * agfit = new TF1("agfit","[0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2])");
	agfit->SetParNames("Area","Mean","Sigma"); //area = sqrt(2pi)*constant*sigma
	TF1 * efit = new TF1("efit","expo");
	TF1 * p0fit = new TF1("p0fit","pol0");
	TF1 * p1fit = new TF1("p1fit","pol1");

	TF1 * m1fit = new TF1("m1fit","([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+[3]");
	m1fit->SetParNames("Area","Mean","Sigma","Flat");

	TF1 * m2fit = new TF1("m2fit","([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+exp([3]+[4]*x)");
	m2fit->SetParNames("Area","Mean","Sigma","expo1","expo2");

	TF1 * m3fit = new TF1("m3fit", "([0]/(sqrt(2*TMath::Pi())*[2])*exp(-((x-[1])*(x-[1]))/2/[2]/[2]))+[3]+[4]*x");
	m3fit->SetParNames("Area", "Mean", "Sigma", "Intercept", "Slope");



	TCanvas * c1 = new TCanvas("c1","fitting",1200,800);

	m1fit->SetLineColor(kBlue);
	m2fit->SetLineColor(kRed);
	m3fit->SetLineColor(kGreen);
	p0fit->SetLineColor(kBlue);
	efit->SetLineColor(kRed);
	p1fit->SetLineColor(kGreen);
	p0fit->SetLineStyle(7);
	efit->SetLineStyle(7);
	p1fit->SetLineStyle(7);
	int maxb, maxh, run;
	double maxx;
	int N = line_max1-1;
	double bg1, bg2, bg;
	double cnt1, cnt2, cnt3;
	double pa1[N][3];//gfit 
	double cal[N][2]; //integral fit func - m1fit,m2fit
	double epa[6]; //expo fit 
	int xbin[N];
	double energy;
	double gfitpa[3], gfiter[3], gfitchi2;
	double agfitpa[3], agfiter[3], agfitchi2;
	double m1fitpa[4], m1fiter[4], m1fitmean, m1fitsig, m1fitchi2;
	double m2fitpa[5], m2fiter[5], m2fitmean, m2fitsig, m2fitchi2;
	double m3fitpa[5], m3fiter[5], m3fitmean, m3fitsig, m3fitchi2;

	double fitpa[3], fiter[3], fitchi2;
	double fixbg1, fixbg2[2], fixbg3[2];//bg line pa. for m1, m2, m3
	double intpa1, intpa2, intpa3, intpa4, intpa5, intpa6, intpa7;
	double m1pcnt, m2pcnt, m3pcnt, intpcnt;
	double m1pcnter, m2pcnter, m3pcnter , intpcnter;
	double m1pcnter1, m1pcnter2, m2pcnter1, m2pcnter2, m3pcnter1, m3pcnter2, intE1, intE2, intE3;

	TTree *gfitpara = new TTree("gfitpara","gaus fitting result");
	gfitpara->Branch("run",&runnum,"run/I");
	gfitpara ->Branch("peak",&energy,"energy/D");
	gfitpara ->Branch("const",&gfitpa[0],"gfitpa0/D");
	gfitpara ->Branch("mean",&gfitpa[1],"gfitpa1/D");
	gfitpara ->Branch("sigma",&gfitpa[2],"gfitpa2/D");
	gfitpara ->Branch("chi2",&gfitchi2,"gfitchi2/D");
	gfitpara ->Branch("const_er",&gfiter[0],"gfiter0/D");
	gfitpara ->Branch("mean_er",&gfiter[1],"gfiter1/D");
	gfitpara ->Branch("sigma_er",&gfiter[2],"gfiter2/D");

	TTree *agfitpara = new TTree("agfitpara","area-gaus fitting result");
	agfitpara->Branch("run",&runnum,"run/I");
	agfitpara ->Branch("peak",&energy,"energy/D");
	agfitpara ->Branch("area",&agfitpa[0],"gfitpa0/D");
	agfitpara ->Branch("mean",&agfitpa[1],"gfitpa1/D");
	agfitpara ->Branch("sigma",&agfitpa[2],"gfitpa2/D");
	agfitpara ->Branch("chi2",&agfitchi2,"gfitchi2/D");
	agfitpara ->Branch("area_er",&agfiter[0],"gfiter0/D");
	agfitpara ->Branch("mean_er",&agfiter[1],"gfiter1/D");
	agfitpara ->Branch("sigma_er",&agfiter[2],"gfiter2/D");

	TTree *m1fitpara = new TTree("m1fitpara","gaus+pol fitting result");
	m1fitpara->Branch("run",&runnum,"run/I");
	m1fitpara ->Branch("peak",&energy,"energy/D");
	m1fitpara ->Branch("area",&m1fitpa[0],"m1fitpa0/D");
	m1fitpara ->Branch("mean",&m1fitpa[1],"m1fitpa1/D");
	m1fitpara ->Branch("sigma",&m1fitpa[2],"m1fitpa2/D");
	m1fitpara ->Branch("pol0",&m1fitpa[3],"m1fitpa3/D");
	m1fitpara ->Branch("chi2",&m1fitchi2,"m1fitchi2/D");
	m1fitpara ->Branch("area_er",&m1fiter[0],"m1fiter0/D");
	m1fitpara ->Branch("mean_er",&m1fiter[1],"m1fiter1/D");
	m1fitpara ->Branch("sigma_er",&m1fiter[2],"m1fiter2/D");
	m1fitpara ->Branch("pol0_er",&m1fiter[3],"m1fitper3/D");

	TTree *m2fitpara = new TTree("m2fitpara","gaus+expo fitting result");
	m2fitpara->Branch("run",&runnum,"run/I");
	m2fitpara ->Branch("peak",&energy,"energy/D");
	m2fitpara ->Branch("area",&m2fitpa[0],"m2fitpa0/D");
	m2fitpara ->Branch("mean",&m2fitpa[1],"m2fitpa1/D");
	m2fitpara ->Branch("sigma",&m2fitpa[2],"m2fitpa2/D");
	m2fitpara ->Branch("exp1",&m2fitpa[3],"m2fitpa3/D");
	m2fitpara ->Branch("exp2",&m2fitpa[4],"m2fitpa4/D");
	m2fitpara ->Branch("chi2",&m2fitchi2,"m2fitchi2/D");
	m2fitpara ->Branch("area_er",&m2fiter[0],"m2fiter0/D");
	m2fitpara ->Branch("mean_er",&m2fiter[1],"m2fiter1/D");
	m2fitpara ->Branch("sigma_er",&m2fiter[2],"m2fiter2/D");
	m2fitpara ->Branch("exp1_er",&m2fiter[3],"m2fitper3/D");
	m2fitpara ->Branch("exp2_er",&m2fiter[4],"m2fitper4/D");


	TTree *m3fitpara = new TTree("m3fitpara", "gaus+pol1 fitting result");
	m3fitpara -> Branch("run", &runnum, "run/I");
	m3fitpara -> Branch("peak", &energy, "energy/D");
	m3fitpara -> Branch("area", &m3fitpa[0], "m3fitpa0/D");
	m3fitpara -> Branch("mean", &m3fitpa[1], "m3fitpa1/D");
	m3fitpara -> Branch("sigma", &m3fitpa[2], "m3fitpa2/D");
	m3fitpara -> Branch("pol0", &m3fitpa[3], "m3fitpa3/D");
	m3fitpara -> Branch("pol1", &m3fitpa[4], "m3fitpa4/D");
	m3fitpara -> Branch("chi2", &m3fitchi2, "m3fitchi2/D");
	m3fitpara -> Branch("area_er", &m3fiter[0], "m3fiter0/D");
	m3fitpara -> Branch("mean_er", &m3fiter[1], "m3fiter1/D");
	m3fitpara -> Branch("sigma_er", &m3fiter[2], "m3fiter2/D");
	m3fitpara -> Branch("pol0_er", &m3fiter[3], "m3fitper3/D");
	m3fitpara -> Branch("pol1_er", &m3fiter[4], "m3fitper4/D");


	int fit_type;
	TTree *fitpara = new TTree("fitpara","best fitting result");
	fitpara->Branch("fit_type",&fit_type,"f_type/I");
	fitpara->Branch("run",&runnum,"run/I");
	fitpara ->Branch("peak",&energy,"energy/D");
	fitpara ->Branch("area",&fitpa[0],"fitpa0/D");
	fitpara ->Branch("mean",&fitpa[1],"fitpa1/D");
	fitpara ->Branch("sigma",&fitpa[2],"fitpa2/D");
	fitpara ->Branch("chi2",&fitchi2,"fitchi2/D");
	fitpara ->Branch("area_er",&fiter[0],"fiter0/D");
	fitpara ->Branch("mean_er",&fiter[1],"fiter1/D");
	fitpara ->Branch("sigma_er",&fiter[2],"fiter2/D");


	TF1 *m1bg[N];
	TF1 *m2bg[N];
	char m1bgname[256];
	char m2bgname[256];

	int bin = binnum;// 1bin = 1kev
	//	int bin = 40000;// 1bin = 100ev
	int mul = bin/4000;
	//	for(int hn = 0; hn <hmax; hn++){

	//		run = hn;
	//
	sprintf(hisfile,"%s/result/RUN%s/Bin%i/his_%s.root", workdir.c_str(), runnumber3, binnum, runnumber6);
	sprintf(res2file,"%s/result/RUN%s/Bin%i/fit_%s.C", workdir.c_str(), runnumber3,binnum,runnumber6);
	sprintf(res3file,"%s/result/RUN%s/Bin%i/fit_%s.root", workdir.c_str(), runnumber3,binnum,runnumber6);
	TFile *hf=new TFile(hisfile);
	TH1D * his_temp = new TH1D("his_temp","",bin,0,4000);
//	his_temp = (TH1D*)hf->Get(chname);
	his_temp = (TH1D*)hf->Get("his_tot");
	his_temp->SetName("his_temp");
	his_temp->SetLineColor(1);

	double sigma_cal;



	for(int i = 0;i<N;i++){
		c1->cd();
		energy = peak[i];
		xbin[i]=his_temp->GetXaxis()->FindBin(peak[i]);
		printf("\n");
		printf("\n=========================\n");
		printf("%.2f keV bin : %i	\n",peak[i],xbin[i]);

		his_temp->GetXaxis()->SetRange(xbin[i]-1*mul,xbin[i]+1*mul);

		sigma_cal = peak[i]*(respa[0] + (respa[1]/peak[i]) + (respa[2]/(sqrt(peak[i]))));
		cout<<"Cal. Sigma : "<<sigma_cal<<endl;

		maxb = his_temp->GetMaximumBin();
		//maxx = maxb/mul;
		maxx = (his_temp->GetBinCenter(maxb))/mul;
		cout<<"maxb(bin) : "<<maxb<<endl;
		cout<<"maxx(kev) : "<<maxx<<endl;
		maxh = his_temp->GetMaximum();

		his_temp->GetXaxis()->SetRange(0,4000*mul);
//		gfit->SetRange((maxx-5),(maxx+5));
//		agfit->SetRange((maxx-5),(maxx+5));
//		m1fit->SetRange((maxx-10),(maxx+10));
//		m2fit->SetRange((maxx-10),(maxx+10));
//		m3fit->SetRange((maxx-10),(maxx+10));
		gfit->SetRange((maxx-3*sigma_cal),(maxx+3*sigma_cal));
		agfit->SetRange((maxx-3*sigma_cal),(maxx+3*sigma_cal));
		m1fit->SetRange((maxx-6*sigma_cal),(maxx+6*sigma_cal));
		m2fit->SetRange((maxx-6*sigma_cal),(maxx+6*sigma_cal));
		m3fit->SetRange((maxx-6*sigma_cal),(maxx+6*sigma_cal));

		efit->SetRange((maxx-30),(maxx+30));
		p0fit->SetRange((maxx-30),(maxx+30));
		p1fit->SetRange((maxx-30),(maxx+30));

		his_temp->Fit("p0fit","RQ0+");
		fixbg1 = p0fit->GetParameter(0);

		his_temp->Fit("efit","RQ0+");
		fixbg2[0] = efit->GetParameter(0);
		fixbg2[1] = efit->GetParameter(1);

		his_temp->Fit("p1fit","RQ0+");
		fixbg3[0] = p1fit->GetParameter(0);
		fixbg3[1] = p1fit->GetParameter(1);
		//gaus fit
		his_temp->Fit("gfit","RQ0+");
		gfitpa[0] = gfit->GetParameter(0);
		gfitpa[1] = gfit->GetParameter(1);
		gfitpa[2] = gfit->GetParameter(2);
		gfiter[0] = gfit->GetParError(0);
		gfiter[1] = gfit->GetParError(1);
		gfiter[2] = gfit->GetParError(2);

		gfitchi2 = (gfit->GetChisquare())/(gfit->GetNDF());

		bg1 = 0; bg2 = 0; 
		int bgbin = (((int) sigma_cal*3)+1)*mul;
		cout<<"bgbin : "<<bgbin<<endl;
		//for(int j1=0;j1<(5*mul);j1++){
		for(int j1=0;j1<bgbin;j1++){
		//	bg1 = bg1 + (his_temp->GetBinContent((maxb)-(gfitpa[2]*3+5)*mul+j1));
		//	bg2 = bg2 + (his_temp->GetBinContent((maxb)+(gfitpa[2]*3+5)*mul-j1));
			bg1 = bg1 + (his_temp->GetBinContent(maxb-bgbin+j1));
			bg2 = bg2 + (his_temp->GetBinContent(maxb+bgbin-j1));
			//		cout<<maxx-10+j<<"kev: "<<bg1<<" / "<<maxx+10-j<<"kev : "<<bg2<<endl;
		}
		//fixbg1=(bg1+bg2)/(10*mul); //bg avr.

		//area gaus fit

		//agfit->SetParameters(2.5*gfitpa[1]*gfitpa[2],gfitpa[1],gfitpa[2]);
		agfit->SetParameters(2.5*gfitpa[1]*gfitpa[2],peak[i],sigma_cal);
		//agfit->SetParameters(2.5*gfitpa2*gfitpa3,maxx,2);
		his_temp->Fit("agfit","R0Q+");
		agfitpa[0] = agfit->GetParameter(0);
		agfitpa[1] = agfit->GetParameter(1);
		agfitpa[2] = agfit->GetParameter(2);
		agfiter[0] = agfit->GetParError(0);
		agfiter[1] = agfit->GetParError(1);
		agfiter[2] = agfit->GetParError(2);
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
		//m1fit->SetParameters(agfitpa1,agfitpa2,agfitpa3,bg);
		m1fit->SetParameters(agfitpa[0],peak[i],sigma_cal,fixbg1);


		//peak check
		his_temp->Fit("m1fit","RQ0+");
		m1fitpa[1] = m1fit->GetParameter(1); //mean
		if(abs(m1fitpa[1]-peak[i])>3*(sigma_cal)){
			cout <<"m1 test fitting mean "<<m1fitpa[1]<<" for peak "<<peak[i]<<endl;
			cout<<">> peak fixed at"<<peak[i]<<endl;
			m1fit->FixParameter(1,maxx);
		}

		//sigma check
		his_temp->Fit("m1fit","RQ0+");
		m1fitpa[2] = m1fit->GetParameter(2); //sigma
		if((m1fitpa[2]<(0.8*sigma_cal))||(m1fitpa[2]>(3*sigma_cal))){
			cout <<"m1 test fitting sigma "<<m1fitpa[2]<<" for peak "<<peak[i]<<endl;
			cout<<">> sigma fixed at "<<sigma_cal<<endl;
			m1fit->FixParameter(2,sigma_cal);
		}

		//re-check peak
		his_temp->Fit("m1fit","RQ0+");
		m1fitpa[1] = m1fit->GetParameter(1); //mean
		if(abs(m1fitpa[1]-peak[i])>3*(sigma_cal)){
			cout <<"m1 test fitting mean "<<m1fitpa[1]<<" for peak "<<peak[i]<<endl;
			cout<<">> peak fixed at"<<peak[i]<<endl;
			m1fit->FixParameter(1,maxx);
		}

		his_temp->Fit("m1fit","R+");
		m1fitpa[0] = m1fit->GetParameter(0);
		m1fitpa[1] = m1fit->GetParameter(1);
		m1fitpa[2] = m1fit->GetParameter(2);
		m1fitpa[3] = m1fit->GetParameter(3);
		m1fiter[0] = m1fit->GetParError(0);
		m1fiter[1] = m1fit->GetParError(1);
		m1fiter[2] = m1fit->GetParError(2);
		m1fiter[3] = m1fit->GetParError(3);
		m1fitchi2 = (m1fit->GetChisquare())/(m1fit->GetNDF());

		m1fitpara->Fill();


		cout<<""<<endl;
		cout<<"[m2 fitting]"<<endl;
		cout<<""<<endl;
		//m2fit

		m2fit->ReleaseParameter(0);
		m2fit->ReleaseParameter(1);
		m2fit->ReleaseParameter(2);
		m2fit->ReleaseParameter(3);
		m2fit->ReleaseParameter(4);

		m2fit->SetParameters(agfitpa[0],maxx,sigma_cal,fixbg2[0],fixbg2[1]);
		m2fit->FixParameter(3,fixbg2[0]);
		m2fit->FixParameter(4,fixbg2[1]);

		//peak check
		his_temp->Fit("m2fit","RQ0+");
		m2fitpa[1] = m2fit->GetParameter(1);
		if(abs(m2fitpa[1]-peak[i])>3*(sigma_cal)){
			cout <<"m2 test fitting mean "<<m2fitpa[1]<<" for peak "<<peak[i]<<endl;
			cout<<">> peak fixed at"<<peak[i]<<endl;
			m2fit->FixParameter(1,maxx);
		}

		//sigma check
		his_temp->Fit("m2fit","RQ0+");
		m2fitpa[2] = m2fit->GetParameter(2);
		if((m2fitpa[2]<(0.8*sigma_cal))||(m2fitpa[2]>(3*sigma_cal))){
			cout <<"m2 test fitting sigma "<<m2fitpa[2]<<" for peak "<<peak[i]<<endl;
			cout<<">> sigma fixed at "<<sigma_cal<<endl;
			m2fit->FixParameter(2,sigma_cal);
		}

		//peak re-check
		his_temp->Fit("m2fit","RQ0+");
		m2fitpa[1] = m2fit->GetParameter(1);
		if(abs(m2fitpa[1]-peak[i])>3*(sigma_cal)){
			cout <<"m2 test fitting mean "<<m2fitpa[1]<<" for peak "<<peak[i]<<endl;
			cout<<">> peak fixed at"<<peak[i]<<endl;
			m2fit->FixParameter(1,maxx);
		}

		his_temp->Fit("m2fit","R+");

		m2fitpa[0] = m2fit->GetParameter(0);
		m2fitpa[1] = m2fit->GetParameter(1);
		m2fitpa[2] = m2fit->GetParameter(2);
		m2fitpa[3] = m2fit->GetParameter(3);
		m2fitpa[4] = m2fit->GetParameter(4);
		m2fiter[0] = m2fit->GetParError(0);
		m2fiter[1] = m2fit->GetParError(1);
		m2fiter[2] = m2fit->GetParError(2);
		m2fiter[3] = m2fit->GetParError(3);
		m2fiter[4] = m2fit->GetParError(4);
		m2fitchi2 = (m2fit->GetChisquare())/(m2fit->GetNDF());

		/*
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
		   */
		m2fitpara->Fill();


		cout<<""<<endl;
		cout<<"[m3 fitting]"<<endl;
		cout<<""<<endl;
		//m3fit

		m3fit->ReleaseParameter(0);
		m3fit->ReleaseParameter(1);
		m3fit->ReleaseParameter(2);
		m3fit->ReleaseParameter(3);
		m3fit->ReleaseParameter(4);

		m3fit->SetParameters(agfitpa[0],maxx,sigma_cal,fixbg3[0],fixbg3[1]);
		m3fit->FixParameter(3,fixbg3[0]);
		m3fit->FixParameter(4,fixbg3[1]);

		//peak check
		his_temp->Fit("m3fit","RQ0+");
		m3fitpa[1] = m3fit->GetParameter(1);
		if(abs(m3fitpa[1]-peak[i])>3*(sigma_cal)){
			cout <<"m3 test fitting mean "<<m3fitpa[1]<<" for peak "<<peak[i]<<endl;
			cout<<">> peak fixed at"<<peak[i]<<endl;
			m3fit->FixParameter(1,maxx);
		}

		//sigma check
		his_temp->Fit("m3fit","RQ0+");
		m3fitpa[2] = m3fit->GetParameter(2);
		if((m3fitpa[2]<(0.8*sigma_cal))||(m3fitpa[2]>(3*sigma_cal))){
			cout <<"m3 test fitting sigma "<<m3fitpa[2]<<" for peak "<<peak[i]<<endl;
			cout<<">> sigma fixed at "<<sigma_cal<<endl;
			m3fit->FixParameter(2,sigma_cal);

		}

		//peak re-check
		his_temp->Fit("m3fit","RQ0+");
		m3fitpa[1] = m3fit->GetParameter(1);
		if(abs(m3fitpa[1]-peak[i])>3*(sigma_cal)){
			cout <<"m3 test fitting mean "<<m3fitpa[1]<<" for peak "<<peak[i]<<endl;
			cout<<">> peak fixed at"<<peak[i]<<endl;
			m3fit->FixParameter(1,maxx);
		}

		his_temp->Fit("m3fit","R+");

		m3fitpa[0] = m3fit->GetParameter(0);
		m3fitpa[1] = m3fit->GetParameter(1);
		m3fitpa[2] = m3fit->GetParameter(2);
		m3fitpa[3] = m3fit->GetParameter(3);
		m3fitpa[4] = m3fit->GetParameter(4);
		m3fiter[0] = m3fit->GetParError(0);
		m3fiter[1] = m3fit->GetParError(1);
		m3fiter[2] = m3fit->GetParError(2);
		m3fiter[3] = m3fit->GetParError(3);
		m3fiter[4] = m3fit->GetParError(4);
		m3fitchi2 = (m3fit->GetChisquare())/(m3fit->GetNDF());

		m3fitpara->Fill();
		//find best fit

		cout<<""<<endl;
		cout<<""<<endl;
		cout<<">>chi squre check<<"<<endl;
		cout<<""<<endl;
		cout<<"m1fit chi^2/NDF = "<<m1fitchi2<<endl;
		cout<<"m2fit chi^2/NDF = "<<m2fitchi2<<endl;
		cout<<"m3fit chi^2/NDF = "<<m3fitchi2<<endl;
		cout<<""<<endl;

		if((abs(1.-m1fitchi2)<abs(1.-m2fitchi2)) && (abs(1.-m1fitchi2)<abs(1.-m3fitchi2))){
			cout<<"bestfit : m1 fit "<<endl;
			fit_type = 1;
			fitpa[0] = m1fitpa[0];
			fitpa[1] = m1fitpa[1];
			fitpa[2] = m1fitpa[2];
			fiter[0] = m1fiter[0];
			fiter[1] = m1fiter[1];
			fiter[2] = m1fiter[2];
			fitchi2 = m1fitchi2;
//			c2->cd();
		}
		else if((abs(1.-m2fitchi2)<abs(1.-m1fitchi2)) && (abs(1.-m2fitchi2)<abs(1.-m3fitchi2))){
			cout<<"bestfit : m2 fit "<<endl;
			fit_type = 2;
			fitpa[0] = m2fitpa[0];
			fitpa[1] = m2fitpa[1];
			fitpa[2] = m2fitpa[2];
			fiter[0] = m2fiter[0];
			fiter[1] = m2fiter[1];
			fiter[2] = m2fiter[2];
			fitchi2 = m2fitchi2;
//			c2->cd();
		}
		else if((abs(1.-m3fitchi2)<abs(1.-m1fitchi2)) && (abs(1.-m3fitchi2)<abs(1.-m2fitchi2))){
			cout<<"bestfit : m3 fit "<<endl;
			fit_type = 3;
			fitpa[0] = m3fitpa[0];
			fitpa[1] = m3fitpa[1];
			fitpa[2] = m3fitpa[2];
			fiter[0] = m3fiter[0];
			fiter[1] = m3fiter[1];
			fiter[2] = m3fiter[2];
			fitchi2 = m3fitchi2;
//			c2->cd();
		}

		cout<<fitpa[0]<<endl;
		cout<<fitpa[1]<<endl;
		cout<<fitpa[2]<<endl;
		fitpara->Fill();
		cout<<""<<endl;
		cout<<""<<endl;

	}

	char histitle[256];
	sprintf(histitle,"source peaks Fitting");
	his_temp->SetTitle(histitle);
	his_temp->SetXTitle("Energy[keV]");
	his_temp->SetYTitle("Counts/kev");

	c1->cd();
	his_temp->Draw();
	TFile r3f (res3file,"RECREATE");
	r3f.cd();
	gfitpara->Write();
	m1fitpara->Write();
	m2fitpara->Write();
	m3fitpara->Write();
	fitpara->Write();
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
	c1->SaveAs(res2file);
}
