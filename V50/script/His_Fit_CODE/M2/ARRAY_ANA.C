//R__LOAD_LIBRARY(libGui)
//R__LOAD_LIBRARY(libFramework)
//R__LOAD_LIBRARY(libRawObjs)
//R__LOAD_LIBRARY(libCupObjs)

{
	char filename[256];
	cout<<"file name[without .root] :";
	//cout<< "K40? U238? Th232? Y88? Co57? Co58? Co60? Zn65? Cs137? Mn54? Mo100?: ";
	scanf("%s",filename);

	double e1,e2,e3;


	// Th coincidence
	e1= 583.187;
	e2= 2614.511;
/*	
	// Y coincidence
	e1= 898.042;
	e2= 1836.063;
	

	// Mo coincidence
	e1=539.512;
	e2=590.782;
*/

	e3= e1+e2;

	char datafile[256];
	sprintf(datafile,"/home/kkw/geant4.9.6.p02/work/SIM_v4/HpgeSim/result/ARRAY/%s.root",filename);

	char outfile_tot[256];
	sprintf(outfile_tot,"/home/kkw/geant4.9.6.p02/work/SIM_v4/HpgeSim/result/ARRAY/his_%s_tot.root",filename);
//	sprintf(outfile_tot, "../../results/cmd108/final/%s_eff_simul_tot.root", filename);

	char outfile1[256];
	sprintf(outfile1,"/home/kkw/geant4.9.6.p02/work/SIM_v4/HpgeSim/result/ARRAY/his_%s_1mul.root",filename);
//	sprintf(outfile1, "../../results/cmd108/final/%s_eff_simul_1mul.root",filename);

	char outfile2[256];
	sprintf(outfile2,"/home/kkw/geant4.9.6.p02/work/SIM_v4/HpgeSim/result/ARRAY/his_%s_2mul.root",filename);
//	sprintf(outfile2, "../../results/cmd108/final/%s_eff_simul_2mul.root",filename);
/*
	char outfile_mul[256];
	sprintf(outfile_mul,"/home/kkw/geant4.9.6.p02/work/SIM_v4/HpgeSim/result/ARRAY/ana_%s_mul2.txt",filename);
//	sprintf(outfile_mul, "../../results/cmd108/final/%s_eff_simul_num.txt",filename);

	FILE *mul_txt = fopen(outfile_mul,"w");
	//////////root file load
*/
	TChain * chain = new TChain("event_tree");
	chain->Add(datafile);

	TLeaf * fNvertex = chain->GetLeaf("fNvertex");
	TLeaf * TotEdep = chain->GetLeaf("TotEdep");
	TLeaf * fCellEdep = chain->GetLeaf("fCell.Edep");
	TLeaf * fCellIndex = chain->GetLeaf("fCell.Index");
	TLeaf * fVertexke = chain->GetLeaf("fVertex.ke");
	TLeaf * fTrackID = chain->GetLeaf("fTrack.TrackID"); //id1: 215&234, id2:332&350


	double MinE = 0;
	double MaxE = 4000;
	//double binN = MaxE*2; //0.5kev
	double binN = MaxE; //1kev
	//	double binN = MaxE*1; //1kev
	char detname[256];
	char mul1detname[256];
	char mul2detname[256];
	char coindetname[256];
	char detname_det[256];
	char mul1detname_det[256];

	//tot graph
	TH1F *tot;
	tot = new TH1F("tot","All event with 1keV bin",binN, MinE, MaxE);
	tot->GetXaxis()->SetTitle("Energy(keV)");
	tot->GetYaxis()->SetTitle("Count");

	TH1F *det[14];
	for(int i = 0; i<14; i++) {
		sprintf(detname, "det%i", i);
		sprintf(detname_det, "All event. Det %i with 1keV bin",i);
		det[i] =new TH1F(detname,detname_det, binN, MinE, MaxE);
		det[i]->GetXaxis()->SetTitle("Energy(keV)");
		det[i]->GetYaxis()->SetTitle("Count");
	}

	//mul1 graph
	TH1F *mul1_tot;
	mul1_tot = new TH1F("mul1_tot", "1MUL event with 1keV bin", binN, MinE, MaxE);
	mul1_tot->GetXaxis()->SetTitle("Energy(keV)");
	mul1_tot->GetYaxis()->SetTitle("Count");
	
	TH1F *mul1_det[14];
	for(int i = 0; i<14; i++) {
		sprintf(mul1detname,"mul1_det%i",i);
		sprintf(mul1detname_det, "1MUL event. Det %i with 1keV bin",i);
		mul1_det[i] =new TH1F(mul1detname,mul1detname_det, binN, MinE, MaxE);
		mul1_det[i]->GetXaxis()->SetTitle("Energy(keV)");
		mul1_det[i]->GetYaxis()->SetTitle("Count");
	}

	//mul2 graph
	TH1F *mul2_1dim;
	mul2_1dim = new TH1F("mul2_1dim", "", binN, MinE, MaxE);
	mul2_1dim->GetXaxis()->SetTitle("Energy(keV)");
	mul2_1dim->GetYaxis()->SetTitle("Count");
	
	TH2F *mul2;
	mul2 = new TH2F("mul2","", binN, MinE, MaxE, binN, MinE, MaxE);
	mul2->GetXaxis()->SetTitle("Energy1(keV)");
	mul2->GetYaxis()->SetTitle("Energy2(kev)");

	// mul2 dec graph
	TH1F *mul2_det[14];
	for(int j = 0; j<14; j++){
		sprintf(mul2detname, "mul2_det%i", j);
		mul2_det[j] =new TH1F(mul2detname, "", 14, 0, 14);
		mul2_det[j]->GetXaxis()->SetTitle("Det.");
		mul2_det[j]->GetYaxis()->SetTitle("Count");
	}

	//mul2 coin graph
	TH1F *coin_det[14];
	for(int k = 0; k<14; k++){
		sprintf(coindetname, "coin_det%i", k);
		coin_det[k] =new TH1F(coindetname, "", 14, 0, 14);
		coin_det[k]->GetXaxis()->SetTitle("Det.");
		coin_det[k]->GetYaxis()->SetTitle("Count");
	}

	// mul2 dec ntuple
	double mul2_e1, mul2_e2, mul2_d1, mul2_d2, P_info;
	TNtuple * mul2_para = new TNtuple("mul2_para", "mul2", "mul2_e1:mul2_e2:mul2_d1:mul2_d2:P_info");

	// mul3 ntuple
	double mul3_e1, mul3_e2, mul3_e3, mul3_d1, mul3_d2, mul3_d3;
	TNtuple * mul3_para = new TNtuple("mul3_para", "mul3", "mul3_e1:mul3_e2:mul3_e3:mul3_d1:mul3_d2:mul3_d3:mul3_tot");


	int nent = event_tree->GetEntries();

	int nvertex, maxcell, cell, proc;
	double totedep, size;
	int cellindex[14]; 
	double celledep[14], ke[14], trackid[14];

	int count1 = 0;
	int count2 = 0;
	int count3 = 0;
	int count4 = 0;

	maxcell = fCellEdep->GetMaximum();

	proc=0;

	for(int i = 0; i < nent; i++){

		if(proc == i){
			cout <<"process...("<<proc<<"/"<<nent<<")"<<endl;
			proc=proc+100000;
		}

		chain->GetEntry(i);
		nvertex = fNvertex->GetValue();
		totedep=TotEdep->GetValue();

		if(totedep>0){

			cell = fCellEdep->GetLen();

			for(int ii=0; ii<cell; ii++) {

				celledep[ii] = (fCellEdep->GetValue(ii))*1000.;
				cellindex[ii] = fCellIndex->GetValue(ii);
				ke[ii] = fVertexke->GetValue(ii)*1000;
				trackid[ii] = fTrackID->GetValue(ii);

				det[cellindex[ii]]->Fill(celledep[ii]); // tot
				tot->Fill(celledep[ii]);
			}

			if(cell==1){ //MUL1

				mul1_det[cellindex[0]]->Fill(celledep[0]);
				mul1_tot->Fill(celledep[0]);
		
				count1 = count1+1;
				//cout<<"cnt1 : "<<count1<<endl;
			}
//MUL 2 part

			if(cell==2) { //MUL2


				if (celledep[0] >= celledep[1]) { mul2->Fill(celledep[0], celledep[1]); }
				if (celledep[0] <  celledep[1]) { mul2->Fill(celledep[1], celledep[0]); }

				mul2_1dim->Fill(celledep[0]);
				mul2_1dim->Fill(celledep[1]);

				//	mul2->Fill(ke2,celledep2);
				mul2_det[cellindex[0]]->Fill(cellindex[1]);
				mul2_det[cellindex[1]]->Fill(cellindex[0]);

				//coincidence check
				if( ( celledep[0]+celledep[1] > e3-2 ) && (celledep[0]+celledep[1] < e3+2 ) ) { // E1+E2 = e1+e2 +/-2kev

					if( ( celledep[0] > e1-1 ) && ( celledep[0] < e1+1) ) {	// E1 = e1 +/- 1kev
						mul2_para->Fill(celledep[0], celledep[1], cellindex[0], cellindex[1], 1);
						coin_det[cellindex[0]]->Fill(cellindex[1]);
						coin_det[cellindex[1]]->Fill(cellindex[0]);
					}

					else if( (celledep[0] > e2-1) && (celledep[0] < e2+1 ) ) {	// E1 = e2 +/- 1kev
						mul2_para->Fill(celledep[0],celledep[1],cellindex[0],cellindex[1],1);
						coin_det[cellindex[0]]->Fill(cellindex[1]);
						coin_det[cellindex[1]]->Fill(cellindex[0]);
					}

					else {
						mul2_para->Fill(celledep[0], celledep[1], cellindex[0], cellindex[1], 2);
					}
				}

				else {
					mul2_para->Fill(celledep[0],celledep[1],cellindex[0],cellindex[1],5);
				}
				count2 = count2+1;
				//cout<<"cnt2 : "<<count2<<endl;
			}


			//	if(cell==3) { mul3_para->Fill(celledep[0],celledep[1],celledep[2],cellindex[0],cellindex[1],cellindex[2],(totedep*1000)); }

			//cout<<size<<endl;
		}
	}


	//output

	//tot
	TFile *write = new TFile(outfile_tot, "RECREATE");

//	can_tot->cd();
//	tot->Draw();
	tot->Write();

	for(int i=0; i<14; i++) {
//		can->cd(i+1);
//		det[i]->Draw();
		det[i]->Write();
//		can->SetGrid();
	}
	write->Close();


	//mul1
//	TCanvas *can1 = new TCanvas("can1", "", 1000, 800);
//	can1->Divide(4, 4);
//	can1->SetGrid();

	TFile *write1 = new TFile(outfile1,"RECREATE");

	for(int i = 0; i<14; i++) {
//		can1->cd(i+1);
//		mul1_det[i]->Draw();
		mul1_det[i]->Write();
//		can1->SetGrid();
	}
	mul1_tot->Write();
	write1->Close();


	//mul2 write

	TFile *write2 = new TFile(outfile2,"RECREATE");

	int p1[14][14], pt[14][14];

	for(int i = 0; i<14; i++) {
//		can2_1->cd(i+1);
//		mul2_det[i]->Draw();
		for(int j = 1; j<15; j++) {
			pt[i][j-1] = mul2_det[i]->GetBinContent(j);
		}
		mul2_det[i]->Write();
//		can2_1->SetGrid();

//		can2_2->cd(i+1);
//		coin_det[i]->Draw();
		for(int k = 1; k<15; k++) {
			p1[i][k-1] = coin_det[i]->GetBinContent(k);
		}
		coin_det[i]->Write();
//		can2_2->SetGrid();
	}

	//mul_txt = fopen(outfile_mul,"w");
/*
	printf("\nMUL2 counts	\n\n"
			"MUL2 T0 T1 T2 T3 T4 T5 T6 T7 T8 T9 T10 T11 T12 T13\n");
	fprintf(mul_txt,"\nMUL2 counts	\n\n"
			"MUL2 T0 T1 T2 T3 T4 T5 T6 T7 T8 T9 T10 T11 T12 T13\n");

	for(int l=0;l<14;l++){
		printf("T%i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n"
				,l,pt[l][0],pt[l][1],pt[l][2],pt[l][3],pt[l][4],pt[l][5],pt[l][6],pt[l][7],pt[l][8],pt[l][9],pt[l][10],pt[l][11],pt[l][12],pt[l][13]);
		fprintf(mul_txt,"T%i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n"
				,l,pt[l][0],pt[l][1],pt[l][2],pt[l][3],pt[l][4],pt[l][5],pt[l][6],pt[l][7],pt[l][8],pt[l][9],pt[l][10],pt[l][11],pt[l][12],pt[l][13]);
	}

	printf("\nMUL2 P1 counts	\n\n"
			"P1 T0 T1 T2 T3 T4 T5 T6 T7 T8 T9 T10 T11 T12 T13\n");
	fprintf(mul_txt, "\nMUL2 P1 counts	\n\n"
			"P1 T0 T1 T2 T3 T4 T5 T6 T7 T8 T9 T10 T11 T12 T13\n");

	for(int l=0;l<14;l++){
		printf("T%i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n"
				,l,p1[l][0],p1[l][1],p1[l][2],p1[l][3],p1[l][4],p1[l][5],p1[l][6],p1[l][7],p1[l][8],p1[l][9],p1[l][10],p1[l][11],p1[l][12],p1[l][13]);
		fprintf(mul_txt, "T%i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n"
				,l,p1[l][0],p1[l][1],p1[l][2],p1[l][3],p1[l][4],p1[l][5],p1[l][6],p1[l][7],p1[l][8],p1[l][9],p1[l][10],p1[l][11],p1[l][12],p1[l][13]);
	}
*/
//	can2->cd();
//	mul2->Draw();
	mul2->Write();
//	can2_1dim->cd();
//	mul2_1dim->Draw();
	mul2_1dim->Write();
	mul2_para->Write();
	write2->Close();

/*
	//mul3
	TFile *write3 = new TFile(outfile3,"RECREATE");
	int mul3_evt;
	double mul3_lnum[7];
	mul3_evt = mul3_para->GetEntries();

	int mul3_count1 = 0;
	int mul3_count2 = 0;
			if((mul3_lnum[0]<=(e1+1)&&mul3_lnum[0]>=(e1-1))||
					(mul3_lnum[1]<=(e1+1)&&mul3_lnum[1]>=(e1-1))||
					(mul3_lnum[2]<=(e1+1)&&mul3_lnum[2]>=(e1-1))){

				mul3_count2 = mul3_count2 + 1;	
				printf("%.0f %.2f %.0f %.2f %.0f %.2f %.2f\n",mul3_lnum[3],mul3_lnum[0],mul3_lnum[4],mul3_lnum[1],mul3_lnum[5],mul3_lnum[2],mul3_lnum[6]);
				fprintf(mul_txt,"%.2f %.2f %.2f %.2f %.2f %.2f %.2f\n",mul3_lnum[3],mul3_lnum[0],mul3_lnum[4],mul3_lnum[1],mul3_lnum[5],mul3_lnum[2],mul3_lnum[6]);
			}
		}
	}



	mul3_para->Write();
	write3->Close();
*/
}


