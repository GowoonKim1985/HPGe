#include "TFile.h"
#include "TH1.h"
#include <fstream>

void hist_to_ascii() {
    TFile *f = TFile::Open("his_000657.v1.root");
    if (!f || f->IsZombie()) return;

    TH1 *h = (TH1*)f->Get("his_tot3");

    char hisname[256];
    char txtname[256];

    //    for(int i=1; i<=14; i++){
    //      if(i!=2&&i!=4&&i!=10){
    //	sprintf(hisname,"his%i",i);
    //	sprintf(txtname,"his_det_%i.txt",i);
    //    TH1 *h = (TH1*)f->Get(hisname);

    if (!h) return;

    std::ofstream out("his_det_all.txt");
    //    std::ofstream out(txtname);    
    for (int i = 1; i <= h->GetNbinsX(); i++) {
        out << i << "\t" << h->GetBinContent(i) << "\n";
    }
    out.close();
    //    }
    //    }
}
