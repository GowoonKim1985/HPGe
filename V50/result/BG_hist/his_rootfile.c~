#include "TFile.h"
#include "TH1.h"
#include <iostream>
#include <vector>
#include <string>

void his_rootfile() {
    // 1️⃣ 원본 ROOT 파일 열기
    TFile *fin = TFile::Open("his_000650.v1.root", "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "입력 ROOT 파일을 열 수 없습니다." << std::endl;
        return;
    }

    // 2️⃣ 새 ROOT 파일 만들기
    TFile *fout = new TFile("his_th1.root", "RECREATE");

    // 3️⃣ 일반 his1~his14 중 필요한 것만 복사
    for (int i = 1; i <= 14; i++) {
        if (i == 2 || i == 4 || i == 10) continue; // 제외 목록

        char hisname[64];
        sprintf(hisname, "his%i", i);

        TH1 *h = (TH1*)fin->Get(hisname);
        if (!h) {
            std::cerr << hisname << " 히스토그램을 찾을 수 없습니다." << std::endl;
            continue;
        }

        // Y축 타이틀 변경
        h->GetYaxis()->SetTitle("Count/(0.25 keV)");
	
        fout->cd();  // 출력 파일로 전환
        h->Write();  // 같은 이름으로 저장
        std::cout << hisname << " 저장 완료" << std::endl;
    }

    // 4️⃣ his_tot3을 his_tot으로 이름 바꿔서 저장
    TH1 *h_tot = (TH1*)fin->Get("his_tot3");
    if (h_tot) {
        fout->cd();
	h_tot->GetYaxis()->SetTitle("Count/(0.25 keV)");

        h_tot->SetName("his_tot");  // 이름 변경
        h_tot->Write();
        std::cout << "his_tot3 → his_tot 이름으로 저장 완료" << std::endl;
    } else {
        std::cerr << "his_tot3 히스토그램을 찾을 수 없습니다." << std::endl;
    }

    // 5️⃣ 마무리
    fout->Close();
    fin->Close();

    std::cout << "\n✅ 선택된 히스토그램이 selected_hists.root 에 저장되었습니다." << std::endl;
}
