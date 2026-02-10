#include <iostream>
#include "TMath.h"
#include "TMinuit.h"


void draw_run650_check() {



  double respa1[3]={0.498182,0.000373447 ,4.3987e-08 };
  double respa2[3]={0.4457, 0.000605946, -6.31339e-08};

  TF1 *res1 = new TF1("res1", "sqrt([0]*[0]+[1]*x+[2]*x*x)", 500, 1800);
  TF1 *res2 = new TF1("res2", "sqrt([0]*[0]+[1]*x+[2]*x*x)", 500, 1800);
  res1->SetParameters(respa1[0], respa1[1], respa1[2]);
  res2->SetParameters(respa2[0], respa2[1], respa2[2]);

  res1->SetLineColor(kRed);
  res2->SetLineColor(kBlue);
  // 점 데이터
  const int dp=3;
  /*
  //old cal set
  double x[dp] = {6.09126e+02, 7.67416e+02,  9.11054e+02, 1.12001e+03};
  double ex[dp] = {2.18868e-02, 1.79930e-01, 9.11295e-02, 5.09417e-02};
  double y[dp] = { 8.52464e-01 ,   2.67469e+00, 1.02395e+00, 9.33861e-01};
  double ey[dp] = {2.35830e-02,2.91398e-01,  1.19865e-01, 5.67681e-02};
  */
  //new cal set v1
  double x[dp] = {6.09416e+02, 9.11297e+02, 1.12039e+03 };
  double ex[dp] = { 1.09171e-02, 6.32965e-02 ,  2.74479e-02 };
  double y[dp] = { 7.17671e-01,  0.856 ,  8.44341e-01};
  double ey[dp] = { 9.70076e-03, 0.0105,  2.66512e-02 };
  

  TGraphErrors *points = new TGraphErrors(dp, x, y, ex, ey);

  points->SetMarkerStyle(21); 
  points->SetMarkerSize(1.5);
  points->SetMinimum(0.5);
  points->SetMaximum(1.6);
  //  points->SetMarkerColor(kRed);

  // 캔버스
  TCanvas *c1 = new TCanvas("c1", "Graph Example", 800, 600);
  points->GetXaxis()->SetTitle("Energy[keV]");
  points->GetYaxis()->SetTitle("Sigma[keV]");

 points->Draw("AP"); // P = point  
  // 먼저 라인 하나 그리고
 //  res1->Draw("SAME");
  cout<<"test1"<<endl;
  // 두 번째 라인 덧그리기
  res2->Draw("SAME");
  cout<<"test2"<<endl;
  // 점 그리기
 
  cout<<"test3"<<endl;


  c1->Modified();
  c1->Update();
}
