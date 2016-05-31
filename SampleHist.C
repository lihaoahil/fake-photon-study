#include<iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"

void SampleHist(){ //main

  TFile *gammaTemp = new TFile("template/gammaTemplateMC.root");
  TFile *hadronTemp = new TFile("template/DYhadronicTemplateMCv2.root");
 

  TH1D *data = (TH1D*)hadronTemp->Get("targetEBSigmaIEtaIEta");
  TH1D *num = (TH1D*)hadronTemp->Get("hEBSigmaIEtaIEta");
  TH1D *denom = (TH1D*)gammaTemp->Get("GammaEBSigmaIEtaIEta");
  data->SetLineColor(kBlack);
  denom->SetLineColor(kRed);
  num->SetLineColor(kBlue);


  TH1D *data_scaled = (TH1D*)hadronTemp->Get("targetEBSigmaIEtaIEta");
  TH1D *num_scaled = (TH1D*)hadronTemp->Get("hEBSigmaIEtaIEta");
  TH1D *denom_scaled = (TH1D*)gammaTemp->Get("GammaEBSigmaIEtaIEta");
  data_scaled->SetLineColor(kBlack);
  denom_scaled->SetLineColor(kRed);
  num_scaled->SetLineColor(kBlue);
  double norm = data -> GetEntries();
      double datascale = norm/data -> Integral();
      double denomscale = norm/denom_scaled -> Integral();
      double numscale = norm/num_scaled -> Integral();

      //data -> TH1::Sumw2();
      data_scaled -> TH1::Scale(datascale);
      denom_scaled -> TH1::Scale(denomscale);
      num_scaled -> TH1::Scale(numscale);
  
  TCanvas *c1 = new TCanvas("c1","",1600,900);
  
  c1->SetGrid();
 
  c1->Divide(2);
  
  data->SetTitle("Samples_All_Pt_range");
  

  c1->cd(1);
  data->Draw();
  denom->Draw("same");
  num->Draw("same");

  TLegend *leg =  new TLegend(0.78,0.65,0.98,0.75);
    leg->AddEntry(data,"Fit Target");
    leg->AddEntry(denom,"Photon Template");
    leg->AddEntry(num,"Hadron Template");
    leg->Draw();

  c1->cd(2);

  data_scaled->SetTitle("Scaled samples_All_Pt_range");
  
  data_scaled->Draw();
  denom_scaled->Draw("same");
  num_scaled->Draw("same");


  c1->SaveAs("test_result_DY_v2/SampleHist.png");
  //delete c1;










}//end of the main function