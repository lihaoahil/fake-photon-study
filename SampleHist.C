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
  double norm = data_scaled -> GetEntries();
      double datascale = norm/data_scaled -> Integral();
      double denomscale = norm/denom_scaled -> Integral();
      double numscale = norm/num_scaled -> Integral();
      std::cout<< norm <<" " << datascale <<" "<<denomscale<<" "<<numscale<<std::endl;
      //data -> TH1::Sumw2();
      data_scaled -> TH1::Scale(datascale);
      denom_scaled -> TH1::Scale(denomscale);
      num_scaled -> TH1::Scale(numscale);
  
  TCanvas *canvas = new TCanvas("canvas1","",1600,900);
  
  canvas->SetGrid();
 
  canvas->Divide(2);
  
  

  canvas->cd(2);
  denom->SetTitle("Samples_All_Pt_range");
  denom->Draw();
  data->Draw("same");
  num->Draw("same");

  TLegend *leg =  new TLegend(0.78,0.65,0.98,0.75);
    leg->AddEntry(data,"Fit Target");
    leg->AddEntry(denom,"Photon Template");
    leg->AddEntry(num,"Hadron Template");
    leg->Draw();



  canvas->cd(1);
  denom_scaled->SetTitle("Scaled samples_All_Pt_range");
  denom_scaled->Draw("");
  data_scaled->Draw("same");
  num_scaled->Draw("same");


  canvas->SaveAs("test_result_DY_v2/SampleHist.png");
  delete canvas;










}//end of the main function