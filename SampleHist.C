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
  denom->SetTitle("Samples_All_Pt_range");


  TH1D *dataScaled = (TH1D*)hadronTemp->Get("targetEBSigmaIEtaIEta");
  TH1D *numScaled = (TH1D*)hadronTemp->Get("hEBSigmaIEtaIEta");
  TH1D *denomScaled = (TH1D*)gammaTemp->Get("GammaEBSigmaIEtaIEta");
  dataScaled->SetLineColor(kBlack);
  denomScaled->SetLineColor(kRed);
  numScaled->SetLineColor(kBlue);
  denomScaled->SetTitle("Scaled samples_All_Pt_range");

  double norm = dataScaled -> GetEntries();
      double datascale = norm/dataScaled -> Integral();
      double denomscale = norm/denomScaled -> Integral();
      double numscale = norm/numScaled -> Integral();
      std::cout<< norm <<" " << datascale <<" "<<denomscale<<" "<<numscale<<std::endl;
      //data -> TH1::Sumw2();
      dataScaled -> TH1::Scale(datascale);
      denomScaled -> TH1::Scale(denomscale);
      numScaled -> TH1::Scale(numscale);
  
  TCanvas *canvas = new TCanvas("canvas1","",1600,900);
  canvas->SetGrid(); 
  canvas->Divide(2);
  
  

  canvas->cd(1);
  denom->Draw("");
  data->Draw("same");
  num->Draw("same");
  TLegend *leg =  new TLegend(0.78,0.65,0.98,0.75);
    leg->AddEntry(data,"Fit Target");
    leg->AddEntry(denom,"Photon Template");
    leg->AddEntry(num,"Hadron Template");
    leg->Draw();



  canvas->cd(2);
  denomScaled->Draw("");
  dataScaled->Draw("same");
  numScaled->Draw("same");


  canvas->SaveAs("test_result_DY_v2/SampleHist.png");
  delete canvas;










}//end of the main function