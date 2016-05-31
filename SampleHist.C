#include<iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TFractionFitter.h"
#include "TObjArray.h"
#include "TGraphErrors.h"
void SampleHist(){ //main

  TFile *gammaTemp = new TFile("template/gammaTemplateMC.root");
  TFile *hadronTemp = new TFile("template/DYhadronicTemplateMCv2.root");
 

  TH1D *data = (TH1D*)hadronTemp->Get(targetEBSigmaIEtaIEta);
  TH1D *num = (TH1D*)hadronTemp->Get(hEBSigmaIEtaIEta);
  TH1D *denom = (TH1D*)gammaTemp->Get(gammaEBSigmaIEtaIEta);
  
  data->SetLineColor(kBlack);
  denom->SetLineColor(kRed);
  num->SetLineColor(kBlue);
  
  TCanvas *c1 = new TCanvas("c1","",1600,900);
  
  c1->SetGrid();
 
  
  data->SetTitle("Samples_All_Pt_range");
  
  data->Draw();
  denom->Draw("same");
  num->Draw("same");

  TLegend *leg =  new TLegend(0.78,0.65,0.98,0.75);
    leg->AddEntry(data,"Fit Target");
    leg->AddEntry(denom,"Photon Template");
    leg->AddEntry(num,"Hadron Template");
    leg->Draw();

  c1->SaveAs("test_result_DY_v2/SampleHist.png");
  //delete c1;










}//end of the main function