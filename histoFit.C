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
void histoFit(){ //main

  TFile *gammaTemp = new TFile("template/gammaTemplateMC.root");
  TFile *hadronTemp = new TFile("DYhadronicTemplateMC.root");
 

  Int_t binrange[] = {20,22,24,26,28,30,32,34,36,38,42,46,50,60,70,90,140};
  
  TH1D *data;
  TH1D *mc0;
  TH1D *mc1;
  TH1D* result;
  TString dataName, mc0Name, mc1Name, latexFitpar,latexStatus, imageName, histTitle1, histTitle2 ;
  TObjArray *mc;
  TFractionFitter* fit;
  Int_t status;
  Double_t p0, p1, errP0, errP1,numerator,denominator;
  

  //hadron fraction result
  Double_t hadronFraction[16];
  Double_t binPos[16];
  Double_t ex[16];
  Double_t ey[16];
  for(unsigned i(0);i<16;i++)
    {
      ex[i] =(binrange[i+1] - binrange[i])/2;
      ey[i] = 0;
      hadronFraction[i] = -1;
      binPos[i] = (binrange[i]+binrange[i+1])/2;
    }

  Double_t norm, datascale,mc0scale,mc1scale;
  for(unsigned ibin(2); ibin<16;ibin++)
    {
     
      dataName. Form("FitTarget_Pt_bin(%d,%d)",binrange[ibin],binrange[ibin+1]);
      mc0Name. Form("PhotonTemplate_Pt_bin(%d,%d)",binrange[ibin],binrange[ibin+1]);
      mc1Name. Form("HadronicTemplate_Pt_bin(%d,%d)",binrange[ibin],binrange[ibin+1]);
      imageName. Form("test_result_DY/Pt_bin(%d,%d).png",binrange[ibin],binrange[ibin+1]);
      histTitle1. Form("SigmaIEtaIEta template fits in Pt bin (%d,%d)",binrange[ibin],binrange[ibin+1]);
      histTitle2. Form("fit result");
      data = (TH1D*)hadronTemp->Get(dataName);//fit target
      mc0 = (TH1D*)gammaTemp->Get(mc0Name);//photonic sample
      mc1 = (TH1D*)hadronTemp->Get(mc1Name);//hadronic sample

      norm = data -> GetEntries();
      datascale = norm/data -> Integral();
      mc0scale = norm/mc0 -> Integral();
      mc1scale = norm/mc1 -> Integral();

      //data -> TH1::Sumw2();
      data -> TH1::Scale(datascale);
      mc0 -> TH1::Scale(mc0scale);
      mc1 -> TH1::Scale(mc1scale);

      mc = new TObjArray(2);//put histograms in this array
      mc->Add(mc0);
      mc->Add(mc1);
  

      fit = new TFractionFitter(data, mc); //initialise
      fit->Constrain(0,0.0,1.0);//constrain fraction 0 to be between 0 and 1
      fit->Constrain(1,0.0,1.0);//constrain fraction 1 to be between 0 and 1
      status = fit->Fit();// perform the fit;
      std::cout<< "--------------bin:  "<< ibin<<"-----------------" << std::endl;
      std::cout<< "fit status:  "<< status << std::endl;

      if (status == 0) //if the fit does not converge, check on fit result
  {
    result = (TH1D*) fit->GetPlot();
    
    fit->GetResult(0,p0,errP0);
    fit->GetResult(1,p1,errP1);
    std::cout<< p0 <<"   "<< errP0 <<std::endl;
    std::cout<< p1 <<"   "<< errP1 <<std::endl;
    
    numerator = mc1->Integral(0,25);
    denominator = mc0->Integral(0,25); 
    std::cout<< numerator <<"   "<< denominator <<std::endl;
    numerator *= p1;
    denominator *= p0;
    std::cout<< numerator <<"   "<< denominator <<std::endl;
    hadronFraction[ibin] = numerator/denominator;
    ey[ibin] = errP1;
    std::cout<< "fraction:  "<<hadronFraction[ibin]<<"+/-"<<errP1 << std::endl;
    std::cout<< "---------------------------------------" << std::endl;

    latexFitpar . Form("Hadron fraction: %.3f #pm %.3f", p1, errP1);
    latexStatus . Form("Status: %d", status);

    //define the canvas
    TCanvas *canvas = new TCanvas("c1","",1600,900);
    canvas->Divide(2);
    
    //draw original histograms of target, pho, hadron on the left
    canvas->cd(1);
    data->SetLineColor(kBlack);
    mc0->SetLineColor(kRed);
    mc1->SetLineColor(kBlue);

    data -> SetOption("EP");
    data -> SetTitle(histTitle1);

    data->Draw();
    mc0->Draw("same");
    mc1->Draw("same");

    TLegend *leg =  new TLegend(0.78,0.65,0.98,0.75);
    leg->AddEntry(data,"Fit Target");
    leg->AddEntry(mc0,"Photon Template");
    leg->AddEntry(mc1,"Hadron Template");
    leg->Draw();



    //Then draw the fit result on the left
    canvas->cd(2);
    result -> SetTitle(histTitle2);
    
    result->Draw("");
    data->Draw("sameEP");

    TLatex *lat = new TLatex();
    lat->SetTextSize(0.025);
    lat->DrawLatexNDC(0.1,0.86,latexFitpar);
    lat->DrawLatexNDC(0.1,0.88,latexStatus);

    canvas->SaveAs(imageName);
    delete canvas;

    
  }
      
    }

  
  
  TCanvas *c1 = new TCanvas("c1","",1600,900);
  
  c1->SetGrid();
  TGraphErrors *gr = new TGraphErrors(16,binPos,hadronFraction,ex,ey);
  
  //delete point that convergence is not concluded
  for(unsigned j(0);j<16;j++)
    {
      if(hadronFraction[j]<0) gr->RemovePoint(j);
    }
  gr->RemovePoint(0);
  gr->SetTitle("Hadron Fraction");
  gr->GetYaxis()->SetRangeUser(0,1.0);
  gr->Draw("AP");

  c1->SaveAs("test_result_DY/HadronFraction.png");
  //delete c1;










}//end of the main function