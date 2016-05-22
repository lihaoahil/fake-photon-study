#include<iostream>
#include<vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TChain.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "photonCut.h"

void hadronTemplate(){//main  
  
  TChain *es =new TChain("ggNtuplizer/EventTree");
  es->Add("/export/cmss2/mengleis/SingleMuon/SingleMuon-Run2015D-PR4-MINIAOD-part1.root");
  es->Add("/export/cmss2/mengleis/SingleMuon/SingleMuon-Run2015D-PR4-MINIAOD-part2.root");

  //TFile *file = new TFile("/export/home/hl2/fake-photon-study/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-part2.root");
  //TTree *es = (TTree*)file->Get("ggNtuplizer/EventTree");


  TFile *hSample = new TFile("hadronicTemplate.root","RECREATE");



    //reco-photon
    int                    nPho=0;
    std::vector<float>    *phoEta =0;
    std::vector<float>    *phoPhi = 0;
    std::vector<float>    *phoEt = 0;
    std::vector<float>    *phoE = 0;
    std::vector<int>      *phoEleVeto = 0;
    std::vector<float>    *phoHoverE = 0;
    std::vector<float>    *phoSigmaIEtaIEtaFull5x5 = 0;
    std::vector<float>    *phoPFChIso = 0;
    std::vector<float>    *phoPFNeuIso = 0;
    std::vector<float>    *phoPFPhoIso = 0;
    std::vector<float>    *phoR9 = 0;
    std::vector<UShort_t> *phoIDbit = 0;

    //reco-muon/antimuon
    int                   nMu = 0; 
    std::vector<float>    *muPt = 0;
    std::vector<float>    *muEn = 0;
    std::vector<float>    *muEta = 0;
    std::vector<float>    *muPhi = 0;
    std::vector<Int_t>    *muCharge =0;
    std::vector<Int_t>    *muType =0;
    std::vector<bool>    *muIsTightID=0;

    float rho = 0;//rho
    ULong64_t  HLTEleMuX = 0; //HLT

    //Missing transverse energy
    float pfMET = 0;


    UShort_t decision = 0;
    unsigned counter = 0;
    bool isBarrel, isEndCap;

    int flagFSR;
    

//Set branches to photon variables


    es->SetBranchAddress("HLTEleMuX",&HLTEleMuX);
    es->SetBranchAddress("pfMET",&pfMET);
    es->SetBranchAddress("nPho",&nPho);
    es->SetBranchAddress("phoEta",&phoEta);  
    es->SetBranchAddress("phoPhi",&phoPhi);
    es->SetBranchAddress("phoEt",&phoEt);  
    es->SetBranchAddress("phoE",&phoE);
    es->SetBranchAddress("phoEleVeto",&phoEleVeto); //electron veto
    es->SetBranchAddress("phoHoverE",&phoHoverE);   //H over E
    es->SetBranchAddress("phoSigmaIEtaIEtaFull5x5",&phoSigmaIEtaIEtaFull5x5);
    es->SetBranchAddress("phoPFChIso",&phoPFChIso); //Charged hadron isolation
    es->SetBranchAddress("phoPFNeuIso",&phoPFNeuIso);//neutral hadron isolation
    es->SetBranchAddress("phoPFPhoIso",&phoPFPhoIso);//photon isolation
    es->SetBranchAddress("phoIDbit",&phoIDbit); 
    es->SetBranchAddress("rho",&rho);
    es->SetBranchAddress("phoR9",&phoR9);

    es->SetBranchAddress("nMu",&nMu);
    es->SetBranchAddress("muPt",&muPt);
    es->SetBranchAddress("muEn",&muEn);
    es->SetBranchAddress("muEta",&muEta);
    es->SetBranchAddress("muPhi",&muPhi);
    es->SetBranchAddress("muCharge",&muCharge);
    es->SetBranchAddress("muType",&muType);
    es->SetBranchAddress("muIsTightID",&muIsTightID);

    //TH1D *MuMuHist = new TH1D("InvariantMassMuMuHist", "Mu+Mu-", 90, 30, 120);
    //TH1D *MuMuGammaHist = new TH1D("InvariantMassMuMuGammaHist", "Mu+Mu-Gamma", 60, 60, 120);
    //TH1D *MuMuplusMuMuGammaHist = new TH1D("MuMu+MuMuGamma", "Mu+Mu-PlusMu+Mu-Gamma", 100, 90, 240);
    TH1D *hEBSigmaIEtaIEta = new TH1D("hEBSigmaIEtaIEta", "hadron_SigmaIEtaIEta_EB_Sample_Total",100, 0, 0.02);
    //TH1D *hEESigmaIEtaIEta = new TH1D("hEESigmaIEtaIEta", "hadron_SigmaIEtaIEta_EE_Sample_Total",100, 0.01, 0.05);
    TH1D *hEBEt = new TH1D("hadronEBEt", "hadron_SigmaIEt_EB_Sample_Full_Energy_Spectral",70, 0, 140);
    //TH1D *hEEEt = new TH1D("hadronEEEt", "hadron_SigmaIEt_EE_Sample_Full_Energy_Spectral",70, 0, 140);

    TH1D *taretEBSigmaIEtaIEta = new TH1D("targetEBSigmaIEtaIEta", "target_SigmaIEtaIEta_EB_Sample_Total",100, 0, 0.02);
    //TH1D *targetEESigmaIEtaIEta = new TH1D("targetEESigmaIEtaIEta", "target_SigmaIEtaIEta_EE_Sample_Total",100, 0.01, 0.05);
    TH1D *targetEBEt = new TH1D("targetEBEt", "target_SigmaIEt_EB_Sample_Full_Energy_Spectral",70, 0, 140);
    //TH1D *targetEEEt = new TH1D("targetEEEt", "target_SigmaIEt_EE_Sample_Full_Energy_Spectral",70, 0, 140);

    //Pt Binning of the photonic Sample's SigmaIEtaiEta histogram
    int binrange[] = {20,22,24,26,28,30,32,34,36,38,42,46,50,60,70,90,140};
    TH1D *hTemplateHist[100];
    TH1D *targetHist[100];
    TString hbintitle, targetbintitle;
    for(unsigned ibin(0); ibin<16;ibin++)
      {
    hbintitle. Form("HadronicTemplate_Pt_bin(%d,%d)",binrange[ibin],binrange[ibin+1]);
    targetbintitle. Form("FitTarget_Pt_bin(%d,%d)",binrange[ibin],binrange[ibin+1]);
    hTemplateHist[ibin] = new TH1D(hbintitle, "pt_bin template for hadronic fraction",50,0,0.02);
    targetHist[ibin] = new TH1D(targetbintitle, "pt_bin for target",50,0,0.02);
      }



    Double_t InvariantMassMuMu, InvariantMassMuMuGamma;
    TLorentzVector Mu4Momentum0, Mu4Momentum1, Gamma4Momentum, MuMu4Momentum, MuMuGamma4Momentum;

    Int_t isGamma(0);// how many qualified gamma is found.
    unsigned gammaID(0);//index of the qualified photon.
    unsigned muonID(0);
    Double_t muonPt(0);

    Double_t dEta0, dEta1, dPhi0, dPhi1, deltaR0, deltaR1;


    const unsigned nEvts = es->GetEntries(); 
    std::cout << " nEvts=" << nEvts << std::endl;
 //--------------------------------------------------------------------------------------
 //loop over all the events
     for (unsigned ievt(0); ievt<nEvts; ++ievt)   
      {

    if (ievt%10000==0) std::cout << " -- Processing event " << ievt << std::endl;
        //std::cout << " -- Processing event " << ievt << std::endl;
    


    es->GetEntry(ievt);
    
    //cut on missing transverse energy
    if(pfMET >70) continue;

    //number of particles
    if( nMu == 0) continue;//should have at least one reco-muons.
    if( nPho == 0 ) continue;//should have at least one reco-photon.


      
    //two muons pass TightIDmuon selection
    muonPt=0;
    
        for(unsigned imu(0);imu < nMu; imu++)
      {
        if (!(*muIsTightID)[imu]) continue;
        if(fabs((*muEta)[imu])>2.4) continue;
        if(muonPt<(*muPt)[imu]) 
          {
        muonPt = (*muPt)[imu];
        muonID = imu;
          }
      }
    if(muonPt == 0) continue;
         
    

    //loop for gamma 

    for(unsigned ipho(0); ipho<nPho; ipho++)
      {
        isBarrel = fabs((*phoEta)[ipho]) < 1.444;
        isEndCap = (fabs((*phoEta)[ipho]) > 1.560 && fabs((*phoEta)[ipho]) < 2.5);
        if(!isBarrel && !isEndCap) continue;
        if((*phoEt)[ipho]<25) continue;
        if((*phoEleVeto)[ipho] != 0 ) continue;
        decision = LooseCut((*phoEta)[ipho], (*phoHoverE)[ipho], (*phoSigmaIEtaIEtaFull5x5)[ipho], (*phoPFChIso)[ipho], (*phoPFNeuIso)[ipho], (*phoPFPhoIso)[ipho], (*phoEt)[ipho], rho);

         //fill the histograms for barrels now.
        if(!isBarrel)continue;

        //FSR cut
        dEta0 = fabs((*phoEta)[ipho]-(*muEta)[muonID]);
        dPhi0 = fabs((*phoPhi)[ipho]-(*muPhi)[muonID]);
        deltaR0 = TMath::Sqrt(dEta0*dEta0+dPhi0*dPhi0);
        if(deltaR0<0.8) continue;
        flagFSR = 0;
        for(unsigned imu(0); imu < nMu; imu++)
          {
        dEta1 = fabs((*phoEta)[ipho]-(*muEta)[imu]);
        dPhi1 = fabs((*phoPhi)[ipho]-(*muPhi)[imu]);
        deltaR1 = TMath::Sqrt(dEta0*dEta0+dPhi0*dPhi0);
        if((deltaR1<0.3) && ((*muPt)[imu]>2))
          {
            flagFSR = 1;
          }
          }
        if(flagFSR == 1) continue;
        
        //target
        if((((decision >> 2) &1) ==1) && (((decision >> 3) &1) ==1) && (((decision >> 4) &1) ==1)  && (((decision >> 5) &1) ==1) )
          {
        targetEBSigmaIEtaIEta -> Fill((*phoSigmaIEtaIEtaFull5x5)[ipho]);
        targetEBEt -> Fill ((*phoEt)[ipho]);
        for(unsigned ibin(0); ibin<16;ibin++)
           {
             if(((*phoEt)[ipho]>binrange[ibin]) && ((*phoEt)[ipho]<binrange[ibin+1]) ) targetHist[ibin]->Fill((*phoSigmaIEtaIEtaFull5x5)[ipho]);
           }
          }
          //hadronic template
        if( (((decision >> 2) &1) ==1) && (((decision >> 3) &1) ==1) && (((decision >> 5) &1) ==1) )
          {
        
        if((*phoPFChIso)[ipho] >15)//select hardronic template by put upper bound of I_ch.
          {
            hEBSigmaIEtaIEta -> Fill((*phoSigmaIEtaIEtaFull5x5)[ipho]);
            hEBEt -> Fill ((*phoEt)[ipho]);
            for(unsigned ibin(0); ibin<16;ibin++)
              {
            if(((*phoEt)[ipho]>binrange[ibin]) && ((*phoEt)[ipho]<binrange[ibin+1]) ) hTemplateHist[ibin]->Fill((*phoSigmaIEtaIEtaFull5x5)[ipho]);
              }
          }

          }

      }

    




        
 



      }//end of event loop
     
     
     hSample->Write();
     std::cout << "All histograms have been saved in hadronicTemplate.root."<<std::endl;
    
     //std::cout<<counter<<std::endl;
       



     /*
    //Draw the histogram
    TCanvas *c1 =new TCanvas("c1","",1600,900);
    c1->Divide(2);
    
    c1->cd(1);
    GammaEBSigmaIEtaIEta -> Draw();

    c1->cd(2);
    GammaEESigmaIEtaIEta -> Draw();

    c1->SaveAs("test.png");
     */
    
    


}//end of main().