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
#include "include/photonCut.h"

void gammaTemplateDY(){//main
  
  TChain *es =new TChain("ggNtuplizer/EventTree");
  es->Add("/export/cmss2/mengleis/SM/DYJetsToLL/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-part1.root");
  es->Add("/export/cmss2/mengleis/SM/DYJetsToLL/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-part2.root");

 

  TFile *gSample = new TFile("gammaTemplateMC.root","RECREATE");



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


    ULong64_t  HLTEleMuX = 0; //HLT
    UShort_t decision = 0;
    float rho = 0;
    bool isBarrel, isEndCap;

 

//Set branches to photon variables


    es->SetBranchAddress("HLTEleMuX",&HLTEleMuX);

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
    TH1D *GammaEBSigmaIEtaIEta = new TH1D("GammaEBSigmaIEtaIEta", "Gamma_SigmaIEtaIEta_EB_Sample_Total",100, 0, 0.02);
    //TH1D *GammaEESigmaIEtaIEta = new TH1D("GammaEESigmaIEtaIEta", "Gamma_SigmaIEtaIEta_EE_Sample_Total",100, 0.01, 0.05);
    TH1D *GammaEBEt = new TH1D("GammaEBEt", "Gamma_SigmaIEt_EB_Sample_Full_Energy_Spectral",70, 0, 140);
    //TH1D *GammaEEEt = new TH1D("GammaEEEt", "Gamma_SigmaIEt_EE_Sample_Full_Energy_Spectral",70, 0, 140);

    //Pt Binning of the photonic Sample's SigmaIEtaiEta histogram
    int binrange[] = {20,22,24,26,28,30,32,34,36,38,42,46,50,60,70,90,140};
    TH1D *phoTemplateHist[100];
    TString bintitle;
    for(unsigned ibin(0); ibin<16;ibin++)
      {
    bintitle. Form("PhotonTemplate_Pt_bin(%d,%d)",binrange[ibin],binrange[ibin+1]);
    phoTemplateHist[ibin] = new TH1D(bintitle, "pt_bin template for photon",50,0,0.02);
      }


    Double_t InvariantMassMuMu, InvariantMassMuMuGamma;
    TLorentzVector Mu4Momentum0, Mu4Momentum1, Gamma4Momentum, MuMu4Momentum, MuMuGamma4Momentum;

    Int_t isGamma(0);// how many qualified gamma is found.
    unsigned gammaID(0);//index of the qualified photon.

    Double_t dEta0, dEta1, dPhi0, dPhi1, deltaR0, deltaR1;


    const unsigned nEvts = es->GetEntries(); 
    std::cout << " nEvts=" << nEvts << std::endl;

     for (unsigned ievt(0); ievt<nEvts; ++ievt)    //loop over all the events
      {

    if (ievt%10000==0) std::cout << " -- Processing event " << ievt << std::endl;
        //std::cout << " -- Processing event " << ievt << std::endl;
    


    es->GetEntry(ievt);
    
        if( ( ((HLTEleMuX >> 20) & 1) == 0 ) && (( (HLTEleMuX >> 21) & 1) == 0 ) ) continue; //not necessary, just in case.

    if( nMu != 2) continue;//should have exactly two reco-muons.(events with more than one qualified di-muon pair are very rare and for such case the energy is spread out.)
    if( nPho == 0 ) continue;//should have at least one reco-photon.
    
    //two muons pass TightID
    if(!(*muIsTightID)[0] || !(*muIsTightID)[1] ) continue;
         
    //each reco-muon should have Pt>10.5 and |eta|<2.4.
        if( (*muPt)[0] < 10.5  || fabs((*muEta)[0])>2.4) continue;
    if( (*muPt)[1] < 10.5  || fabs((*muEta)[1])>2.4) continue;
    

    ////need two reco-muons with opposite charges.
    if( (*muCharge)[0] + (*muCharge)[1] != 0 ) continue;
    //counter++;

    //loop for gamma 
    isGamma=0;
    gammaID=0;
    for(unsigned ipho(0); ipho<nPho; ipho++)
      {
        isBarrel = fabs((*phoEta)[ipho]) < 1.444;
        isEndCap = (fabs((*phoEta)[ipho]) > 1.560 && fabs((*phoEta)[ipho]) < 2.5);
        if(!isBarrel && !isEndCap) continue;
        if((*phoEt)[ipho]<20) continue;
        
        dEta0 = fabs((*phoEta)[ipho]-(*muEta)[0]);
        dEta1 = fabs((*phoEta)[ipho]-(*muEta)[1]);
        dPhi0 = fabs((*phoPhi)[ipho]-(*muPhi)[0]);
        dPhi1 = fabs((*phoPhi)[ipho]-(*muPhi)[1]);
        deltaR0 = TMath::Sqrt(dEta0*dEta0+dPhi0*dPhi0);
        deltaR1 = TMath::Sqrt(dEta1*dEta1+dPhi1*dPhi1);

        if(deltaR0>0.8 && deltaR1>0.8) continue;
        if( deltaR0<deltaR1  && (*muPt)[1]>21 ) 
          { 
          isGamma++;
          gammaID=ipho;
          }
        if( deltaR0>deltaR1  && (*muPt)[0]>21 )
          {
        isGamma++;
        gammaID=ipho;
          }

      }

    



    if(isGamma == 1)
      {
        Mu4Momentum0.SetPtEtaPhiE((*muPt)[0],(*muEta)[0],(*muPhi)[0],(*muEn)[0]);
        Mu4Momentum1.SetPtEtaPhiE((*muPt)[1],(*muEta)[1],(*muPhi)[1],(*muEn)[1]);
        Gamma4Momentum.SetPtEtaPhiE((*phoEt)[gammaID],(*phoEta)[gammaID],(*phoPhi)[gammaID],(*phoE)[gammaID]);

        MuMu4Momentum = Mu4Momentum0 + Mu4Momentum1;
        MuMuGamma4Momentum = MuMu4Momentum + Gamma4Momentum;

        InvariantMassMuMu = MuMu4Momentum.M();
        InvariantMassMuMuGamma = MuMuGamma4Momentum.M();

        if(InvariantMassMuMu<35 || InvariantMassMuMuGamma<60 || InvariantMassMuMuGamma>120 ) continue;
        if((InvariantMassMuMu+InvariantMassMuMuGamma)>180)continue;
        
        decision = LooseCut((*phoEta)[gammaID], (*phoHoverE)[gammaID], (*phoSigmaIEtaIEtaFull5x5)[gammaID], (*phoPFChIso)[gammaID], (*phoPFNeuIso)[gammaID], (*phoPFPhoIso)[gammaID], (*phoEt)[gammaID], rho);

        
        //if(((decision >> 1) &1) !=1) continue;
        if(((decision >> 2) &1) !=1) continue;
        if(((decision >> 3) &1) !=1) continue;
        if(((decision >> 4) &1) !=1) continue;
        if(((decision >> 5) &1) !=1) continue;
        //isBarrel = fabs((*phoEta)[gammaID]) < 1.444;
        //isEndCap = (fabs((*phoEta)[gammaID]) > 1.560 && fabs((*phoEta)[gammaID]) < 2.5);MuMuHist->Fill(InvariantMassMuMu);
        //MuMuHist->Fill(InvariantMassMuMu);
        //MuMuGammaHist->Fill(InvariantMassMuMuGamma);
        //MuMuplusMuMuGammaHist->Fill(InvariantMassMuMu+InvariantMassMuMuGamma);

          
        if(isBarrel)
          {
        GammaEBSigmaIEtaIEta -> Fill((*phoSigmaIEtaIEtaFull5x5)[gammaID]);
        GammaEBEt -> Fill ((*phoEt)[gammaID]);

         for(unsigned ibin(0); ibin<16;ibin++)
           {
             if(((*phoEt)[gammaID]>binrange[ibin]) && ((*phoEt)[gammaID]<binrange[ibin+1]) ) phoTemplateHist[ibin]->Fill((*phoSigmaIEtaIEtaFull5x5)[gammaID]);
           }
          }
        /*if(isEndCap) 
          {
        GammaEESigmaIEtaIEta -> Fill((*phoSigmaIEtaIEtaFull5x5)[gammaID]);
        GammaEEEt -> Fill ((*phoEt)[gammaID]);
        }*/
       
          
      }



      }//end of event loop
     
     
     gSample->Write();
     std::cout << "All histograms have been saved in gammaTemplateDY.root."<<std::endl;
    
  


}//end of main().