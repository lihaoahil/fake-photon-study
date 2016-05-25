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
#include "include/printProgBar.h"

void ProgBarTest(){//main

    int nEvts = 34556728;
    int percent;

    for(unsigned ievt(0);ievt<nEvts;ievt++)
    {
        if (ievt%10000==0)
        {
          percent = 100*ievt/nEvts;
          printProgBar(percent, ievt);
        }
    }  
  

    


}//end of main().
