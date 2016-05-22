#include "TMath.h"

typedef struct Photon
{
  public:
    float Et;
    float eta;
    float phi;
    int EleVeto;
    float HoverE;
    float Sigmaietaieta;
    float ChIso;
    float NeuIso;
    float PhoIso;
    UShort_t phoIDbit;

  Photon(float Et_, float eta_, float phi_, int EleVeto_, float HoverE_, float SigmaIEtaIEta_, float ChIso_, float NeuIso_, float PhoIso_, UShort_t phoIDbit_)
  {
    this->Et = Et_;
    this->eta = eta_;
    this->phi = phi_;
    this->EleVeto = EleVeto_;
    this->HoverE = HoverE_;
    this->Sigmaietaieta = SigmaIEtaIEta_;
    this->ChIso = ChIso_;
    this->NeuIso = NeuIso_;
    this->PhoIso = PhoIso_;
    this->phoIDbit = phoIDbit_;
   }
};
   
float EAch[] = {0.0157, 0.0143, 0.0115, 0.0094, 0.0095, 0.0068, 0.0053};
float EAneu[]= {0.0143, 0.0210, 0.0148, 0.0082, 0.0124, 0.0186, 0.0320};
float EApho[]= {0.0725, 0.0604, 0.0320, 0.0512, 0.0766, 0.0949, 0.1160};
float eRange[]={0, 1.0,  1.479, 2.0, 2.2, 2.3, 2.4, 2.5};

int RhoCorrection(float eta, float rho)
{
     int EAindex(-1);
     for(unsigned i(0); i< 7; i++)
        if(eRange[i] < fabs(eta) && eRange[i+1] > fabs(eta))EAindex = i;
     return EAindex;
} 


UShort_t LooseCut(float eta, float HoverE, float sigma, float chIso, float neuIso, float phoIso, float pt, float rho)
{
         bool isBarrel = fabs(eta) < 1.444;
         bool isEndCap = (fabs(eta) > 1.560 && fabs(eta) < 2.5);
         bool passHE(true),passSIGMA(true),passCH(true),passNEU(true),passPHO(true);
         UShort_t  decision(0);
         int  iEA = RhoCorrection(eta, rho);
         if(iEA >= 0)
         {
            chIso = chIso-rho*EAch[iEA];
            neuIso = neuIso-rho*EAneu[iEA];
            phoIso = phoIso-rho*EApho[iEA];
         }
         if(isBarrel)
         {
           if(HoverE > 0.05) passHE = false;
           if(sigma > 0.0103) passSIGMA = false;
           if(chIso > 2.44) passCH = false;
           if(neuIso > 2.57 + TMath::Exp(pt*0.0044 + 0.5809))passNEU = false;
           if(phoIso > 1.92 + 0.0043*pt)passPHO = false;
         }
         else if(isEndCap)
         {
           if(HoverE > 0.05) passHE = false;
           if(sigma > 0.0277) passSIGMA = false;
           if(chIso > 1.84) passCH = false;
           if(neuIso > 4.00 + TMath::Exp(pt*0.0044 + 0.9402))passNEU = false;
           if(phoIso > 2.15 + 0.0041*pt)passPHO = false;
         }
         if(passSIGMA)decision |= 1 << 1;
         if(passHE)decision |= 1 << 2;
         if(passNEU)decision |= 1 << 3;
         if(passCH)decision |= 1 << 4;
         if(passPHO)decision |= 1 << 5;
         return decision;
}


// TO USE LOOSECUT IN MAIN
{

      std::vector<Photon> phovec;

      phovec.clear();
      for(unsigned iPho(0); iPho < (*ggphoET).size(); iPho++){
         phovec.push_back(Photon((*ggphoET)[iPho], (*ggphoEta)[iPho], (*ggphoPhi)[iPho],(*ggphoEleVeto)[iPho], (*ggphoHoverE)[iPho] , (*ggphoSigmaIEtaIEta)[iPho], (*ggphoPFChIso)[iPho], (*ggphoPFNeuIso)[iPho], (*ggphoPFPhoIso)[iPho], (*phoIDbit)[iPho]));
      }


   UShort_t decision = LooseCut(phovec[iPho].eta, phovec[iPho].HoverE, phovec[iPho].Sigmaietaieta, phovec[iPho].ChIso, phovec[iPho].NeuIso, phovec[iPho].PhoIso, phovec[iPho].Et, rho);

   if(((decision >> 1) &1) ==1)std::cout << "pass sigmaietaieta";
   if(((decision >> 2) &1) ==1)std::cout << "pass H/E";
   if(((decision >> 3) &1) ==1)std::cout << "pass rho-corrected Neutral hadron isolation";
   if(((decision >> 4) &1) ==1)std::cout << "pass rho-corrected Charged hadron isolation";
   if(((decision >> 5) &1) ==1)std::cout << "pass rho-corrected photon isolation";
}

