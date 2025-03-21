
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "initOniaTree.C"
#include <TH3.h>

void readList(string filename, std::vector<string> &dataMap) {
  std::ifstream inputFile(filename);
  if (!inputFile) {
    std::cerr << "Error opening file." << std::endl;
    return;
  }
  std::string line;
  while (std::getline(inputFile, line)) {
    dataMap.push_back(line);
  }
  inputFile.close();
}

// void tree2DataSet()
RooDataSet *tree2DataSet(string inputName, string TreeName, RooRealVar *ctau,
                         RooRealVar *mass) {
  vector<string> InputFileNames;
  // readList("list.txt",InputFileNames);
  readList(inputName, InputFileNames);
  // string TreeName("O2rtdimuonall");

  TChain *theTree = new TChain(TreeName.c_str(), "");
  cout << "[INFO] Creating TChain" << endl;
  getTChain(theTree, InputFileNames, TreeName); // Import files to TChain
  initOniaTree(theTree, TreeName);              // Initialize the Onia Tree
  iniBranch(theTree, TreeName);                 // Initialize the Branches

  // RooRealVar* mass    = new RooRealVar("invMass","#mu#mu mass", 1.0, 6.0,
  // "GeV/c^{2}"); RooRealVar* ctau    = new RooRealVar("ctau","c_{#tau}",
  // -100000.0, 100000.0, "mm"); RooRealVar* ctauErr = new
  // RooRealVar("ctauErr","#sigma_{c#tau}", -100000.0, 100000.0, "mm");
  // RooRealVar* ptQQ    = new RooRealVar("pt","#mu#mu p_{T}", -1.0, 10000.0,
  // "GeV/c"); RooRealVar* rapQQ   = new RooRealVar("rap","#mu#mu y", -5., 5.,
  // ""); RooRealVar* chi21   = new RooRealVar("chi21","#chi^2_{MFT-MCH}",
  // -1000,1000000, ""); RooRealVar* chi22   = new
  // RooRealVar("chi22","#chi^2_{MFT-MCH}", -1000,1000000, "");

  // RooArgSet* cols = new RooArgSet(*mass, *ctau, *ctauErr, *ptQQ, *rapQQ,
  // *chi21, *chi22);
  RooArgSet *cols = new RooArgSet(*ctau, *mass); // PENSER A AJOUTER MASSE

  RooDataSet *data = new RooDataSet("dOS", "dOS", *cols);

  Long64_t nentries = theTree->GetEntries();

  cout << "[INFO] Starting to process " << nentries << " nentries" << endl;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    theTree->GetEntry(jentry);

    /*ROOT::Math::PtEtaPhiMVector v1MC(fPtMC1, fEtaMC1, fPhiMC1, 0.105658);
    ROOT::Math::PtEtaPhiMVector v2MC(fPtMC2, fEtaMC2, fPhiMC2, 0.105658);
    ROOT::Math::PtEtaPhiMVector v12MC = v1MC + v2MC;
    float Tauz1MC = (fMCPosZ - fVz1) * v12MC.M() / TMath::Abs(v12MC.Pz());
    float Tauz2MC = (fMCPosZ - fVz2) * v12MC.M() / TMath::Abs(v12MC.Pz());
    float fTauzMC = (Tauz1MC + Tauz2MC) / 60.;*/ //facteur /60. pour MC car quelqu'un a décidé de mettre c en cm/ns

    // cout << "Vertex MC : " << fVz1 << endl;

    mass->setVal(fMass);
    // float ct= fTauz*299792458.e-7*10.;
    // float ctErr= fTauzErr*299792458.e-7*10.;

    /*float tauz_mPDG = fTauz * 3.096916 / fMass;
    float tauz_mPDG_Err = fTauzErr * 3.096916 / fMass;

    float ct = tauz_mPDG * 300;
    float ctErr = tauz_mPDG_Err * 300;*/

    /*****************************/
    float ct = fTauz * 300.;
    float ctErr = fTauzErr * 300.;
    /*****************************/

    //float ct = fTauz * 300. * 3.097 / fMass;

    //  float ct = fTauz;
    //  float ctErr = fTauzErr;
    // float fTauzbis = fTauz/fMass*3.097;
    // ctau->setVal(fTauz); //rajouter putain de facteur 30
    ctau->setVal(ct);
    //ctauErr->setVal(ctErr);
    // ROOT::Math::PtEtaPhiMVector v12(fPt, fEta, fPhi, fMass);
    // ptQQ->setVal(fPt);
    // chi21->setVal(fChi2MatchMCHMFT1);
    // chi22->setVal(fChi2MatchMCHMFT2);
    // rapQQ->setVal(std::abs(v12.Rapidity()));
    // if (fPt1 > 0.7 && fPt2 > 0.7 && fEta1 < -2.5 && fEta1 > - 3.6 && fEta2 <
    // -2.5 && fEta2 > - 3.6 && fChi2pca < 8 && fNumContrib > 5 && fSVertex
    // >-10. && fSVertex < 10.){

    float rap = -ROOT::Math::log(
        (ROOT::Math::sqrt(fMass * fMass + fPt * fPt * ROOT::Math::cosh(fEta) *
                                              ROOT::Math::cosh(fEta)) +
         fPt * ROOT::Math::sinh(fEta)) /
        (ROOT::Math::sqrt(fMass * fMass + fPt * fPt)));
    if (fEta1 < -2.5 && fEta1 > -3.6 && fEta2 < -2.5 && fEta2 > -3.6) {
      if (fChi2MatchMCHMFT1 <= 30. && fChi2MatchMCHMFT2 <= 30.) {
        if (fPt1 > 0.7 && fPt2 > 0.7) {
          // if (fChi2MatchMCHMFT1 <= 35. && fChi2MatchMCHMFT2 <= 35. &&
          // fChi2MatchMCHMFT1 >= 30. && fChi2MatchMCHMFT2 >= 30.) {
          //  if (fMcMask1 < 1. && fMcMask1 < 1.) {
          if (fIsAmbig1 == 0 && fIsAmbig2 == 0) {
            if (fSign == 0) {
              if (rap > 2.5 && rap < 3.6) {
                if (fPt > 0. && fPt < 1.) {
                  if (fMass > 2.5 && fMass < 3.5) {
                    // if (fTauz > -3.e-2 && fTauz < 3.e-2) {
                    // if (fTauz > -1e-4 && fTauz < 6e-3) {
                    if (ct > -10. && ct < 10.) {
                      // if (ct > -1.e-1 && ct < 1.5) {
                      // if(ct > -1.e-1 && ct < 1.5){
                      //  if (fPdgCode1 != 100443 && fPdgCode2 != 100443) { //
                      //  ON GARDE LE && pour une fois que
                      //   ça fonctionne
                      // if (fPdgCode1 == 443 && fPdgCode2 == 443) {
                      // if (fMcMask1 < 1. && fMcMask1 < 1.) {
                      data->add(*cols, 1.0); // Signal and background dimuons
                      // std::cout << "mass = " << fMass << std::endl;
                      //}
                      //}
                    } // ctau
                  } // mass
                } //pT
              } // rap
              //}
            } // sign
          } // ambig
        } // single pT
      } // chi2MCHMFT
    } // eta
  } // for loop
  theTree->Reset();
  delete theTree;
  return data;
};
