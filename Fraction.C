#include "RooPlot.h"
#include "RooRealVar.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "THnSparse.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TString.h"
#include "TSystem.h"
#include "initOniaTree.C"
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>

void Fraction() {
  double bins[] = {0., 1., 2., 3., 4., 5., 6., 8., 10., 15., 30.};
  double binsbis[] = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14.};
  TH1F *Fraction = new TH1F("fB", "", 10, bins);
  TH1F *Fraction_New = new TH1F("fB_New", "", 9, bins);

  /*TFile *fB_pT01 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_Full_01_Classic.root");
  TH1F *fB_01 = static_cast<TH1F *>(fB_pT01->Get("fB_PDG_01"));
  Fraction->SetBinContent(1, fB_01->GetBinContent(1));
  Fraction->SetBinError(1, fB_01->GetBinContent(2));
  TFile *fB_pT12 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_Full_12_Classic.root");
  TH1F *fB_12 = static_cast<TH1F *>(fB_pT12->Get("fB_PDG_12"));
  Fraction->SetBinContent(2, fB_12->GetBinContent(1));
  Fraction->SetBinError(2, fB_12->GetBinContent(2));
  TFile *fB_pT23 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_Full_23_Classic.root");
  TH1F *fB_23 = static_cast<TH1F *>(fB_pT23->Get("fB_PDG_23"));
  Fraction->SetBinContent(3, fB_23->GetBinContent(1));
  Fraction->SetBinError(3, fB_23->GetBinContent(2));
  TFile *fB_pT34 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_Full_34_Classic.root");
  TH1F *fB_34 = static_cast<TH1F *>(fB_pT34->Get("fB_PDG_34"));
  Fraction->SetBinContent(4, fB_34->GetBinContent(1));
  Fraction->SetBinError(4, fB_34->GetBinContent(2));
  TFile *fB_pT45 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_Full_45_Classic.root");
  TH1F *fB_45 = static_cast<TH1F *>(fB_pT45->Get("fB_PDG_45"));
  Fraction->SetBinContent(5, fB_45->GetBinContent(1));
  Fraction->SetBinError(5, fB_45->GetBinContent(2));
  TFile *fB_pT56 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_Full_56_Classic.root");
  TH1F *fB_56 = static_cast<TH1F *>(fB_pT56->Get("fB_PDG_56"));
  Fraction->SetBinContent(6, fB_56->GetBinContent(1));
  Fraction->SetBinError(6, fB_56->GetBinContent(2));
  TFile *fB_pT68 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_Full_68_Classic.root");
  TH1F *fB_68 = static_cast<TH1F *>(fB_pT68->Get("fB_PDG_68"));
  Fraction->SetBinContent(7, fB_68->GetBinContent(1));
  Fraction->SetBinError(7, fB_68->GetBinContent(2));
  TFile *fB_pT810 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_Full_810_Classic.root");
  TH1F *fB_810 = static_cast<TH1F *>(fB_pT810->Get("fB_PDG_810"));
  Fraction->SetBinContent(8, fB_810->GetBinContent(1));
  Fraction->SetBinError(8, fB_810->GetBinContent(2));
  TFile *fB_pT1030 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_Full_1030_Classic.root");
  TH1F *fB_1030 = static_cast<TH1F *>(fB_pT1030->Get("fB_PDG_1030"));
  Fraction->SetBinContent(9, fB_1030->GetBinContent(1));
  Fraction->SetBinError(9, fB_1030->GetBinContent(2));

  Fraction->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  Fraction->GetYaxis()->SetTitle("f_{B}");
  Fraction->SetMinimum(0.);
  Fraction->SetMaximum(0.5);
  Fraction->SetMarkerStyle(34);
  Fraction->SetMarkerSize(2);
  Fraction->SetMarkerColor(4);*/

  /*TFile *fB_PDG_pT01 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_01.root");
  TH1F *fB_PDG_01 = static_cast<TH1F *>(fB_PDG_pT01->Get("fB_PDG_01"));
  Fraction_PDG->SetBinContent(1, fB_PDG_01->GetBinContent(1));
  Fraction_PDG->SetBinError(1, fB_PDG_01->GetBinContent(2));
  TFile *fB_PDG_pT12 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_12.root");
  TH1F *fB_PDG_12 = static_cast<TH1F *>(fB_PDG_pT12->Get("fB_PDG_12"));
  Fraction_PDG->SetBinContent(2, fB_PDG_12->GetBinContent(1));
  Fraction_PDG->SetBinError(2, fB_PDG_12->GetBinContent(2));
  TFile *fB_PDG_pT23 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_23.root");
  TH1F *fB_PDG_23 = static_cast<TH1F *>(fB_PDG_pT23->Get("fB_PDG_23"));
  Fraction_PDG->SetBinContent(3, fB_PDG_23->GetBinContent(1));
  Fraction_PDG->SetBinError(3, fB_PDG_23->GetBinContent(2));
  TFile *fB_PDG_pT34 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_34.root");
  TH1F *fB_PDG_34 = static_cast<TH1F *>(fB_PDG_pT34->Get("fB_PDG_34"));
  Fraction_PDG->SetBinContent(4, fB_PDG_34->GetBinContent(1));
  Fraction_PDG->SetBinError(4, fB_PDG_34->GetBinContent(2));
  TFile *fB_PDG_pT45 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_45.root");
  TH1F *fB_PDG_45 = static_cast<TH1F *>(fB_PDG_pT45->Get("fB_PDG_45"));
  Fraction_PDG->SetBinContent(5, fB_PDG_45->GetBinContent(1));
  Fraction_PDG->SetBinError(5, fB_PDG_45->GetBinContent(2));
  TFile *fB_PDG_pT530 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_530.root");
  TH1F *fB_PDG_530 = static_cast<TH1F *>(fB_PDG_pT530->Get("fB_PDG_530"));
  Fraction_PDG->SetBinContent(6, fB_PDG_530->GetBinContent(1));
  Fraction_PDG->SetBinError(6, fB_PDG_530->GetBinContent(2));

  Fraction_PDG->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  Fraction_PDG->GetYaxis()->SetTitle("f_{B}");
  Fraction_PDG->SetMinimum(0.);
  Fraction_PDG->SetMaximum(0.5);
  Fraction_PDG->SetMarkerStyle(34);
  Fraction_PDG->SetMarkerSize(2);
  Fraction_PDG->SetMarkerColor(2);*/

  TH1F *LHCb = new TH1F("fB_LHCb", "", 13, binsbis); //3 3.5
  //01
  LHCb->SetBinContent(1, 0.092); //t'as mis les valeurs en % patate
  LHCb->SetBinError(1, 0.003);
  //12
  LHCb->SetBinContent(2, 0.109); //t'as mis les valeurs en % patate
  LHCb->SetBinError(2, 0.002);
  //23
  LHCb->SetBinContent(3, 0.125); //t'as mis les valeurs en % patate
  LHCb->SetBinError(3, 0.002);
  //34
  LHCb->SetBinContent(4, 0.139);
  LHCb->SetBinError(4, 0.002);
  //45
  LHCb->SetBinContent(5, 0.153);
  LHCb->SetBinError(5, 0.003);
  //56
  LHCb->SetBinContent(6, 0.169);
  LHCb->SetBinError(6, 0.003);
  //67
  LHCb->SetBinContent(7, 0.197);
  LHCb->SetBinError(7, 0.005);
  //78
  LHCb->SetBinContent(8, 0.213);
  LHCb->SetBinError(8, 0.006);
  //89
  LHCb->SetBinContent(9, 0.237);
  LHCb->SetBinError(9, 0.008);
  //910
  LHCb->SetBinContent(10, 0.272);
  LHCb->SetBinError(10, 0.011);
  //1011
  LHCb->SetBinContent(11, 0.309);
  LHCb->SetBinError(11, 0.014);
  //1112
  LHCb->SetBinContent(12, 0.281);
  LHCb->SetBinError(12, 0.018);
  //1213
  LHCb->SetBinContent(13, 0.333);
  LHCb->SetBinError(13, 0.022);
  //1314
  LHCb->SetBinContent(14, 0.334);
  LHCb->SetBinError(14, 0.028);

  LHCb->SetMinimum(0.);
  LHCb->SetMaximum(0.5);
  LHCb->SetMarkerColor(kGreen+1);
  LHCb->SetMarkerStyle(22);
  LHCb->SetMarkerSize(2);

  TH1F *LHCbis = new TH1F("fB_LHCbis", "", 13, binsbis); //2.5 3
  //01
  LHCbis->SetBinContent(1, 0.092); //t'as mis les valeurs en % patate
  LHCbis->SetBinError(1, 0.003);
  //12
  LHCbis->SetBinContent(2, 0.112); //t'as mis les valeurs en % patate
  LHCbis->SetBinError(2, 0.002);
  //23
  LHCbis->SetBinContent(3, 0.126); //t'as mis les valeurs en % patate
  LHCbis->SetBinError(3, 0.002);
  //34
  LHCbis->SetBinContent(4, 0.142);
  LHCbis->SetBinError(4, 0.002);
  //45
  LHCbis->SetBinContent(5, 0.17);
  LHCbis->SetBinError(5, 0.003);
  //56
  LHCbis->SetBinContent(6, 0.182);
  LHCbis->SetBinError(6, 0.003);
  //67
  LHCbis->SetBinContent(7, 0.211);
  LHCbis->SetBinError(7, 0.004);
  //78
  LHCbis->SetBinContent(8, 0.23);
  LHCbis->SetBinError(8, 0.006);
  //89
  LHCbis->SetBinContent(9, 0.256);
  LHCbis->SetBinError(9, 0.008);
  //910
  LHCbis->SetBinContent(10, 0.263);
  LHCbis->SetBinError(10, 0.01);
  //1011
  LHCbis->SetBinContent(11, 0.315);
  LHCbis->SetBinError(11, 0.013);
  //1112
  LHCbis->SetBinContent(12, 0.333);
  LHCbis->SetBinError(12, 0.016);
  //1213
  LHCbis->SetBinContent(13, 0.365);
  LHCbis->SetBinError(13, 0.021);
  //1314
  LHCbis->SetBinContent(14, 0.373);
  LHCbis->SetBinError(14, 0.023);

  LHCbis->SetMinimum(0.);
  LHCbis->SetMaximum(0.5);
  LHCbis->SetMarkerColor(kRed+1);
  LHCbis->SetMarkerStyle(22);
  LHCbis->SetMarkerSize(2);

  TFile *AccxEff = TFile::Open("/Users/emiliebarreau/alice/Macros/RatioOutput_newbinning.root");
  TH1F *AccEff = static_cast<TH1F *>(AccxEff->Get("h_ratio")); 

  std::cout << " TEST" << std::endl;
  std::cout << AccEff->GetBinContent(1) << std::endl;
  std::cout << AccEff->GetBinContent(2) << std::endl;
  std::cout << AccEff->GetBinContent(3) << std::endl;
  std::cout << AccEff->GetBinContent(4) << std::endl;
  std::cout << AccEff->GetBinContent(5) << std::endl;
  std::cout << AccEff->GetBinContent(6) << std::endl;
  std::cout << AccEff->GetBinContent(7) << std::endl;
  std::cout << AccEff->GetBinContent(8) << std::endl;
  std::cout << AccEff->GetBinContent(9) << std::endl;
  std::cout << AccEff->GetBinContent(10) << std::endl;

  TFile *fB_pT01 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_Full_01_tzJPsi_free_MEH.root"); 
  TH1F *fB_01 = static_cast<TH1F *>(fB_pT01->Get("fB_PDG_01"));
  double fB01 = 1/(1 + (1 - fB_01->GetBinContent(1))/fB_01->GetBinContent(1) * AccEff->GetBinContent(1));
  double fB01_E = 1/(1 + (1 - fB_01->GetBinContent(2))/fB_01->GetBinContent(2) * AccEff->GetBinContent(1));
  Fraction->SetBinContent(1, fB01);
  Fraction->SetBinError(1, fB01_E);

  TFile *fB_pT12 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_Full_12_tzJPsi_free_MEH.root");
  TH1F *fB_12 = static_cast<TH1F *>(fB_pT12->Get("fB_PDG_12"));
  double fB12 = 1/(1 + (1 - fB_12->GetBinContent(1))/fB_12->GetBinContent(1) * AccEff->GetBinContent(2));
  double fB12_E = 1/(1 + (1 - fB_12->GetBinContent(2))/fB_12->GetBinContent(2) * AccEff->GetBinContent(2));
  Fraction->SetBinContent(2, fB12);
  Fraction->SetBinError(2, fB12_E);

  TFile *fB_pT23 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_Full_23_tzJPsi_free_MEH.root");
  TH1F *fB_23 = static_cast<TH1F *>(fB_pT23->Get("fB_PDG_23"));
  double fB23 = 1/(1 + (1 - fB_23->GetBinContent(1))/fB_23->GetBinContent(1) * AccEff->GetBinContent(3));
  double fB23_E = 1/(1 + (1 - fB_23->GetBinContent(2))/fB_23->GetBinContent(2) * AccEff->GetBinContent(3));
  Fraction->SetBinContent(3, fB23);
  Fraction->SetBinError(3, fB23_E);

  TFile *fB_pT34 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_Full_34_tzJPsi_free_MEH.root");
  TH1F *fB_34 = static_cast<TH1F *>(fB_pT34->Get("fB_PDG_34"));
  double fB34 = 1/(1 + (1 - fB_34->GetBinContent(1))/fB_34->GetBinContent(1) * AccEff->GetBinContent(4));
  double fB34_E = 1/(1 + (1 - fB_34->GetBinContent(2))/fB_34->GetBinContent(2) * AccEff->GetBinContent(4));
  Fraction->SetBinContent(4, fB34);
  Fraction->SetBinError(4, fB34_E);

  TFile *fB_pT45 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_Full_45_tzJPsi_free_MEH.root");
  TH1F *fB_45 = static_cast<TH1F *>(fB_pT45->Get("fB_PDG_45"));
  double fB45 = 1/(1 + (1 - fB_45->GetBinContent(1))/fB_45->GetBinContent(1) * AccEff->GetBinContent(5));
  double fB45_E = 1/(1 + (1 - fB_45->GetBinContent(2))/fB_45->GetBinContent(2) * AccEff->GetBinContent(5));
  Fraction->SetBinContent(5, fB45);
  Fraction->SetBinError(5, fB45_E);

  TFile *fB_pT56 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_Full_56_tzJPsi_free_MEH.root");
  TH1F *fB_56 = static_cast<TH1F *>(fB_pT56->Get("fB_PDG_56"));
  double fB56 = 1/(1 + (1 - fB_56->GetBinContent(1))/fB_56->GetBinContent(1) * AccEff->GetBinContent(6));
  double fB56_E = 1/(1 + (1 - fB_56->GetBinContent(2))/fB_56->GetBinContent(2) * AccEff->GetBinContent(6));
  Fraction->SetBinContent(6, fB56);
  Fraction->SetBinError(6, fB56_E);

  TFile *fB_pT68 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_Full_68_tzJPsi_free_MEH.root");
  TH1F *fB_68 = static_cast<TH1F *>(fB_pT68->Get("fB_PDG_68"));
  double fB68 = 1/(1 + (1 - fB_68->GetBinContent(1))/fB_68->GetBinContent(1) * AccEff->GetBinContent(7));
  double fB68_E = 1/(1 + (1 - fB_68->GetBinContent(2))/fB_68->GetBinContent(2) * AccEff->GetBinContent(7));
  Fraction->SetBinContent(7, fB68);
  Fraction->SetBinError(7, fB68_E);

  TFile *fB_pT810 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_Full_810_tzJPsi_free_MEH.root");
  TH1F *fB_810 = static_cast<TH1F *>(fB_pT810->Get("fB_PDG_810"));
  double fB810 = 1/(1 + (1 - fB_810->GetBinContent(1))/fB_810->GetBinContent(1) * AccEff->GetBinContent(8));
  double fB810_E = 1/(1 + (1 - fB_810->GetBinContent(2))/fB_810->GetBinContent(2) * AccEff->GetBinContent(8));
  Fraction->SetBinContent(8, fB810);
  Fraction->SetBinError(8, fB810_E);

  TFile *fB_pT1015 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_Full_1015_tzJPsi_free_MEH.root");
  TH1F *fB_1015 = static_cast<TH1F *>(fB_pT1015->Get("fB_PDG_1015"));
  double fB1015 = 1/(1 + (1 - fB_1015->GetBinContent(1))/fB_1015->GetBinContent(1) * AccEff->GetBinContent(9));
  double fB1015_E = 1/(1 + (1 - fB_1015->GetBinContent(2))/fB_1015->GetBinContent(2) * AccEff->GetBinContent(9));
  Fraction->SetBinContent(9, fB_1015->GetBinContent(1));
  Fraction->SetBinError(9, fB_1015->GetBinContent(2));

  TFile *fB_pT1530 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_Full_1530_tzJPsi_free_MEH.root");
  TH1F *fB_1530 = static_cast<TH1F *>(fB_pT1530->Get("fB_PDG_1530"));
  double fB1530 = 1/(1 + (1 - fB_1530->GetBinContent(1))/fB_1530->GetBinContent(1) * AccEff->GetBinContent(10));
  double fB1530_E = 1/(1 + (1 - fB_1530->GetBinContent(2))/fB_1530->GetBinContent(2) * AccEff->GetBinContent(10));
  Fraction->SetBinContent(10, fB_1530->GetBinContent(1));
  Fraction->SetBinError(10, fB_1530->GetBinContent(2));

  /*TFile *fB_pT1030 = TFile::Open("/Users/emiliebarreau/alice/Macros/fB_PDG_Full_1030_tzJPsi_free_MEH.root");
  TH1F *fB_1030 = static_cast<TH1F *>(fB_pT1030->Get("fB_PDG_1030"));
  double fB1030 = 1/(1 + (1 - fB_1030->GetBinContent(1))/fB_1030->GetBinContent(1) * 1.073322);
  double fB1030_E = 1/(1 + (1 - fB_1030->GetBinContent(2))/fB_1030->GetBinContent(2) * 1.073322);
  Fraction->SetBinContent(9, fB1030);
  Fraction->SetBinError(9, fB1030_E);*/

  Fraction->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  Fraction->GetYaxis()->SetTitle("f_{B}");
  Fraction->SetMinimum(0.);
  Fraction->SetMaximum(0.5);
  Fraction->SetMarkerStyle(34);
  Fraction->SetMarkerSize(2);
  Fraction->SetMarkerColor(4);

  TCanvas *c1 = new TCanvas();
  c1->SetWindowSize(2000, 1600);
  c1->SetTitle("Non-prompt fraction");

  c1->cd();
  Fraction->GetYaxis()->SetTitle("f_{B}");
  //Fraction->GetYaxis()->SetTitleOffset(0.8);
  Fraction->Draw("P E1");
  LHCbis->SetLineColor(kRed+1);
  LHCbis->Draw("SAME P E1");
  LHCb->SetLineColor(kGreen+1);
  LHCb->Draw("SAME P E1");

  TLegend *legend = new TLegend(0.53, 0.16, 0.78, 0.36);
  legend->SetTextSize(0.03);
  legend->SetLineColor(0);
  legend->AddEntry(Fraction, "ALICE 2023pp -3.6 < y < -2.5", "lp");
  legend->AddEntry(LHCbis, "2.5 < y_{LHCb} < 3", "lp");
  legend->AddEntry(LHCb, "3 < y_{LHCb} < 3.5", "lp");
  legend->Draw();

  std::cout << fB01 << " " << fB01_E << std::endl;
  std::cout << fB12 << " " << fB12_E << std::endl;
  std::cout << fB23 << " " << fB23_E << std::endl;
  std::cout << fB34 << " " << fB34_E << std::endl;
  std::cout << fB45 << " " << fB45_E << std::endl;
  std::cout << fB56 << " " << fB56_E << std::endl;
  std::cout << fB68 << " " << fB68_E << std::endl;
  std::cout << fB810 << " " << fB810_E << std::endl;
  //std::cout << fB1030 << " " << fB1030_E << std::endl;
}