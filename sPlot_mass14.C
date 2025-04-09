#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooAddModel.h"
#include "RooGaussModel.h"
#include "RooResolutionModel.h"
#include "RooNumConvolution.h"
#include "RooAddition.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooChi2Var.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooDecay.h"
#include "RooDoubleCB.h"
#include "RooExponential.h"
#include "RooExtendPdf.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooNA60.h"
#include "RooPlot.h"
#include "RooProdPdf.h"
#include "RooProduct.h"
#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooTruthModel.h"
#include "RooVWG.h"
#include "RooVWG_C.h"
#include "RooWorkspace.h"
#include "RooBumpExp.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLegend.h"
#include "tree2DataSet.C"
#include "RooExp_Affine.h"
#include "RooNA60.h"
#include "RooKeysPdf.h"
#include "TH1D.h"
#include <TF2.h>
#include <TH1.h>
#include <TH2.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <TMath.h>
#include "TLine.h"

using namespace RooFit;
using namespace RooStats;

void FitMassModel(RooWorkspace &);
void DoSPlot(RooWorkspace &);
void MakePlots(RooWorkspace &);
void AddTauzSignalModel(RooWorkspace &);
void AddTauzBkgModel(RooWorkspace &);
void TemplateFit(RooWorkspace &);
void MakeSignalPlot(RooWorkspace &);
void MakeBkgPlot(RooWorkspace &);
void MakeTotalPlot(RooWorkspace &);
void MakeFit2D(RooWorkspace &);
void Make2DPlot(RooWorkspace &);

void sPlot_mass14() {
  // create a workspace to manage the project
  RooWorkspace wspace{"myWSbisbisbisbisbisbisbisbisbis"};

  // add the mass signal and background models
  FitMassModel(wspace);

  // new dataset with sWeights added for every event
  DoSPlot(wspace);

  // make some plots showing discriminating variable & control variable after
  MakePlots(wspace);

  //TemplateFit(wspace);

  // add the tauz signal model
  //AddTauzSignalModel(wspace);

  // add the tauz bkg model
  //AddTauzBkgModel(wspace);

  // make tauz signal plot
  //MakeSignalPlot(wspace);

  // make tauz bkg plot
  //MakeBkgPlot(wspace);

  //MakeTotalPlot(wspace);

  //MakeFit2D(wspace);
}

TH1* rebinctauBkghist(RooWorkspace& ws, TH1 *hist, double xmin, double xmax)
{
  TH1 *hcopy = (TH1*) hist->Clone("hcopy");

  // range of the new hist
  int imin = hcopy->FindBin(xmin);
  if (imin>=hcopy->GetNbinsX()) imin=1;
  int imax = hcopy->FindBin(0.999999*xmax);
  if (imax<=1) imax=hcopy->GetNbinsX();

  vector<double> newbins;
  newbins.push_back(hcopy->GetBinLowEdge(imin));
  for (int i=imin; i<=imax; i++) {
    if (hcopy->GetBinContent(i)>0.1) {
      newbins.push_back(hcopy->GetBinLowEdge(i)+hcopy->GetBinWidth(i));
    } else {
      int nrebin=2;
      for (i++; i<=imax; i++) {
        if (hcopy->GetBinContent(i)>0.00000001) {
          newbins.push_back(hcopy->GetBinLowEdge(i)+hcopy->GetBinWidth(i));
          hcopy->SetBinContent(i,hcopy->GetBinContent(i)/nrebin);
          break;
        }
        nrebin++;
      }
    }
  }

  if (xmin < newbins[1]) newbins[0] = xmin;
  if (xmax > newbins[newbins.size()-2]) newbins[newbins.size()-1] = xmax;

  TH1 *ans = hcopy->Rebin(newbins.size()-1,"hnew",newbins.data());

  delete hcopy;
  return ans;
};

//____________________________________
void FitMassModel(RooWorkspace &myws) {
  gROOT->ProcessLine(".L ./RooDoubleCB.cxx+");
  gROOT->ProcessLine(".L ./RooExp_Affine.cxx+");
  gROOT->ProcessLine(".L ./RooBumpExp.cxx+");
  gROOT->ProcessLine(".L ./RooVWG.cxx+");
  gROOT->ProcessLine(".L ./RooBumpVWG.cxx+");
  gROOT->ProcessLine(".L ./RooNA60.cxx+");
  //gROOT->ProcessLine(".L ./RooVWG.cxx+");

  double massMin = 2.5;
  double massMax = 3.5;
  double tau_min = -10.;
  double tau_max = 10.;
  double ptmin = 2;
  double ptmax = 3;
  double chi2opti = 30;
  bool fixtails = false;

  RooRealVar *mass =
      new RooRealVar("mass", "M_{#mu#mu}", massMin, massMax, "GeV/c^{2}");
  RooRealVar *tau = new RooRealVar("tau", "l_{J/#Psi}", tau_min, tau_max, "mm");

  mass->setRange("Mass", massMin, massMax);
  
  RooDataSet *data = tree2DataSet("AO2D_NEW_DCA.txt","O2rtdimuonall", tau, mass, ptmin, ptmax, chi2opti);
  //RooDataSet *data = tree2DataSet("/Users/emiliebarreau/alice/Python_Scripts/LHC23_pass4_skimmed_FULL/AO2D_NEW_KF.txt","O2rtdimuonall", tau, mass);
  //RooDataSet *data = tree2DataSet("/Users/emiliebarreau/alice/Python_Scripts/LHC23_pass4_skimmed_FULL/AO2D_NEW_DCA.txt","O2rtdimuonall", tau, mass);
  //RooDataSet *data = tree2DataSet("/Users/emiliebarreau/alice/Python_Scripts/LHC23_pass4/LHC23_pass4_DCA.txt","O2rtdimuonall", tau, mass);
  //RooDataSet *data = tree2DataSet("/Users/emiliebarreau/alice/Python_Scripts/LHC23_pass4_skimmed_FULL/LHC23_pass4_skimmed_DCA.txt", "O2rtdimuonall", tau, mass);
  //RooDataSet *data = tree2DataSet("/Users/emiliebarreau/alice/Python_Scripts/LHC22o_pass7/LHC22o_pass7_DCA.txt","O2rtdimuonall", tau, mass);
  
  //double numEntries = data->sumEntries(); 
  double numEntries = data->numEntries();

  myws.import(*mass);
  myws.import(*tau);
  myws.import(*data, Rename("data"));

  myws.factory(Form("nsig[%.0f,%.0f]", 0., numEntries));
  myws.factory(Form("nbkg[%.0f,%.0f]", 0., numEntries));

  /*************************************/
  myws.factory(Form("pT1[%.0f]", ptmin));
  myws.factory(Form("pT2[%.0f]", ptmax));
  /*************************************/

  stringstream ss;
  ss << myws.var("pT1")->getValV() ;
  string test = ss.str();
  stringstream sss;
  sss << myws.var("pT2")->getValV() ;
  string testbis = sss.str();

  TFile *Param_01 = TFile::Open(Form("corr_frac_%s%s_EM_newChi2.root", test.c_str(), testbis.c_str()));
  TH1F *Histo_01 = static_cast<TH1F *>(Param_01->Get(Form("corr_frac_%s%s", test.c_str(), testbis.c_str())));

  //if bkg is exp
  //myws.factory(Form("A[%.4f]", Histo_01->GetBinContent(1)));
  //myws.factory(Form("B[%.4f]", Histo_01->GetBinContent(2)));
  //myws.factory(Form("C[%.4f]", Histo_01->GetBinContent(3)));
  //myws.factory(Form("cste[%.4f]", Histo_01->GetBinContent(4)));
  //if bkg is VWG
  myws.factory(Form("A1[%.4f]", Histo_01->GetBinContent(1)));
  myws.factory(Form("B1[%.4f]", Histo_01->GetBinContent(2)));
  myws.factory(Form("C1[%.4f]", Histo_01->GetBinContent(3)));
  myws.factory(Form("cste1[%.4f]", Histo_01->GetBinContent(4)));

  //if sig is CB2
  //myws.factory(Form("mJPsi[%.4f,%.4f,%.4f]",  3.0980  , 3.0 , 3.2));
  //myws.factory(Form("sigma[%.4f,%.4f,%.4f]",  0.0487  , 0.01, 0.08));
  //myws.factory(Form("n1[%.4f,%.4f,%.4f]",     3.5150  , 2.99, 3.98));
  //myws.factory(Form("alpha1[%.4f,%.4f,%.4f]", 0.7078  , 0.53, 0.86));
  //myws.factory(Form("n2[%.4f,%.4f,%.4f]",     1.0158  , 0.16, 0.75));
  //myws.factory(Form("alpha2[%.4f,%.4f,%.4f]", 2.8946  , 2.05, 3.51));

  //if sig is NA60
  myws.factory(Form("mJPsi[%.4f,%.4f,%.4f]",  3.09954   , 3.0 , 3.2));
  myws.factory(Form("sigma[%.4f,%.4f,%.4f]", 0.0471576  , 0.01, 0.08));
  myws.factory(Form("alpha1[%.4f,%.4f,%.4f]",-0.366456  , -3.0, 3.00));
  myws.factory(Form("p1[%.4f,%.4f,%.4f]",    0.146592   , -2.0, 2.00));
  myws.factory(Form("p2[%.4f,%.4f,%.4f]",     1.34788   , -2.0, 3.00));
  myws.factory(Form("p3[%.4f,%.4f,%.4f]",    0.146592   , -2.0, 2.0));
  myws.factory(Form("alpha2[%.4f,%.4f,%.4f]", 1.19215   ,  0.0, 4.0));
  myws.factory(Form("p12[%.4f,%.4f,%.4f]",   0.294032   , -1.0, 1.0));
  myws.factory(Form("p22[%.4f,%.4f,%.4f]",   2.000000   , -5.0, 5.0));
  myws.factory(Form("p32[%.4f,%.4f,%.4f]",   0.294032   , -2.0, 2.0));
                                              
  //For Exp
  //myws.factory(Form("lambdas[%.4f,%.4f,%.4f]", -2.6155, -4., 1.2));

  //For VWG
  myws.factory(Form("A[%.4f,%.4f,%.4f]", 1.88, -2., 2.));
  myws.factory(Form("B[%.4f,%.4f,%.4f]", 0.94, -2., 2.));
  myws.factory(Form("C[%.4f,%.4f,%.4f]", -0.01, -1., 1.));
  myws.factory(Form("cste[%.4f]", 0.0));
  //myws.var("cste")->setConstant(kTRUE);

  //0-1 values DCA Chi2 30
  /*//bkg Exp sig CB2 data tails
  myws.factory(Form("mJPsi[%.4f]", 3.0980));
  myws.factory(Form("sigma[%.4f]", 0.0487));
  myws.factory(Form("n1[%.4f]", 5.8494));
  myws.factory(Form("alpha1[%.4f]", 0.6926));
  myws.factory(Form("n2[%.4f]", 0.1531));
  myws.factory(Form("alpha2[%.4f]", 2.2491));
  myws.factory(Form("lambdas[%.4f]", -2.6155));*/
  /*bkg Exp sig CB2 MC tails
  myws.factory(Form("mJPsi[%.4f]",  3.09777  ));
  myws.factory(Form("sigma[%.4f]",  0.0488182));
  myws.factory(Form("n1[%.4f]",     5.85     ));
  myws.factory(Form("alpha1[%.4f]", 0.69     ));
  myws.factory(Form("n2[%.4f]",     0.15     ));
  myws.factory(Form("alpha2[%.4f]", 2.25     ));
  myws.factory(Form("lambdas[%.4f]",-2.61884 ));*/
  /*bkg VWG sig CB2 data tails
  myws.factory(Form("mJPsi[%.4f]", 3.09688));
  myws.factory(Form("sigma[%.4f]", 0.04831));
  myws.factory(Form("n1[%.4f]", 1.93564));
  myws.factory(Form("alpha1[%.4f]",  0.751267));
  myws.factory(Form("n2[%.4f]", 1.27298));
  myws.factory(Form("alpha2[%.4f]", 1.88225));
  myws.factory(Form("A[%.4f]", 0.342657));
  myws.factory(Form("B[%.4f]", 0.946094));
  myws.factory(Form("C[%.4f]", -0.0177515));
 *bkg VWG sig NA60 data tails
  myws.factory(Form("mJPsi[%.4f]",  3.09814));
  myws.factory(Form("sigma[%.4f]",  0.0476203));
  myws.factory(Form("alpha1[%.4f]", -0.843727));
  myws.factory(Form("p1[%.4f]",     0.277282));
  myws.factory(Form("p2[%.4f]",     1       ));
  myws.factory(Form("p3[%.4f]",     0.133241));
  myws.factory(Form("alpha2[%.4f]", 1.40556));
  myws.factory(Form("p12[%.4f]",    0.293451 ));
  myws.factory(Form("p22[%.4f]",    2.04106  ));
  myws.factory(Form("p32[%.4f]",    0.304531 ));
  myws.factory(Form("A[%.4f]",      0.462995  ));
  myws.factory(Form("B[%.4f]",      1.01737   ));
  myws.factory(Form("C[%.4f]",      -0.0104918));
  myw myws.factory(Form("cste[%.4f]", 0));*/
  /*bkg VWG sig CB2 MC tails
  myws.factory(Form("mJPsi[%.4f]",  ));
  myws.factory(Form("sigma[%.4f]",  ));
  myws.factory(Form("n1[%.4f]",     ));
  myws.factory(Form("alpha1[%.4f]", ));
  myws.factory(Form("n2[%.4f]",     ));
  myws.factory(Form("alpha2[%.4f]", ));
  myws.factory(Form("A[%.4f]",      ));
  myws.factory(Form("B[%.4f]",      ));
  myws.factory(Form("C[%.4f]",      ));
  myws.factory(Form("cste[%.4f]",   ));*/

  //1-2 values DCA Chi2 30
  /*bkg Exp sig CB2 data tails
  myws.factory(Form("mJPsi[%.4f]", 3.0959));
  myws.factory(Form("sigma[%.4f]", 0.0461));
  myws.factory(Form("n1[%.4f]", 6.3014));
  myws.factory(Form("alpha1[%.4f]", 0.7097));
  myws.factory(Form("n2[%.4f]", 58.8212));
  myws.factory(Form("alpha2[%.4f]", 1.8013));
  myws.factory(Form("lambdas[%.4f]", -1.5753));*/
  /*bkg Exp sig CB2 MC tail >> does not work. chi2/ndf is about 9
  myws.factory(Form("mJPsi[%.4f]",   3.10177  ));
  myws.factory(Form("sigma[%.4f]",   0.0418373));
  myws.factory(Form("n1[%.4f]",      3.93     ));
  myws.factory(Form("alpha1[%.4f]",  0.56     ));
  myws.factory(Form("n2[%.4f]",      0.880033 ));
  myws.factory(Form("alpha2[%.4f]",  2.5      ));
  myws.factory(Form("lambdas[%.4f]", -1.54918 ));*/
  /*bkg VWG sig CB2 data tails
  myws.factory(Form("mJPsi[%.4f]",  3.09501));
  myws.factory(Form("sigma[%.4f]",  0.0478742));
  myws.factory(Form("n1[%.4f]",     4.46303));
  myws.factory(Form("alpha1[%.4f]", 0.763432));
  myws.factory(Form("n2[%.4f]",     1.3603));
  myws.factory(Form("alpha2[%.4f]", 2.49067));
  myws.factory(Form("A[%.4f]",      0.217098));
  myws.factory(Form("B[%.4f]",      1.12318));
  myws.factory(Form("C[%.4f]",      0.000562958));
  myws.factory(Form("cste[%.4f]",   0.0156945));*/
 /*bkg VWG sig NA60 data tails
  myws.factory(Form("mJPsi[%.4f]",   3.0954     ));
  myws.factory(Form("sigma[%.4f]",   0.0479411  ));
  myws.factory(Form("alpha1[%.4f]",  -0.969597  ));
  myws.factory(Form("p1[%.4f]",      0.280235   ));
  myws.factory(Form("p2[%.4f]",      0.980661   ));
  myws.factory(Form("p3[%.4f]",      0.0715572  ));
  myws.factory(Form("alpha2[%.4f]",  1.70148    ));
  myws.factory(Form("p12[%.4f]",     0.26947    ));
  myws.factory(Form("p22[%.4f]",     2.02555    ));
  myws.factory(Form("p32[%.4f]",     0.290416   ));
  myws.factory(Form("A[%.4f]",       1.07803    ));
  myws.factory(Form("B[%.4f]",       1.00568    ));
  myws.factory(Form("C[%.4f]",       -0.00311195));
  myw myws.factory(Form("cste[%.4f]", 0));*/
  /*bkg VWG sig CB2 MC tails
  myws.factory(Form("mJPsi[%.4f]",  ));
  myws.factory(Form("sigma[%.4f]",  ));
  myws.factory(Form("n1[%.4f]",     ));
  myws.factory(Form("alpha1[%.4f]", ));
  myws.factory(Form("n2[%.4f]",     ));
  myws.factory(Form("alpha2[%.4f]", ));
  myws.factory(Form("A[%.4f]",      ));
  myws.factory(Form("B[%.4f]",      ));
  myws.factory(Form("C[%.4f]",      ));
  myws.factory(Form("cste[%.4f]",   ));*/

  //2-3 values DCA Chi2 35
  /*bkg Exp sig CB2 data tails
  myws.factory(Form("mJPsi[%.4f]", 3.0951));
  myws.factory(Form("sigma[%.4f]", 0.0500));
  myws.factory(Form("n1[%.4f]", 23.6946));
  myws.factory(Form("alpha1[%.4f]", 0.8144));
  myws.factory(Form("n2[%.4f]", 0.0));
  myws.factory(Form("alpha2[%.4f]", 3.0575));
  myws.factory(Form("lambdas[%.4f]", -2.1092));*/
  /*bkg Exp sig CB2 MC tail >> does not work chi2/ndf is about 11s
  myws.factory(Form("mJPsi[%.4f]",   3.1001   ));
  myws.factory(Form("sigma[%.4f]",   0.0451444));
  myws.factory(Form("n1[%.4f]",      3.47     ));
  myws.factory(Form("alpha1[%.4f]",  0.7      ));
  myws.factory(Form("n2[%.4f]",      0.83     ));
  myws.factory(Form("alpha2[%.4f]",  2.60007  ));
  myws.factory(Form("lambdas[%.4f]", -1.90223 ));*/
  /*bkg VWG sig CB2 data tails
  myws.factory(Form("mJPsi[%.4f]",  3.09502));
  myws.factory(Form("sigma[%.4f]",  0.0500455));
  myws.factory(Form("n1[%.4f]",     10.7335));
  myws.factory(Form("alpha1[%.4f]", 0.815036));
  myws.factory(Form("n2[%.4f]",     5.77415));
  myws.factory(Form("alpha2[%.4f]", 5.9344));
  myws.factory(Form("A[%.4f]",      0.2519));
  myws.factory(Form("B[%.4f]",      0.71853));
  myws.factory(Form("C[%.4f]",      0.0160211));
  myws.factory(Form("cste[%.4f]",   0.00306005));*/
  /*bkg VWG sig CB2 MC tails
  myws.factory(Form("mJPsi[%.4f]",  ));
  myws.factory(Form("sigma[%.4f]",  ));
  myws.factory(Form("n1[%.4f]",     ));
  myws.factory(Form("alpha1[%.4f]", ));
  myws.factory(Form("n2[%.4f]",     ));
  myws.factory(Form("alpha2[%.4f]", ));
  myws.factory(Form("A[%.4f]",      ));
  myws.factory(Form("B[%.4f]",      ));
  myws.factory(Form("C[%.4f]",      ));
  myws.factory(Form("cste[%.4f]",   ));*/

  //3-4 values DCA Chi2 40
  /*bkg Exp sig CB2 data tails
  myws.factory(Form("mJPsi[%.4f]", 3.0936));
  myws.factory(Form("sigma[%.4f]", 0.0541));
  myws.factory(Form("n1[%.4f]", 3.3052));
  myws.factory(Form("alpha1[%.4f]", 0.9609));
  myws.factory(Form("n2[%.4f]", 1.9236));
  myws.factory(Form("alpha2[%.4f]", 2.2969));
  myws.factory(Form("lambdas[%.4f]", -1.6397));*/
  /*bkg Exp sig CB2 MC tail >> does not work chi2 is about 5s
  myws.factory(Form("mJPsi[%.4f]",   3.09792  ));
  myws.factory(Form("sigma[%.4f]",   0.0507406));
  myws.factory(Form("n1[%.4f]",      3.49     ));
  myws.factory(Form("alpha1[%.4f]",  0.75     ));
  myws.factory(Form("n2[%.4f]",      0.990003 ));
  myws.factory(Form("alpha2[%.4f]",  2.57     ));
  myws.factory(Form("lambdas[%.4f]", -1.79626 ));*/
  /*bkg VWG sig CB2 data tails
  myws.factory(Form("mJPsi[%.4f]",  3.09408   ));
  myws.factory(Form("sigma[%.4f]",  0.0529031 ));
  myws.factory(Form("n1[%.4f]",     9.25478   ));
  myws.factory(Form("alpha1[%.4f]", 0.896013  ));
  myws.factory(Form("n2[%.4f]",     8.21295   ));
  myws.factory(Form("alpha2[%.4f]", 2.40457   );
  myws.factory(Form("A[%.4f]",      2.3492    );
  myws.factory(Form("B[%.4f]",      0.530676  ));
  myws.factory(Form("C[%.4f]",      -0.0162912));
  myws.factory(Form("cste[%.4f]",   0.0707107 ));*/
  /*bkg VWG sig CB2 MC tails
  myws.factory(Form("mJPsi[%.4f]",  ));
  myws.factory(Form("sigma[%.4f]",  ));
  myws.factory(Form("n1[%.4f]",     ));
  myws.factory(Form("alpha1[%.4f]", ));
  myws.factory(Form("n2[%.4f]",     ));
  myws.factory(Form("alpha2[%.4f]", ));
  myws.factory(Form("A[%.4f]",      ));
  myws.factory(Form("B[%.4f]",      ));
  myws.factory(Form("C[%.4f]",      ));
  myws.factory(Form("cste[%.4f]",   ));*/

  //4-5 values DCA Chi2 45
  /*bkg Exp sig CB2 data tails
  myws.factory(Form("mJPsi[%.4f]", 3.0923));
  myws.factory(Form("sigma[%.4f]", 0.0589));
  myws.factory(Form("n1[%.4f]", 3.2672));
  myws.factory(Form("alpha1[%.4f]", 1.0243));
  myws.factory(Form("n2[%.4f]", 2.4629));
  myws.factory(Form("alpha2[%.4f]", 2.1946));
  myws.factory(Form("lambdas[%.4f]", -1.2638));*/
  /*bkg Exp sig CB2 MC tails >> does not work chi2/ndf is about 6
  myws.factory(Form("mJPsi[%.4f]",   3.09835 ));
  myws.factory(Form("sigma[%.4f]",   0.053679));
  myws.factory(Form("n1[%.4f]",      3.82    ));
  myws.factory(Form("alpha1[%.4f]",  0.75    ));
  myws.factory(Form("n2[%.4f]",      1.31013 ));
  myws.factory(Form("alpha2[%.4f]",  3.09936 ));
  myws.factory(Form("lambdas[%.4f]", -1.04197));*/
  /*bkg VWG sig CB2 data tails
  myws.factory(Form("mJPsi[%.4f]",  3.09485  ));
  myws.factory(Form("sigma[%.4f]",  0.0563413));
  myws.factory(Form("n1[%.4f]",     70.3444  ));
  myws.factory(Form("alpha1[%.4f]", 0.79876  ));
  myws.factory(Form("n2[%.4f]",     59.7456  ));
  myws.factory(Form("alpha2[%.4f]", 1.8853   );
  myws.factory(Form("A[%.4f]",      1.56326  );
  myws.factory(Form("B[%.4f]",      1.19978  ));
  myws.factory(Form("C[%.4f]",      0.0280159));
  myws.factory(Form("cste[%.4f]",   -0.132399));*/
  /*bkg VWG sig CB2 MC tails
  myws.factory(Form("mJPsi[%.4f]",  ));
  myws.factory(Form("sigma[%.4f]",  ));
  myws.factory(Form("n1[%.4f]",     ));
  myws.factory(Form("alpha1[%.4f]", ));
  myws.factory(Form("n2[%.4f]",     ));
  myws.factory(Form("alpha2[%.4f]", ));
  myws.factory(Form("A[%.4f]",      ));
  myws.factory(Form("B[%.4f]",      ));
  myws.factory(Form("C[%.4f]",      ));
  myws.factory(Form("cste[%.4f]",   ));*/

  //5-6 values DCA Chi2 50
  /*bkg Exp sig CB2 data tails
  myws.factory(Form("mJPsi[%.4f]", 3.0961));
  myws.factory(Form("sigma[%.4f]", 0.0584));
  myws.factory(Form("n1[%.4f]", 30.7687));
  myws.factory(Form("alpha1[%.4f]", 0.8261));
  myws.factory(Form("n2[%.4f]", 32.8423));
  myws.factory(Form("alpha2[%.4f]", 4.8511));
  myws.factory(Form("lambdas[%.4f]", -0.9498));*/
  /*bkg Exp sig CB2 MC tail >> works but allowing tail parameters to fluctuate within MC uncertaintiess
  myws.factory(Form("mJPsi[%.4f]",   3.09451  ));
  myws.factory(Form("sigma[%.4f]",   0.059862 ));
  myws.factory(Form("n1[%.4f]",      3.69999  ));
  myws.factory(Form("alpha1[%.4f]",  1.00925  ));
  myws.factory(Form("n2[%.4f]",      1.05924  ));
  myws.factory(Form("alpha2[%.4f]",  2.91     ));
  myws.factory(Form("lambdas[%.4f]", -0.893192));*/
  /*bkg VWG sig CB2 data tails
  myws.factory(Form("mJPsi[%.4f]", 3.09623     ));
  myws.factory(Form("sigma[%.4f]", 0.0576572   ));
  myws.factory(Form("n1[%.4f]",    42.3767     ));
  myws.factory(Form("alpha1[%.4f]",0.805869    ));
  myws.factory(Form("n2[%.4f]",    49.7249     ));
  myws.factory(Form("alpha2[%.4f]",1.84187     );
  myws.factory(Form("A[%.4f]",     0.619382    );
  myws.factory(Form("B[%.4f]",     1.36987     ));
  myws.factory(Form("C[%.4f]",     0.0300386   ));
  myws.factory(Form("cste[%.4f]",  -0.00249156 ));*/
  /*bkg VWG sig CB2 MC tails
  myws.factory(Form("mJPsi[%.4f]",  ));
  myws.factory(Form("sigma[%.4f]",  ));
  myws.factory(Form("n1[%.4f]",     ));
  myws.factory(Form("alpha1[%.4f]", ));
  myws.factory(Form("n2[%.4f]",     ));
  myws.factory(Form("alpha2[%.4f]", ));
  myws.factory(Form("A[%.4f]",      ));
  myws.factory(Form("B[%.4f]",      ));
  myws.factory(Form("C[%.4f]",      ));
  myws.factory(Form("cste[%.4f]",   ));*/

  //6-8 values DCA Chi2 60
  /*bkg Exp sig CB2 data tails
  myws.factory(Form("mJPsi[%.4f]", 3.0991));
  myws.factory(Form("sigma[%.4f]", 0.0596));
  myws.factory(Form("n1[%.4f]", 23.1637));
  myws.factory(Form("alpha1[%.4f]", 0.8285));
  myws.factory(Form("n2[%.4f]", 7.2662));
  myws.factory(Form("alpha2[%.4f]", 4.3479));
  myws.factory(Form("lambdas[%.4f]", -0.7814));*/
  /*bkg Exp sig CB2 MC tails >> works but allowing tail parameters to fluctuate within MC uncertaintiess
  myws.factory(Form("mJPsi[%.4f]",   3.09928  ));
  myws.factory(Form("sigma[%.4f]",   0.0597928));
  myws.factory(Form("n1[%.4f]",      2.72     ));
  myws.factory(Form("alpha1[%.4f]",  0.959998 ));
  myws.factory(Form("n2[%.4f]",      0.770274 ));
  myws.factory(Form("alpha2[%.4f]",  2.47     ));
  myws.factory(Form("lambdas[%.4f]", -0.781627));*/
  /*bkg VWG sig CB2 data tails
  myws.factory(Form("mJPsi[%.4f]",  3.10048   ));
  myws.factory(Form("sigma[%.4f]",  0.0576714 ));
  myws.factory(Form("n1[%.4f]",     69.3282   ));
  myws.factory(Form("alpha1[%.4f]", 0.772313  ));
  myws.factory(Form("n2[%.4f]",     7.71559   ));
  myws.factory(Form("alpha2[%.4f]", 1.87408   );
  myws.factory(Form("A[%.4f]",      0.80793   );
  myws.factory(Form("B[%.4f]",      1.54734   ));
  myws.factory(Form("C[%.4f]",      0.034879  ));
  myws.factory(Form("cste[%.4f]",   -0.0551348));*/
  /*bkg VWG sig CB2 MC tails
  myws.factory(Form("mJPsi[%.4f]",  ));
  myws.factory(Form("sigma[%.4f]",  ));
  myws.factory(Form("n1[%.4f]",     ));
  myws.factory(Form("alpha1[%.4f]", ));
  myws.factory(Form("n2[%.4f]",     ));
  myws.factory(Form("alpha2[%.4f]", ));
  myws.factory(Form("A[%.4f]",      ));
  myws.factory(Form("B[%.4f]",      ));
  myws.factory(Form("C[%.4f]",      ));
  myws.factory(Form("cste[%.4f]",   ));*/

  //8-10 values DCA Chi2 65
  /*bkg Exp sig CB2 data tails
  myws.factory(Form("mJPsi[%.4f]", 3.0997));
  myws.factory(Form("sigma[%.4f]", 0.0607));
  myws.factory(Form("n1[%.4f]", 85.2370));
  myws.factory(Form("alpha1[%.4f]", 0.8919));
  myws.factory(Form("n2[%.4f]", 59.9946));
  myws.factory(Form("alpha2[%.4f]", 5.6002));
  myws.factory(Form("lambdas[%.4f]", -0.7075));*/
  /*bkg Exp sig CB2 MC tails >> works but allowing tail parameters to fluctuate within MC uncertaintiess
  myws.factory(Form("mJPsi[%.4f]",   3.0987   ));
  myws.factory(Form("sigma[%.4f]",   0.0614964));
  myws.factory(Form("n1[%.4f]",      2.95989  ));
  myws.factory(Form("alpha1[%.4f]",  1.09999  ));
  myws.factory(Form("n2[%.4f]",      2.24931  ));
  myws.factory(Form("alpha2[%.4f]",  2.44     ));
  myws.factory(Form("lambdas[%.4f]", -0.637181));*/
  /*bkg VWG sig CB2 data tails
  myws.factory(Form("mJPsi[%.4f]",  3.10022  ));
  myws.factory(Form("sigma[%.4f]",  0.0588446));
  myws.factory(Form("n1[%.4f]",     9.95232  ));
  myws.factory(Form("alpha1[%.4f]", 0.866689 ));
  myws.factory(Form("n2[%.4f]",     4.89691  ));
  myws.factory(Form("alpha2[%.4f]", 1.63752  );
  myws.factory(Form("A[%.4f]",      1.12131  );
  myws.factory(Form("B[%.4f]",      0.828839 ));
  myws.factory(Form("C[%.4f]",      0.0212713));
  myws.factory(Form("cste[%.4f]",   0.139733 ));*/
  /*bkg VWG sig CB2 MC tails
  myws.factory(Form("mJPsi[%.4f]",  ));
  myws.factory(Form("sigma[%.4f]",  ));
  myws.factory(Form("n1[%.4f]",     ));
  myws.factory(Form("alpha1[%.4f]", ));
  myws.factory(Form("n2[%.4f]",     ));
  myws.factory(Form("alpha2[%.4f]", ));
  myws.factory(Form("A[%.4f]",      ));
  myws.factory(Form("B[%.4f]",      ));
  myws.factory(Form("C[%.4f]",      ));
  myws.factory(Form("cste[%.4f]",   ));*/

  //10-30 values DCA Chi2 90
  /*bkg Exp sig CB2 data tails
  myws.factory(Form("mJPsi[%.4f]", 3.1066));
  myws.factory(Form("sigma[%.4f]", 0.0621));
  myws.factory(Form("n1[%.4f]", 9.3165));
  myws.factory(Form("alpha1[%.4f]", 0.7615));
  myws.factory(Form("n2[%.4f]", 3.1386));
  myws.factory(Form("alpha2[%.4f]", 1.3483));
  myws.factory(Form("lambdas[%.4f]", -1.0664));*/
  /*bkg Exp sig CB2 MC tails
  myws.factory(Form("mJPsi[%.4f]",   ));
  myws.factory(Form("sigma[%.4f]",   ));
  myws.factory(Form("n1[%.4f]",      ));
  myws.factory(Form("alpha1[%.4f]",  ));
  myws.factory(Form("n2[%.4f]",      ));
  myws.factory(Form("alpha2[%.4f]",  ));
  myws.factory(Form("lambdas[%.4f]", ));*/
  /*bkg VWG sig CB2 data tails
  myws.factory(Form("mJPsi[%.4f]",  3.1067   ));
  myws.factory(Form("sigma[%.4f]",  0.0589177));
  myws.factory(Form("n1[%.4f]",     83.6916  ));
  myws.factory(Form("alpha1[%.4f]", 0.703115 ));
  myws.factory(Form("n2[%.4f]",     7.69012  ));
  myws.factory(Form("alpha2[%.4f]", 1.17105  );
  myws.factory(Form("A[%.4f]",      0.921675 );
  myws.factory(Form("B[%.4f]",      1.86516  ));
  myws.factory(Form("C[%.4f]",      0.0413855));
  myws.factory(Form("cste[%.4f]",   -0.24591 ));*/
  /*bkg VWG sig CB2 MC tails
  myws.factory(Form("mJPsi[%.4f]",  ));
  myws.factory(Form("sigma[%.4f]",  ));
  myws.factory(Form("n1[%.4f]",     ));
  myws.factory(Form("alpha1[%.4f]", ));
  myws.factory(Form("n2[%.4f]",     ));
  myws.factory(Form("alpha2[%.4f]", ));
  myws.factory(Form("A[%.4f]",      ));
  myws.factory(Form("B[%.4f]",      ));
  myws.factory(Form("C[%.4f]",      ));
  myws.factory(Form("cste[%.4f]",   ));*/

  //10-15 values DCA Chi2 90
  /*bkg Exp sig CB2 data tails
  myws.factory(Form("mJPsi[%.4f]", 3.1049));
  myws.factory(Form("sigma[%.4f]", 0.0595));
  myws.factory(Form("n1[%.4f]", 7.0747));
  myws.factory(Form("alpha1[%.4f]", 0.8011));
  myws.factory(Form("n2[%.4f]", 6.0106));
  myws.factory(Form("alpha2[%.4f]", 1.2714));
  myws.factory(Form("lambdas[%.4f]", -0.9021));*/
  /*bkg Exp sig CB2 MC tails
  myws.factory(Form("mJPsi[%.4f]",   3.10316  ));
  myws.factory(Form("sigma[%.4f]",   0.0655498));
  myws.factory(Form("n1[%.4f]",      2.6608   ));
  myws.factory(Form("alpha1[%.4f]",  1.0681   ));
  myws.factory(Form("n2[%.4f]",      2        ));
  myws.factory(Form("alpha2[%.4f]",  2.329    ));
  myws.factory(Form("lambdas[%.4f]", -0.623135));*/
  /*bkg VWG sig CB2 data tails
  myws.factory(Form("mJPsi[%.4f]",  3.1119   ));
  myws.factory(Form("sigma[%.4f]",  0.0420464));
  myws.factory(Form("n1[%.4f]",     36.6923  ));
  myws.factory(Form("alpha1[%.4f]", 0.459533 ));
  myws.factory(Form("n2[%.4f]",     41.7978  ));
  myws.factory(Form("alpha2[%.4f]", 0.685303 );
  myws.factory(Form("A[%.4f]",      0.173982 );
  myws.factory(Form("B[%.4f]",      1.23614  ));
  myws.factory(Form("C[%.4f]",      0.015476 ));*/
  /*bkg VWG sig CB2 MC tails
  myws.factory(Form("mJPsi[%.4f]",  ));
  myws.factory(Form("sigma[%.4f]",  ));
  myws.factory(Form("n1[%.4f]",     ));
  myws.factory(Form("alpha1[%.4f]", ));
  myws.factory(Form("n2[%.4f]",     ));
  myws.factory(Form("alpha2[%.4f]", ));
  myws.factory(Form("A[%.4f]",      ));
  myws.factory(Form("B[%.4f]",      ));
  myws.factory(Form("C[%.4f]",      ));*/

  //15-30 values DCA Chi2 90
  /*bkg Exp sig CB2 data tails
  myws.factory(Form("mJPsi[%.4f]", 3.1117));
  myws.factory(Form("sigma[%.4f]", 0.0622));
  myws.factory(Form("n1[%.4f]", 99.9159));
  myws.factory(Form("alpha1[%.4f]", 0.5817));
  myws.factory(Form("n2[%.4f]", 59.8058));
  myws.factory(Form("alpha2[%.4f]", 0.8516));
  myws.factory(Form("lambdas[%.4f]", -1.0701));*/
  /*bkg Exp sig CB2 MC tails
  myws.factory(Form("mJPsi[%.4f]",   3.10313  ));
  myws.factory(Form("sigma[%.4f]",   0.08     ));
  myws.factory(Form("n1[%.4f]",      2.5243   ));
  myws.factory(Form("alpha1[%.4f]",  1.3211   ));
  myws.factory(Form("n2[%.4f]",      2        ));
  myws.factory(Form("alpha2[%.4f]",  2.0936   ));
  myws.factory(Form("lambdas[%.4f]", -0.702181));*/
  /*bkg VWG sig CB2 data tails
  myws.factory(Form("mJPsi[%.4f]",  3.11241   ));
  myws.factory(Form("sigma[%.4f]",  0.0619954 ));
  myws.factory(Form("n1[%.4f]",     99.957    ));
  myws.factory(Form("alpha1[%.4f]", 0.541307  ));
  myws.factory(Form("n2[%.4f]",     59.9901   ));
  myws.factory(Form("alpha2[%.4f]", 0.843528  );
  myws.factory(Form("A[%.4f]",      1.72926   );
  myws.factory(Form("B[%.4f]",      0.00673999));
  myws.factory(Form("C[%.4f]",      0.273343  ));*/
  /*bkg VWG sig CB2 MC tails
  myws.factory(Form("mJPsi[%.4f]",  ));
  myws.factory(Form("sigma[%.4f]",  ));
  myws.factory(Form("n1[%.4f]",     ));
  myws.factory(Form("alpha1[%.4f]", ));
  myws.factory(Form("n2[%.4f]",     ));
  myws.factory(Form("alpha2[%.4f]", ));
  myws.factory(Form("A[%.4f]",      ));
  myws.factory(Form("B[%.4f]",      ));
  myws.factory(Form("C[%.4f]",      ));*/

  //if sig is CB2
  //myws.factory("RooDoubleCB::m_signal(mass, mJPsi, sigma, alpha1, n1, alpha2, n2)");

  //if sig is NA60
  myws.factory("RooNA60::m_signal(mass, mJPsi, sigma, alpha1, p1, p2, p3, alpha2, p12, p22, p32)");

  //For Exp bkg
  //myws.factory("RooExponential::m_bkg(mass, lambdas, true)");
  //myws.factory("RooBumpExp::m_bkg(mass, lambdas, A, B, C, cste)");

  //For VWG background
  //myws.factory("RooVWG::m_bkg(mass, A, B, C)");
  myws.factory("RooBumpVWG::m_bkg(mass, A, B, C, cste, A1, B1, C1, cste1)");

  RooAddPdf massModel("massModel", "massModel",
                      {*myws.pdf("m_bkg"), *myws.pdf("m_signal")},
                      {*myws.var("nbkg"), *myws.var("nsig")});
  //massModel.setNormRange("Mass");

  //For MC tail fix parmeters
  if(fixtails) {
     myws.var("n1")->setConstant(kTRUE);
     myws.var("alpha1")->setConstant(kTRUE);
     myws.var("n2")->setConstant(kTRUE);
     myws.var("alpha2")->setConstant(kTRUE);
  }

  RooFitResult *fitResult(0x0);
  fitResult = massModel.fitTo(*data, PrintLevel(-1), Save(), Range(massMin, massMax)); //, Range(massMin, massMax), Save()); // voir Range()... pour non prompt
  //fitResult = massModel.fitTo(*data, RooFit::Extended, PrintLevel(-1), Save(), Range(massMin, massMax)); //, Range(massMin, massMax), Save()); // voir Range()... pour non prompt
  fitResult->Print("v");
  myws.import(massModel);

  myws.var("nsig")->setConstant(kFALSE);
  myws.var("nbkg")->setConstant(kFALSE);
}

//____________________________________
void DoSPlot(RooWorkspace &myws) {
  std::cout << "Calculate sWeights" << std::endl;

  RooAbsPdf *massModel = myws.pdf("massModel");
  RooRealVar *nsig = myws.var("nsig");
  RooRealVar *nbkg = myws.var("nbkg");
  RooDataSet &data = static_cast<RooDataSet &>(*myws.data("data"));

  RooMsgService::instance().setSilentMode(true);

  std::cout << "\n\n------------------------------------------\nThe dataset "
               "before creating sWeights:\n";
  data.Print();

  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

  myws.pdf("massModel")
      ->getParameters(RooArgSet(*myws.var("mass")))
      ->setAttribAll("Constant", kTRUE);

  massModel->Print("t");

  RooStats::SPlot sData{"sData", "An SPlot", data, massModel,
                        RooArgList(*nsig, *nbkg)};

  std::cout << "\n\nThe dataset after creating sWeights:\n";
  data.Print();

  // Check that our weights have the desired properties

  std::cout
      << "\n\n------------------------------------------\n\nCheck SWeights:"
      << std::endl;

  std::cout << std::endl
            << "Yield of signal is\t" << nsig->getVal()
            << ".  From sWeights it is " << sData.GetYieldFromSWeight("nsig")
            << std::endl;

  std::cout << "Yield of bkg is\t" << nbkg->getVal()
            << ".  From sWeights it is " << sData.GetYieldFromSWeight("nbkg")
            << std::endl
            << std::endl;

  for (Int_t i = 0; i < 10; i++) {
    std::cout << "sig Weight for event " << i << std::right << std::setw(12)
              << sData.GetSWeight(i, "nsig") << "  bkg Weight" << std::setw(12)
              << sData.GetSWeight(i, "nbkg") << "  Total Weight"
              << std::setw(12) << sData.GetSumOfEventSWeight(i) << std::endl;
  }

  std::cout << std::endl;

  // import this new dataset with sWeights
  std::cout << "import new dataset with sWeights" << std::endl;
  myws.import(data, Rename("dataWithSWeights"));

  RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);
}

//____________________________________
void MakePlots(RooWorkspace &myws) {

  TCanvas *cdata = new TCanvas("sPlot", "sPlot demo", 1000, 800);
  cdata->Divide(1, 2);

  //For sig
  RooAbsPdf *massModel = myws.pdf("massModel");
  RooAbsPdf *m_signal = myws.pdf("m_signal");
  RooAbsPdf *m_bkg = myws.pdf("m_bkg");
  RooRealVar *tau = myws.var("tau");
  RooRealVar *mass = myws.var("mass");
  RooRealVar *sigma_signal = myws.var("sigma");
  RooRealVar *mJPsi = myws.var("mJPsi");

  ////For CB2 sig
  //RooRealVar *n1 = myws.var("n1");
  //RooRealVar *n2 = myws.var("n2");
  //RooRealVar *alpha1 = myws.var("alpha1");
  //RooRealVar *alpha2 = myws.var("alpha2");
  //RooRealVar *mJPsi = myws.var("mJPsi");

  //For NA60 sig
  RooRealVar *alpha1 = myws.var("alpha1");
  RooRealVar *p1 = myws.var("p1");
  RooRealVar *p2 = myws.var("p2");
  RooRealVar *p3 = myws.var("p3");
  RooRealVar *alpha2 = myws.var("alpha2");
  RooRealVar *p12 = myws.var("p12");
  RooRealVar *p22 = myws.var("p22");
  RooRealVar *p32 = myws.var("p32");

  //For Exp background
  //RooRealVar *lambda = myws.var("lambdas");

  //For VWG background
  RooRealVar *A = myws.var("A");
  RooRealVar *B = myws.var("B");
  RooRealVar *C = myws.var("C");
  RooRealVar *cste = myws.var("cste");

  RooRealVar *nsig = myws.var("nsig");
  RooRealVar *nbkg = myws.var("nbkg");

  RooRealVar *pT1 = myws.var("pT1");
  RooRealVar *pT2 = myws.var("pT2");

  std::cout << "mean : " << mJPsi->getValV() << std::endl;
  std::cout << "sigma : " << sigma_signal->getValV() << std::endl;

  ////if sig is CB2
  //std::cout << "n1 : " << n1->getValV() << std::endl;
  //std::cout << "alpha1 : " << alpha1->getValV() << std::endl;
  //std::cout << "n2 : " << n2->getValV() << std::endl;
  //std::cout << "alpha2 : " << alpha2->getValV() << std::endl;

  //if sig is NA60
  std::cout << "alpha1 : " << alpha1->getValV() << std::endl;
  std::cout << "p1 : " << p1->getValV() << std::endl;
  std::cout << "p2 : " << p2->getValV() << std::endl;
  std::cout << "p3 : " << p3->getValV() << std::endl;
  std::cout << "alpha2 : " << alpha2->getValV() << std::endl;
  std::cout << "p12 : " << p12->getValV() << std::endl;
  std::cout << "p22 : " << p22->getValV() << std::endl;
  std::cout << "p32 : " << p32->getValV() << std::endl;

  //if bkg is Exp
  //std::cout << "lambda : " << lambda->getValV() << std::endl;

  //if bkg is VWG
  std::cout << "A : " << A->getValV() << std::endl;
  std::cout << "B : " << B->getValV() << std::endl;
  std::cout << "C : " << C->getValV() << std::endl;
  std::cout << "cste : " << cste->getValV() << std::endl;

  int nBins = min(int(round((3.5 - 2.5) / 0.025)), 1000); // 40
  //int nBins = 10;

  std::cout << "NB BINS TEST1 : " << nBins << std::endl;

  auto &data = static_cast<RooDataSet &>(*myws.data("dataWithSWeights"));

  Double_t numEntries = data.sumEntries();
  //Double_t numEntries = data.numEntries();

  double massMin = 2.5;
  double massMax = 3.5;

  std::cout << "make weighted data sets" << std::endl;

  RooDataSet dataw_sig{data.GetName(), data.GetTitle(), &data,
                       *data.get(),    nullptr,         "nsig_sw"};
  RooDataSet dataw_bkg{data.GetName(), data.GetTitle(), &data,
                       *data.get(),    nullptr,         "nbkg_sw"};

  //RooPlot *frame = mass->frame(Title("Fit of model to discriminating variable"), Bins(nBins));

  RooPlot *frame = mass->frame(Bins(nBins));
  frame->SetTitle(" ");

  data.plotOn(frame);
  massModel->plotOn(
      frame, Components(*m_signal), LineStyle(kDashed), LineColor(kRed+1),
      Normalization(numEntries, RooAbsReal::NumEvent), Range(massMin, massMax),
      NormRange("Mass"), Name("m_signal"));
  massModel->plotOn(frame, Components(*m_bkg), LineStyle(kDashed),
                    LineColor(kGreen+1),
                    Normalization(numEntries, RooAbsReal::NumEvent), //RooAbsReal::isConstant
                    Range(massMin, massMax), NormRange("Mass"), Name("m_bkg"));
  massModel->plotOn(frame, Name("massModel"),
                    Normalization(numEntries, RooAbsReal::NumEvent),
                    Range(massMin, massMax), NormRange("Mass"));

  RooPlot *frameTMP = (RooPlot *)frame->Clone("TMP");
  RooHist *hpull = frameTMP->pullHist(nullptr, "massModel", true);
  //hpull->SetName("hpull");
  //RooPlot *frame4 = mass->frame(Title("Pull Distribution (#sigma)"));
  //gStyle->SetTitleFontSize(0.06);
  RooPlot *frame4 = mass->frame();
  frame4->SetTitle(" ");

  TLine* l_middle = new TLine(frame4->GetXaxis()->GetXmin(), 0, frame4->GetXaxis()->GetXmax(), 0);
  l_middle->SetLineColor(kBlack);
  l_middle->SetLineWidth(1);
  TLine* l_up = new TLine(frame4->GetXaxis()->GetXmin(), 2, frame4->GetXaxis()->GetXmax(), 2);
  l_up->SetLineColor(kBlack);
  l_up->SetLineWidth(1);
  l_up->SetLineStyle(kDashed);
  TLine* l_down = new TLine(frame4->GetXaxis()->GetXmin(), -2, frame4->GetXaxis()->GetXmax(), -2);
  l_down->SetLineColor(kBlack);
  l_down->SetLineWidth(1);
  l_down->SetLineStyle(kDashed);

  frame4->GetYaxis()->SetTitle("#frac{|hist - curve|}{#delta(hist)}");
  frame4->GetXaxis()->SetTitleSize(0.08);
  frame4->GetXaxis()->SetLabelSize(0.08);
  frame4->GetXaxis()->SetTitleOffset(0.95);
  frame4->GetYaxis()->SetTitleSize(0.08);
  frame4->GetYaxis()->SetLabelSize(0.08);
  frame4->GetYaxis()->SetTitleOffset(0.35);
  frame4->SetNdivisions(6, "Y");
  frame4->addPlotable(hpull, "PX");
  TCanvas *cFig = new TCanvas("cMassFig_PP", "cMassFig", 1000, 800);
  TPad *pad1 = new TPad("pad1_PP", "", 0, 0.24, 1, 1);
  TPad *pad2 = new TPad("pad2_PP", "", 0, 0, 1, .238);

  float titlesize = 0.035;
  float titleoffset = 0.5;
  float labelsize = 0.035;

  //pad2->SetBottomMargin(0.15);
  pad2->SetBottomMargin(0.2);

  frame->GetXaxis()->SetTitleSize(titlesize);
  //frame->GetXaxis()->SetTitleOffset(titleoffset);
  frame->GetXaxis()->SetLabelSize(labelsize);
  frame->GetYaxis()->SetTitleSize(titlesize);
  //frame->GetYaxis()->SetTitleOffset(titleoffset);
  frame->GetYaxis()->SetLabelSize(labelsize);

  cFig->cd();
  pad1->Draw();
  pad1->cd();
  frame->Draw();

  TLatex *t00 = new TLatex();
  t00->SetNDC();
  t00->SetTextSize(0.05);
  t00->DrawLatex(0.259519, 0.926995, "Fit of model to discriminating variable");

  TLatex *t0 = new TLatex();
  t0->SetNDC();
  t0->SetTextSize(0.035);
  t0->DrawLatex(0.410822, 0.00679117, "Pull Distribution (#sigma)");
  //t0->DrawLatex(0.410822, 0.959117, "Pull Distribution (#sigma)");

  TLatex *t = new TLatex();
  t->SetNDC();
  //t->SetTextSize(0.032);
  t->SetTextSize(0.04);
  float dy = 0;

  t->SetTextSize(0.03);
  t->DrawLatex(0.7175, 0.69 - dy, Form("%.1f < |y^{#mu#mu}| < %.1f", 2.5, 3.6));
  dy += 0.045;
  t->DrawLatex(0.6763, 0.69-dy, Form("%.0f GeV/c < p_{T}^{#mu#mu} < %.0f GeV/c", pT1->getValV(), pT2->getValV()));
  // t->DrawLatex(0.5175, 0.69-dy, Form("%g < p_{T}^{#mu#mu} < %g
  // GeV/c",cut.dMuon.Pt.Min,cut.dMuon.Pt.Max)); dy+=0.045;
  t->DrawLatex(0.2, 0.89 - dy, Form("Fit parameters : "));
  dy += 0.045;
  t->DrawLatex(0.2, 0.89 - dy, Form("#sigma = %.4f#pm%.4f", sigma_signal->getValV(), sigma_signal->getError()));
  dy += 0.045;

  ////if sig is CB2
  //t->DrawLatex(0.2, 0.89 - dy, Form("n_{1} = %.4f#pm%.4f", n1->getValV(), n1->getError()));
  //dy += 0.045;
  //t->DrawLatex(0.2, 0.89 - dy, Form("#alpha_{1} = %.4f#pm%.4f", alpha1->getValV(), alpha1->getError()));
  //dy += 0.045;
  //t->DrawLatex(0.2, 0.89 - dy, Form("n_{2} = %.4f#pm%.4f", n2->getValV(), n2->getError()));
  //dy += 0.045;
  //t->DrawLatex(0.2, 0.89 - dy, Form("#alpha_{2} = %.4f#pm%.4f", alpha2->getValV(), alpha2->getError()));
  //dy += 0.045;

  ////if sig is NA60
  t->DrawLatex(0.2, 0.89 - dy, Form("#alpha_{1} = %.4f#pm%.4f", alpha1->getValV(), alpha1->getError()));
  dy += 0.045;
  t->DrawLatex(0.2, 0.89 - dy, Form("p_{1} = %.4f#pm%.4f", p1->getValV(), p1->getError()));
  dy += 0.045;
  t->DrawLatex(0.2, 0.89 - dy, Form("p_{2} = %.4f#pm%.4f", p2->getValV(), p2->getError()));
  dy += 0.045;
  t->DrawLatex(0.2, 0.89 - dy, Form("p_{3} = %.4f#pm%.4f", p3->getValV(), p3->getError()));
  dy += 0.045;
  t->DrawLatex(0.2, 0.89 - dy, Form("#alpha_{2} = %.4f#pm%.4f", alpha2->getValV(), alpha2->getError()));
  dy += 0.045;
  t->DrawLatex(0.2, 0.89 - dy, Form("p_{12} = %.4f#pm%.4f", p12->getValV(), p12->getError()));
  dy += 0.045;
  t->DrawLatex(0.2, 0.89 - dy, Form("p_{22} = %.4f#pm%.4f", p22->getValV(), p22->getError()));
  dy += 0.045;
  t->DrawLatex(0.2, 0.89 - dy, Form("p_{32} = %.4f#pm%.4f", p32->getValV(), p32->getError()));
  dy += 0.045;

  ////if bkg is Exp
  //t->DrawLatex(0.2, 0.89 - dy, Form("#lambda = %.4f#pm%.4f", lambda->getValV(), lambda->getError()));
  //dy += 0.045;

  //if bkg is VWG
  t->DrawLatex(0.2, 0.89 - dy, Form("#A = %.4f#pm%.4f", A->getValV(), A->getError()));
  dy += 0.045;
  t->DrawLatex(0.2, 0.89 - dy, Form("#B = %.4f#pm%.4f", B->getValV(), B->getError()));
  dy += 0.045;
  t->DrawLatex(0.2, 0.89 - dy, Form("#C = %.4f#pm%.4f", C->getValV(), C->getError()));
  dy += 0.045;
  t->DrawLatex(0.2, 0.89 - dy, Form("#cste = %.4f#pm%.4f", cste->getValV(), cste->getError()));
  dy += 0.045;

  t->DrawLatex(0.2, 0.89 - dy, Form("n_{J/#Psi} = %.0f#pm%.0f", nsig->getValV(), nsig->getError()));
  dy += 0.045;
  t->DrawLatex(0.2, 0.89 - dy, Form("n_{bkg} = %.0f#pm%.0f", nbkg->getValV(), nbkg->getError()));
  dy += 0.045;
  t->DrawLatex(0.2, 0.89 - dy, Form("m_{J/#Psi} = %.4f#pm%.4f", mJPsi->getValV(), mJPsi->getError()));
  dy += 0.045;

  double ymin = 0.7802;
  TLegend *leg = new TLegend(0.6275, ymin, 0.7880, 0.8909);
  leg->SetTextSize(0.035);
  leg->SetLineColor(0);
  TLegendEntry *l1 = leg->AddEntry("m_signal", "Signal", "l");
  l1->SetLineColor(kRed+1);
  TLegendEntry *l2 = leg->AddEntry("m_bkg", "Background", "l");
  l2->SetLineColor(kGreen+1);
  TLegendEntry *l3 = leg->AddEntry("massModel", "Total fit", "l");
  l3->SetLineColor(4);
  leg->Draw("same");

  pad1->Update();
  cFig->cd();

  pad2->Draw();
  pad2->cd();
  frame4->Draw("E");

  l_middle->Draw("same");
  l_up->Draw("same");
  l_down->Draw("same");

  map<string, double> binWidth;
  binWidth["MASS"] = 0.025;

  double chi2 = 0;
  unsigned int ndof = 0;
  pad2->cd();
  TLatex *t1 = new TLatex();
  t1->SetNDC();
  t1->SetTextSize(0.08);
  unsigned int nFitPar = massModel->getParameters(data)
                             ->selectByAttrib("Constant", kFALSE)
                             ->getSize();
  TH1 *hdatact = data.createHistogram("hdatact", *mass,
                                      Binning(frame->GetNbinsX(),
                                              frame->GetXaxis()->GetXmin(),
                                              frame->GetXaxis()->GetXmax()));
  RooHist *hpull1 = frame->pullHist(0, 0, true);
  double *ypulls = hpull1->GetY();
  unsigned int nFullBins = 0;
  for (int i = 0; i < nBins; i++) {
    if (hdatact->GetBinContent(i + 1) > 1.0) {
      chi2 += ypulls[i] * ypulls[i];
      nFullBins++;
    }
  }
  ndof = nFullBins - nFitPar;

  //t1->DrawLatex(0.13, 0.8, Form("#chi^{2}/ndof = %.0f / %d ", chi2, ndof));
  pad2->Update();
  pad1->cd();
  t->DrawLatex(0.7263, 0.6, Form("#chi^{2}/ndof = %.0f / %d ", chi2, ndof));
  pad1->Update();

  cdata->cd(1);

  RooPlot *frame2 = tau->frame(Title("Lz distribution with s weights to project out signal"), Bins(100));
  dataw_sig.plotOn(frame2, DataError(RooAbsData::SumW2), Normalization(numEntries, RooAbsReal::NumEvent));
  frame2->GetXaxis()->SetTitleSize(0.045);
  //frame->GetXaxis()->SetTitleOffset(titleoffset);
  frame2->GetXaxis()->SetLabelSize(0.05);
  frame2->GetYaxis()->SetTitleSize(0.045);
  //frame->GetYaxis()->SetTitleOffset(titleoffset);
  frame2->GetYaxis()->SetLabelSize(0.05);
  frame2->Draw();

  cdata->cd(2);
  RooPlot *frame3 = tau->frame(Title("Lz distribution with s weights to project out bkg"), Bins(100));
  // dataw_bkg.plotOn(frame3, DataError(RooAbsData::SumW2));
  dataw_bkg.plotOn(frame3, DataError(RooAbsData::SumW2));
  frame3->GetXaxis()->SetTitleSize(0.045);
  //frame->GetXaxis()->SetTitleOffset(titleoffset);
  frame3->GetXaxis()->SetLabelSize(0.05);
  frame3->GetYaxis()->SetTitleSize(0.045);
  //frame->GetYaxis()->SetTitleOffset(titleoffset);
  frame3->GetYaxis()->SetLabelSize(0.05);
  frame3->Draw();

  myws.import(dataw_sig, Rename("dataw_sig"));
  myws.import(dataw_bkg, Rename("dataw_bkg"));
}

void TemplateFit(RooWorkspace &myws) {
    
    RooRealVar *tau = myws.var("tau");
    //RooRealVar ctau("ctau", "ctau", -10, 10);
    //myws.import(ctau);

    auto &dataSet = static_cast<RooDataSet &>(*myws.data("dataw_bkg"));
    myws.import(dataSet, Rename("dataSet"));
    
    /*RooDataSet dataSet("dataSet", "dataSet", RooArgSet(*tau));

    for (int i = 0; i < 10000; i++) {
        //ctau.setVal(-log(gRandom->Rndm()) * 2.0);  // Exponential decay
        tau->setVal(gRandom->Gaus(0., 1.)); 
        if (tau->getVal() > -10 && tau->getVal() < 10)
            dataSet.add(RooArgSet(*tau));
    }
    
    myws.import(dataSet);*/

    double ctauMin = -10;
    double ctauMax = 10;
    int nBins = 200;

    TH1D* hTot = (TH1D*)myws.data("dataSet")->createHistogram("CTAU_Bkg", *myws.var("tau"), Binning(nBins, ctauMin, ctauMax));
    TCanvas c1("c1", "Canvas", 800, 600);
    hTot->Draw();
    c1.SaveAs("template_histo.png");

    hTot->ClearUnderflowAndOverflow();
    for (int i=0; i<=hTot->GetNbinsX(); i++) { if (hTot->GetBinContent(i)<0) { hTot->SetBinContent(i, 0.0000000001); } }
    TH1* hClean = rebinctauBkghist(myws, hTot, ctauMin, ctauMax);
    string dataName = "dh_bkg";
    RooDataHist* dataHist = new RooDataHist(dataName.c_str(), "", *myws.var("tau"), hClean);

    if (dataHist==NULL) { cout << "[ERROR] DataHist used failed!" << endl;}
    if (dataHist->sumEntries()==0) { cout << "[ERROR] DataHist used is empty!" << endl; }
    if (fabs(dataHist->sumEntries()-hClean->GetSumOfWeights())>0.001) { cout << "[ERROR] DataHist used is invalid!  " << endl;}
    myws.import(*dataHist);
    
//    RooKeysPdf pdf("pdf", "pdf", *ws.var("ctau"), *((RooDataSet*)ws.data(dataName.c_str())), RooKeysPdf::MirrorAsymBoth, 2);
    RooKeysPdf pdf("pdf", "pdf", *myws.var("tau"), *((RooDataSet*)myws.data(dataName.c_str())), RooKeysPdf::NoMirror, 1.5);
    myws.import(pdf);
    
    TCanvas c("c", "Canvas", 800, 600);
    RooPlot* frame = tau->frame();
    dataSet.plotOn(frame);
    pdf.plotOn(frame);
    frame->Draw();
    c.SaveAs("template_plot.png");
}

void AddTauzSignalModel(RooWorkspace &myws) {
  //bool full_model = false;
  RooRealVar *tau = myws.var("tau");

  auto &data_sig_cut = static_cast<RooDataSet &>(*myws.data("dataw_sig"));
  //auto &data_sig_cut = static_cast<RooDataSet &>(*myws.data("dataw_sig")->reduce(Cut("tau < 0.")));

  double numEntries = data_sig_cut.sumEntries();

  //double tau_max = 0.;
  double tau_max = 10.;
  double tau_min = -10.; //-11

  tau->setRange("Tau", tau_min, tau_max);

  /*myws.factory(Form("tzJPsi[%.6f,%.6f,%.6f]", 0., -1.e-2, 1.e-2));
  myws.factory(Form("sigma1[%.6f,%.6f,%.6f]", 2.e-2, 1.e-2, 10.)); //10 0.5
  myws.factory(Form("sigma2[%.6f,%.6f,%.6f]", 6.e-2, 1.e-2, 10.));
  myws.factory(Form("lambdat[%.6f,%.6f,%.6f]", 1.e-1, 1.e-2, 10.));
  myws.factory(Form("c[%.6f,%.6f,%.6f]", 0., 0., 300.));
  myws.factory(Form("d[%.6f,%.6f,%.6f]", 0., 0., 0.));
  myws.factory(Form("f2res[%.4f,%.4f,%.4f]", 0.5, 0., 1.)); // 0.3
  myws.factory(Form("fres[%.4f,%.4f,%.4f]", 0.3, 0., 1.));*/ // 0.52

  /*myws.factory(Form("tzJPsi[%.6f]", 0.));
  myws.factory(Form("sigma1[%.6f,%.6f,%.6f]", 2.e-2, 1.e-2, 0.5)); //10 0.5
  myws.factory(Form("sigma2[%.6f,%.6f,%.6f]", 6.e-2, 1.e-2, 0.5));
  myws.factory(Form("n11[%.6f,%.6f,%.6f]", 12., 0., 80.));
  myws.factory(Form("alpha11[%.6f,%.6f,%.6f]", 1., 0.1, 5.));
  myws.factory(Form("fres[%.4f,%.4f,%.4f]", 0.3, 0., 1.));*/

  /*myws.factory(Form("tzJPsi[%.6f]", 0.));
  myws.factory(Form("sigma1[%.6f,%.6f,%.6f]", 2.e-2, 1.e-2, 9.)); //10 0.5
  myws.factory(Form("sigma2[%.6f,%.6f,%.6f]", 6.e-2, 1.e-2, 9.));
  myws.factory(Form("n11[%.6f,%.6f,%.6f]", 12., 0., 80.));
  myws.factory(Form("alpha11[%.6f,%.6f,%.6f]", 1., 0.1, 10.));
  myws.factory(Form("fres[%.4f,%.4f,%.4f]", 0.3, 0., 1.));*/

  //0-1 values DCA Chi2 30
  //myws.factory(Form("tzJPsi[%.6f]", 0.));
  /*myws.factory(Form("sigma1[%.6f]", 0.225724)); //10 0.5
  myws.factory(Form("sigma2[%.6f]", 0.098771));
  myws.factory(Form("n11[%.6f]", 1.34093));
  myws.factory(Form("alpha11[%.6f]", 1.828183));
  myws.factory(Form("fres[%.4f]", 0.259245));
  myws.factory(Form("tzJPsi[%.6f]", 0.005333));*/

  //1-2 values DCA Chi2 30
  //myws.factory(Form("tzJPsi[%.6f]", 0.));
  /*myws.factory(Form("sigma1[%.6f]", 0.161464)); //10 0.5
  myws.factory(Form("sigma2[%.6f]", 0.078631));
  myws.factory(Form("n11[%.6f]", 1.521372));
  myws.factory(Form("alpha11[%.6f]", 1.590860));
  myws.factory(Form("fres[%.4f]", 0.327360));
  myws.factory(Form("tzJPsi[%.6f]", 0.001962));*/

  //2-3 values DCA Chi2 35
  //myws.factory(Form("tzJPsi[%.6f]", 0.));
  /*myws.factory(Form("sigma1[%.6f]", 0.143049)); //10 0.5
  myws.factory(Form("sigma2[%.6f]", 0.074224));
  myws.factory(Form("n11[%.6f]", 1.480944));
  myws.factory(Form("alpha11[%.6f]", 1.667820));
  myws.factory(Form("fres[%.4f]", 0.340513));
  myws.factory(Form("tzJPsi[%.6f]", 0.001552));*/

  //3-4 values DCA Chi2 40
  //myws.factory(Form("tzJPsi[%.6f]", 0.));
  /*myws.factory(Form("sigma1[%.6f]", 0.114288)); //10 0.5
  myws.factory(Form("sigma2[%.6f]", 0.060975));
  myws.factory(Form("n11[%.6f]", 1.411540));
  myws.factory(Form("alpha11[%.6f]", 1.667450));
  myws.factory(Form("fres[%.4f]", 0.432404));
  myws.factory(Form("tzJPsi[%.6f]", 0.003240));*/

  //4-5 values DCA Chi2 45
  //myws.factory(Form("tzJPsi[%.6f]", 0.));
  /*myws.factory(Form("sigma1[%.6f]", 0.100990)); //10 0.5
  myws.factory(Form("sigma2[%.6f]", 0.065271));
  myws.factory(Form("n11[%.6f]", 1.521497));
  myws.factory(Form("alpha11[%.6f]", 1.710038));
  myws.factory(Form("fres[%.4f]", 0.346746));
  myws.factory(Form("tzJPsi[%.6f]", 0.003599));*/

  //5-6 values DCA Chi2 50
  //myws.factory(Form("tzJPsi[%.6f]", 0.));
  /*myws.factory(Form("sigma1[%.6f]", 0.085472)); //10 0.5
  myws.factory(Form("sigma2[%.6f]", 0.051492));
  myws.factory(Form("n11[%.6f]", 1.488315));
  myws.factory(Form("alpha11[%.6f]", 1.619788));
  myws.factory(Form("fres[%.4f]", 0.482244));
  myws.factory(Form("tzJPsi[%.6f]", 0.002315));*/

  //6-8 values DCA Chi2 60
  //myws.factory(Form("tzJPsi[%.6f]", 0.));
  /*myws.factory(Form("sigma1[%.6f]", 0.075992)); //10 0.5
  myws.factory(Form("sigma2[%.6f]", 0.045563));
  myws.factory(Form("n11[%.6f]", 1.368173));
  myws.factory(Form("alpha11[%.6f]", 1.846924));
  myws.factory(Form("fres[%.4f]", 0.428932));
  myws.factory(Form("tzJPsi[%.6f]", 0.000330));*/

  //8-10 values DCA Chi2 65
  //myws.factory(Form("tzJPsi[%.6f]", 0.));
  /*myws.factory(Form("sigma1[%.6f]", 0.066366)); //10 0.5
  myws.factory(Form("sigma2[%.6f]", 0.037309));
  myws.factory(Form("n11[%.6f]", 1.612392));
  myws.factory(Form("alpha11[%.6f]", 1.666926));
  myws.factory(Form("fres[%.4f]", 0.447506));
  myws.factory(Form("tzJPsi[%.6f]", 0.001158));*/

  //10-30 values DCA Chi2 90
  //myws.factory(Form("tzJPsi[%.6f]", 0.));
  /*myws.factory(Form("sigma1[%.6f]", 0.073610)); //10 0.5
  myws.factory(Form("sigma2[%.6f]", 0.035784));
  myws.factory(Form("n11[%.6f]", 1.316378));
  myws.factory(Form("alpha11[%.6f]", 2.028249));
  myws.factory(Form("fres[%.4f]", 0.265550));
  myws.factory(Form("tzJPsi[%.6f]", 0.000715));*/

  //10-15 values DCA Chi2 90 tzJPsi free
  //myws.factory(Form("tzJPsi[%.6f]", 0.));
  /*myws.factory(Form("sigma1[%.6f]", 0.078267)); //10 0.5
  myws.factory(Form("sigma2[%.6f]", 0.037413));
  myws.factory(Form("n11[%.6f]", 1.235676));
  myws.factory(Form("alpha11[%.6f]", 2.166185));
  myws.factory(Form("fres[%.4f]", 0.264027));
  myws.factory(Form("tzJPsi[%.6f]", -0.000068));*/

  //15-30 values DCA Chi2 90 tzJPsi free
  //myws.factory(Form("tzJPsi[%.6f]", 0.));
  /*myws.factory(Form("sigma1[%.6f]", 0.056386)); //10 0.5
  myws.factory(Form("sigma2[%.6f]", 0.025935));
  myws.factory(Form("n11[%.6f]", 1.468908));
  myws.factory(Form("alpha11[%.6f]", 1.724188));
  myws.factory(Form("fres[%.4f]", 0.258328));
  myws.factory(Form("tzJPsi[%.6f]", 0.003632));*/

  //myws.factory(Form("tzJPsi[%.6f,%.6f,%.6f]", 0., -0.05, 0.05));

  //myws.factory("RooTruthModel::delta_exp(tau)");
  // myws.factory("RooGaussModel::delta_exp(tau, 0., 4.e-4)");
  //myws.factory("RooGaussian::gauss1(tau, tzJPsi, sigma1)");
  //myws.factory("RooGaussian::gauss2(tau, tzJPsi, sigma2)");
  //myws.factory("RooDecay::double_exp(tau, lambdat, delta_exp, RooDecay::DoubleSided)");
  //myws.factory("RooExp_Affine::double_exp(tau, lambdat, c, d)");

  myws.factory("RooGaussModel::gauss1(tau, tzJPsi, sigma1)");
  myws.factory("RooGaussModel::gauss2(tau, tzJPsi, sigma2)");
  myws.factory("RooDoubleCB::DCB(tau, tzJPsi, sigma2, alpha11, n11, alpha11, n11)");

  /*RooAddPdf model("model", "",
                  RooArgList(*myws.pdf("gauss1"), *myws.pdf("gauss2"),
                             *myws.pdf("double_exp")),
                  RooArgList(*myws.var("fres"), *myws.var("f2res")));

  model.setNormRange("Tau");
  myws.import(model);
  
  model.fitTo(data_sig_cut); //, Range(tau_min, tau_max)); //, SumW2Error(data_sig_cut.isWeighted())); // Range(-2.e-3, -1.e-8)
  model.Print("t");
  myws.import(data_sig_cut, Rename("data_sig_cut"));
  myws.import(model, RecycleConflictNodes());*/

  RooAddModel resolution("resolution", "", RooArgList(*myws.pdf("gauss1"), *myws.pdf("gauss2")), *myws.var("fres"));
  myws.import(resolution);

    /*RooAddPdf model("model", "",
                  RooArgList(*myws.pdf("gauss1"), *myws.pdf("DCB")),
                  RooArgList(*myws.var("fres")));

    model.setNormRange("Tau");
    myws.import(model);
    //myws.factory("RooExtendPdf::model(model_noext,nsig)");

    model.fitTo(data_sig_cut, SumW2Error(kTRUE)); //, Range(tau_min, tau_max)); //, SumW2Error(data_sig_cut.isWeighted())); // Range(-2.e-3, -1.e-8)
    model.Print("t");
    myws.import(data_sig_cut, Rename("data_sig_cut"));
    myws.import(model, RecycleConflictNodes());*/


  /*RooAddPdf tausigModel("tausigModel", "", RooArgList(*myws.pdf("gauss1"), *myws.pdf("gauss2"), *myws.pdf("double_exp")), RooArgList(*myws.var("fres"), *myws.var("f2res")));
  myws.import(tausigModel);
  myws.factory(Form("lambda_np[%.6f,%.6f,%.6f]", 0.6, 0.3, 0.9));
  myws.factory(Form("mean_np[%.6f]", 0.));
  myws.factory(Form("sigma_np[%.6f]", 0.1));

  //myws.var("lambda_np")->setConstant(kTRUE);

  myws.factory(Form("fB[%.4f,%.4f,%.4f]", 0.1, 0., 1.));

  myws.factory("RooGaussModel::delta_np(tau, mean_np, sigma_np)");
  myws.factory("Decay::exp_np(tau, lambda_np, delta_np, RooDecay::SingleSided)");

  RooAddPdf model("model", "", *myws.pdf("exp_np"), *myws.pdf("tausigModel"), *myws.var("fB"));
  model.setNormRange("Tau");
  myws.import(model);*/

    RooAddPdf tausigModel("tausigModel", "", RooArgList(*myws.pdf("gauss1"), *myws.pdf("DCB")), RooArgList(*myws.var("fres")));
    myws.import(tausigModel);
    myws.factory(Form("lambda_np[%.6f,%.6f,%.6f]", 0.7, 0.3, 0.9));
    //myws.factory(Form("mean_np[%.6f]", 0.));
    //myws.factory(Form("sigma_np[%.6f]", 0.1));

    myws.factory("RooTruthModel::delta_expbis(tau)");

    myws.factory(Form("fB[%.4f,%.4f,%.4f]", 0.1, 0., 1.));

    //myws.factory("RooGaussModel::delta_np(tau, mean_np, sigma_np)");

    myws.factory("Decay::exp_np(tau, lambda_np, resolution, RooDecay::SingleSided)");

    //myws.factory("Decay::exp_npbis(tau, lambda_np, delta_expbis, RooDecay::SingleSided)");
    //myws.factory("RooFFTConvPdf::exp_np(tau, resolution, exp_npbis)");

    RooAddPdf model("model", "", *myws.pdf("exp_np"), *myws.pdf("tausigModel"), *myws.var("fB"));
    model.setNormRange("Tau");
    myws.import(model);

    model.fitTo(data_sig_cut, SumW2Error(kTRUE)); //, Range(tau_min, tau_max)); //, SumW2Error(data_sig_cut.isWeighted())); // Range(-2.e-3, -1.e-8)
    model.Print("t");
    myws.import(data_sig_cut, Rename("data_sig_cut"));
    myws.import(model, RecycleConflictNodes());


    /*RooAddPdf model_noext("model_noext", "", *myws.pdf("exp_np"), *myws.pdf("tausigModel"), *myws.var("fB"));
    model_noext.setNormRange("Tau");
    myws.import(model_noext);

    myws.factory("RooExtendPdf::model(model_noext,nsig)");
    auto *model = myws.pdf("model");

    model->fitTo(data_sig_cut, SumW2Error(kTRUE), Extended(kTRUE)); //, Range(tau_min, tau_max)); //, SumW2Error(data_sig_cut.isWeighted())); // Range(-2.e-3, -1.e-8)
    model->Print("t");
    myws.import(data_sig_cut, Rename("data_sig_cut"));
    myws.import(*model, RecycleConflictNodes());*/
}

//=====================================================//
void AddTauzBkgModel(RooWorkspace &myws) {
  RooRealVar *tau = myws.var("tau");
  //auto &data_bkg = static_cast<RooDataSet &>(*myws.data("dataw_bkg"));
  auto &data_bkg = static_cast<RooDataSet &>(*myws.data("dataw_bkg")->reduce(Cut("tau > -10 && tau < 10")));

  double tau_max = 10.;
  double tau_min = -10.;

  tau->setRange("Tau_bis", tau_min, tau_max);

  /*myws.factory(Form("b[%.4f,%.4f,%.4f]", 0.8, 0.5, 0.9));
  myws.factory(Form("f_DLIV[%.4f,%.4f,%.4f]", 0.4, 0.3, 1.)); //0.9, 0.7, 0.9
  myws.factory(Form("f_DFSS[%.4f,%.4f,%.4f]", 0.7, 0.3, 1.)); //0.6, 0.4, 0.6
  myws.factory(Form("lambda_DSS[%.4f,%.4f,%.4f]", 0.3, 0.1, 0.9)); //Red
  myws.factory(Form("lambda_DF[%.4f,%.4f,%.4f]", 0.2, 0., 0.9)); //Green 2.5 6.e-4
  myws.factory(Form("lambda_DDS[%.4f,%.4f,%.4f]", 1.5, 0., 3.));
  myws.factory(Form("tzJPsi_bkg[%.6f,%.6f,%.6f]", 0., -1., 1.));*/

  //0-1 values DCA Chi2 30 tzJPsi free NEW
  /*myws.factory(Form("b[%.4f]", 0.677929));
  myws.factory(Form("f_DLIV[%.4f]", 0.708500)); 
  myws.factory(Form("f_DFSS[%.4f]", 0.594487)); 
  myws.factory(Form("lambda_DSS[%.4f]", 0.376400));
  myws.factory(Form("lambda_DF[%.4f]", 0.226335));
  myws.factory(Form("lambda_DDS[%.4f]", 2.450728));
  myws.factory(Form("tzJPsi_bkg[%.4f]", 0.034425));*/

  //1-2 values DCA Chi2 30 tzJPsi free NEW
  /*myws.factory(Form("b[%.4f]", 0.732702));
  myws.factory(Form("f_DLIV[%.4f]", 0.741752)); 
  myws.factory(Form("f_DFSS[%.4f]", 0.610291)); 
  myws.factory(Form("lambda_DSS[%.4f]", 0.347962));
  myws.factory(Form("lambda_DF[%.4f]", 0.202082));
  myws.factory(Form("lambda_DDS[%.4f]", 2.271579));
  myws.factory(Form("tzJPsi_bkg[%.4f]", 0.031188));*/

  //2-3 values DCA Chi2 35 tzJPsi free NEW
  /*myws.factory(Form("b[%.4f]", 0.717970));
  myws.factory(Form("f_DLIV[%.4f]", 0.771084)); 
  myws.factory(Form("f_DFSS[%.4f]", 0.681752)); 
  myws.factory(Form("lambda_DSS[%.4f]", 0.332718));
  myws.factory(Form("lambda_DF[%.4f]", 0.199907));
  myws.factory(Form("lambda_DDS[%.4f]", 2.111054));
  myws.factory(Form("tzJPsi_bkg[%.4f]", 0.021992));*/

  //3-4 values DCA Chi2 40 tzJPsi free NEW
  /*myws.factory(Form("b[%.4f]", 0.825606));
  myws.factory(Form("f_DLIV[%.4f]", 0.797110)); 
  myws.factory(Form("f_DFSS[%.4f]", 0.605539)); 
  myws.factory(Form("lambda_DSS[%.4f]", 0.316187));
  myws.factory(Form("lambda_DF[%.4f]", 0.141081));
  myws.factory(Form("lambda_DDS[%.4f]", 1.806764));
  myws.factory(Form("tzJPsi_bkg[%.4f]", 0.044440));*/

  //4-5 values DCA Chi2 45 tzJPsi free NEW
  /*myws.factory(Form("b[%.4f]", 0.818906));
  myws.factory(Form("f_DLIV[%.4f]", 0.822345)); 
  myws.factory(Form("f_DFSS[%.4f]", 0.668166)); 
  myws.factory(Form("lambda_DSS[%.4f]", 0.325406));
  myws.factory(Form("lambda_DF[%.4f]", 0.148014));
  myws.factory(Form("lambda_DDS[%.4f]", 1.689273));
  myws.factory(Form("tzJPsi_bkg[%.4f]", 0.035436));*/

  //5-6 values DCA Chi2 50 tzJPsi free NEW
  /*myws.factory(Form("b[%.4f]", 0.796209));
  myws.factory(Form("f_DLIV[%.4f]", 0.841410)); 
  myws.factory(Form("f_DFSS[%.4f]", 0.688501)); 
  myws.factory(Form("lambda_DSS[%.4f]", 0.366528));
  myws.factory(Form("lambda_DF[%.4f]", 0.137910));
  myws.factory(Form("lambda_DDS[%.4f]", 1.630785));
  myws.factory(Form("tzJPsi_bkg[%.4f]", 0.026687));*/

  //6-8 values DCA Chi2 60 tzJPsi free NEW
  /*myws.factory(Form("b[%.4f]", 0.99991));
  myws.factory(Form("f_DLIV[%.4f]", 0.849761)); 
  myws.factory(Form("f_DFSS[%.4f]", 0.466887)); 
  myws.factory(Form("lambda_DSS[%.4f]", 0.350077));
  myws.factory(Form("lambda_DF[%.4f]", 0.094008));
  myws.factory(Form("lambda_DDS[%.4f]", 1.322568));
  myws.factory(Form("tzJPsi_bkg[%.4f]", 0.077613));*/

  //6-8 values DCA Chi2 60 tzJPsi free NEW NEW
  /*myws.factory(Form("b[%.4f]", 0.85));
  myws.factory(Form("f_DLIV[%.4f]", 0.836706)); 
  myws.factory(Form("f_DFSS[%.4f]", 0.646222)); 
  myws.factory(Form("lambda_DSS[%.4f]", 0.350981));
  myws.factory(Form("lambda_DF[%.4f]", 0.110253));
  myws.factory(Form("lambda_DDS[%.4f]", 1.387303));
  myws.factory(Form("tzJPsi_bkg[%.4f]", 0.037230));*/

  //8-10 values DCA Chi2 65 tzJPsi free NEW
  /*myws.factory(Form("b[%.4f]", 0.896247));
  myws.factory(Form("f_DLIV[%.4f]", 0.833252)); 
  myws.factory(Form("f_DFSS[%.4f]", 0.766203)); 
  myws.factory(Form("lambda_DSS[%.4f]", 0.350221));
  myws.factory(Form("lambda_DF[%.4f]", 0.083046));
  myws.factory(Form("lambda_DDS[%.4f]", 1.236831));
  myws.factory(Form("tzJPsi_bkg[%.4f]", 0.009110));*/

  //10-30 values DCA Chi2 90 tzJPsi free NEW
  /*myws.factory(Form("b[%.4f]", 0.964557));
  myws.factory(Form("f_DLIV[%.4f]", 0.856408)); 
  myws.factory(Form("f_DFSS[%.4f]", 0.734517)); 
  myws.factory(Form("lambda_DSS[%.4f]", 0.379131));
  myws.factory(Form("lambda_DF[%.4f]", 0.089167));
  myws.factory(Form("lambda_DDS[%.4f]", 0.970180));
  myws.factory(Form("tzJPsi_bkg[%.4f]", 0.032484));*/

  //10-15 values DCA Chi2 90 tzJPsi free NEW
  /*myws.factory(Form("b[%.4f]", 0.85));
  myws.factory(Form("f_DLIV[%.4f]", 0.856529)); 
  myws.factory(Form("f_DFSS[%.4f]", 0.802960)); 
  myws.factory(Form("lambda_DSS[%.4f]", 0.401445));
  myws.factory(Form("lambda_DF[%.4f]", 0.131094));
  myws.factory(Form("lambda_DDS[%.4f]", 1.036481));
  myws.factory(Form("tzJPsi_bkg[%.4f]", 0.020180));*/

  //15-30 values DCA Chi2 90 tzJPsi free NEW
  /*myws.factory(Form("b[%.4f]", 0.84999));
  myws.factory(Form("f_DLIV[%.4f]", 0.896680)); 
  myws.factory(Form("f_DFSS[%.4f]", 0.897764)); 
  myws.factory(Form("lambda_DSS[%.4f]", 0.405195));
  myws.factory(Form("lambda_DF[%.4f]", 0.129420));
  myws.factory(Form("lambda_DDS[%.4f]", 0.982057));
  myws.factory(Form("tzJPsi_bkg[%.4f]", -0.0187884));*/

  //myws.factory("RooGaussModel::delta(tau, 0., 0.07)"); // 0.07
  //myws.factory("RooGaussModel::delta_test(tau, 0.01, 0.1)"); //0.01, 0.16 for integrated -20 20 //0., 0.07 for integrated -10 10 //0.005, 0.1 for 0-1 //0.01, 0.1 for 1-2
  //myws.factory("RooGaussModel::delta_test(tau, 0., 0.15)");

  /*myws.factory(Form("b[%.4f,%.4f,%.4f]", 0.7, 0., 1.));
  myws.factory(Form("f_DLIV[%.4f,%.4f,%.4f]", 0.4, 0., 1.)); //0.9, 0.7, 0.9
  myws.factory(Form("f_DFSS[%.4f,%.4f,%.4f]", 0.6, 0., 1.)); //0.6, 0.4, 0.6
  myws.factory(Form("lambda_DSS[%.4f,%.4f,%.4f]", 2.5e-1, 2.e-1, 1.)); //Red
  myws.factory(Form("lambda_DF[%.4f,%.4f,%.4f]", 2.e-2, 1.e-1, 1.));
  myws.factory(Form("n12[%.6f,%.6f,%.6f]", 12., 0., 80.));
  myws.factory(Form("alpha12[%.6f,%.6f,%.6f]", 1., 0.1, 5.));
  myws.factory(Form("n22[%.6f,%.6f,%.6f]", 12., 0., 80.));
  myws.factory(Form("alpha22[%.6f,%.6f,%.6f]", 1., 0.1, 5.));
  myws.factory(Form("mid[%.4f,%.4f,%.4f]", 0., -0.5, 0.5));*/

  //myws.factory("RooDoubleCB::DCB_bkg(tau, mid, sigma2, alpha12, n12, alpha12, n12)");

  myws.factory("RooTruthModel::delta(tau)");

  /*myws.factory("Decay::exp_DSS(tau, lambda_DSS, delta, RooDecay::SingleSided)");
  myws.factory("Decay::exp_DF(tau, lambda_DF, delta, RooDecay::Flipped)");
  myws.factory("Decay::exp_DDS(tau, lambda_DDS, delta, RooDecay::DoubleSided)");*/

  //myws.factory("RooTruthModel::delta_bis(tau)");
  //myws.factory("RooFFTConvPdf::conv(tau, delta_bis, delta)");

  /*myws.factory("RooFFTConvPdf::conv1(tau, resolution, exp_DSS)");
  myws.factory("RooFFTConvPdf::conv2(tau, resolution, exp_DF)");
  myws.factory("RooFFTConvPdf::conv3(tau, resolution, exp_DDS)");*/

  /*myws.factory("RooFFTConvPdf::conv1(tau, model, exp_DSS)");
  myws.factory("RooFFTConvPdf::conv2(tau, model, exp_DF)");
  myws.factory("RooFFTConvPdf::conv3(tau, model, exp_DDS)");*/

  /*myws.factory("SUM::model1(f_DFSS*exp_DSS, exp_DF)");
  myws.factory("SUM::model3(f_DLIV*model1, conv3)");
  myws.factory("SUM::model2(b*model3, resolution)");*/

  /*RooAddPdf model1("model1", "", RooArgList(*myws.pdf("conv1"), *myws.pdf("conv2")), *myws.var("f_DFSS"));
  myws.import(model1);
  RooAddPdf model3("model3", "", RooArgList(*myws.pdf("model1"), *myws.pdf("conv3")), *myws.var("f_DLIV"));
  myws.import(model3);
  RooAddPdf model2("model2", "", RooArgList(*myws.pdf("model3"), *myws.pdf("resolution")), *myws.var("b"));
  myws.import(model2);*/

  myws.factory("RooGaussModel::gauss1_bkg(tau, tzJPsi_bkg, sigma1)");
  myws.factory("RooGaussModel::gauss2_bkg(tau, tzJPsi_bkg, sigma2)");

  RooAddModel resolution_bkg("resolution_bkg", "", RooArgList(*myws.pdf("gauss1_bkg"), *myws.pdf("gauss2_bkg")), *myws.var("fres"));
  myws.import(resolution_bkg);

  myws.factory("Decay::exp_DSS(tau, lambda_DSS, resolution_bkg, RooDecay::SingleSided)");
  myws.factory("Decay::exp_DF(tau, lambda_DF, resolution_bkg, RooDecay::Flipped)");
  myws.factory("Decay::exp_DDS(tau, lambda_DDS, resolution_bkg, RooDecay::DoubleSided)");

  myws.factory("SUM::model1(f_DFSS*exp_DSS, exp_DF)");
  myws.factory("SUM::model3(f_DLIV*model1, exp_DDS)");
  myws.factory("SUM::model2(b*model3, resolution_bkg)");

  //myws.factory("SUM::model2_noext(b*model3, resolution_bkg)");
  //myws.factory("RooExtendPdf::model2(model2_noext,nbkg)");

  myws.pdf("model2")->setNormRange("Tau_bis");

  myws.pdf("model2")->fitTo(data_bkg, Range(tau_min, tau_max));
  //myws.pdf("model2")->fitTo(data_bkg, Range(tau_min, tau_max), Extended(kTRUE)); //, Extended(kTRUE)
  myws.pdf("model2")->Print("t");
  myws.import(data_bkg, Rename("data_bkg"));
}

//=====================================================//
void MakeSignalPlot(RooWorkspace &myws) {
  RooRealVar *tau = myws.var("tau");
  RooRealVar *tzJPsi = myws.var("tzJPsi");
  RooRealVar *sigma1 = myws.var("sigma1");
  RooRealVar *sigma2 = myws.var("sigma2");
  RooRealVar *n11 = myws.var("n11");
  RooRealVar *alpha11 = myws.var("alpha11");
  RooRealVar *lambdat = myws.var("lambdat");
  RooRealVar *c = myws.var("c");
  RooRealVar *d = myws.var("d");
  RooRealVar *fres = myws.var("fres");
  RooRealVar *f2res = myws.var("f2res");
  RooAbsPdf *gauss1 = myws.pdf("gauss1");
  RooAbsPdf *gauss2 = myws.pdf("gauss2");
  RooAbsPdf *double_exp = myws.pdf("double_exp");
  RooAbsPdf *model = myws.pdf("model");
  RooAbsPdf *DCB = myws.pdf("DCB");
  RooRealVar *fB = myws.var("fB");

  RooRealVar *pT1 = myws.var("pT1");
  RooRealVar *pT2 = myws.var("pT2");

  std::cout << "sigma1 : " << sigma1->getValV() << std::endl;
  std::cout << "sigma2 : " << sigma2->getValV() << std::endl;
  std::cout << "fres : " << fres->getValV() << std::endl;
  std::cout << "alpha1 : " << alpha11->getValV() << std::endl;
  std::cout << "n1 : " << n11->getValV() << std::endl;

  double tau_max = 0.;
  //double tau_max = 10.;
  double tau_min = -10.;

  auto &data_sig = static_cast<RooDataSet &>(*myws.data("data_sig_cut"));

  //Double_t numEntries = data_sig.sumEntries();
  Double_t numEntries = data_sig.sumEntries();

  //int nBins = min(int(round((3.5 - 2.5) / 0.025)), 1000); //mettre même valeur que plot
  int nBins = 80;

  //RooPlot *frame = tau->frame(Title("Fit of l_{z} signal"), Bins(100));
  RooPlot *frame = tau->frame(Bins(nBins));
  frame->SetTitle(" ");

  data_sig.plotOn(frame); //, Normalization(numEntries, RooAbsReal::NumEvent), Range(tau_min, tau_max), NormRange("Tau"));
  model->plotOn(frame, Components(*gauss1), LineStyle(kDashed), LineColor(kRed+1),
                Name("gauss1")); //, Normalization(numEntries, RooAbsReal::NumEvent), Range(tau_min, tau_max), NormRange("Tau"));
  model->plotOn(frame, Components(*DCB), LineStyle(kDashed), LineColor(kGreen+1), Name("DCB"));
  //model->plotOn(frame, Components(*gauss2), LineStyle(kDashed), LineColor(kGreen+1), Name("gauss2")); //, Normalization(numEntries, RooAbsReal::NumEvent), Range(tau_min, tau_max), NormRange("Tau"));
  //model->plotOn(frame, Components(*double_exp), LineStyle(kDashed), LineColor(kBlack), Name("double_exp")); //, Normalization(numEntries, RooAbsReal::NumEvent), Range(tau_min, tau_max), NormRange("Tau"));
  model->plotOn(frame, Name("model")); //, Normalization(numEntries, RooAbsReal::NumEvent),Range(tau_min, tau_max), NormRange("Tau"));

  RooPlot *frameTMP = (RooPlot *)frame->Clone("TMP");
  RooHist *hpull = frameTMP->pullHist(0, "model", true);
  hpull->SetName("hpull");
  //RooPlot *frame2 = tau->frame(Title("Pull Distribution (#sigma)"));
  RooPlot *frame2 = tau->frame();
  frame2->SetTitle(" ");
  frame2->SetTitleSize(0.07);

  TLine* l_middle = new TLine(frame2->GetXaxis()->GetXmin(), 0, frame2->GetXaxis()->GetXmax(), 0);
  l_middle->SetLineColor(kBlack);
  l_middle->SetLineWidth(1);
  TLine* l_up = new TLine(frame2->GetXaxis()->GetXmin(), 2, frame2->GetXaxis()->GetXmax(), 2);
  l_up->SetLineColor(kBlack);
  l_up->SetLineWidth(1);
  l_up->SetLineStyle(kDashed);
  TLine* l_down = new TLine(frame2->GetXaxis()->GetXmin(), -2, frame2->GetXaxis()->GetXmax(), -2);
  l_down->SetLineColor(kBlack);
  l_down->SetLineWidth(1);
  l_down->SetLineStyle(kDashed);

  frame2->GetYaxis()->SetTitle("#frac{|hist - curve|}{#delta(hist)}");
  frame2->GetXaxis()->SetTitleSize(0.08);
  frame2->GetXaxis()->SetLabelSize(0.08);
  frame2->GetXaxis()->SetTitleOffset(0.95);
  frame2->GetYaxis()->SetTitleSize(0.08);
  frame2->GetYaxis()->SetLabelSize(0.08);
  frame2->GetYaxis()->SetTitleOffset(0.35);
  //frame2->GetYaxis()->SetTitleSize(0.07);
  //frame2->GetXaxis()->SetLabelSize(0.07);
  //frame2->GetYaxis()->SetLabelSize(0.07);
  frame2->SetNdivisions(6, "Y");
  frame2->addPlotable(hpull, "PX");

  TCanvas *cFig = new TCanvas("cTauzFig_PP", "cTauzFig", 1000, 800);
  TPad *pad1 = new TPad("pad1_TauzPP", "", 0, 0.23, 1, 1);
  TPad *pad2 = new TPad("pad2_TauzPP", "", 0, 0, 1, .228);

  float titlesize = 0.035;
  float titleoffset = 0.5;
  float labelsize = 0.035;

  //pad2->SetBottomMargin(0.15);
  pad2->SetBottomMargin(0.2);

  frame->GetXaxis()->SetTitleSize(titlesize);
  //frame->GetXaxis()->SetTitleOffset(titleoffset);
  frame->GetXaxis()->SetLabelSize(labelsize);
  frame->GetYaxis()->SetTitleSize(titlesize);
  //frame->GetYaxis()->SetTitleOffset(titleoffset);
  frame->GetYaxis()->SetLabelSize(labelsize);

  // void setMassRange3
  TH1 *h = data_sig.createHistogram("hist", *tau,
                                    Binning(frame->GetNbinsX(),
                                            frame->GetXaxis()->GetXmin(),
                                            frame->GetXaxis()->GetXmax()));
  Double_t YMax = h->GetBinContent(h->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i = 1; i <= h->GetNbinsX(); i++)
    if (h->GetBinContent(i) > 0)
      YMin = min(YMin, h->GetBinContent(i));
  Double_t Yup(0.), Ydown(0.);
  Ydown = max(YMin - (YMax - YMin) * (0.1 / (1.0 - 0.1 - 0.4)), 0.0);
  Yup = YMax + (YMax - YMin) * (0.4 / (1.0 - 0.1 - 0.4));
  frame->GetYaxis()->SetRangeUser(Ydown, Yup);
  delete h;

  cFig->cd();
  pad1->Draw();
  pad1->cd();
  frame->Draw();

  TLatex *t00 = new TLatex();
  t00->SetNDC();
  t00->SetTextSize(0.05);
  t00->DrawLatex(0.410822, 0.926995, "Fit of l_{z} signal");

  TLatex *t0 = new TLatex();
  t0->SetNDC();
  t0->SetTextSize(0.035);
  t0->DrawLatex(0.410822, 0.00679117, "Pull Distribution (#sigma)");

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextSize(0.032);
  float dy = 0;

  t->SetTextSize(0.03);
  t->DrawLatex(0.7175, 0.69 - dy, Form("%.1f < |y^{#mu#mu}| < %.1f", 2.5, 3.6));
  dy += 0.045;
  t->DrawLatex(0.6763, 0.69-dy, Form("%.0f GeV/c < p_{T}^{#mu#mu} < %.0f GeV/c", pT1->getValV(), pT2->getValV()));
  // t->DrawLatex(0.5175, 0.69-dy, Form("%g < p_{T}^{#mu#mu} < %g
  // GeV/c",cut.dMuon.Pt.Min,cut.dMuon.Pt.Max)); dy+=0.045;
  t->DrawLatex(0.15, 0.9 - dy, Form("Fit parameters : "));
  dy += 0.045;
  t->DrawLatex(
      0.15, 0.9 - dy,
      Form("#sigma_{1} = %.6f#pm%.6f", sigma1->getValV(), sigma1->getError()));
  dy += 0.045;
  t->DrawLatex(
      0.15, 0.9 - dy,
      Form("#sigma_{2} = %.6f#pm%.6f", sigma2->getValV(), sigma2->getError()));
  dy += 0.045;
  //t->DrawLatex(0.15, 0.9 - dy, Form("#lambda = %.6f#pm%.6f", lambdat->getValV(), lambdat->getError()));
  //dy += 0.045;
  //t->DrawLatex(0.15, 0.9 - dy, Form("c = %.6f#pm%.6f", c->getValV(), c->getError()));
  //dy += 0.045;
  //t->DrawLatex(0.15, 0.9 - dy, Form("d = %.6f#pm%.6f", d->getValV(), d->getError()));
  //dy += 0.045;
  t->DrawLatex(0.15, 0.9 - dy, Form("f_{res} = %.6f#pm%.6f", fres->getValV(), fres->getError()));
  dy += 0.045;
  t->DrawLatex(0.15, 0.9 - dy, Form("n_{1} = %.6f#pm%.6f", n11->getValV(), n11->getError()));
  dy += 0.045;
  t->DrawLatex(0.15, 0.9 - dy, Form("#alpha_{1} = %.6f#pm%.6f", alpha11->getValV(), alpha11->getError()));
  dy += 0.045;
  //t->DrawLatex(0.15, 0.9 - dy, Form("fB = %.6f#pm%.6f", myws.var("fB")->getValV(), myws.var("fB")->getError()));
  //dy += 0.045;
  //t->DrawLatex(0.15, 0.9 - dy, Form("#lambda_{np} = %.6f#pm%.6f", myws.var("lambda_np")->getValV(), myws.var("lambda_np")->getError()));
  //dy += 0.045;
  //t->DrawLatex(0.15, 0.9 - dy, Form("f_{2res} = %.6f#pm%.6f", f2res->getValV(), f2res->getError()));
  //dy += 0.045;
  t->DrawLatex(0.15, 0.9 - dy, Form("#tau_{z}^{J/#Psi} = %.6f#pm%.6f", tzJPsi->getValV(), tzJPsi->getError()));
  dy += 0.045;
  //t->DrawLatex(0.15, 0.9 - dy, Form("f_{B} = %.6f#pm%.6f", fB->getValV(), fB->getError()));
  //dy += 0.045;

  double ymin = 0.7802;
  TLegend *leg = new TLegend(0.707415, 0.423544, 0.877756, 0.552576);
  //leg->SetTextSize(0.03);
  leg->SetTextSize(0.035);
  leg->SetLineColor(0);
  TLegendEntry *l1 = leg->AddEntry("gauss1", "gauss1", "l");
  l1->SetLineColor(kRed+1);
  /*TLegendEntry *l2 = leg->AddEntry("gauss2", "gauss2", "l");
  l2->SetLineColor(kGreen+1);
  TLegendEntry* l3 = leg->AddEntry("double_exp","exp","l");
  l3->SetLineColor(1);*/
  TLegendEntry *l2 = leg->AddEntry("DCB", "DCB", "l");
  l2->SetLineColor(kGreen+1);
  TLegendEntry *l4 = leg->AddEntry("model", "Total fit", "l");
  l4->SetLineColor(4);
  leg->Draw("same");

  // pad1->Update();
  cFig->cd();

  pad2->Draw();
  pad2->cd();
  frame2->Draw();

  l_middle->Draw("same");
  l_up->Draw("same");
  l_down->Draw("same");

  map<string, double> binWidth;
  binWidth["MASS"] = 0.025;

  double chi2 = 0;
  unsigned int ndof = 0;
  pad2->cd();
  TLatex *t1 = new TLatex();
  t1->SetNDC();
  t1->SetTextSize(0.08);
  unsigned int nFitPar = model->getParameters(data_sig)
                             ->selectByAttrib("Constant", kFALSE)
                             ->getSize();
  TH1 *hdatact = data_sig.createHistogram(
      "hdatact", *tau,
      Binning(frame->GetNbinsX(), frame->GetXaxis()->GetXmin(),
              frame->GetXaxis()->GetXmax()));
  RooHist *hpull1 = frame->pullHist(0, 0, true);
  double *ypulls = hpull1->GetY();
  unsigned int nFullBins = 0;
  //int nBins = 110;
  for (int i = 0; i < nBins; i++) {
    if (hdatact->GetBinContent(i + 1) > 1.0) {
      chi2 += ypulls[i] * ypulls[i];
      nFullBins++;
    }
  }
  ndof = nFullBins - nFitPar;
  //t1->DrawLatex(0.13, 0.8, Form("#chi^{2}/ndof = %.0f / %d ", chi2, ndof));
  pad2->Update();

  pad1->cd();
  t->DrawLatex(0.7263, 0.6, Form("#chi^{2}/ndof = %.0f / %d ", chi2, ndof));
  pad1->Update();
}

//=====================================================//
void MakeBkgPlot(RooWorkspace &myws) {
  RooRealVar *tau = myws.var("tau");
  RooRealVar *b = myws.var("b");
  RooRealVar *tzJPsi_bkg = myws.var("tzJPsi_bkg");
  RooRealVar *f_DLIV = myws.var("f_DLIV");
  RooRealVar *f_DFSS = myws.var("f_DFSS");
  RooRealVar *lambda_DSS = myws.var("lambda_DSS");
  RooRealVar *lambda_DF = myws.var("lambda_DF");
  RooRealVar *lambda_DDS = myws.var("lambda_DDS");
  RooAbsPdf *delta = myws.pdf("delta");
  RooAbsPdf *exp_DSS = myws.pdf("exp_DSS");
  RooAbsPdf *exp_DF = myws.pdf("exp_DF");
  RooAbsPdf *exp_DDS = myws.pdf("exp_DDS");
  RooAbsPdf *model1 = myws.pdf("model1");
  RooAbsPdf *model2 = myws.pdf("model2");
  RooAbsPdf *resolution = myws.pdf("resolution");
  RooAbsPdf *resolution_bkg = myws.pdf("resolution_bkg");
  RooAbsPdf *DCB = myws.pdf("DCB");
  RooAbsPdf *conv3 = myws.pdf("conv3");
  RooAbsPdf *conv2 = myws.pdf("conv2");

  RooRealVar *pT1 = myws.var("pT1");
  RooRealVar *pT2 = myws.var("pT2");

  std::cout << "f_DLIV : " << f_DLIV->getValV() << std::endl;
  std::cout << "b : " << b->getValV() << std::endl;
  std::cout << "f_DFSS : " << f_DFSS->getValV() << std::endl;
  std::cout << "lambda_DSS : " << lambda_DSS->getValV() << std::endl;
  std::cout << "lambda_DF : " << lambda_DF->getValV() << std::endl;
  std::cout << "lambda_DDS : " << lambda_DDS->getValV() << std::endl;
  std::cout << "tzJPsi : " << tzJPsi_bkg->getValV() << std::endl;

  auto &data_bkg = static_cast<RooDataSet &>(*myws.data("dataw_bkg"));

  double tau_max = 10.;
  double tau_min = -10.;

  Double_t numEntries = data_bkg.sumEntries();

  int nBins = 80;
  //int nBins = min(int(round((3.5 - 2.5) / 0.025)), 1000);

  //RooPlot *frame = tau->frame(Title("Fit of l_{z} bkg"), Bins(100));
  RooPlot *frame = tau->frame(Bins(nBins));
  frame->SetTitle(" ");

  data_bkg.plotOn(frame);
  //model1->plotOn(frame, Components(*exp_DSS), LineStyle(kDashed), LineColor(kRed+1), Name("exp_DSS"), Normalization(0.72)); // 0.6
  //model1->plotOn(frame, Components(*exp_DF), LineStyle(kDashed), LineColor(kGreen+1), Name("exp_DF"), Normalization(0.28)); // Normalization(norm2, RooAbsReal::NumEvent)); 0.4
  //model1->plotOn(frame, Components(*conv2), LineStyle(kDashed), LineColor(kGreen+1), Name("conv2"), Normalization(0.28));
  model2->plotOn(frame, Components(*exp_DDS), LineStyle(kDashed),LineColor(kBlack), Name("exp_DDS"), Normalization(numEntries, RooAbsReal::NumEvent), Range(tau_min, tau_max), NormRange("Tau_bis"));
  model2->plotOn(frame, Components(*model1), LineStyle(kDashed),LineColor(kRed+1), Name("model1"), Normalization(numEntries, RooAbsReal::NumEvent), Range(tau_min, tau_max), NormRange("Tau_bis"));
  //model2->plotOn(frame, Components(*conv3), LineStyle(kDashed),LineColor(kBlack), Name("conv3"), Normalization(numEntries, RooAbsReal::NumEvent), Range(tau_min, tau_max), NormRange("Tau_bis"));
  //model2->plotOn(frame, Components(*delta), LineStyle(kDashed),LineColor(kCyan+1), Name("delta"), Normalization(numEntries, RooAbsReal::NumEvent), Range(tau_min, tau_max), NormRange("Tau_bis"));
  //model2->plotOn(frame, Components(*resolution), LineStyle(kDashed),LineColor(kCyan+1), Name("resolution"), Normalization(numEntries, RooAbsReal::NumEvent), Range(tau_min, tau_max), NormRange("Tau_bis"));
  model2->plotOn(frame, Components(*resolution_bkg), LineStyle(kDashed),LineColor(kCyan+1), Name("resolution_bkg"), Normalization(numEntries, RooAbsReal::NumEvent), Range(tau_min, tau_max), NormRange("Tau_bis"));
  //model2->plotOn(frame, Components(*DCB), LineStyle(kDashed),LineColor(kCyan+1), Name("DCB"), Normalization(numEntries, RooAbsReal::NumEvent), Range(tau_min, tau_max), NormRange("Tau_bis"));
  model2->plotOn(frame, Name("model2"), Normalization(numEntries, RooAbsReal::NumEvent), Range(tau_min, tau_max), NormRange("Tau_bis"));

  RooPlot *frameTMP = (RooPlot *)frame->Clone("TMP");
  RooHist *hpull = frameTMP->pullHist(0, "model2", true);
  hpull->SetName("hpull");
  //RooPlot *frame2 = tau->frame(Title("Pull Distribution (#sigma)"));
  RooPlot *frame2 = tau->frame();
  frame2->SetTitle(" ");

  frame2->SetTitleSize(0.07);
  
  TLine* l_middle = new TLine(frame2->GetXaxis()->GetXmin(), 0, frame2->GetXaxis()->GetXmax(), 0);
  l_middle->SetLineColor(kBlack);
  l_middle->SetLineWidth(1);
  TLine* l_up = new TLine(frame2->GetXaxis()->GetXmin(), 2, frame2->GetXaxis()->GetXmax(), 2);
  l_up->SetLineColor(kBlack);
  l_up->SetLineWidth(1);
  l_up->SetLineStyle(kDashed);
  TLine* l_down = new TLine(frame2->GetXaxis()->GetXmin(), -2, frame2->GetXaxis()->GetXmax(), -2);
  l_down->SetLineColor(kBlack);
  l_down->SetLineWidth(1);
  l_down->SetLineStyle(kDashed);

  frame2->GetYaxis()->SetTitle("#frac{|hist - curve|}{#delta(hist)}");
  frame2->GetXaxis()->SetTitleSize(0.08);
  frame2->GetXaxis()->SetLabelSize(0.08);
  frame2->GetXaxis()->SetTitleOffset(0.95);
  frame2->GetYaxis()->SetTitleSize(0.08);
  frame2->GetYaxis()->SetLabelSize(0.08);
  frame2->GetYaxis()->SetTitleOffset(0.35);
  //frame2->GetYaxis()->SetTitleSize(0.07);
  //frame2->GetXaxis()->SetLabelSize(0.07);
  //frame2->GetYaxis()->SetLabelSize(0.07);
  frame2->SetNdivisions(6, "Y");
  frame2->addPlotable(hpull, "PX");

  TCanvas *cFig = new TCanvas("cTauz_PP", "cTauz", 1000, 800);
  TPad *pad1 = new TPad("pad1_Tauz", "", 0, 0.23, 1, 1);
  TPad *pad2 = new TPad("pad2_Tauz", "", 0, 0, 1, .228);

  float titlesize = 0.035;
  float titleoffset = 0.5;
  float labelsize = 0.035;

  //pad2->SetBottomMargin(0.15);
  pad2->SetBottomMargin(0.2);

  frame->GetXaxis()->SetTitleSize(titlesize);
  //frame->GetXaxis()->SetTitleOffset(titleoffset);
  frame->GetXaxis()->SetLabelSize(labelsize);
  frame->GetYaxis()->SetTitleSize(titlesize);
  //frame->GetYaxis()->SetTitleOffset(titleoffset);
  frame->GetYaxis()->SetLabelSize(labelsize);

  cFig->cd();
  pad1->Draw();
  pad1->cd();
  frame->Draw();

  TLatex *t00 = new TLatex();
  t00->SetNDC();
  t00->SetTextSize(0.05);
  t00->DrawLatex(0.371743, 0.926995, "Fit of l_{z} background");

  TLatex *t0 = new TLatex();
  t0->SetNDC();
  t0->SetTextSize(0.035);
  t0->DrawLatex(0.410822, 0.00679117, "Pull Distribution (#sigma)");

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextSize(0.032);

  float dy = 0;

  t->SetTextSize(0.03);
  t->DrawLatex(0.7175, 0.69 - dy, Form("%.1f < |y^{#mu#mu}| < %.1f", 2.5, 3.6));
  dy += 0.045;
  t->DrawLatex(0.6763, 0.69-dy, Form("%.0f GeV/c < p_{T}^{#mu#mu} < %.0f GeV/c", pT1->getValV(), pT2->getValV()));
  // t->DrawLatex(0.5175, 0.69-dy, Form("%g < p_{T}^{#mu#mu} < %g
  // GeV/c",cut.dMuon.Pt.Min,cut.dMuon.Pt.Max)); dy+=0.045;
  t->DrawLatex(0.2, 0.9 - dy, Form("Fit parameters : "));
  dy += 0.045;
  t->DrawLatex(0.2, 0.9 - dy, Form("f_{DLIV} = %.6f#pm%.6f", f_DLIV->getValV(), f_DLIV->getError()));
  dy += 0.045;
  t->DrawLatex(0.2, 0.9 - dy, Form("b = %.6f#pm%.6f", b->getValV(), b->getError()));
  dy += 0.045;
  t->DrawLatex(
      0.2, 0.9 - dy,
      Form("f_{DFSS} = %.6f#pm%.6f", f_DFSS->getValV(), f_DFSS->getError()));
  dy += 0.045;
  t->DrawLatex(0.2, 0.9 - dy,
               Form("#lambda_{DSS} = %.6f#pm%.6f", lambda_DSS->getValV(),
                    lambda_DSS->getError()));
  dy += 0.045;
  t->DrawLatex(0.2, 0.9 - dy,
               Form("#lambda_{DF} = %.6f#pm%.6f", lambda_DF->getValV(),
                    lambda_DF->getError()));
  dy += 0.045;
  t->DrawLatex(0.2, 0.9 - dy, Form("#lambda_{DDS} = %.6f#pm%.6f", lambda_DDS->getValV(), lambda_DDS->getError()));
  dy += 0.045;
  t->DrawLatex(0.2, 0.9 - dy, Form("#tau_{z}^{J/#psi} = %.6f#pm%.6f", tzJPsi_bkg->getValV(), tzJPsi_bkg->getError()));
  dy += 0.045;

  TLegend *leg = new TLegend(0.707415, 0.333054, 0.877756, 0.534143);
  //leg->SetTextSize(0.03);
  leg->SetTextSize(0.035);
  leg->SetLineColor(0);
  TLegendEntry *l1 = leg->AddEntry("exp_DSS", "exp_{DSS} + exp_{DF}", "l");
  l1->SetLineColor(kRed+1);
  //TLegendEntry *l2 = leg->AddEntry("exp_DF", "exp_{DF}", "l");
  //l2->SetLineColor(kGreen+1);
  TLegendEntry *l3 = leg->AddEntry("exp_DDS", "exp_{DDS}", "l");
  l3->SetLineColor(1);
  TLegendEntry *l4 = leg->AddEntry("resolution_bkg", "delta", "l");
  l4->SetLineColor(kCyan+1);
  TLegendEntry *l5 = leg->AddEntry("model3", "Total fit", "l");
  l5->SetLineColor(4);
  leg->Draw("same");

  cFig->cd();

  pad2->Draw();
  pad2->cd();
  frame2->Draw();

  l_middle->Draw("same");
  l_up->Draw("same");
  l_down->Draw("same");

  map<string, double> binWidth;
  binWidth["MASS"] = 0.025;

  double chi2 = 0;
  unsigned int ndof = 0;
  pad2->cd();
  TLatex *t1 = new TLatex();
  t1->SetNDC();
  t1->SetTextSize(0.08);
  unsigned int nFitPar = model2->getParameters(data_bkg)
                             ->selectByAttrib("Constant", kFALSE)
                             ->getSize();
  TH1 *hdatact = data_bkg.createHistogram(
      "hdatact", *tau,
      Binning(frame->GetNbinsX(), frame->GetXaxis()->GetXmin(),
              frame->GetXaxis()->GetXmax()));
  RooHist *hpull1 = frame->pullHist(0, 0, true);
  double *ypulls = hpull1->GetY();
  unsigned int nFullBins = 0;
  std::cout << "NBINS : " << nBins << std::endl;
  for (int i = 0; i < nBins; i++) {
    if (hdatact->GetBinContent(i + 1) > 1.0) {
      //std::cout << "OOOOOOOOOH" << std::endl;
      //std::cout << "ypull i : " << ypulls[i] << std::endl;
      chi2 += ypulls[i] * ypulls[i];
      nFullBins++;
    }
  }
  ndof = nFullBins - nFitPar;
  //t1->DrawLatex(0.13, 0.8, Form("#chi^{2}/ndof = %.0f / %d ", chi2, ndof));
  pad2->Update();

  pad1->cd();
  t->DrawLatex(0.7263, 0.6, Form("#chi^{2}/ndof = %.0f / %d ", chi2, ndof));
  pad1->Update();
}


//=====================================================//
void MakeTotalPlot(RooWorkspace &myws){
  RooRealVar *tau = myws.var("tau");
  RooRealVar *fB = myws.var("fB");
  RooRealVar *lambda_np = myws.var("lambda_np");
  RooAbsPdf *model = myws.pdf("model");
  RooAbsPdf *exp_np = myws.pdf("exp_np");
  RooAbsPdf *tausigModel = myws.pdf("tausigModel");
  RooAbsPdf *model2 = myws.pdf("model2");

  double tau_max = 10.;
  double tau_min = -10.;

  tau->setRange("Lz", tau_min, tau_max);

  auto &data = static_cast<RooDataSet &>(*myws.data("data"));

  //myws.var("nsig")->setConstant(kTRUE);
  //myws.var("nbkg")->setConstant(kTRUE);

  //RooAddPdf tau_Model("tau_Model", "", {*myws.pdf("model2"), *myws.pdf("model")}, {*myws.var("nbkg"), *myws.var("nsig")});

  RooAddPdf tau_Model("tau_Model", "", RooArgList(*myws.pdf("model2"), *myws.pdf("model")), RooArgList(*myws.var("nbkg"), *myws.var("nsig")));

  myws.import(tau_Model);

  tau_Model.setNormRange("Lz");

  tau_Model.fitTo(data, Range(tau_min, tau_max)); 

  Double_t numEntries = data.sumEntries();

  RooPlot *frame = tau->frame(Title("Fit of l_{z} total"), Bins(100));
  data.plotOn(frame);
  model->plotOn(frame, Components(*exp_np), LineStyle(kDashed), LineColor(kGreen), Name("exp_np"), Normalization(numEntries, RooAbsReal::NumEvent), Range(tau_min, tau_max), NormRange("Lz"));
  model->plotOn(frame, Components(*tausigModel), LineStyle(kDashed), LineColor(kRed), Name("tausigModel"), Normalization(numEntries, RooAbsReal::NumEvent), Range(tau_min, tau_max), NormRange("Lz"));
  //tau_Model.plotOn(frame, Components(*myws.pdf("model")), Name("model"), LineStyle(kDashed),LineColor(kRed), Normalization(numEntries, RooAbsReal::NumEvent), Range(tau_min, tau_max), NormRange("Lz"));
  tau_Model.plotOn(frame, Components(*myws.pdf("model2")), Name("model2"), LineStyle(kDashed),LineColor(kBlack), Normalization(numEntries, RooAbsReal::NumEvent), Range(tau_min, tau_max), NormRange("Lz"));
  tau_Model.plotOn(frame, Name("tau_Model"), LineColor(kBlue), Normalization(numEntries, RooAbsReal::NumEvent), Range(tau_min, tau_max), NormRange("Lz"));

  RooPlot *frameTMP = (RooPlot *)frame->Clone("TMP");
  RooHist *hpull = frameTMP->pullHist(0, "tau_Model", true);
  hpull->SetName("hpull");
  RooPlot *frame2 = tau->frame(Title("Pull Distribution (#sigma)"));
  frame2->SetTitleSize(0.07);
  frame2->GetYaxis()->SetTitle("#frac{|hist - curve|}{#delta(hist)}");
  frame2->GetYaxis()->SetTitleSize(0.07);
  frame2->GetXaxis()->SetLabelSize(0.07);
  frame2->GetYaxis()->SetLabelSize(0.07);
  frame2->SetNdivisions(6, "Y");
  frame2->addPlotable(hpull, "PX");

  TCanvas *cFig = new TCanvas("cTauz_PP", "cTauz", 1000, 800);
  TPad *pad1 = new TPad("pad1_Tauz", "", 0, 0.23, 1, 1);
  TPad *pad2 = new TPad("pad2_Tauz", "", 0, 0, 1, .228);

  cFig->cd();
  pad1->Draw();
  pad1->cd();
  frame->Draw();

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextSize(0.032);

  float dy = 0;

  t->SetTextSize(0.03);
  t->DrawLatex(0.7175, 0.69 - dy, Form("%.1f < |y^{#mu#mu}| < %.1f", 2.5, 3.6));
  dy += 0.045;
  t->DrawLatex(0.2, 0.9 - dy, Form("Fit parameters : "));
  dy += 0.045;
  t->DrawLatex(0.2, 0.9 - dy, Form("#lambda_{np} = %.6f#pm%.6f", lambda_np->getValV(), lambda_np->getError()));
  dy += 0.045;
  t->DrawLatex(0.2, 0.9 - dy, Form("f_{B} = %.6f#pm%.6f", fB->getValV(), fB->getError()));
  dy += 0.045;
  t->DrawLatex(0.2, 0.9 - dy, Form("n_{sig} = %.0f#pm%.0f", myws.var("nsig")->getValV(), myws.var("nsig")->getError()));
  dy += 0.045;
  t->DrawLatex(0.2, 0.9 - dy, Form("n_{bkg} = %.0f#pm%.0f", myws.var("nbkg")->getValV(), myws.var("nbkg")->getError()));
  dy += 0.045;

  TLegend *leg = new TLegend(0.6175, 0.75, 0.7880, 0.8809);
  leg->SetTextSize(0.03);
  TLegendEntry *l1 = leg->AddEntry("tausigModel", "prompt", "l");
  l1->SetLineColor(2);
  TLegendEntry *l2 = leg->AddEntry("exp_np", "non-prompt", "l");
  l2->SetLineColor(3);
  TLegendEntry *l3 = leg->AddEntry("model2", "background", "l");
  l3->SetLineColor(1);
  TLegendEntry *l4 = leg->AddEntry("tau_Model", "Total fit", "l");
  l4->SetLineColor(4);
  leg->Draw("same");

  cFig->cd();

  pad2->Draw();
  pad2->cd();
  frame2->Draw();

  map<string, double> binWidth;
  binWidth["MASS"] = 0.025;

  double chi2 = 0;
  unsigned int ndof = 0;
  pad2->cd();
  TLatex *t1 = new TLatex();
  t1->SetNDC();
  t1->SetTextSize(0.08);
  unsigned int nFitPar = tau_Model.getParameters(data)
                             ->selectByAttrib("Constant", kFALSE)
                             ->getSize();
  TH1 *hdatact = data.createHistogram(
      "hdatact", *tau,
      Binning(frame->GetNbinsX(), frame->GetXaxis()->GetXmin(),
              frame->GetXaxis()->GetXmax()));
  RooHist *hpull1 = frame->pullHist(0, 0, true);
  double *ypulls = hpull1->GetY();
  unsigned int nFullBins = 0;
  int nBins = min(int(round((3.5 - 2.5) / binWidth["MASS"])), 1000);
  std::cout << "NBINS : " << nBins << std::endl;
  for (int i = 0; i < nBins; i++) {
    if (hdatact->GetBinContent(i + 1) > 1.0) {
      chi2 += ypulls[i] * ypulls[i];
      nFullBins++;
    }
  }
  ndof = nFullBins - nFitPar;
  t1->DrawLatex(0.13, 0.8, Form("#chi^{2}/ndof = %.0f / %d ", chi2, ndof));
  pad2->Update();
}

//=====================================================//
void MakeFit2D(RooWorkspace &myws) {

  RooRealVar *tau = myws.var("tau");
  RooRealVar *mass = myws.var("mass");

  RooRealVar *fB = myws.var("fB");
  RooRealVar *b = myws.var("b");

  RooRealVar *nsig = myws.var("nsig");
  RooRealVar *nbkg = myws.var("nbkg"); 

  RooAbsPdf *m_signal = myws.pdf("m_signal");
  RooAbsPdf *m_bkg = myws.pdf("m_bkg");
  RooAbsPdf *massModel = myws.pdf("massModel");

  RooAbsPdf *model = myws.pdf("model");
  RooAbsPdf *model2 = myws.pdf("model2");
  RooAbsPdf *model3 = myws.pdf("model3");

  RooAbsPdf *tau_Model = myws.pdf("tau_Model");
  RooAbsPdf *tausig_Model = myws.pdf("tausigModel");

  RooAbsPdf *Key_bkg = myws.pdf("pdf");

  RooRealVar *pT1 = myws.var("pT1");
  RooRealVar *pT2 = myws.var("pT2");

  RooDataSet &data = static_cast<RooDataSet &>(*myws.data("data")); //essayer avec weighted data 

  //auto &data_sig = static_cast<RooDataSet &>(*myws.data("dataw_sig"));

  Double_t numEntries = data.numEntries();

  //Double_t numEntriesbis = data_sig.numEntries();

  double tau_max = 10.;
  double tau_min = -10.;

  tau->setRange("Tauz", tau_min, tau_max);
  mass->setRange("Mass", 2.5, 3.5);

  myws.factory("PROD::Mass_tau_p(tausigModel, m_signal)");
  myws.factory("PROD::Mass_tau_np(exp_np, m_signal)");
  myws.factory("SUM::Mass_tau_signal(fB*Mass_tau_np, Mass_tau_p)");

  //myws.factory("PROD::Mass_tau_bkg_p(delta, m_bkg)");
  //myws.factory("PROD::Mass_tau_bkg_np(model3, m_bkg)");
  //myws.factory("SUM::Mass_tau_bkg(b*Mass_tau_bkg_np, Mass_tau_bkg_p)");

  myws.factory("PROD::Mass_tau_bkg(model2, m_bkg)");
  //myws.factory("PROD::Mass_tau_bkg(pdf, m_bkg)");

  RooAbsPdf *Fit_Total = new RooAddPdf("Fit_Total", "Fit_Total", RooArgList(*myws.pdf("Mass_tau_signal"), *myws.pdf("Mass_tau_bkg")), RooArgList(*myws.var("nsig"), *myws.var("nbkg")));

  Fit_Total->fitTo(data, Save(), PrintLevel(2)); //, Range("Mass", "Tauz"));

  myws.import(*Fit_Total);

  //==============================================================
  map<string, double> binWidth;
  binWidth["MASS"] = 0.025;
  //int nBins = min(int(round((3.5 - 2.5) / binWidth["MASS"])), 1000);
  int nBins = 80;

  /*RooPlot *frame = mass->frame(Title("Mass Projection"), Bins(nBins));
  data.plotOn(frame);
  Fit_Total->plotOn(frame);
  Fit_Total->plotOn(frame, Components(*myws.pdf("m_signal")), LineStyle(kDashed),
                   LineColor(kRed), Name("m_sig"),
                   Normalization(numEntries, RooAbsReal::NumEvent),
                   Range(2.6, 3.5), NormRange("Mass"));
  Fit_Total->plotOn(frame, Components(*myws.pdf("m_bkg")), LineStyle(kDashed),
                   LineColor(kGreen), Name("m_bkg"),
                   Normalization(numEntries, RooAbsReal::NumEvent),
                   Range(2.6, 3.5), NormRange("Mass"));*/

  /*RooPlot *frameTMP = (RooPlot *)frame->Clone("TMP");
  RooHist *hpull = frameTMP->pullHist(0, "Fit_Total", true);
  hpull->SetName("hpull");
  RooPlot *frame2 = mass->frame(Title("Pull Distribution (%)"));
  frame2->SetTitleSize(0.07);
  frame2->GetXaxis()->SetLabelSize(0.07);
  frame2->GetYaxis()->SetLabelSize(0.07);
  frame2->SetNdivisions(6, "Y");
  frame2->addPlotable(hpull, "PX");*/

  /*TCanvas *cFig = new TCanvas("2DFit", "2DFit", 1000, 800);
  TPad *pad1 = new TPad("pad1_2DFit", "", 0, 0.23, 1, 1);
  TPad *pad2 = new TPad("pad2_2DFit", "", 0, 0, 1, .228);

  pad2->SetBottomMargin(0.15);

  cFig->cd();
  pad1->Draw();
  pad1->cd();
  frame->Draw();

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextSize(0.032);

  float dy = 0;

  t->SetTextSize(0.03);
  t->DrawLatex(0.7175, 0.69 - dy, Form("%.1f < |y^{#mu#mu}| < %.1f", 2.5, 3.6));
  dy += 0.045;
  t->DrawLatex(0.6763, 0.69-dy, Form("%.0f GeV/c < p_{T}^{#mu#mu} < %.0f GeV/c", 1., 2.));
  // t->DrawLatex(0.5175, 0.69-dy, Form("%g < p_{T}^{#mu#mu} < %g
  // GeV/c",cut.dMuon.Pt.Min,cut.dMuon.Pt.Max)); dy+=0.045;
  t->DrawLatex(0.15, 0.9 - dy, Form("Fit parameters : "));
  dy += 0.045;
  t->DrawLatex(0.15, 0.9 - dy, Form("N_{signal} = %.0f#pm%.0f", nsig->getValV(), nsig->getError()));
  dy += 0.045;
  t->DrawLatex(
      0.15, 0.9 - dy,
      Form("N_{bkg} = %.0f#pm%.0f", nbkg->getValV(), nbkg->getError()));
  dy += 0.045;
  t->DrawLatex(0.15, 0.9 - dy,
               Form("f_{B} = %.6f#pm%.6f", myws.var("fB")->getValV(),
                    myws.var("fB")->getError()));
  dy += 0.045;*/

  //==============================================================
  //RooPlot *framebis = tau->frame(Title("l_{z} Projection"), Bins(100));
  RooPlot *framebis = tau->frame(Bins(nBins));
  framebis->SetTitle(" ");

  data.plotOn(framebis);
  Fit_Total->plotOn(framebis, Components(*myws.pdf("model2")), LineStyle(kDashed), LineColor(kBlack), Name("model2"), Normalization(numEntries, RooAbsReal::NumEvent), Range(tau_min, tau_max), NormRange("Tauz"));
  //Fit_Total->plotOn(framebis, Components(*myws.pdf("pdf")), LineStyle(kDashed), LineColor(kBlack), Name("pdf"), Normalization(numEntries, RooAbsReal::NumEvent), Range(tau_min, tau_max), NormRange("Tauz"));
  Fit_Total->plotOn(framebis, Components(*myws.pdf("tausigModel")), LineStyle(kDashed), LineColor(kRed+1), Name("tausigModel"), Normalization(numEntries, RooAbsReal::NumEvent), Range(tau_min, tau_max), NormRange("Tauz"));
  Fit_Total->plotOn(framebis, Components(*myws.pdf("exp_np")), LineStyle(kDashed), LineColor(kGreen+1), Name("exp_np"), Normalization(numEntries, RooAbsReal::NumEvent), Range(tau_min, tau_max), NormRange("Tauz"));
  Fit_Total->plotOn(framebis, Name("Fit_Total"), Normalization(numEntries, RooAbsReal::NumEvent), Range(tau_min, tau_max), NormRange("Tauz"));

  RooPlot *frameTMPbis = (RooPlot *)framebis->Clone("TMP");
  RooHist *hpullbis = frameTMPbis->pullHist(0, "Fit_Total", true);
  hpullbis->SetName("hpullbis");
  //RooPlot *frame2bis = tau->frame(Title("Pull Distribution (#sigma)"));
  RooPlot *frame2bis = tau->frame(Bins(25));
  frame2bis->SetTitle(" ");
  frame2bis->SetTitleSize(0.07);

  TLine* l_middle = new TLine(frame2bis->GetXaxis()->GetXmin(), 0, frame2bis->GetXaxis()->GetXmax(), 0);
  l_middle->SetLineColor(kBlack);
  l_middle->SetLineWidth(1);
  TLine* l_up = new TLine(frame2bis->GetXaxis()->GetXmin(), 2, frame2bis->GetXaxis()->GetXmax(), 2);
  l_up->SetLineColor(kBlack);
  l_up->SetLineWidth(1);
  l_up->SetLineStyle(kDashed);
  TLine* l_down = new TLine(frame2bis->GetXaxis()->GetXmin(), -2, frame2bis->GetXaxis()->GetXmax(), -2);
  l_down->SetLineColor(kBlack);
  l_down->SetLineWidth(1);
  l_down->SetLineStyle(kDashed);

  frame2bis->GetYaxis()->SetTitle("#frac{|hist - curve|}{#delta(hist)}");
  frame2bis->GetXaxis()->SetTitleSize(0.08);
  frame2bis->GetXaxis()->SetLabelSize(0.08);
  frame2bis->GetXaxis()->SetTitleOffset(0.95);
  frame2bis->GetYaxis()->SetTitleSize(0.08);
  frame2bis->GetYaxis()->SetLabelSize(0.08);
  frame2bis->GetYaxis()->SetTitleOffset(0.35);
  //frame2bis->GetYaxis()->SetTitleSize(0.07);
  //frame2bis->GetXaxis()->SetLabelSize(0.07);
  //frame2bis->GetYaxis()->SetLabelSize(0.07);
  frame2bis->SetNdivisions(6, "Y");
  frame2bis->addPlotable(hpullbis, "PX");

  TCanvas *cFigbis = new TCanvas("2DFitbis", "2DFitbis", 1000, 800);
  TPad *pad1bis = new TPad("pad1_2DFitbis", "", 0, 0.23, 1, 1); // xlow, ylow, xup, yup
  TPad *pad2bis = new TPad("pad2_2DFitbis", "", 0, 0, 1, .228); // 0.228

  float titlesize = 0.035;
  float titleoffset = 0.5;
  float labelsize = 0.035;

  //pad2bis->SetBottomMargin(0.15);
  pad2bis->SetBottomMargin(0.2);

  framebis->GetXaxis()->SetTitleSize(titlesize);
  //frame->GetXaxis()->SetTitleOffset(titleoffset);
  framebis->GetXaxis()->SetLabelSize(labelsize);
  framebis->GetYaxis()->SetTitleSize(titlesize);
  //frame->GetYaxis()->SetTitleOffset(titleoffset);
  framebis->GetYaxis()->SetLabelSize(labelsize);

  cFigbis->cd();
  pad1bis->Draw();
  pad1bis->cd();
  framebis->Draw();

  TLatex *t00 = new TLatex();
  t00->SetNDC();
  t00->SetTextSize(0.05);
  t00->DrawLatex(0.41984, 0.926995, "l_{z} Projection");

  TLatex *t0 = new TLatex();
  t0->SetNDC();
  t0->SetTextSize(0.035);
  t0->DrawLatex(0.410822, 0.00679117, "Pull Distribution (#sigma)");

  TLatex *t1 = new TLatex();
  t1->SetNDC();
  t1->SetTextSize(0.032);

  float dy1 = 0;

  t1->SetTextSize(0.03);
  t1->DrawLatex(0.7175, 0.69 - dy1,
                Form("%.1f < |y^{#mu#mu}| < %.1f", 2.5, 3.6));
  dy1 += 0.045;
  t1->DrawLatex(0.6763, 0.69-dy1, Form("%.0f GeV/c < p_{T}^{#mu#mu} < %.0f GeV/c", pT1->getValV(), pT2->getValV()));
  // t->DrawLatex(0.5175, 0.69-dy, Form("%g < p_{T}^{#mu#mu} < %g
  // GeV/c",cut.dMuon.Pt.Min,cut.dMuon.Pt.Max)); dy+=0.045;
  t1->DrawLatex(0.15, 0.9 - dy1, Form("Fit parameters : "));
  dy1 += 0.045;
  t1->DrawLatex(
      0.15, 0.9 - dy1,
      Form("N_{signal} = %.0f#pm%.0f", nsig->getValV(), nsig->getError()));
  dy1 += 0.045;
  t1->DrawLatex(
      0.15, 0.9 - dy1,
      Form("N_{bkg} = %.0f#pm%.0f", nbkg->getValV(), nbkg->getError()));
  dy1 += 0.045;
  t1->DrawLatex(0.15, 0.9 - dy1,
                Form("f_{B} = %.6f#pm%.6f", myws.var("fB")->getValV(),
                     myws.var("fB")->getError()));
  dy1 += 0.045;
  t1->DrawLatex(0.15, 0.9 - dy1,
                Form("#lambda_{np} = %.6f#pm%.6f",
                     myws.var("lambda_np")->getValV(),
                     myws.var("lambda_np")->getError()));
  dy1 += 0.045;
  //t1->DrawLatex(0.15, 0.9 - dy1, Form("x_{np} = %.6f#pm%.6f", myws.var("mean_np")->getValV(), myws.var("mean_np")->getError()));
  //dy1 += 0.045;
  //t1->DrawLatex(0.15, 0.9 - dy1, Form("#sigma_{np} = %.6f#pm%.6f", myws.var("sigma_np")->getValV(), myws.var("sigma_np")->getError()));
  //dy1 += 0.045;

  TLegend *leg_tau = new TLegend(0.6175, 0.75, 0.7880, 0.8809);
  //leg_tau->SetTextSize(0.03);
  leg_tau->SetTextSize(0.035);
  leg_tau->SetLineColor(0);
  TLegendEntry *l1 = leg_tau->AddEntry("Fit_Total", "Total fit", "l");
  l1->SetLineColor(4); // Blue
  TLegendEntry *l2 = leg_tau->AddEntry("tausigModel", "prompt", "l");
  l2->SetLineColor(kRed+1); // Red
  //TLegendEntry *l3 = leg_tau->AddEntry("model2", "background", "l");
  TLegendEntry *l3 = leg_tau->AddEntry("pdf", "background", "l");
  l3->SetLineColor(1); // Black
  TLegendEntry *l4 = leg_tau->AddEntry("exp_np", "non-prompt", "l");
  l4->SetLineColor(kGreen+1); // Green
  leg_tau->Draw("same");

  cFigbis->cd();

  pad2bis->Draw();
  pad2bis->cd();
  frame2bis->Draw();

  l_middle->Draw("same");
  l_up->Draw("same");
  l_down->Draw("same");

  double chi2 = 0;
  unsigned int ndof = 0;
  pad2bis->cd();
  TLatex *t1bis = new TLatex();
  t1bis->SetNDC();
  t1bis->SetTextSize(0.08);
  unsigned int nFitPar = Fit_Total->getParameters(data)
                             ->selectByAttrib("Constant", kFALSE)
                             ->getSize();
  TH1 *hdatact = data.createHistogram("hdatact", *tau,
                                      Binning(framebis->GetNbinsX(),
                                              framebis->GetXaxis()->GetXmin(),
                                              framebis->GetXaxis()->GetXmax()));
  RooHist *hpull1 = framebis->pullHist(0, 0, true);
  double *ypulls = hpull1->GetY();
  unsigned int nFullBins = 0;
  std::cout << "NBINS : " << nBins << std::endl;
  for (int i = 0; i < nBins; i++) {
    if (hdatact->GetBinContent(i + 1) > 1.0) {
      chi2 += ypulls[i] * ypulls[i];
      nFullBins++;
    }
  }
  ndof = nFullBins - nFitPar;
  //t1bis->DrawLatex(0.13, 0.8, Form("#chi^{2}/ndof = %.0f / %d ", chi2, ndof));
  pad2bis->Update();

  pad1bis->cd();
  t1->DrawLatex(0.7263, 0.6, Form("#chi^{2}/ndof = %.0f / %d ", chi2, ndof));
  pad1bis->Update();

  /*TH1F *frac = new TH1F(Form("fB_PDG_%.0f%.0f", pT1->getValV(), pT2->getValV()), "", 3, 0, 3);
  frac->GetXaxis()->SetBinLabel(1, "fB");
  frac->SetBinContent(1, fB->getValV());
  frac->GetXaxis()->SetBinLabel(2, "fB_Err");
  frac->SetBinContent(2, fB->getError());

  auto res = TFile::Open(Form("fB_PDG_Full_%.0f%.0f_tzJPsi_free_KF.root", pT1->getValV(), pT2->getValV()), "RECREATE");
  frac->Write(Form("fB_PDG_%.0f%.0f", pT1->getValV(), pT2->getValV()));*/

  //==============================================================
  /*TH2F *ph2 = tau->createHistogram("test", *mass);
  data.fillHistogram(ph2, RooArgList(*tau, *mass));
  Fit_Total->fillHistogram(ph2, RooArgList(*tau, *mass));

  TCanvas *cFigtest = new TCanvas("2DFittest", "2DFittest", 1000, 800);
  TPad *pad1test =
      new TPad("pad1_2DFittest", "", 0, 0, 1, 1); // xlow, ylow, xup, yup
  cFigtest->cd();
  pad1test->Draw();
  pad1test->cd();
  ph2->Draw("surf1");*/
}
