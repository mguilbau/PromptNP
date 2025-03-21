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
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

void Incertitude2() {

  std::vector<float> range = {0., 1., 2., 3., 4., 5., 6., 8., 10., 15., 30.};

  for (int i = 0; i < range.size() - 1; i++) {
    //std::cout << "test : " << range[i] << " " << range[i + 1] << std::endl;

    float pT1 = range[i];
    float pT2 = range[i + 1];

    TFile *fB_pT01 = TFile::Open(Form("/Users/emiliebarreau/alice/Macros/fB_PDG_Full_%.0f%.0f_tzJPsi_free_MEH.root", pT1, pT2));
    TH1F *fB_01 = static_cast<TH1F *>(fB_pT01->Get(Form("fB_PDG_%.0f%.0f", pT1, pT2)));

    TFile *fB_pT01_ambig = TFile::Open(Form("/Users/emiliebarreau/alice/Macros/fB_PDG_Full_%.0f%.0f_chi2_MIN.root", pT1, pT2));
    TH1F *fB_01_ambig = static_cast<TH1F *>(fB_pT01_ambig->Get(Form("fB_PDG_%.0f%.0f", pT1, pT2)));

    TFile *fB_pT01_ambigbis = TFile::Open(Form("/Users/emiliebarreau/alice/Macros/fB_PDG_Full_%.0f%.0f_chi2_MAX.root", pT1, pT2));
    TH1F *fB_01_ambigbis = static_cast<TH1F *>(fB_pT01_ambigbis->Get(Form("fB_PDG_%.0f%.0f", pT1, pT2)));

    double fB = fB_01->GetBinContent(1);
    double fB_E = fB_01->GetBinContent(2);

    double fB_ambig = fB_01_ambig->GetBinContent(1);
    double fB_E_ambig = fB_01_ambig->GetBinContent(2);

    double fB_ambigbis = fB_01_ambigbis->GetBinContent(1);
    double fB_E_ambigbis = fB_01_ambigbis->GetBinContent(2);

    double mean = (fB + fB_ambig + fB_ambigbis) / 3.;

    //std::cout << fB << std::endl;
    //std::cout << fB_E << std::endl;
    //std::cout << fB_ambig << std::endl;
    //std::cout << fB_E_ambig << std::endl;

    // Example set of fit results: { (mean, stat_uncertainty), ... }
    std::vector<std::pair<double, double>> fit_results = {
        {fB, fB_E}, {fB_ambig, fB_E_ambig}, {fB_ambigbis, fB_E_ambigbis}};

    int N = fit_results.size();

    if (N == 0) {
      std::cerr << "No data provided!" << std::endl;
      return;
    }

    // Compute mean of means
    double sum_means = 0;
    double sum_stat_uncertainties = 0;

    for (const auto &result : fit_results) {
      sum_means += result.first;
      sum_stat_uncertainties += result.second;
    }
    double mean_total = sum_means / N;

    // Compute systematic uncertainty (Standard deviation of the means)
    double sum_sq_diff = 0;
    for (const auto &result : fit_results) {
      sum_sq_diff += std::pow(result.first - mean_total, 2);
    }
    double sys_uncertainty = (N > 1) ? std::sqrt(sum_sq_diff / (N - 1)) : 0.0;

    // Compute statistical uncertainty as the mean of the statistical
    // uncertainties
    double stat_uncertainty = sum_stat_uncertainties / N;

    // Compute total uncertainty (quadrature sum of systematic and statistical)
    double total_uncertainty =
        std::sqrt(std::pow(sys_uncertainty, 2) + std::pow(stat_uncertainty, 2));

    // Print results
    //std::cout << "Mean Value: " << mean_total << std::endl;
    //std::cout << "Systematic Uncertainty: " << sys_uncertainty << std::endl;
    //std::cout << "Statistical Uncertainty: " << stat_uncertainty << std::endl;
    //std::cout << "Total Uncertainty: " << total_uncertainty << std::endl;
    std::cout << std::fixed;
    std::cout << std::setprecision(3);
    std::cout << pT1 << "-" << pT2  << " & " << sys_uncertainty << " & " << stat_uncertainty << std::endl;
    //std::cout << pT1 << "-" << pT2 << " & " << mean << " & " << sys_uncertainty / mean << " & " << stat_uncertainty / mean << std::endl;
  }
}