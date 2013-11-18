#pragma once

#include <memory>
#include <map>
#include <tuple>

#include <TGraphAsymmErrors.h>

#include "Extractors/PatExtractor/interface/ScaleFactorService.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class BTaggingEfficiencyProvider {
  public:
    BTaggingEfficiencyProvider(const edm::ParameterSet& settings) {
      TH1::AddDirectory(false);

      const edm::ParameterSet& pset = settings.getParameterSet("b_tagging_efficiency");
      std::string filename = edm::FileInPath(pset.getParameter<std::string>("filename")).fullPath();

      TFile* f = TFile::Open(filename.c_str());

      TH1* binning = static_cast<TH1*>(f->Get("binning"));

      // Get histo from file
      loadEfficiency(0, binning, pset.getParameter<std::string>("b_eff_histo_name"), f);
      loadEfficiency(1, binning, pset.getParameter<std::string>("cjets_fakerate_histo_name"), f);
      loadEfficiency(2, binning, pset.getParameter<std::string>("lightjets_fakerate_histo_name"), f);

      f->Close();
      delete f;
    }

    std::tuple<double, double, double> getEfficiency(ScaleFactorService::Flavor flavor, float pt, float eta) {
      if (pt > 800)
        pt = 800;

      int f = 0;
      switch (flavor) {
        case ScaleFactorService::B:
          f = 0;
          break;

        case ScaleFactorService::C:
          f = 1;
          break;

        case ScaleFactorService::LIGHT:
          f = 2;
          break;
      }

      eta = fabs(eta);
      for (const auto& etaBin: m_efficiency.at(f)) {
        if (eta >= etaBin.first.first && eta < etaBin.first.second) {
          // Look for pt bin
          for (const auto& ptBin: etaBin.second) {
            if (pt >= ptBin.first.first && pt < ptBin.first.second) {
              return ptBin.second;
            }
          }
        }
      }

      return std::make_tuple(1., 0., 0.);
    }

  private:
    void loadEfficiency(int flavor, TH1* binning, const std::string& name, TFile* file) {

      for (int etaBin = 1; etaBin <= binning->GetNbinsY(); etaBin++) {
        double eta_low = binning->GetYaxis()->GetBinLowEdge(etaBin);
        double eta_high = binning->GetYaxis()->GetBinUpEdge(etaBin);

        std::pair<double, double> eta = std::make_pair(eta_low, eta_high);

        TString histo_name = TString::Format("%s_%.02f_%.02f", name.c_str(), eta_low, eta_high);
        TGraphAsymmErrors* h = static_cast<TGraphAsymmErrors*>(file->Get(histo_name));

        for (int ptBin = 1; ptBin <= binning->GetNbinsX(); ptBin++) {
          double pt_low = binning->GetXaxis()->GetBinLowEdge(ptBin);
          double pt_high = binning->GetXaxis()->GetBinUpEdge(ptBin);

          std::pair<double, double> pt = std::make_pair(pt_low, pt_high);

          double dummy, eff, error_up, error_low;

          if (ptBin >= h->GetN()) {
            eff = error_up = error_low = 0;
            int bin = h->GetN() - 1;
            do {
              h->GetPoint(bin, dummy, eff);
              bin--;
            } while (eff == 0);
            error_up = 2 * h->GetErrorYhigh(bin + 1);
            error_low = 2 * h->GetErrorYlow(bin + 1);
          } else {

            int bin = ptBin;
            do {
              h->GetPoint(bin, dummy, eff);
              bin--;
            } while (eff == 0);

            error_up = h->GetErrorYhigh(bin + 1);
            error_low = h->GetErrorYlow(bin + 1);
          }

          //std::cout << "Flavor: " << flavor << " ; eta: [" << eta_low << "; " << eta_high << "] ; pt: [" << pt_low << "; " << pt_high << "] ; " << eff << " +" << error_up << " -" << error_low << std::endl;

          m_efficiency[flavor][eta][pt] = std::make_tuple(eff, error_up, error_low);
        }
      }
    }

    typedef std::map<
      int, // Flavor: 0 = B, 1 = C, 2 = light
      std::map<
        std::pair<double, double>, // Eta binning
        std::map<
          std::pair<double, double>, // Pt binning
          std::tuple<double, double, double>
        >
      >
    > EfficiencyMap;

    EfficiencyMap m_efficiency;
};
