#pragma once

#include <vector>
#include <iostream>
#include <sstream>
#include <memory>

#include <TH1.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TDirectory.h>

class GaussianProfile {

  public:
    GaussianProfile(const std::string& name, int nBinsX, const double* binsX, bool doGraph = true):
      m_name(name), m_prefix("pt"), m_autoBinning(true), m_autoBinningLowPercent(0.4), m_autoBinningHighPercent(0.4), m_nXBins(nBinsX), m_dirty(true), m_doGraph(doGraph) {
        m_XBins.assign(binsX, binsX + nBinsX + 1);
      }

    GaussianProfile(const std::string& name, int nBinsX, const double* binsX, int nBinsY, double yMin, double yMax, bool doGraph = true):
      m_name(name), m_prefix("pt"), m_autoBinning(false), m_autoBinningLowPercent(0), m_autoBinningHighPercent(0), m_nXBins(nBinsX),
      m_nYBins(nBinsY), m_YMin(yMin), m_YMax(yMax), m_dirty(true), m_doGraph(doGraph) {
        m_XBins.assign(binsX, binsX + nBinsX + 1);
      }

    virtual ~GaussianProfile() {
      for (TH1* h: m_profiles) {
        delete h;
      }
    }

    void fill(double x, double y, double weight = 1) {
      if (m_profiles.size() == 0) {
        createProfiles();
      }

      int bin = findBin(x);
      if (bin < 0)
        return;

      m_profiles[bin]->Fill(y, weight);
      m_dirty = true;
    }

    void drawBin(int bin, Option_t* options) {
      if (bin < 0 || bin >= m_nXBins || m_profiles.size() == 0)
        return;

      m_profiles[bin]->Draw(options);
    }

    void drawResponse(Option_t* options) {
      if (! m_doGraph)
        return;

      if (m_dirty) {
        createGraph();
      }

      m_response->Draw(options);
    }

    void drawResolution(Option_t* options) {
      if (! m_doGraph)
        return;

      if (m_dirty) {
        createGraph();
      }

      m_resolution->Draw(options);
    }

    TGraphErrors* getResolutionGraph() {
      if (! m_doGraph)
        return nullptr;

      if (m_dirty) {
        createGraph();
      }

      return m_resolution.get();
    }

    void drawRMS(Option_t* options) {
      if (! m_doGraph)
        return;

      if (m_dirty) {
        createGraph();
      }

      m_rms->Draw(options);
    }

    TGraphErrors* getRMSGraph() {
      if (! m_doGraph)
        return nullptr;

      if (m_dirty) {
        createGraph();
      }

      return m_rms.get();
    }

    void setAutoBinningPercent(double lowPercent, double highPercent) {
      m_autoBinningLowPercent = lowPercent;
      m_autoBinningHighPercent = highPercent;
    }

    void setPrefix(const std::string prefix) {
      m_prefix = prefix;
    }

    void write(TFile* f = NULL) {
      if (m_doGraph && m_dirty) {
        createGraph();
      }

      TDirectory* currentDir = NULL;
      if (f) {
        currentDir = f;
        f->mkdir(m_name.c_str());
        f->cd(m_name.c_str());
      } else {
        currentDir = gDirectory;
        gDirectory->mkdir(m_name.c_str());
        gDirectory->cd(m_name.c_str());      
      }

      for (TH1* h: m_profiles) {
        h->Write();
      }

      if (m_doGraph) {
        if (m_response.get())
          m_response->Write();

        if (m_resolution.get())
          m_resolution->Write();

        if (m_rms.get())
          m_rms->Write();
      }

      currentDir->cd();
    }

  private:

    void createProfiles();
    void createGraph();

    int findBin(double value) {
      for (int i = 0; i < m_nXBins; i++) {
        if (value >= m_XBins[i] && value < m_XBins[i + 1])
          return i;
      }

      return -1;
    }

    std::string m_name;
    std::string m_prefix;

    bool m_autoBinning;
    double m_autoBinningLowPercent;
    double m_autoBinningHighPercent;

    int m_nXBins;
    std::vector<double> m_XBins;

    int m_nYBins;
    double m_YMin;
    double m_YMax;

    bool m_dirty;
    std::vector<TH1*> m_profiles;
    std::shared_ptr<TGraphErrors> m_response;
    std::shared_ptr<TGraphErrors> m_resolution;
    std::shared_ptr<TGraphErrors> m_mean;
    std::shared_ptr<TGraphErrors> m_rms;

    bool m_doGraph;
};
