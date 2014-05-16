#include "GaussianProfile.h"

#include <sstream>
#include <TH1D.h>
#include <TF1.h>


void GaussianProfile::createProfiles() {
  for (int i = 0; i < m_nXBins; i++) {

    std::stringstream ss;
    ss << m_name << "_" << m_prefix << "_" << (int) m_XBins[i] << "_" << (int) m_XBins[i + 1];

    int nBins = m_nYBins;
    double min = m_YMin, max = m_YMax;

    if (m_autoBinning) {
      nBins = 100;
      min = m_XBins[i] - m_autoBinningLowPercent * m_XBins[i];
      max = m_XBins[i + 1] + m_autoBinningHighPercent * m_XBins[i + 1];
    }

    TH1* object = new TH1D(ss.str().c_str(), ss.str().c_str(), nBins, min, max);
    object->SetDirectory(nullptr);
    m_profiles.push_back(object);
  }
}

void GaussianProfile::createGraph() {

  if (m_profiles.size() == 0 || (m_response.get() && !m_dirty))
    return;

  std::stringstream ss;
  ss << m_name << "_resolution";
  std::string resolutionName = ss.str();

  ss.str("");
  ss << m_name << "_response";
  std::string responseName = ss.str();

  ss.str("");
  ss << m_name << "_rms";
  std::string rmsName = ss.str();

  m_response.reset(new TGraphErrors(m_nXBins));
  m_response->SetName(responseName.c_str());

  m_resolution.reset(new TGraphErrors(m_nXBins));
  m_resolution->SetName(resolutionName.c_str());

  m_rms.reset(new TGraphErrors(m_nXBins));
  m_rms->SetName(rmsName.c_str());
  
  // Create gaussian for fitting
  TF1* gauss = new TF1("g", "gaus");

  for (int i = 0; i < m_nXBins; i++) {
    double mean = (m_XBins[i] + m_XBins[i + 1]) / 2.;

    TH1* hist = m_profiles[i];

    //double min = hist->GetXaxis()->GetBinLowEdge(1);
    //double max = hist->GetXaxis()->GetBinUpEdge(hist->GetXaxis()->GetLast());
    
    double min = hist->GetXaxis()->GetXmin();
    double max = hist->GetXaxis()->GetXmax();
    if (m_autoBinning) {
      min = m_XBins[i] * 0.90;
      max = m_XBins[i + 1] * 1.10;
    }

    gauss->SetRange(min, max);
    gauss->SetParameter(1, (min + max) / 2.);
    gauss->SetParLimits(1, min, max);
    gauss->SetParameter(2, 20);

    hist->Fit(gauss, "QR");

    if (! m_autoBinning) {
      //double mean = hist->GetBinCenter(hist->GetMaximumBin());
      double sigma = gauss->GetParameter(2);
      gauss->SetRange(-1 * sigma, 1 * sigma);
      //gauss->FixParameter(1, mean);
      gauss->SetParameter(1, 0);
      gauss->SetParameter(2, 20);

      hist->Fit(gauss, "QR");
    }

    m_response->SetPoint(i, mean, gauss->GetParameter(1));
    m_response->SetPointError(i, 0, gauss->GetParError(1));

    m_rms->SetPoint(i, mean, hist->GetRMS());
    m_rms->SetPointError(i, 0, hist->GetRMSError());

    m_resolution->SetPoint(i, mean, gauss->GetParameter(2));
    m_resolution->SetPointError(i, 0, gauss->GetParError(2));
  }

  delete gauss;
  m_dirty = false;
}
