#include "SortingAlgorithm.h"

void SortingAlgorithm::setObjects(const std::vector<Jet>& jets, const Lepton& lepton, const LorentzVector& met) {
  m_jets = jets;
  m_lepton = lepton;
  m_met = met;
}

bool SortingAlgorithm::computeNeutrinoPz(const LorentzVector& bJet, bool* no_real_sol/* = nullptr*/, bool* pt_corrected/* = nullptr*/) {
  if (m_lepton.p.E() == 0) {
    return false;
  }

  if (no_real_sol)
    *no_real_sol = false;

  if (pt_corrected)
    *pt_corrected = false;
  
  m_neutrino = m_met;
  m_neutrino.SetM(0.);

  const double m_w = 8.04190000E+01; // Value used in Madgraph generation
  const double m_top = 172.5; // Value used in Madgraph generation

  double x = (m_w * m_w - m_lepton.p.M() * m_lepton.p.M() + 2. * (m_neutrino.Px() * m_lepton.p.Px() + m_neutrino.Py() * m_lepton.p.Py())) / (2 * m_lepton.p.E());
  double a = 1 - (m_lepton.p.Pz() * m_lepton.p.Pz()) / (m_lepton.p.E() * m_lepton.p.E());
  double b = -2. * (m_lepton.p.Pz() / m_lepton.p.E()) * x;
  double c = m_neutrino.Pt() * m_neutrino.Pt() - x * x;

  if (!a && !b) {
    return false;
  }

  if (!a) {     
    m_neutrino.SetPz(-1 * c / b);

    return true;
  }

  double delta = b * b - 4 * a *c;

  if (delta < 0) {   // No solution, try to correct MET
    if (no_real_sol)
      *no_real_sol = true;

    double rat = m_neutrino.Py() / m_neutrino.Px();
    double gamma_x = m_lepton.p.Px() * 1. / std::sqrt(1. + std::pow(rat, 2));
    double gamma_y = m_lepton.p.Py() * 1. / std::sqrt(1. + std::pow(1. / rat, 2));
    if (m_neutrino.Px() < 0.) {
      gamma_x *= -1.;
    }
    if (m_neutrino.Py() < 0.) {
      gamma_y *= -1.;
    }
    double u = std::pow(m_lepton.p.Pz() / m_lepton.p.E(), 2) + std::pow((gamma_x + gamma_y) / m_lepton.p.E(), 2) - 1.;
    double v = std::pow(1. / m_lepton.p.E(), 2) * (gamma_x + gamma_y) * (std::pow(m_w, 2) - std::pow(m_lepton.p.M(), 2));
    double w = std::pow((std::pow(m_w, 2) - std::pow(m_lepton.p.M(), 2)) / (2 * m_lepton.p.E()), 2);

    double deltan = v * v - 4 * u * w;

    if (deltan < 0) {
      return false; // Hopeless, MET can't be corrected
    }

    double pt      = 0.;
    double corfact = 0.;

    if (u == 0) {
      pt = -w / v;
      if (pt <= 0) {
        return false; // There is no way out...
      }

      corfact = pt / m_neutrino.Pt();
    } else { // Deltan>=0 and u!=0
      double pt2 = (-v - (sqrt(deltan))) / (2 * u);
      double pt1 = (-v + (sqrt(deltan))) / (2 * u);

      // Pas de correction car negative
      if (pt1 <= 0 && pt2 <= 0) {
        return false;
      }

      if (pt1 > 0 && pt2 < 0) pt = pt1;
      if (pt2 > 0 && pt1 < 0) pt = pt2;
      if (pt1 > 0 && pt2 > 0) {
        (fabs(pt1 - m_neutrino.Pt()) <= fabs(pt2 - m_neutrino.Pt()))
          ? pt = pt1
          : pt = pt2;
      }

      corfact = pt / m_neutrino.Pt();
    }

    // Now we have the correction factor

    m_neutrino.SetPx(corfact * m_neutrino.Px());
    m_neutrino.SetPy(corfact * m_neutrino.Py());

    // Recompute the new parameters

    x = (m_w * m_w - m_lepton.p.M() * m_lepton.p.M() + 2. * (m_neutrino.Px() * m_lepton.p.Px() + m_neutrino.Py() * m_lepton.p.Py())) / (2 * m_lepton.p.E());
    a = 1 - (m_lepton.p.Pz() * m_lepton.p.Pz()) / (m_lepton.p.E() * m_lepton.p.E());
    b = -2. * (m_lepton.p.Pz() / m_lepton.p.E()) * x;
    c = m_neutrino.Pt() * m_neutrino.Pt() - x * x;

    delta = b * b - 4 * a * c;

    if (fabs(delta) < 0.000001) delta = 0.;

    if (delta != 0) {
      std::cout << "---> !!!!! Delta after neutrino pt correction is non-zero: THIS SHOULD NOT HAPPEN !!!!!" << std::endl;
      return false; // This should not happen, but who knows...
    }

    if (pt_corrected)
      *pt_corrected = true;
  }

  // We can go back to the normal path: 

  LorentzVector TopCand1 = m_lepton.p + bJet;
  LorentzVector TopCand2 = m_lepton.p + bJet;

  m_neutrino.SetPz((-b - (sqrt(delta))) / (2 * a));
  TopCand1 += m_neutrino;

  m_neutrino.SetPz((-b + (sqrt(delta))) / (2 * a));
  TopCand2 += m_neutrino;

  double mtt_1 = sqrt(std::max(0., TopCand1.M2()));
  double mtt_2 = sqrt(std::max(0., TopCand2.M2()));

  if (fabs(mtt_1 - m_top) <= fabs(mtt_2 - m_top)) { // Otherwise it's OK
    m_neutrino.SetPz((-b - (sqrt(delta))) / (2 * a));
  }

  return true;
}

