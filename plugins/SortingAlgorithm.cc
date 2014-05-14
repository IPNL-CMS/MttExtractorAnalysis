#include "SortingAlgorithm.h"

void SortingAlgorithm::setObjects(const std::vector<Jet>& jets, const Lepton& lepton, const LorentzVector& met) {
  m_jets = jets;
  m_lepton = lepton;
  m_met = met;
}

bool SortingAlgorithm::computeNeutrinoPz(const LorentzVector& bJet, bool* no_real_sol/* = nullptr*/) {
  if (m_lepton.p.E() == 0)
    return false;

  if (no_real_sol)
    *no_real_sol = false;
  
  m_neutrino = m_met;
  m_neutrino.SetM(0.);

  const double m_w = 8.04190000E+01; // Value used in Madgraph generation
  const double m_top = 172.5; // Value used in Madgraph generation

  double x = (m_w * m_w - m_lepton.p.M() * m_lepton.p.M() + 2. * (m_neutrino.Px() * m_lepton.p.Px() + m_neutrino.Py() * m_lepton.p.Py())) / (2 * m_lepton.p.E());
  double a = 1 - (m_lepton.p.Pz() * m_lepton.p.Pz()) / (m_lepton.p.E() * m_lepton.p.E());
  double b = -2. * (m_lepton.p.Pz() / m_lepton.p.E()) * x;
  double c = m_neutrino.Pt() * m_neutrino.Pt() - x * x;

  if (!a && !b)
    return false;

  if (!a) {     
    m_neutrino.SetPz(-1 * c / b);

    return true;
  }

  double delta = b * b - 4 * a *c;

  if (delta < 0) {   // No solution, try to correct MET
    double rat = m_neutrino.Py() / m_neutrino.Px();

    double u = 4. / (m_lepton.p.E() * m_lepton.p.E()) * ((m_lepton.p.Px() + rat * m_lepton.p.Py()) * (m_lepton.p.Px() + rat * m_lepton.p.Py()) / (1 + rat * rat)
        - (m_lepton.p.E() * m_lepton.p.E()) + (m_lepton.p.Pz() * m_lepton.p.Pz()));

    double v = 4. / (m_lepton.p.E() * m_lepton.p.E()) * (m_w * m_w - m_lepton.p.M() * m_lepton.p.M())
      * (m_lepton.p.Px() + rat * m_lepton.p.Py()) / sqrt(1 + rat * rat);

    double w = (m_w * m_w - m_lepton.p.M() * m_lepton.p.M()) * (m_w * m_w - m_lepton.p.M() * m_lepton.p.M()) / (m_lepton.p.E() * m_lepton.p.E());

    double deltan = v * v - 4 * u * w;

    if (deltan < 0)
      return false; // Hopeless, MET can't be corrected

    double pt      = 0.;
    double corfact = 0.;

    if (u == 0) {
      pt = -w / v;
      if (pt <= 0)
        return false; // There is no way out...

      corfact = pt / m_neutrino.Pt();
    } else { // Deltan>=0 and u!=0
      double pt2 = (v - (sqrt(deltan))) / (2 * u);
      double pt1 = (v + (sqrt(deltan))) / (2 * u);

      // Pas de correction car negative
      if (pt1 <= 0 && pt2 <= 0) return 0;

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

    x = (m_w * m_w - m_lepton.p.M() * m_lepton.p.M() + 2.*(m_neutrino.Px() * m_lepton.p.Px() + m_neutrino.Py() * m_lepton.p.Py())) / (2 * m_lepton.p.E());
    a = 1 - (m_lepton.p.Pz() * m_lepton.p.Pz()) / (m_lepton.p.E() * m_lepton.p.E());
    b = -2.*(m_lepton.p.Pz() / m_lepton.p.E()) * x;
    c = m_neutrino.Px() * m_neutrino.Px() + m_neutrino.Py() * m_neutrino.Py() - x * x;

    delta = b * b - 4 * a * c;

    if (fabs(delta) < 0.000001) delta = 0.;

    if (delta != 0)
      return false; // This should not happen, but who knows...

    if (no_real_sol)
      *no_real_sol = true;
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

