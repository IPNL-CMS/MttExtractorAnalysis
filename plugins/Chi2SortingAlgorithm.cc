#include "Chi2SortingAlgorithm.h"

#include <TTree.h>

Chi2SortingAlgorithm::Chi2SortingAlgorithm(const edm::ParameterSet& cfg, bool isSemiMu) {
  m_isSemiMu = isSemiMu;
  m_useBTagInCombinatorics = cfg.getParameter<bool>("use_btag_in_combinatorics");

  m_w = cfg.getParameter<double>("w_mass");
  m_top = cfg.getParameter<double>("top_mass");
  m_b = cfg.getParameter<double>("b_mass");

  // Chi2 optimal values and errors
  chi2_hadronic_top_mass = cfg.getParameter<double>("hadronic_top_mass");

  chi2_hadronic_top_mass = cfg.getParameter<double>("hadronic_top_mass");
  chi2_leptonic_top_mass_semimu = cfg.getParameter<double>("leptonic_top_mass_semimu");
  chi2_leptonic_top_mass_semie = cfg.getParameter<double>("leptonic_top_mass_semie");
  chi2_hadronic_w_mass = cfg.getParameter<double>("hadronic_w_mass");
  chi2_pt_ttbar_system = cfg.getParameter<double>("pt_ttbar_system");
  chi2_ht_frac = cfg.getParameter<double>("ht_frac");

  chi2_sigma_hadronic_top_mass = cfg.getParameter<double>("sigma_hadronic_top_mass");
  chi2_sigma_leptonic_top_mass_semimu = cfg.getParameter<double>("sigma_leptonic_top_mass_semimu");
  chi2_sigma_leptonic_top_mass_semie = cfg.getParameter<double>("sigma_leptonic_top_mass_semie");
  chi2_sigma_hadronic_w_mass = cfg.getParameter<double>("sigma_hadronic_w_mass");
  chi2_sigma_pt_ttbar_system = cfg.getParameter<double>("sigma_pt_ttbar_system");
  chi2_sigma_ht_frac = cfg.getParameter<double>("sigma_ht_frac");

  m_usePtSystInChi2 = cfg.getParameter<bool>("use_pt_syst");
  m_useHtFracInChi2 = cfg.getParameter<bool>("use_ht_frac");

  std::cout << std::endl << "Chi^2 sorting algorithm: Using" << std::endl;
  std::cout << "\t- Hadronic top mass: " << chi2_hadronic_top_mass << " +/- " << chi2_sigma_hadronic_top_mass << std::endl;
  std::cout << "\t- Leptonic top mass (semi-mu): " << chi2_leptonic_top_mass_semimu << " +/- " << chi2_sigma_leptonic_top_mass_semimu << std::endl;
  std::cout << "\t- Leptonic top mass (semi-e): " << chi2_leptonic_top_mass_semie << " +/- " << chi2_sigma_leptonic_top_mass_semie << std::endl;
  std::cout << "\t- Hadronic W mass: " << chi2_hadronic_w_mass << " +/- " << chi2_sigma_hadronic_w_mass << std::endl;

  if (m_usePtSystInChi2)
    std::cout << "\t- pT tt system: " << chi2_pt_ttbar_system << " +/- " << chi2_sigma_pt_ttbar_system << std::endl;

  if (m_useHtFracInChi2)
    std::cout << "\t- HT fraction: " << chi2_ht_frac << " +/- " << chi2_sigma_ht_frac << std::endl;

  chi2_sigma_hadronic_top_mass_square = chi2_sigma_hadronic_top_mass * chi2_sigma_hadronic_top_mass;
  chi2_sigma_leptonic_top_mass_semimu_square = chi2_sigma_leptonic_top_mass_semimu * chi2_sigma_leptonic_top_mass_semimu;
  chi2_sigma_leptonic_top_mass_semie_square = chi2_sigma_leptonic_top_mass_semie * chi2_sigma_leptonic_top_mass_semie;
  chi2_sigma_hadronic_w_mass_square = chi2_sigma_hadronic_w_mass * chi2_sigma_hadronic_w_mass;
  chi2_sigma_pt_ttbar_system_square = chi2_sigma_pt_ttbar_system * chi2_sigma_pt_ttbar_system;
  chi2_sigma_ht_frac_square = chi2_sigma_ht_frac * chi2_sigma_ht_frac;

  m_selectedNeutrinoP4_AfterChi2 = new LorentzVector();
  m_selectedLeptonicBP4_AfterChi2 = new LorentzVector();
  m_selectedHadronicBP4_AfterChi2 = new LorentzVector();
  m_selectedFirstJetP4_AfterChi2 = new LorentzVector();
  m_selectedSecondJetP4_AfterChi2 = new LorentzVector();
  m_hadTopP4_AfterChi2 = new LorentzVector();
  m_lepTopP4_AfterChi2 = new LorentzVector();

}

Chi2SortingAlgorithm::~Chi2SortingAlgorithm() {
  delete m_selectedNeutrinoP4_AfterChi2;
  delete m_selectedLeptonicBP4_AfterChi2;
  delete m_selectedHadronicBP4_AfterChi2;
  delete m_selectedFirstJetP4_AfterChi2;
  delete m_selectedSecondJetP4_AfterChi2;
  delete m_hadTopP4_AfterChi2;
  delete m_lepTopP4_AfterChi2;
}

void Chi2SortingAlgorithm::addBranches(TTree& tree) {
  tree.Branch("numComb_chi2"       , &m_mtt_NumComb_chi2           , "numComb_chi2/I");
  tree.Branch("solChi2"            , &m_mtt_SolChi2                , "solChi2[numComb_chi2]/F");
  tree.Branch("bestSolChi2"        , &m_mtt_BestSolChi2           , "bestSolChi2/F");
  // If true, it means that we have selected the correct four reco jets (ie, reco jets coming from tt decay)
  tree.Branch("recoJetsAssociatedWithChi2" , &m_mtt_recoJetsAssociatedWithChi2    , "recoJetsAssociatedWithChi2/O");
  // If true, it means that we have selected the correct four reco jets (ie, reco jets coming from tt decay)
  // and each jets is correctly positionned.
  tree.Branch("recoJetsAssociatedWellPlacedWithChi2", &m_mtt_recoJetsAssociatedWellPlacedWithChi2, "recoJetsAssociatedWellPlacedWithChi2/O");
  // Selected object with Chi2
  tree.Branch("mLepW_AfterChi2"        , &m_mLepW_AfterChi2       , "mLepW_AfterChi2/F");
  tree.Branch("mHadW_AfterChi2"        , &m_mHadW_AfterChi2       , "mHadW_AfterChi2/F");

  tree.Branch("mLepTop_AfterChi2"      , &m_mLepTop_AfterChi2     , "mLepTop_AfterChi2/F");
  tree.Branch("lepTopPt_AfterChi2"     , &m_lepTopPt_AfterChi2    , "lepTopPt_AfterChi2/F");
  tree.Branch("lepTopEta_AfterChi2"    , &m_lepTopEta_AfterChi2   , "lepTopEta_AfterChi2/F");
  tree.Branch("lepTopP4_AfterChi2"    , &m_lepTopP4_AfterChi2);

  tree.Branch("mHadTop_AfterChi2"      , &m_mHadTop_AfterChi2     , "mHadTop_AfterChi2/F");
  tree.Branch("hadTopPt_AfterChi2"     , &m_hadTopPt_AfterChi2    , "hadTopPt_AfterChi2/F");
  tree.Branch("hadTopEta_AfterChi2"    , &m_hadTopEta_AfterChi2   , "hadTopEta_AfterChi2/F");
  tree.Branch("hadTopP4_AfterChi2"    , &m_hadTopP4_AfterChi2);

  tree.Branch("pt_tt_AfterChi2"        , &m_pt_tt_AfterChi2       , "pt_tt_AfterChi2/F");
  tree.Branch("eta_tt_AfterChi2"       , &m_eta_tt_AfterChi2      , "eta_tt_AfterChi2/F");
  tree.Branch("beta_tt_AfterChi2"      , &m_beta_tt_AfterChi2     , "eta_tt_AfterChi2/F");
  tree.Branch("mtt_AfterChi2"          , &m_mtt_AfterChi2         , "mtt_AfterChi2/F");

  tree.Branch("neutrino_no_real_solution_AfterChi2", &m_neutrino_no_real_solution_AfterChi2, "neutrino_no_real_solution_AfterChi2/I");

  tree.Branch("selectedLeptonicBIndex_AfterChi2"     , &m_selectedLeptonicBIndex_AfterChi2    , "selectedLeptonicBIndex_AfterChi2/I");
  tree.Branch("selectedHadronicBIndex_AfterChi2"     , &m_selectedHadronicBIndex_AfterChi2    , "selectedHadronicBIndex_AfterChi2/I");
  tree.Branch("selectedHadronicFirstJetIndex_AfterChi2"  , &m_selectedHadronicFirstJetIndex_AfterChi2  , "selectedHadronicFirstJetIndex_AfterChi2/I");
  tree.Branch("selectedHadronicSecondJetIndex_AfterChi2" , &m_selectedHadronicSecondJetIndex_AfterChi2 , "selectedHadronicSecondJetIndex_AfterChi2/I");

  tree.Branch("selectedNeutrinoP4_AfterChi2"    , &m_selectedNeutrinoP4_AfterChi2);
  tree.Branch("selectedLeptonicBP4_AfterChi2"   , &m_selectedLeptonicBP4_AfterChi2);
  tree.Branch("selectedHadronicBP4_AfterChi2"   , &m_selectedHadronicBP4_AfterChi2);
  tree.Branch("selectedFirstJetP4_AfterChi2"    , &m_selectedFirstJetP4_AfterChi2);
  tree.Branch("selectedSecondJetP4_AfterChi2"   , &m_selectedSecondJetP4_AfterChi2);
}

void Chi2SortingAlgorithm::work() {
  size_t n_jets = m_jets.size();
  int n_btaggedjets = 0;

  double allJetsPt = 0.;
  for (const Jet& i : m_jets) {
    allJetsPt += i.p.pt();
    if (m_useBTagInCombinatorics && i.isBTagged) {
      n_btaggedjets++;
    }
  }

  int numberoflightjets = n_jets - n_btaggedjets;

  if (numberoflightjets < 2)
    return; // if we dont have at least 2 non b-tagged jets, chi2 is -1

  double minChi2 = std::numeric_limits<double>::infinity();

  Jet leptonicBJet;
  Jet hadronicBJet;
  Jet hadronicLightJet1;
  Jet hadronicLightJet2;

  /*
   * Indices:
   * Jet 1: First b-jet (leptonic)
   * Jet 2: Second b-jet (hadronic)
   * Jet 3: First light jet
   * Jet 4: Second light jet
   */

  for (size_t bj1 = 0; bj1 < n_jets; ++bj1)
  {
    for (size_t bj2 = 0; bj2 < n_jets; ++bj2)
    {
      if (bj2 == bj1)
        continue; //dont pick the one you already used

      for (size_t j3 = 0; j3 < n_jets; ++j3)
      {
        // dont pick the two jets you used, or btagged jets
        if (j3 == bj1 || j3 == bj2 || (m_useBTagInCombinatorics && m_jets[j3].isBTagged))
          continue;

        for (size_t j4 = j3 + 1; j4 < n_jets; ++j4)
        {
          //dont pick the two jets you used, or btagged jets
          if (j4 == bj1 || j4 == bj2 || (m_useBTagInCombinatorics && m_jets[j4].isBTagged))
            continue;

          bool res = computeNeutrinoPz(m_jets[bj1].p);

          if (!res)
            return; // We will never get anything with this event

          double currentChi2 = chi2(m_jets[bj1].p, m_jets[bj2].p, m_jets[j3].p, m_jets[j4].p, allJetsPt);
          if (currentChi2 < minChi2) {
            minChi2 = currentChi2;

            m_mtt_SolChi2[m_mtt_NumComb_chi2] = currentChi2;
            m_mtt_NumComb_chi2++;

            leptonicBJet = m_jets[bj1];
            hadronicBJet = m_jets[bj2];
            hadronicLightJet1 = m_jets[j3];
            hadronicLightJet2 = m_jets[j4];
          }
        } // j4 close
      } // j3 close
    } // j2 close
  } // we out of the loop over jet pairings

  if (m_mtt_NumComb_chi2 > 0) {
    m_selectedHadronicBIndex_AfterChi2 = hadronicBJet.index;
    m_selectedLeptonicBIndex_AfterChi2 = leptonicBJet.index;
    m_selectedHadronicFirstJetIndex_AfterChi2 = hadronicLightJet1.index;
    m_selectedHadronicSecondJetIndex_AfterChi2 = hadronicLightJet2.index;

    // Correct MET again
    bool no_real_sol = false;
    computeNeutrinoPz(leptonicBJet.p, &no_real_sol);

    m_neutrino_no_real_solution_AfterChi2 = (no_real_sol) ? 1 : 0;

    // Copy P4 for Tree
    *m_selectedNeutrinoP4_AfterChi2 = m_neutrino;
    *m_selectedLeptonicBP4_AfterChi2 = leptonicBJet.p;
    *m_selectedHadronicBP4_AfterChi2 = hadronicBJet.p;
    *m_selectedFirstJetP4_AfterChi2 = hadronicLightJet1.p;
    *m_selectedSecondJetP4_AfterChi2 = hadronicLightJet2.p;

    /**
     * Compute Mtt before doing KinFit
     */
    m_mLepW_AfterChi2   = (m_neutrino + m_lepton).M();
    m_mHadW_AfterChi2   = (hadronicLightJet1.p + hadronicLightJet2.p).M();

    LorentzVector lepTopP4_AfterChi2 = (m_lepton + m_neutrino + leptonicBJet.p);
    *m_lepTopP4_AfterChi2 = lepTopP4_AfterChi2;
    m_mLepTop_AfterChi2 = lepTopP4_AfterChi2.M();
    m_lepTopPt_AfterChi2 = lepTopP4_AfterChi2.Pt();
    m_lepTopEta_AfterChi2 = lepTopP4_AfterChi2.Eta();

    LorentzVector hadTopP4_AfterChi2 = (hadronicLightJet1.p + hadronicLightJet2.p + hadronicBJet.p);
    *m_hadTopP4_AfterChi2 = (hadronicLightJet1.p + hadronicLightJet2.p + hadronicBJet.p);
    m_mHadTop_AfterChi2 = hadTopP4_AfterChi2.M();
    m_hadTopPt_AfterChi2 = hadTopP4_AfterChi2.Pt();
    m_hadTopEta_AfterChi2 = hadTopP4_AfterChi2.Eta();

    LorentzVector res = (lepTopP4_AfterChi2 + hadTopP4_AfterChi2);
    m_mtt_AfterChi2     = res.M();
    m_pt_tt_AfterChi2   = res.Pt();
    m_eta_tt_AfterChi2   = res.Eta();
    m_beta_tt_AfterChi2  = fabs(res.Pz() / res.E());

    m_mtt_BestSolChi2 = minChi2;
  }

}

double Chi2SortingAlgorithm::chi2(const LorentzVector& leptonicBJet, const LorentzVector& hadronicBJet, const LorentzVector& hadronicLightJet1, const LorentzVector& hadronicLightJet2, double allJetsPt) {
  float MW          = sqrt(std::max(0., (hadronicLightJet1 + hadronicLightJet2).M2()));
  float MtopH       = sqrt(std::max(0., (hadronicLightJet1 + hadronicLightJet2 + hadronicBJet).M2()));
  float MtopL       = sqrt(std::max(0., (m_neutrino + m_lepton + leptonicBJet).M2()));
  float SolPtSystem = (hadronicLightJet1.Pt() + hadronicLightJet2.Pt() + leptonicBJet.Pt() + hadronicBJet.Pt()) / allJetsPt;
  float TTbarSystemPt = ((hadronicLightJet1 + hadronicLightJet2 + leptonicBJet + hadronicBJet + m_neutrino + m_lepton).Pt());

  float chi2 = ((MtopH - chi2_hadronic_top_mass) * (MtopH - chi2_hadronic_top_mass) / (chi2_sigma_hadronic_top_mass_square)) + 
    ((MW - chi2_hadronic_w_mass) * (MW - chi2_hadronic_w_mass) / (chi2_sigma_hadronic_w_mass_square));

  if (m_useHtFracInChi2) {
    chi2 += ((SolPtSystem - chi2_ht_frac) * (SolPtSystem - chi2_ht_frac) / (chi2_sigma_ht_frac_square));
  }

  if (m_usePtSystInChi2) {
    chi2 += ((TTbarSystemPt - chi2_pt_ttbar_system) * (TTbarSystemPt - chi2_pt_ttbar_system) / chi2_sigma_pt_ttbar_system_square);
  }

  /// lepton dependant part
  (m_isSemiMu)
    ? chi2 += (MtopL - chi2_leptonic_top_mass_semimu) * (MtopL - chi2_leptonic_top_mass_semimu) / (chi2_sigma_leptonic_top_mass_semimu_square)
    : chi2 += (MtopL - chi2_leptonic_top_mass_semie) * (MtopL - chi2_leptonic_top_mass_semie) / (chi2_sigma_leptonic_top_mass_semie_square);

  return  chi2;
}

void Chi2SortingAlgorithm::reset() {
  m_neutrino_no_real_solution_AfterChi2 = -1;

  m_mtt_recoJetsAssociatedWithChi2 = false;
  m_mtt_recoJetsAssociatedWellPlacedWithChi2 = false;
  m_mtt_BestSolChi2 = -1.;
  m_mtt_NumComb_chi2 = 0;

  m_mtt_AfterChi2            = -1.;
  m_pt_tt_AfterChi2          = -1.;
  m_eta_tt_AfterChi2         = -1.;
  m_beta_tt_AfterChi2        = -1.;
  m_mLepTop_AfterChi2        = -1.;
  m_mHadTop_AfterChi2        = -1.;
  m_mHadW_AfterChi2          = -1.;

  m_lepTopPt_AfterChi2  = -1;
  m_lepTopEta_AfterChi2 = -1;
  m_hadTopPt_AfterChi2  = -1;
  m_hadTopEta_AfterChi2 = -1;

  m_selectedLeptonicBIndex_AfterChi2            = -1;
  m_selectedHadronicBIndex_AfterChi2            = -1;
  m_selectedHadronicFirstJetIndex_AfterChi2     = -1;
  m_selectedHadronicSecondJetIndex_AfterChi2    = -1;

  m_selectedNeutrinoP4_AfterChi2->SetCoordinates(0, 0, 0, 0);
  m_selectedLeptonicBP4_AfterChi2->SetCoordinates(0, 0, 0, 0);
  m_selectedHadronicBP4_AfterChi2->SetCoordinates(0, 0, 0, 0);
  m_selectedFirstJetP4_AfterChi2->SetCoordinates(0, 0, 0, 0);
  m_selectedSecondJetP4_AfterChi2->SetCoordinates(0, 0, 0, 0);
  m_lepTopP4_AfterChi2->SetCoordinates(0, 0, 0, 0);
  m_hadTopP4_AfterChi2->SetCoordinates(0, 0, 0, 0);
}
