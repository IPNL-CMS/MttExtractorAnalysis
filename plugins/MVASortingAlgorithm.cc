#include "MVASortingAlgorithm.h"
#include "TMVA/Reader.h"
#include <Math/GenVector/VectorUtil.h>

#include <TTree.h>

MVASortingAlgorithm::MVASortingAlgorithm(const edm::ParameterSet& cfg) {
  m_useBTagInCombinatorics = cfg.getParameter<bool>("use_btag_in_combinatorics");

  m_MVAWeightFilename = edm::FileInPath(cfg.getParameter<std::string>("weights")).fullPath();
  m_MVAMethodName = cfg.getParameter<std::string>("name");
  m_MVACut = cfg.getParameter<bool>("cut");
  if (m_MVACut) {
    m_MVACutValue = cfg.getParameter<double>("cut_value");
  }

  m_MVAReader.reset(new TMVA::Reader("V"));

  m_MVAReader->AddVariable("lightJet1p2_Pt", &m_mva_lightJet1p2_Pt);
  m_MVAReader->AddVariable("leptonic_B_Pt", &m_mva_leptonic_B_Pt);
  m_MVAReader->AddVariable("leptonic_W_Pt", &m_mva_leptonic_W_Pt);
  m_MVAReader->AddVariable("leptonic_Top_Pt", &m_mva_leptonic_Top_Pt);
  m_MVAReader->AddVariable("leptonic_Top_M", &m_mva_leptonic_Top_M);
  m_MVAReader->AddVariable("hadronic_B_Pt", &m_mva_hadronic_B_Pt);
  m_MVAReader->AddVariable("hadronic_W_Pt", &m_mva_hadronic_W_Pt);
  m_MVAReader->AddVariable("hadronic_W_M", &m_mva_hadronic_W_M);
  m_MVAReader->AddVariable("hadronic_Top_Pt", &m_mva_hadronic_Top_Pt);
  m_MVAReader->AddVariable("hadronic_Top_M", &m_mva_hadronic_Top_M);

  m_MVAReader->AddVariable("delta_phi_tops", &m_mva_delta_phi_tops);
  m_MVAReader->AddVariable("delta_phi_lightjets", &m_mva_delta_phi_lightjets);
  m_MVAReader->AddVariable("delta_phi_W", &m_mva_delta_phi_W);
  m_MVAReader->AddVariable("delta_R_tops", &m_mva_delta_R_tops);
  m_MVAReader->AddVariable("delta_R_lightjets", &m_mva_delta_R_lightjets);
  m_MVAReader->AddVariable("delta_R_W", &m_mva_delta_R_W);

  m_MVAReader->BookMVA(m_MVAMethodName, m_MVAWeightFilename.c_str());

  m_selectedNeutrinoP4_AfterMVA = new LorentzVector();
  m_selectedLeptonicBP4_AfterMVA = new LorentzVector();
  m_selectedHadronicBP4_AfterMVA = new LorentzVector();
  m_selectedFirstJetP4_AfterMVA = new LorentzVector();
  m_selectedSecondJetP4_AfterMVA = new LorentzVector();
  m_hadTopP4_AfterMVA = new LorentzVector();
  m_lepTopP4_AfterMVA = new LorentzVector();
}

MVASortingAlgorithm::~MVASortingAlgorithm() {
  delete m_selectedNeutrinoP4_AfterMVA;
  delete m_selectedLeptonicBP4_AfterMVA;
  delete m_selectedHadronicBP4_AfterMVA;
  delete m_selectedFirstJetP4_AfterMVA;
  delete m_selectedSecondJetP4_AfterMVA;
  delete m_hadTopP4_AfterMVA;
  delete m_lepTopP4_AfterMVA;
}

void MVASortingAlgorithm::work() {
  double maxMVAValue = std::numeric_limits<double>::lowest();

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

          float mvaValue = 0;

          Jet currentLeptonicBJet = m_jets[bj1];
          Jet currentHadronicBJet = m_jets[bj2];
          Jet currentHadronicLightJet1 = m_jets[j3];
          Jet currentHadronicLightJet2 = m_jets[j4];

          LorentzVector leptonic_W = m_lepton.p + m_neutrino;
          LorentzVector hadronic_W = currentHadronicLightJet1.p + currentHadronicLightJet2.p;

          LorentzVector leptonic_Top = leptonic_W + currentLeptonicBJet.p;
          LorentzVector hadronic_Top = hadronic_W + currentHadronicBJet.p;

          m_mva_lightJet1p2_Pt = currentHadronicLightJet1.p.pt() + currentHadronicLightJet2.p.pt();

          m_mva_leptonic_B_Pt = currentLeptonicBJet.p.pt();
          m_mva_hadronic_B_Pt = currentHadronicBJet.p.pt();

          m_mva_leptonic_W_Pt = leptonic_W.pt();
          //m_mva_leptonic_W_M = leptonic_W.M();
          m_mva_leptonic_Top_Pt = leptonic_Top.pt();
          m_mva_leptonic_Top_M = leptonic_Top.M();

          m_mva_hadronic_W_Pt = hadronic_W.pt();
          m_mva_hadronic_W_M = hadronic_W.M();
          m_mva_hadronic_Top_Pt = hadronic_Top.pt();
          m_mva_hadronic_Top_M = hadronic_Top.M();

          m_mva_delta_phi_tops = ROOT::Math::VectorUtil::DeltaPhi(hadronic_Top, leptonic_Top);
          m_mva_delta_phi_W = ROOT::Math::VectorUtil::DeltaPhi(hadronic_W, leptonic_W);
          m_mva_delta_phi_lightjets = ROOT::Math::VectorUtil::DeltaPhi(currentHadronicLightJet2.p, currentHadronicLightJet1.p);

          m_mva_delta_R_tops = ROOT::Math::VectorUtil::DeltaR(hadronic_Top, leptonic_Top);
          m_mva_delta_R_W = ROOT::Math::VectorUtil::DeltaR(hadronic_W, leptonic_W);
          m_mva_delta_R_lightjets = ROOT::Math::VectorUtil::DeltaR(currentHadronicLightJet2.p, currentHadronicLightJet1.p);

          m_mva_ht_fraction = (currentLeptonicBJet.p.pt() + currentHadronicBJet.p.pt() + currentHadronicLightJet1.p.pt() + currentHadronicLightJet2.p.pt()) / (allJetsPt);

          mvaValue = m_MVAReader->EvaluateMVA(m_MVAMethodName);

          if (mvaValue > maxMVAValue) {
            maxMVAValue = mvaValue;

            m_mtt_SolMVA[m_mtt_NumComb_MVA] = mvaValue;
            m_mtt_NumComb_MVA++;

            leptonicBJet = m_jets[bj1];
            hadronicBJet = m_jets[bj2];
            hadronicLightJet1 = m_jets[j3];
            hadronicLightJet2 = m_jets[j4];
          }
        } // j4 close
      } // j3 close
    } // j2 close
  } // we out of the loop over jet pairings

  if (m_mtt_NumComb_MVA > 0) {

    if (m_MVACut) {
      m_mtt_isSelMVA = (maxMVAValue < m_MVACut) ? 0 : 1;
    } else {
      m_mtt_isSelMVA = 1;
    }

    m_selectedHadronicBIndex_AfterMVA = hadronicBJet.index;
    m_selectedLeptonicBIndex_AfterMVA = leptonicBJet.index;
    m_selectedHadronicFirstJetIndex_AfterMVA = hadronicLightJet1.index;
    m_selectedHadronicSecondJetIndex_AfterMVA = hadronicLightJet2.index;

    // Correct MET again
    bool no_real_sol = false;
    computeNeutrinoPz(leptonicBJet.p, &no_real_sol);

    m_neutrino_no_real_solution_AfterMVA = (no_real_sol) ? 1 : 0;

    // Copy P4 for Tree
    *m_selectedNeutrinoP4_AfterMVA = m_neutrino;
    *m_selectedLeptonicBP4_AfterMVA = leptonicBJet.p;
    *m_selectedHadronicBP4_AfterMVA = hadronicBJet.p;
    *m_selectedFirstJetP4_AfterMVA = hadronicLightJet1.p;
    *m_selectedSecondJetP4_AfterMVA = hadronicLightJet2.p;

    /**
     * Compute Mtt before doing KinFit
     */
    m_mLepW_AfterMVA = (m_neutrino + m_lepton.p).M();
    m_mHadW_AfterMVA = (hadronicLightJet1.p + hadronicLightJet2.p).M();

    LorentzVector lepTopP4_AfterMVA = (m_lepton.p + m_neutrino + leptonicBJet.p);
    m_mLepTop_AfterMVA = lepTopP4_AfterMVA.M();
    m_lepTopPt_AfterMVA = lepTopP4_AfterMVA.pt();
    m_lepTopEta_AfterMVA = lepTopP4_AfterMVA.Eta();

    LorentzVector hadTopP4_AfterMVA = (hadronicLightJet1.p + hadronicLightJet2.p + hadronicBJet.p);
    m_mHadTop_AfterMVA = hadTopP4_AfterMVA.M();
    m_hadTopPt_AfterMVA = hadTopP4_AfterMVA.pt();
    m_hadTopEta_AfterMVA = hadTopP4_AfterMVA.Eta();

    LorentzVector res = (lepTopP4_AfterMVA + hadTopP4_AfterMVA);
    m_mtt_AfterMVA = res.M();
    m_pt_tt_AfterMVA = res.pt();
    m_eta_tt_AfterMVA = res.Eta();
    m_beta_tt_AfterMVA = fabs(res.Pz() / res.E());

    m_mtt_BestSolMVA = maxMVAValue;
  }
}

void MVASortingAlgorithm::reset() {
  m_mtt_isSelMVA = 0;
  m_mtt_recoJetsAssociatedWithMVA = false;
  m_mtt_recoJetsAssociatedWellPlacedWithMVA = false;

  m_mtt_BestSolMVA = -1.;
  m_mtt_NumComb_MVA = 0;

  m_mtt_AfterMVA = -1.;
  m_mtt_resolution_AfterMVA = -1.;
  m_pt_tt_AfterMVA = -1.;
  m_eta_tt_AfterMVA = -1.;
  m_beta_tt_AfterMVA = -1.;
  m_mLepTop_AfterMVA = -1.;
  m_mHadTop_AfterMVA = -1.;
  m_mHadW_AfterMVA = -1.;

  m_lepTopPt_AfterMVA = -1;
  m_lepTopEta_AfterMVA = -1;
  m_hadTopPt_AfterMVA = -1;
  m_hadTopEta_AfterMVA = -1;

  m_selectedLeptonicBIndex_AfterMVA = -1;
  m_selectedHadronicBIndex_AfterMVA = -1;
  m_selectedHadronicFirstJetIndex_AfterMVA = -1;
  m_selectedHadronicSecondJetIndex_AfterMVA = -1;

  m_selectedNeutrinoP4_AfterMVA->SetCoordinates(0, 0, 0, 0);
  m_selectedLeptonicBP4_AfterMVA->SetCoordinates(0, 0, 0, 0);
  m_selectedHadronicBP4_AfterMVA->SetCoordinates(0, 0, 0, 0);
  m_selectedFirstJetP4_AfterMVA->SetCoordinates(0, 0, 0, 0);
  m_selectedSecondJetP4_AfterMVA->SetCoordinates(0, 0, 0, 0);
  m_lepTopP4_AfterMVA->SetCoordinates(0, 0, 0, 0);
  m_hadTopP4_AfterMVA->SetCoordinates(0, 0, 0, 0);
}

void MVASortingAlgorithm::addBranches(TTree& tree) {
  tree.Branch("numComb_MVA" , &m_mtt_NumComb_MVA , "numComb_MVA/I");
  tree.Branch("solMVA" , &m_mtt_SolMVA , "solMVA[numComb_MVA]/F");
  tree.Branch("bestSolMVA" , &m_mtt_BestSolMVA , "bestSolMVA/F");
  tree.Branch("isSelMVA" , &m_mtt_isSelMVA , "isSelMVA/I");
  // If true, it means that we have selected the correct four reco jets (ie, reco jets coming from tt decay)
  tree.Branch("recoJetsAssociatedWithMVA" , &m_mtt_recoJetsAssociatedWithMVA , "recoJetsAssociatedWithMVA/O");
  // If true, it means that we have selected the correct four reco jets (ie, reco jets coming from tt decay)
  // and each jets is correctly positionned.
  tree.Branch("recoJetsAssociatedWellPlacedWithMVA", &m_mtt_recoJetsAssociatedWellPlacedWithMVA, "recoJetsAssociatedWellPlacedWithMVA/O");
  // Selected object with MVA
  tree.Branch("mLepW_AfterMVA" , &m_mLepW_AfterMVA , "mLepW_AfterMVA/F");
  tree.Branch("mHadW_AfterMVA" , &m_mHadW_AfterMVA , "mHadW_AfterMVA/F");

  tree.Branch("mLepTop_AfterMVA" , &m_mLepTop_AfterMVA , "mLepTop_AfterMVA/F");
  tree.Branch("lepTopPt_AfterMVA" , &m_lepTopPt_AfterMVA , "lepTopPt_AfterMVA/F");
  tree.Branch("lepTopEta_AfterMVA" , &m_lepTopEta_AfterMVA , "lepTopEta_AfterMVA/F");
  tree.Branch("lepTopP4_AfterMVA" , &m_lepTopP4_AfterMVA);

  tree.Branch("mHadTop_AfterMVA" , &m_mHadTop_AfterMVA , "mHadTop_AfterMVA/F");
  tree.Branch("hadTopPt_AfterMVA" , &m_hadTopPt_AfterMVA , "hadTopPt_AfterMVA/F");
  tree.Branch("hadTopEta_AfterMVA" , &m_hadTopEta_AfterMVA , "hadTopEta_AfterMVA/F");
  tree.Branch("hadTopP4_AfterMVA" , &m_hadTopP4_AfterMVA);

  tree.Branch("pt_tt_AfterMVA" , &m_pt_tt_AfterMVA , "pt_tt_AfterMVA/F");
  tree.Branch("eta_tt_AfterMVA" , &m_eta_tt_AfterMVA , "eta_tt_AfterMVA/F");
  tree.Branch("beta_tt_AfterMVA" , &m_beta_tt_AfterMVA , "eta_tt_AfterMVA/F");
  tree.Branch("mtt_AfterMVA" , &m_mtt_AfterMVA , "mtt_AfterMVA/F");
  tree.Branch("mttResolution_AfterMVA" , &m_mtt_resolution_AfterMVA , "mttResolution_AfterMVA/F");

  tree.Branch("neutrino_no_real_solution_AfterMVA", &m_neutrino_no_real_solution_AfterMVA, "neutrino_no_real_solution_AfterMVA/I");

  // Index of selected particles inside respective collection for mtt computation
  tree.Branch("selectedLeptonicBIndex_AfterMVA" , &m_selectedLeptonicBIndex_AfterMVA , "selectedLeptonicBIndex_AfterMVA/I");
  tree.Branch("selectedHadronicBIndex_AfterMVA" , &m_selectedHadronicBIndex_AfterMVA , "selectedHadronicBIndex_AfterMVA/I");
  tree.Branch("selectedHadronicFirstJetIndex_AfterMVA" , &m_selectedHadronicFirstJetIndex_AfterMVA , "selectedHadronicFirstJetIndex_AfterMVA/I");
  tree.Branch("selectedHadronicSecondJetIndex_AfterMVA" , &m_selectedHadronicSecondJetIndex_AfterMVA , "selectedHadronicSecondJetIndex_AfterMVA/I");

  tree.Branch("selectedNeutrinoP4_AfterMVA" , &m_selectedNeutrinoP4_AfterMVA);
  tree.Branch("selectedLeptonicBP4_AfterMVA" , &m_selectedLeptonicBP4_AfterMVA);
  tree.Branch("selectedHadronicBP4_AfterMVA" , &m_selectedHadronicBP4_AfterMVA);
  tree.Branch("selectedFirstJetP4_AfterMVA" , &m_selectedFirstJetP4_AfterMVA);
  tree.Branch("selectedSecondJetP4_AfterMVA" , &m_selectedSecondJetP4_AfterMVA);
}

