#include "KinFitSortingAlgorithm.h"
#include <TopQuarkAnalysis/TopKinFitter/interface/CovarianceMatrix.h>

#include <TTree.h>

KinFitSortingAlgorithm::KinFitSortingAlgorithm(const edm::ParameterSet& cfg, bool isSemiMu) :
    m_maxNrIter               (cfg.getParameter<unsigned>     ("maxNrIter"           )),
    m_maxDeltaS               (cfg.getParameter<double>       ("maxDeltaS"           )),
    m_maxF                    (cfg.getParameter<double>       ("maxF"                )),
    m_jetParam                (cfg.getParameter<unsigned>     ("jetParametrisation"  )),
    m_lepParam                (cfg.getParameter<unsigned>     ("lepParametrisation"  )),
    m_metParam                (cfg.getParameter<unsigned>     ("metParametrisation"  )),
    m_constraints             (cfg.getParameter<std::vector<unsigned> >("constraints")),
    m_mW                      (cfg.getParameter<double>       ("mW"                  )),
    m_mTop                    (cfg.getParameter<double>       ("mTop"                )),
    m_jetEnergyResolutionScaleFactors(cfg.getParameter<std::vector<double> >("jetEnergyResolutionScaleFactors")),
    m_jetEnergyResolutionEtaBinning(cfg.getParameter<std::vector<double> >("jetEnergyResolutionEtaBinning")),
    m_udscResolutions(0), m_bResolutions(0), m_lepResolutions(0), m_metResolutions(0),
    m_mtt_gen_vs_mtt_reco_linearity("mtt_gen_vs_mtt_reco_linearity_AfterKF", nBins_for_gaussian_profile, bins_for_gaussian_profile),
    m_mtt_gen_vs_mtt_reco_resolution("mtt_gen_vs_mtt_reco_resolution_AfterKF", nBins_for_gaussian_profile, bins_for_gaussian_profile, 100, -1.2, 1.2)
{
  m_isSemiMu = isSemiMu;

  m_useBTagInCombinatorics = cfg.getParameter<bool>("use_btag_in_combinatorics");

  if (cfg.exists("udscResolutions") && cfg.exists("bResolutions") && cfg.exists("lepResolutions") && cfg.exists("metResolutions")) {
    m_udscResolutions = cfg.getParameter<std::vector<edm::ParameterSet> >("udscResolutions");
    m_bResolutions    = cfg.getParameter<std::vector<edm::ParameterSet> >("bResolutions"   );
    m_lepResolutions  = cfg.getParameter<std::vector<edm::ParameterSet> >("lepResolutions" );
    m_metResolutions  = cfg.getParameter<std::vector<edm::ParameterSet> >("metResolutions" );
  }
  else if (cfg.exists("udscResolutions") || cfg.exists("bResolutions") || cfg.exists("lepResolutions") || cfg.exists("metResolutions")) {
    throw cms::Exception("Configuration") << "Parameters 'udscResolutions', 'bResolutions', 'lepResolutions', 'metResolutions' should be used together.\n";
  }
 
  m_fitter.reset(new TtSemiLepKinFitter(param(m_jetParam), param(m_lepParam), param(m_metParam), m_maxNrIter, m_maxDeltaS, m_maxF,
              constraints(m_constraints), m_mW, m_mTop, &m_udscResolutions, &m_bResolutions, &m_lepResolutions, &m_metResolutions,
              &m_jetEnergyResolutionScaleFactors, &m_jetEnergyResolutionEtaBinning));

  m_selectedNeutrinoP4_AfterKF = new LorentzVector();
  m_selectedLeptonP4_AfterKF = new LorentzVector();
  m_selectedLeptonicBP4_AfterKF = new LorentzVector();
  m_selectedHadronicBP4_AfterKF = new LorentzVector();
  m_selectedFirstJetP4_AfterKF = new LorentzVector();
  m_selectedSecondJetP4_AfterKF = new LorentzVector();
  m_hadTopP4_AfterKF = new LorentzVector();
  m_lepTopP4_AfterKF = new LorentzVector();
}

KinFitSortingAlgorithm::~KinFitSortingAlgorithm() {
  delete m_selectedNeutrinoP4_AfterKF;
  delete m_selectedLeptonP4_AfterKF;
  delete m_selectedLeptonicBP4_AfterKF;
  delete m_selectedHadronicBP4_AfterKF;
  delete m_selectedFirstJetP4_AfterKF;
  delete m_selectedSecondJetP4_AfterKF;
  delete m_hadTopP4_AfterKF;
  delete m_lepTopP4_AfterKF;
}

TtSemiLepKinFitter::Param KinFitSortingAlgorithm::param(unsigned val) const {
  TtSemiLepKinFitter::Param result;
  switch (val) {
    case TtSemiLepKinFitter::kEMom:
      result = TtSemiLepKinFitter::kEMom;
      break;
    case TtSemiLepKinFitter::kEtEtaPhi:
      result = TtSemiLepKinFitter::kEtEtaPhi;
      break;
    case TtSemiLepKinFitter::kEtThetaPhi:
      result = TtSemiLepKinFitter::kEtThetaPhi;
      break;
    default:
      throw cms::Exception("Configuration")  << "Chosen jet parametrization is not supported: " << val << std::endl;
      break;
  }

  return result;
}

TtSemiLepKinFitter::Constraint KinFitSortingAlgorithm::constraint(unsigned val) const {
  TtSemiLepKinFitter::Constraint result;
  switch (val) {
    case TtSemiLepKinFitter::kWHadMass:
      result = TtSemiLepKinFitter::kWHadMass;
      break;
    case TtSemiLepKinFitter::kWLepMass:
      result = TtSemiLepKinFitter::kWLepMass;
      break;
    case TtSemiLepKinFitter::kTopHadMass:
      result = TtSemiLepKinFitter::kTopHadMass;
      break;
    case TtSemiLepKinFitter::kTopLepMass:
      result = TtSemiLepKinFitter::kTopLepMass;
      break;
    case TtSemiLepKinFitter::kNeutrinoMass:
      result = TtSemiLepKinFitter::kNeutrinoMass;
      break;
    case TtSemiLepKinFitter::kEqualTopMasses:
      result = TtSemiLepKinFitter::kEqualTopMasses;
      break;
    case TtSemiLepKinFitter::kSumPt:
      result = TtSemiLepKinFitter::kSumPt;
      break;

    default:
      throw cms::Exception("Configuration")  << "Chosen fit constraint is not supported: " << val << std::endl;
      break;
  }

  return result;
}

std::vector<TtSemiLepKinFitter::Constraint> KinFitSortingAlgorithm::constraints(std::vector<unsigned>& val) const
{
  std::vector<TtSemiLepKinFitter::Constraint> result;
  for (unsigned i = 0; i < val.size(); ++i) {
    result.push_back(constraint(val[i]));
  }

  return result; 
}

void KinFitSortingAlgorithm::work() {
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

  TtSemiLepKinFitter& fitter = *m_fitter;
  fitter.setVerbosity(0);

  double minChiSquare = std::numeric_limits<double>::infinity();

  TLorentzVector leptonP4 = toTLorentzVector(m_lepton.p);
  int leptonCharge = m_lepton.charge;

  Jet leptonicBJet;
  Jet hadronicBJet;
  Jet hadronicLightJet1;
  Jet hadronicLightJet2;

  LorentzVector lepton;
  LorentzVector neutrino;

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

          if (!computeNeutrinoPz(m_jets[bj1].p))
            m_neutrino.SetPz(0.);

          TLorentzVector neutrinoP4 = toTLorentzVector(m_neutrino);
          TLorentzVector currentLeptonicBJetP4 = toTLorentzVector(m_jets[bj1].p);
          TLorentzVector currentHadronicBJetP4 = toTLorentzVector(m_jets[bj2].p);
          TLorentzVector currentHadronicLightJet1P4 = toTLorentzVector(m_jets[j3].p);
          TLorentzVector currentHadronicLightJet2P4 = toTLorentzVector(m_jets[j4].p);

          if (fitter.fit(currentHadronicLightJet1P4, currentHadronicLightJet2P4, currentHadronicBJetP4, currentLeptonicBJetP4, leptonP4, neutrinoP4, leptonCharge, (! m_isSemiMu) ? CovarianceMatrix::kElectron : CovarianceMatrix::kMuon) != 0)
            continue;

          double chi2 = fitter.fitS();

          if (chi2 < minChiSquare) {
            minChiSquare = chi2;
            m_kf_converged = true;
            m_mtt_NumComb_kf++;

            leptonicBJet = m_jets[bj1];
            hadronicBJet = m_jets[bj2];
            hadronicLightJet1 = m_jets[j3];
            hadronicLightJet2 = m_jets[j4];

            leptonicBJet.p = fitter.fittedLepB().p4();
            hadronicBJet.p = fitter.fittedHadB().p4();
            hadronicLightJet1.p =fitter.fittedHadP().p4();
            hadronicLightJet2.p = fitter.fittedHadQ().p4();

            lepton = fitter.fittedLepton().p4();
            neutrino = fitter.fittedNeutrino().p4();

            m_kf_proba = fitter.fitProb();
            m_kf_chisquare = chi2;
          }
        } // j4 close
      } // j3 close
    } // j2 close
  } // we out of the loop over jet pairings

  if (m_mtt_NumComb_kf > 0) {
    m_selectedHadronicBIndex_AfterKF = hadronicBJet.index;
    m_selectedLeptonicBIndex_AfterKF = leptonicBJet.index;
    m_selectedHadronicFirstJetIndex_AfterKF = hadronicLightJet1.index;
    m_selectedHadronicSecondJetIndex_AfterKF = hadronicLightJet2.index;

    // Copy P4 for Tree
    *m_selectedNeutrinoP4_AfterKF = neutrino;
    *m_selectedLeptonP4_AfterKF = lepton;
    *m_selectedLeptonicBP4_AfterKF = leptonicBJet.p;
    *m_selectedHadronicBP4_AfterKF = hadronicBJet.p;
    *m_selectedFirstJetP4_AfterKF = hadronicLightJet1.p;
    *m_selectedSecondJetP4_AfterKF = hadronicLightJet2.p;

    /**
     * Compute Mtt before doing KinFit
     */
    m_mLepW_AfterKF   = (neutrino + lepton).M();
    m_mHadW_AfterKF   = (hadronicLightJet1.p + hadronicLightJet2.p).M();

    LorentzVector lepTopP4_AfterKF = (lepton + neutrino + leptonicBJet.p);
    *m_lepTopP4_AfterKF = lepTopP4_AfterKF;
    m_mLepTop_AfterKF = lepTopP4_AfterKF.M();
    m_lepTopPt_AfterKF = lepTopP4_AfterKF.Pt();
    m_lepTopEta_AfterKF = lepTopP4_AfterKF.Eta();

    LorentzVector hadTopP4_AfterKF = (hadronicLightJet1.p + hadronicLightJet2.p + hadronicBJet.p);
    *m_hadTopP4_AfterKF = (hadronicLightJet1.p + hadronicLightJet2.p + hadronicBJet.p);
    m_mHadTop_AfterKF = hadTopP4_AfterKF.M();
    m_hadTopPt_AfterKF = hadTopP4_AfterKF.Pt();
    m_hadTopEta_AfterKF = hadTopP4_AfterKF.Eta();

    LorentzVector res = (lepTopP4_AfterKF + hadTopP4_AfterKF);
    m_mtt_AfterKF     = res.M();
    m_pt_tt_AfterKF   = res.Pt();
    m_eta_tt_AfterKF   = res.Eta();
    m_beta_tt_AfterKF  = fabs(res.Pz() / res.E());

    m_mtt_BestSolKF = minChiSquare;
  }

  m_mtt_gen_vs_mtt_reco_linearity.fill(m_mtt_gen, m_mtt_AfterKF);
  m_mtt_gen_vs_mtt_reco_resolution.fill(m_mtt_gen, (m_mtt_AfterKF - m_mtt_gen) / m_mtt_gen);
}

void KinFitSortingAlgorithm::reset() {
  m_mtt_recoJetsAssociatedWithKF = false;
  m_mtt_recoJetsAssociatedWellPlacedWithKF = false;
  m_mtt_BestSolKF = -1.;
  m_mtt_NumComb_kf = 0;

  m_mtt_AfterKF            = -1.;
  m_pt_tt_AfterKF          = -1.;
  m_eta_tt_AfterKF         = -1.;
  m_beta_tt_AfterKF        = -1.;
  m_mLepTop_AfterKF        = -1.;
  m_mHadTop_AfterKF        = -1.;
  m_mHadW_AfterKF          = -1.;

  m_lepTopPt_AfterKF  = -1;
  m_lepTopEta_AfterKF = -1;
  m_hadTopPt_AfterKF  = -1;
  m_hadTopEta_AfterKF = -1;

  m_kf_converged = false;
  m_kf_chisquare = -1;
  m_kf_proba = -1;

  m_selectedLeptonicBIndex_AfterKF            = -1;
  m_selectedHadronicBIndex_AfterKF            = -1;
  m_selectedHadronicFirstJetIndex_AfterKF     = -1;
  m_selectedHadronicSecondJetIndex_AfterKF    = -1;

  m_selectedNeutrinoP4_AfterKF->SetCoordinates(0, 0, 0, 0);
  m_selectedLeptonP4_AfterKF->SetCoordinates(0, 0, 0, 0);
  m_selectedLeptonicBP4_AfterKF->SetCoordinates(0, 0, 0, 0);
  m_selectedHadronicBP4_AfterKF->SetCoordinates(0, 0, 0, 0);
  m_selectedFirstJetP4_AfterKF->SetCoordinates(0, 0, 0, 0);
  m_selectedSecondJetP4_AfterKF->SetCoordinates(0, 0, 0, 0);
  m_lepTopP4_AfterKF->SetCoordinates(0, 0, 0, 0);
  m_hadTopP4_AfterKF->SetCoordinates(0, 0, 0, 0);
}

void KinFitSortingAlgorithm::addBranches(TTree& tree) {
  tree.Branch("numComb_kf"       , &m_mtt_NumComb_kf           , "numComb_kf/I");
  tree.Branch("bestSolKF"        , &m_mtt_BestSolKF           , "bestSolKF/F");
  // If true, it means that we have selected the correct four reco jets (ie, reco jets coming from tt decay)
  tree.Branch("recoJetsAssociatedWithKF" , &m_mtt_recoJetsAssociatedWithKF    , "recoJetsAssociatedWithKF/O");
  // If true, it means that we have selected the correct four reco jets (ie, reco jets coming from tt decay)
  // and each jets is correctly positionned.
  tree.Branch("recoJetsAssociatedWellPlacedWithKF", &m_mtt_recoJetsAssociatedWellPlacedWithKF, "recoJetsAssociatedWellPlacedWithKF/O");
  // Selected object with KF
  tree.Branch("mLepW_AfterKF"        , &m_mLepW_AfterKF       , "mLepW_AfterKF/F");
  tree.Branch("mHadW_AfterKF"        , &m_mHadW_AfterKF       , "mHadW_AfterKF/F");

  tree.Branch("mLepTop_AfterKF"      , &m_mLepTop_AfterKF     , "mLepTop_AfterKF/F");
  tree.Branch("lepTopPt_AfterKF"     , &m_lepTopPt_AfterKF    , "lepTopPt_AfterKF/F");
  tree.Branch("lepTopEta_AfterKF"    , &m_lepTopEta_AfterKF   , "lepTopEta_AfterKF/F");
  tree.Branch("lepTopP4_AfterKF"    , &m_lepTopP4_AfterKF);

  tree.Branch("mHadTop_AfterKF"      , &m_mHadTop_AfterKF     , "mHadTop_AfterKF/F");
  tree.Branch("hadTopPt_AfterKF"     , &m_hadTopPt_AfterKF    , "hadTopPt_AfterKF/F");
  tree.Branch("hadTopEta_AfterKF"    , &m_hadTopEta_AfterKF   , "hadTopEta_AfterKF/F");
  tree.Branch("hadTopP4_AfterKF"    , &m_hadTopP4_AfterKF);

  tree.Branch("pt_tt_AfterKF"        , &m_pt_tt_AfterKF       , "pt_tt_AfterKF/F");
  tree.Branch("eta_tt_AfterKF"       , &m_eta_tt_AfterKF      , "eta_tt_AfterKF/F");
  tree.Branch("beta_tt_AfterKF"      , &m_beta_tt_AfterKF     , "eta_tt_AfterKF/F");
  tree.Branch("mtt_AfterKF"          , &m_mtt_AfterKF         , "mtt_AfterKF/F");

  tree.Branch("selectedLeptonicBIndex_AfterKF"     , &m_selectedLeptonicBIndex_AfterKF    , "selectedLeptonicBIndex_AfterKF/I");
  tree.Branch("selectedHadronicBIndex_AfterKF"     , &m_selectedHadronicBIndex_AfterKF    , "selectedHadronicBIndex_AfterKF/I");
  tree.Branch("selectedHadronicFirstJetIndex_AfterKF"  , &m_selectedHadronicFirstJetIndex_AfterKF  , "selectedHadronicFirstJetIndex_AfterKF/I");
  tree.Branch("selectedHadronicSecondJetIndex_AfterKF" , &m_selectedHadronicSecondJetIndex_AfterKF , "selectedHadronicSecondJetIndex_AfterKF/I");

  tree.Branch("selectedNeutrinoP4_AfterKF"    , &m_selectedNeutrinoP4_AfterKF);
  tree.Branch("selectedLeptonP4_AfterKF"      , &m_selectedLeptonP4_AfterKF);
  tree.Branch("selectedLeptonicBP4_AfterKF"   , &m_selectedLeptonicBP4_AfterKF);
  tree.Branch("selectedHadronicBP4_AfterKF"   , &m_selectedHadronicBP4_AfterKF);
  tree.Branch("selectedFirstJetP4_AfterKF"    , &m_selectedFirstJetP4_AfterKF);
  tree.Branch("selectedSecondJetP4_AfterKF"   , &m_selectedSecondJetP4_AfterKF);

  tree.Branch("kf_chisquare", &m_kf_chisquare, "kf_chisquare/F");
  tree.Branch("kf_proba", &m_kf_proba, "kf_proba/F");
  tree.Branch("kf_converged", &m_kf_converged, "kf_converged/O");

}


