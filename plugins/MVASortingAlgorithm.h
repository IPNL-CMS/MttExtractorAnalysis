#pragma once

#include "SortingAlgorithm.h"
#include <FWCore/ParameterSet/interface/ParameterSet.h>

namespace TMVA {
  class Reader;
}

class MVASortingAlgorithm: public SortingAlgorithm {
  public:
    MVASortingAlgorithm(const edm::ParameterSet& cfg);
    ~MVASortingAlgorithm();
    virtual void work();
    virtual void reset();
    virtual void addBranches(TTree& tree);

  private:
    bool m_useBTagInCombinatorics;

    std::string m_MVAWeightFilename;
    std::string m_MVAMethodName;
    bool m_MVACut;
    double m_MVACutValue;
    std::shared_ptr<TMVA::Reader> m_MVAReader;

    float m_mva_lightJet1p2_Pt;
    float m_mva_leptonic_B_Pt;
    float m_mva_leptonic_W_Pt;
    float m_mva_leptonic_Top_Pt;
    float m_mva_leptonic_Top_M;
    float m_mva_hadronic_B_Pt;
    float m_mva_hadronic_W_Pt;
    float m_mva_hadronic_W_M;
    float m_mva_hadronic_Top_Pt;
    float m_mva_hadronic_Top_M;

    float m_mva_delta_phi_tops;
    float m_mva_delta_phi_lightjets;
    float m_mva_delta_phi_W;
    float m_mva_delta_R_tops;
    float m_mva_delta_R_lightjets;
    float m_mva_delta_R_W;

    float m_mva_ht_fraction;

    int m_mtt_NumComb_MVA;
    float m_mtt_SolMVA[1000];
    float m_mtt_BestSolMVA;

    // MVA values
    int m_mtt_isSelMVA;
    bool m_mtt_recoJetsAssociatedWithMVA;
    bool m_mtt_recoJetsAssociatedWellPlacedWithMVA;

    float m_mtt_AfterMVA;
    float m_mtt_resolution_AfterMVA;
    float m_eta_tt_AfterMVA;
    float m_pt_tt_AfterMVA;
    float m_beta_tt_AfterMVA;
    float m_mLepTop_AfterMVA;
    float m_mHadTop_AfterMVA;
    float m_mHadW_AfterMVA;
    float m_mLepW_AfterMVA;

    int m_neutrino_no_real_solution_AfterMVA;

    float m_lepTopPt_AfterMVA;
    float m_lepTopEta_AfterMVA;
    LorentzVector* m_lepTopP4_AfterMVA;

    float m_hadTopPt_AfterMVA;
    float m_hadTopEta_AfterMVA;
    LorentzVector* m_hadTopP4_AfterMVA;

    // Indexes of selected particles for mtt computation, after MVA selection
    int m_selectedLeptonicBIndex_AfterMVA;
    int m_selectedHadronicBIndex_AfterMVA;
    int m_selectedHadronicFirstJetIndex_AfterMVA;
    int m_selectedHadronicSecondJetIndex_AfterMVA;

    LorentzVector* m_selectedNeutrinoP4_AfterMVA;
    LorentzVector* m_selectedLeptonicBP4_AfterMVA;
    LorentzVector* m_selectedHadronicBP4_AfterMVA;
    LorentzVector* m_selectedFirstJetP4_AfterMVA;
    LorentzVector* m_selectedSecondJetP4_AfterMVA;
};
