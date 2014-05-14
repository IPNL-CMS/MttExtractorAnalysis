#pragma once

#include "SortingAlgorithm.h"
#include <FWCore/ParameterSet/interface/ParameterSet.h>

class Chi2SortingAlgorithm: public SortingAlgorithm {
  public:
    Chi2SortingAlgorithm(const edm::ParameterSet& cfg, bool isSemiMu);
    ~Chi2SortingAlgorithm();
    virtual void work();
    virtual void reset();
    virtual void addBranches(TTree& tree);
    virtual int getSelectedLeptonicBIndex() {
        return m_selectedLeptonicBIndex_AfterChi2;
    };
    virtual int getSelectedHadronicBIndex() {
        return m_selectedHadronicBIndex_AfterChi2;
    };
    virtual int getSelectedHadronicFirstJetIndex() {
        return m_selectedHadronicFirstJetIndex_AfterChi2;
    };
    virtual int getSelectedHadronicSecondJetIndex() {
        return m_selectedHadronicSecondJetIndex_AfterChi2;
    };
    virtual void setRecoJetsAssociated(const bool& recoJetsAssociated) {
        m_mtt_recoJetsAssociatedWithChi2 = recoJetsAssociated;
    };
    virtual void setRecoJetsAssociatedWellPlaced(const bool& recoJetsAssociatedWellPlaced) {
        m_mtt_recoJetsAssociatedWellPlacedWithChi2 = recoJetsAssociatedWellPlaced;
    };

  private:
    double chi2(const LorentzVector& leptonicBJet, const LorentzVector& hadronicBJet, const LorentzVector& hadronicLightJet1, const LorentzVector& hadronicLightJet2, double allJetsPt);

    bool m_useBTagInCombinatorics;
    bool m_isSemiMu;

    bool m_mtt_recoJetsAssociatedWithChi2;
    bool m_mtt_recoJetsAssociatedWellPlacedWithChi2;

    int m_mtt_NumComb_chi2;
    float m_mtt_SolChi2[1000];
    float m_mtt_BestSolChi2;

    // Values after Chi2 selection
    float m_mtt_AfterChi2;
    float m_eta_tt_AfterChi2;
    float m_pt_tt_AfterChi2;
    float m_beta_tt_AfterChi2;
    float m_mLepTop_AfterChi2;
    float m_mHadTop_AfterChi2;
    float m_mHadW_AfterChi2;
    float m_mLepW_AfterChi2;

    int   m_neutrino_no_real_solution_AfterChi2;

    float m_lepTopPt_AfterChi2;
    float m_lepTopEta_AfterChi2;
    LorentzVector* m_lepTopP4_AfterChi2;

    float m_hadTopPt_AfterChi2;
    float m_hadTopEta_AfterChi2;
    LorentzVector* m_hadTopP4_AfterChi2;

    int m_selectedLeptonicBIndex_AfterChi2;
    int m_selectedHadronicBIndex_AfterChi2;
    int m_selectedHadronicFirstJetIndex_AfterChi2;
    int m_selectedHadronicSecondJetIndex_AfterChi2;

    LorentzVector* m_selectedNeutrinoP4_AfterChi2;
    LorentzVector* m_selectedLeptonicBP4_AfterChi2;
    LorentzVector* m_selectedHadronicBP4_AfterChi2;
    LorentzVector* m_selectedFirstJetP4_AfterChi2;
    LorentzVector* m_selectedSecondJetP4_AfterChi2;

    /// constant masses
    double m_w;  
    double m_b;  
    double m_top;

    // Chi2
    double chi2_hadronic_top_mass;
    double chi2_leptonic_top_mass_semimu;
    double chi2_leptonic_top_mass_semie;
    double chi2_hadronic_w_mass;
    double chi2_pt_ttbar_system;
    double chi2_ht_frac;

    double chi2_sigma_hadronic_top_mass;
    double chi2_sigma_leptonic_top_mass_semimu;
    double chi2_sigma_leptonic_top_mass_semie;
    double chi2_sigma_hadronic_w_mass;
    double chi2_sigma_pt_ttbar_system;
    double chi2_sigma_ht_frac;

    double chi2_sigma_hadronic_top_mass_square;
    double chi2_sigma_leptonic_top_mass_semimu_square;
    double chi2_sigma_leptonic_top_mass_semie_square;
    double chi2_sigma_hadronic_w_mass_square;
    double chi2_sigma_pt_ttbar_system_square;
    double chi2_sigma_ht_frac_square;

    bool m_usePtSystInChi2;
    bool m_useHtFracInChi2;

};
