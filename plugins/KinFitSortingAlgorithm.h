#pragma once

#include "SortingAlgorithm.h"
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include "TopQuarkAnalysis/TopKinFitter/interface/TtSemiLepKinFitter.h"
#include "Extractors/MttExtractorAnalysis/plugins/GaussianProfile.h"

class KinFitSortingAlgorithm: public SortingAlgorithm {
  public:
    KinFitSortingAlgorithm(const edm::ParameterSet& cfg, bool isSemiMu);
    ~KinFitSortingAlgorithm();

    virtual void work();
    virtual void reset();
    virtual void addBranches(TTree& tree);
    virtual void end() {
      m_mtt_gen_vs_mtt_reco_linearity.write(); 
      m_mtt_gen_vs_mtt_reco_resolution.write(); 
    };

    virtual int getSelectedLeptonicBIndex() {
      return m_selectedLeptonicBIndex_AfterKF;
    };

    virtual int getSelectedHadronicBIndex() {
      return m_selectedHadronicBIndex_AfterKF;
    };

    virtual int getSelectedHadronicFirstJetIndex() {
      return m_selectedHadronicFirstJetIndex_AfterKF;
    };

    virtual int getSelectedHadronicSecondJetIndex() {
      return m_selectedHadronicSecondJetIndex_AfterKF;
    };

    virtual void setRecoJetsAssociated(const bool& recoJetsAssociated) {
      m_mtt_recoJetsAssociatedWithKF = recoJetsAssociated;
    };

    virtual void setRecoJetsAssociatedWellPlaced(const bool& recoJetsAssociatedWellPlaced) {
      m_mtt_recoJetsAssociatedWellPlacedWithKF = recoJetsAssociatedWellPlaced;
    };

  private:

    template<class T> static inline TLorentzVector toTLorentzVector(const T& v) {
      return TLorentzVector(v.Px(), v.Py(), v.Pz(), v.E());
    }

    GaussianProfile m_mtt_gen_vs_mtt_reco_linearity;
    GaussianProfile m_mtt_gen_vs_mtt_reco_resolution;

    TtSemiLepKinFitter::Param param(unsigned val) const;
    TtSemiLepKinFitter::Constraint constraint(unsigned val) const;
    std::vector<TtSemiLepKinFitter::Constraint> constraints(std::vector<unsigned>& val) const;

    /// maximal number of iterations to be performed for the fit
    unsigned int m_maxNrIter;
    /// maximal chi2 equivalent
    double m_maxDeltaS;
    /// maximal deviation for contstraints
    double m_maxF;
    unsigned int m_jetParam;
    unsigned int m_lepParam;
    unsigned int m_metParam;
    /// constrains
    std::vector<unsigned> m_constraints;
    double m_mW;
    double m_mTop;
    /// scale factors for jet energy resolution
    std::vector<double> m_jetEnergyResolutionScaleFactors;
    std::vector<double> m_jetEnergyResolutionEtaBinning;
    /// config-file-based object resolutions
    std::vector<edm::ParameterSet> m_udscResolutions;
    std::vector<edm::ParameterSet> m_bResolutions;
    std::vector<edm::ParameterSet> m_lepResolutions;
    std::vector<edm::ParameterSet> m_metResolutions;

    std::shared_ptr<TtSemiLepKinFitter> m_fitter;

    bool m_useBTagInCombinatorics;
    bool m_isSemiMu;

    bool m_mtt_recoJetsAssociatedWithKF;
    bool m_mtt_recoJetsAssociatedWellPlacedWithKF;

    int m_mtt_NumComb_kf;
    float m_mtt_BestSolKF;

    // Values after KF selection
    float m_mtt_AfterKF;
    float m_eta_tt_AfterKF;
    float m_pt_tt_AfterKF;
    float m_beta_tt_AfterKF;
    float m_mLepTop_AfterKF;
    float m_mHadTop_AfterKF;
    float m_mHadW_AfterKF;
    float m_mLepW_AfterKF;

    float m_lepTopPt_AfterKF;
    float m_lepTopEta_AfterKF;
    LorentzVector* m_lepTopP4_AfterKF;

    float m_hadTopPt_AfterKF;
    float m_hadTopEta_AfterKF;
    LorentzVector* m_hadTopP4_AfterKF;

    float m_kf_chisquare;
    float m_kf_proba;
    bool m_kf_converged;

    int m_selectedLeptonicBIndex_AfterKF;
    int m_selectedHadronicBIndex_AfterKF;
    int m_selectedHadronicFirstJetIndex_AfterKF;
    int m_selectedHadronicSecondJetIndex_AfterKF;
    
    LorentzVector* m_selectedNeutrinoP4_AfterKF;
    LorentzVector* m_selectedLeptonP4_AfterKF;
    LorentzVector* m_selectedLeptonicBP4_AfterKF;
    LorentzVector* m_selectedHadronicBP4_AfterKF;
    LorentzVector* m_selectedFirstJetP4_AfterKF;
    LorentzVector* m_selectedSecondJetP4_AfterKF;

};

