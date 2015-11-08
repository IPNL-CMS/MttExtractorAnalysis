#pragma once

/*****************************
 * chanel :
 *  0  -> no ttbar decays
 *  1  -> isSemiElectronic
 *  2  -> isSemiMuonic
 *  3  -> isSemiTauic
 *  4  -> isFullHadronic
 *  5  -> isDiElectronic
 *  6  -> isDiMuonic
 *  7  -> isDiTauic
 *  8  -> isElectroMuonic
 *  9  -> isElectroTauic
 *  10 -> isMuoTauic
 ******************************/

#include <vector>
#include <boost/regex.hpp>

#include <Extractors/PatExtractor/interface/ExtractorPlugin.h>
#include "Extractors/MttExtractorAnalysis/plugins/SortingAlgorithmFactory.h"

#include <TTree.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TRef.h>
#include <TH2.h>

class AnalysisSettings;
class EventExtractor;
class MuonExtractor;
class ElectronExtractor;
class METExtractor;
class VertexExtractor;
class KinFit;
class HLTExtractor;
class PatExtractor;
class MCExtractor;
class JetMETExtractor;
class BTaggingEfficiencyProvider;

namespace edm {
  class EventSetup;
}

namespace TMVA {
  class Reader;
}

namespace mtt {

  class mtt_analysis: public patextractor::Plugin {
    public:
      mtt_analysis(const edm::ParameterSet& cmsswSettings);
      ~mtt_analysis();

      // TTbar selection
      virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup, PatExtractor& extractor);
      virtual void analyze(const edm::EventSetup& iSetup, PatExtractor& extractor);
      void computeEventShapeVariables();
      virtual void endJob() {
        m_reco_bjets->Write(); // Number of tagged B in reco which are true B jets
        m_reco_fake_bjets_among_cjets->Write(); // Number of tagged B in reco which are in reality C jets
        m_reco_fake_bjets_among_lightjets->Write(); // Number of tagged B in reco which are in reality light jets

        m_gen_bjets->Write();
        m_gen_cjets->Write();
        m_gen_lightjets->Write();

        m_selected_jets_flavor->Write();
        m_selected_jets_flavor_btagged->Write();

        m_selected_light_jets_sf->Write();
        m_selected_b_jets_sf->Write();
        m_selected_c_jets_sf->Write();

        m_1b_sf_flavor->Write();
        m_1b_sf->Write();

        m_2b_sf_flavor->Write();
        m_2b_sf->Write();

        for (auto& algo: m_sorting_algortihms) {
          algo->end();
        }
      }

    private:

      //Selections
      int MuonSel();
      int ElectronSel();
      int JetSel();
      int VertexSel();
      int METSel();

      void loopOverCombinations();
      bool hasRecoPartner(int mcIndex) const;
      bool jetComesFromTTDecay(int mcIndex) const;
      void checkIfSolutionIsCorrect();
      bool isSolutionMatched(uint32_t leptonicBIndex, uint32_t hadronicBIndex, uint32_t hadronicFirstJetIndex, uint32_t hadronicSecondJetIndex);
      void MCidentification();
      void reset();
      void fillTree();

      bool isBJet(unsigned int index);
      float getBTagProbabilityMC();
      float getBTagProbabilityData();

      bool   m_MAIN_doSemiMu;

      TTree*  m_tree_Mtt;
      KinFit* m_KinFit;

      std::shared_ptr<EventExtractor> m_event;
      std::shared_ptr<MCExtractor>    m_MC;

      std::shared_ptr<VertexExtractor> m_vertex;

      float m_MET_Pt_Min;

      std::shared_ptr<MuonExtractor> m_muon;
      float m_MU_Iso_max;
      float m_MU_Pt_min;
      float m_MU_Eta_max;

      std::shared_ptr<MuonExtractor> m_muon_loose;
      float m_MU_Pt_min_loose;
      float m_MU_Eta_max_loose;
      float m_MU_Iso_max_loose;

      std::shared_ptr<ElectronExtractor> m_electron;
      float m_ELE_Pt_min;
      float m_ELE_Eta_max;
      float m_ELE_Iso_max;

      std::shared_ptr<ElectronExtractor> m_electron_loose;
      float m_ELE_Pt_min_loose;
      float m_ELE_Eta_max_loose;
      float m_ELE_Iso_max_loose;

      std::shared_ptr<JetMETExtractor> m_jetMet;
      uint32_t m_JET_N_sel_max;
      float m_JET_Pt_min;
      float m_JET_Eta_max;
      float m_JET_btag_TCHPT;
      float m_JET_btag_CSVL;
      float m_JET_btag_CSVM;
      float m_JET_btag_CSVT;

      TLorentzVector* m_refLept;
      int m_refLeptCharge;

      // Triggers
      bool m_trigger_passed;

      //MC stuff
      int m_MC_channel;
      float m_MC_mtt;
      int m_nPU; // number of interactions

      // Indexes of gen particle in a semi-lept event
      int m_leptonIndex;
      int m_neutrinoIndex;

      int m_leptonicBIndex;
      int m_hadronicBIndex;
      int m_leptonicTopIndex;

      int m_firstJetIndex;
      int m_secondJetIndex;

      float m_MC_hadronicWMass;
      float m_MC_leptonicWMass;
      float m_MC_hadronicTopMass;
      float m_MC_leptonicTopMass;
      float m_MC_pt_tt;
      float m_MC_eta_tt;
      float m_MC_beta_tt;

      TClonesArray* m_MC_lepton_p4;
      TClonesArray* m_MC_neutrino_p4;

      TClonesArray* m_MC_leptonic_B_p4;
      TClonesArray* m_MC_hadronic_B_p4;

      TClonesArray* m_MC_lightJet1_p4;
      TClonesArray* m_MC_lightJet2_p4;

      TClonesArray* m_MC_Top1_p4;
      TClonesArray* m_MC_Top2_p4;

      /// Number of lepton/neutrino from Top->W and quark b from Top
      int nEle;
      int nMu;
      int nTau;
      int nNuEle;
      int nNuMu;
      int nNuTau;
      int nQuarkb;
      int nTop;
      std::vector<TLorentzVector> Top;

      std::vector<int> m_selJetsIds;
      float AllJetsPt;

      int SelLeptIdx;
      //Reco stuff

      int m_mtt_isSel;
      bool m_mtt_eventIsAssociable; // If true, each parton from the event has a reco object associated.
      bool m_mtt_recoJetsAssociatedWithChi2;
      bool m_mtt_recoJetsAssociatedWellPlacedWithChi2;

      // Indexes of selected particles for mtt computation
      int m_selectedLeptonIndex_in_loose_collection; // Index of the selected lepton inside the Extractor "loose muon" collection
      int m_selectedLeptonIndex_in_array; // Index of the selected lepton inside the m_mtt_Muon* arrays
      TClonesArray* m_selectedLeptonP4;

      int m_mtt_NGoodMuons;
      int m_mtt_NLooseGoodMuons;
      float m_mtt_MuonPt[20];
      float m_mtt_MuonEta[20];
      float m_mtt_MuRelIso[20];

      int m_mtt_NGoodElectrons;
      float m_mtt_ElectronPt[20];
      float m_mtt_ElectronEta[20];
      float m_mtt_ElRelIso[20];

      float m_mtt_1stjetpt;
      float m_mtt_2ndjetpt;
      float m_mtt_3rdjetpt;
      float m_mtt_4thjetpt;
      float m_mtt_MET;

      int m_mtt_NJets;
      int m_mtt_NGoodJets;
      int m_mtt_NBtaggedJets_TCHPT;
      int m_mtt_NBtaggedJets_CSVL;
      int m_mtt_NBtaggedJets_CSVM;
      int m_mtt_NBtaggedJets_CSVT;

      float m_b_tagging_efficiency;

      float m_mtt_GoodJetEta[1000];
      float m_mtt_JetEta[1000];
      float m_mtt_JetPt[1000];
      int 	m_mtt_JetPuId[1000];

      //generic variables for 2D cut
      int pass2Dcut;

      TVector3 jet3P2D;
      float minjetpt2D;
      float DrMin;
      float pTRel;
      float costheta;

      int nGoodElectrons_veto;
      int Elepass2Dcut_veto;

      //variables for semie selection

      bool itsaZ;
      TVector3 el3P;

      int nGoodMuons_veto;
      //variables for jet selection
      int  nGoodJets;

      float m_lepton_weight;
      float m_lepton_weight_error_low;
      float m_lepton_weight_error_high;

      float m_btag_weight;
      float m_btag_weight_error_low;
      float m_btag_weight_error_high;

      // Event shape variables
      float m_sphericity;
      float m_aplanarity;
      float m_circularity;
      float m_isotropy;
      float m_D;
      float m_C;

      // Cut ; -1 event drop before arriving to this cut ; 0 cut failed, 1 cut passed
      int m_pass_vertex_cut;
      int m_pass_met_cut;
      int m_pass_lepton_cut;
      unsigned int m_lepton_sel_flag;
      int m_pass_jet_cut;
      bool m_do_mtt_reco;

      // b-tagging efficiencies
      TH2* m_gen_bjets; // Number of b-jets at generator level for an event
      TH2* m_gen_cjets; // Number of c-jets at generator level for an event
      TH2* m_gen_lightjets; // Number of light-jets at generator level for an event

      TH2* m_reco_bjets; // Number of b-jets at reco level correctly tagged
      TH2* m_reco_fake_bjets_among_cjets; // Number of tagged b-jets at reco level which are C jets
      TH2* m_reco_fake_bjets_among_lightjets; // Number of tagged b-jets at reco level which are light jets

      TH1* m_selected_jets_flavor;
      TH1* m_selected_jets_flavor_btagged;

      TH1* m_selected_b_jets_sf;
      TH1* m_selected_c_jets_sf;
      TH1* m_selected_light_jets_sf;

      TH1* m_1b_sf_flavor;
      TH1* m_1b_sf;
      TH1* m_2b_sf_flavor;
      TH1* m_2b_sf;

      std::shared_ptr<BTaggingEfficiencyProvider> m_b_tagging_efficiency_provider;

      std::vector<std::shared_ptr<SortingAlgorithm>> m_sorting_algortihms;
  };

}
