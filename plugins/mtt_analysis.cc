#include <Extractors/MttExtractorAnalysis/plugins/mtt_analysis.h>

#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>

#include "TH2.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "TMVA/Reader.h"

#include "Extractors/PatExtractor/interface/AnalysisSettings.h"
#include "Extractors/PatExtractor/interface/MCExtractor.h"
#include "Extractors/PatExtractor/interface/HLTExtractor.h"
#include "Extractors/PatExtractor/interface/MuonExtractor.h"
#include "Extractors/PatExtractor/interface/ElectronExtractor.h"
#include "Extractors/PatExtractor/interface/METExtractor.h"
#include "Extractors/PatExtractor/interface/VertexExtractor.h"
#include "Extractors/PatExtractor/interface/EventExtractor.h"
#include "Extractors/PatExtractor/interface/PatExtractor.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include <FWCore/Framework/interface/EventSetup.h>
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <FWCore/ParameterSet/interface/FileInPath.h>

#include "Extractors/MttExtractorAnalysis/plugins/SortingAlgorithmFactory.h"
#include "Extractors/MttExtractorAnalysis/plugins/KinFit.h"
#include "Extractors/MttExtractorAnalysis/plugins/BTaggingEfficiencyProvider.h"

using namespace std;

namespace mtt {

mtt_analysis::mtt_analysis(const edm::ParameterSet& cmsswSettings):
  Plugin(cmsswSettings)
{

  m_MAIN_doSemiMu = cmsswSettings.getParameter<bool>("do_semimu");

  m_MC_lepton_p4 = new TClonesArray("TLorentzVector");
  m_MC_neutrino_p4 = new TClonesArray("TLorentzVector");

  m_MC_leptonic_B_p4 = new TClonesArray("TLorentzVector");
  m_MC_hadronic_B_p4 = new TClonesArray("TLorentzVector");

  m_MC_lightJet1_p4 = new TClonesArray("TLorentzVector");
  m_MC_lightJet2_p4 = new TClonesArray("TLorentzVector");

  m_MC_Top1_p4 = new TClonesArray("TLorentzVector");
  m_MC_Top2_p4 = new TClonesArray("TLorentzVector");

  m_selectedLeptonP4 = new TClonesArray("TLorentzVector");

  edm::VParameterSet algos = cmsswSettings.getParameterSetVector("sorting_algortihms");
  for (const edm::ParameterSet& algo: algos) {
    bool enable = algo.getParameter<bool>("enable");
    if (! enable)
      continue;

    std::string name = algo.getParameter<std::string>("name");
    const edm::ParameterSet cfg = algo.getParameterSet("configuration");

    m_sorting_algortihms.push_back(SortingAlgorithmFactory::create(name, cfg, m_MAIN_doSemiMu));
  }

  reset();

  /// Tree definition
  m_tree_Mtt = new TTree("Mtt", "Analysis info");

  /// Branches definition

  m_tree_Mtt->Branch("MC_channel"         , &m_MC_channel             , "MC_channel/I");
  m_tree_Mtt->Branch("MC_mtt"             , &m_MC_mtt                 , "MC_mtt/F");
  m_tree_Mtt->Branch("MC_nPU"             , &m_nPU                    , "m_nPU/I");

  // Indexes of gen particles inside the MC collection. Only valid for semi-lept events
  m_tree_Mtt->Branch("MC_leptonIndex"     , &m_leptonIndex            , "MC_leptonIndex/I");
  m_tree_Mtt->Branch("MC_neutrinoIndex"   , &m_neutrinoIndex          , "MC_neutrinoIndex/I");
  m_tree_Mtt->Branch("MC_leptonicTopIndex", &m_leptonicTopIndex       , "MC_leptonicTopIndex/I");
  m_tree_Mtt->Branch("MC_leptonicBIndex"  , &m_leptonicBIndex         , "MC_leptonicBIndex/I");

  m_tree_Mtt->Branch("MC_hadronicBIndex"  , &m_hadronicBIndex         , "MC_hadronicBIndex/I");
  m_tree_Mtt->Branch("MC_hadronicFirstJetIndex" , &m_firstJetIndex    , "MC_hadronicFirstJetIndex/I");
  m_tree_Mtt->Branch("MC_hadronicSecondJetIndex", &m_secondJetIndex   , "MC_hadronicSecondJetIndex/I");

  m_tree_Mtt->Branch("MC_hadronicWMass"   , &m_MC_hadronicWMass       , "MC_hadronicWMass/F");
  m_tree_Mtt->Branch("MC_leptonicWMass"   , &m_MC_leptonicWMass       , "MC_leptonicWMass/F");
  m_tree_Mtt->Branch("MC_hadronicTopMass" , &m_MC_hadronicTopMass     , "MC_hadronicTopMass/F");
  m_tree_Mtt->Branch("MC_leptonicTopMass" , &m_MC_leptonicTopMass     , "MC_leptonicTopMass/F");
  m_tree_Mtt->Branch("MC_pt_tt"           , &m_MC_pt_tt               , "MC_pt_tt/F");
  m_tree_Mtt->Branch("MC_eta_tt"          , &m_MC_eta_tt              , "MC_eta_tt/F");
  m_tree_Mtt->Branch("MC_beta_tt"         , &m_MC_beta_tt             , "MC_beta_tt/F");

  // Easy access
  m_tree_Mtt->Branch("MC_lepton_p4"      , &m_MC_lepton_p4, 32000, -1);
  m_tree_Mtt->Branch("MC_neutrino_p4"    , &m_MC_neutrino_p4, 32000, -1);
  m_tree_Mtt->Branch("MC_leptonic_B_p4"  , &m_MC_leptonic_B_p4, 32000, -1);
  m_tree_Mtt->Branch("MC_hadronic_B_p4"  , &m_MC_hadronic_B_p4, 32000, -1);
  m_tree_Mtt->Branch("MC_lightJet1_B_p4" , &m_MC_lightJet1_p4, 32000, -1);
  m_tree_Mtt->Branch("MC_lightJet2_B_p4" , &m_MC_lightJet2_p4, 32000, -1);
  m_tree_Mtt->Branch("MC_Top1_p4"        , &m_MC_Top1_p4, 32000, -1);
  m_tree_Mtt->Branch("MC_Top2_p4"        , &m_MC_Top2_p4, 32000, -1);

  m_tree_Mtt->Branch("nGoodMuons"         , &m_mtt_NGoodMuons         , "nGoodMuons/I");
  m_tree_Mtt->Branch("muonPt"             , &m_mtt_MuonPt             , "muonPt[nGoodMuons]/F");
  m_tree_Mtt->Branch("muonEta"            , &m_mtt_MuonEta            , "muonEta[nGoodMuons]/F");
  m_tree_Mtt->Branch("muRelIso"           , &m_mtt_MuRelIso           , "muRelIso[nGoodMuons]/F");

  m_tree_Mtt->Branch("nGoodElectrons"     , &m_mtt_NGoodElectrons     , "nGoodElectrons/I");
  m_tree_Mtt->Branch("electronPt"         , &m_mtt_ElectronPt         , "electronPt[nGoodElectrons]/F");
  m_tree_Mtt->Branch("electronEta"        , &m_mtt_ElectronEta        , "electronEta[nGoodElectrons]/F");
  m_tree_Mtt->Branch("elRelIso"           , &m_mtt_ElRelIso           , "elRelIso[nGoodElectrons]/F");

  m_tree_Mtt->Branch("1stjetpt"           , &m_mtt_1stjetpt           , "1stjetpt/F");
  m_tree_Mtt->Branch("2ndjetpt"           , &m_mtt_2ndjetpt           , "2ndjetpt/F");
  m_tree_Mtt->Branch("3rdjetpt"           , &m_mtt_3rdjetpt           , "3rdjetpt/F");
  m_tree_Mtt->Branch("4thjetpt"           , &m_mtt_4thjetpt           , "4thjetpt/F");

  m_tree_Mtt->Branch("nJets"              , &m_mtt_NJets              , "nJets/I");
  m_tree_Mtt->Branch("jetEta"             , &m_mtt_JetEta             , "jetEta[nJets]/F");
  m_tree_Mtt->Branch("jetPt"              , &m_mtt_JetPt              , "jetPt[nJets]/F");

  m_tree_Mtt->Branch("nBtaggedJets_TCHPT" , &m_mtt_NBtaggedJets_TCHPT    , "nBtaggedJets_TCHPT/I");
  m_tree_Mtt->Branch("nBtaggedJets_CSVL" ,  &m_mtt_NBtaggedJets_CSVL    , "nBtaggedJets_CSVL/I");
  m_tree_Mtt->Branch("nBtaggedJets_CSVM" ,  &m_mtt_NBtaggedJets_CSVM    , "nBtaggedJets_CSVM/I");
  m_tree_Mtt->Branch("nBtaggedJets_CSVT" ,  &m_mtt_NBtaggedJets_CSVT    , "nBtaggedJets_CSVT/I");

  m_tree_Mtt->Branch("MET"                , &m_mtt_MET                   , "MET/F");

  m_tree_Mtt->Branch("isSel"              , &m_mtt_isSel                 , "isSel/I");

  m_tree_Mtt->Branch("eventIsAssociable"  , &m_mtt_eventIsAssociable     , "eventIsAssociable/O");

  m_tree_Mtt->Branch("trigger_passed", &m_trigger_passed, "trigger_passed/O");

  // Weights and errors from differents scale factors
  m_tree_Mtt->Branch("lepton_weight", &m_lepton_weight, "lepton_weight/F");
  m_tree_Mtt->Branch("lepton_weight_error_low", &m_lepton_weight_error_low, "lepton_weight_error_low/F");
  m_tree_Mtt->Branch("lepton_weight_error_high", &m_lepton_weight_error_high, "lepton_weight_error_high/F");

  m_tree_Mtt->Branch("btag_weight", &m_btag_weight, "btag_weight/F");
  m_tree_Mtt->Branch("btag_weight_error_low", &m_btag_weight_error_low, "btag_weight_error_low/F");
  m_tree_Mtt->Branch("btag_weight_error_high", &m_btag_weight_error_high, "btag_weight_error_high/F");

  // cuts
  m_tree_Mtt->Branch("pass_vertex_cut", &m_pass_vertex_cut, "pass_vertex_cut/I");
  m_tree_Mtt->Branch("pass_met_cut", &m_pass_met_cut, "pass_met_cut/I");
  m_tree_Mtt->Branch("pass_lepton_cut", &m_pass_lepton_cut, "pass_lepton_cut/I");
  m_tree_Mtt->Branch("lepton_cut_flag", &m_lepton_sel_flag, "lepton_cut_flag/i");
  m_tree_Mtt->Branch("pass_jet_cut", &m_pass_jet_cut, "pass_jet_cut/I");
  m_tree_Mtt->Branch("do_mtt_reco", &m_do_mtt_reco, "do_mtt_reco/O");

  m_tree_Mtt->Branch("selectedLeptonIndex_in_loose_collection", &m_selectedLeptonIndex_in_loose_collection, "selectedLeptonIndex_in_loose_collection/I");
  m_tree_Mtt->Branch("selectedLeptonIndex_in_array", &m_selectedLeptonIndex_in_array, "selectedLeptonIndex_in_array/I");
  m_tree_Mtt->Branch("selectedLeptonP4", &m_selectedLeptonP4, 32000, -1);

  // Create branches associated to sorting algos
  for (auto& algo: m_sorting_algortihms)
    algo->addBranches(*m_tree_Mtt);

  // METSel()
  m_MET_Pt_Min = cmsswSettings.getParameter<edm::ParameterSet>("met").getParameter<double>("pt_min");

  // MuonSel()
  m_MU_Pt_min_loose = cmsswSettings.getParameter<edm::ParameterSet>("muons_loose").getParameter<double>("pt_min");
  m_MU_Eta_max_loose = cmsswSettings.getParameter<edm::ParameterSet>("muons_loose").getParameter<double>("eta_max");
  m_MU_Iso_max_loose = cmsswSettings.getParameter<edm::ParameterSet>("muons_loose").getParameter<double>("isolation_max");

  m_MU_Pt_min = cmsswSettings.getParameter<edm::ParameterSet>("muons_tight").getParameter<double>("pt_min");
  m_MU_Eta_max = cmsswSettings.getParameter<edm::ParameterSet>("muons_tight").getParameter<double>("eta_max");
  m_MU_Iso_max = cmsswSettings.getParameter<edm::ParameterSet>("muons_tight").getParameter<double>("isolation_max");

  // ElectronSel()
  m_ELE_Pt_min_loose = cmsswSettings.getParameter<edm::ParameterSet>("electrons_loose").getParameter<double>("pt_min");
  m_ELE_Eta_max_loose = cmsswSettings.getParameter<edm::ParameterSet>("electrons_loose").getParameter<double>("eta_max");
  m_ELE_Iso_max_loose = cmsswSettings.getParameter<edm::ParameterSet>("electrons_loose").getParameter<double>("isolation_max");

  m_ELE_Pt_min = cmsswSettings.getParameter<edm::ParameterSet>("electrons_tight").getParameter<double>("pt_min");
  m_ELE_Eta_max = cmsswSettings.getParameter<edm::ParameterSet>("electrons_tight").getParameter<double>("eta_max");
  m_ELE_Iso_max = cmsswSettings.getParameter<edm::ParameterSet>("electrons_tight").getParameter<double>("isolation_max");

  // JetSel()
  m_JET_Pt_min  = cmsswSettings.getParameter<edm::ParameterSet>("jets").getParameter<double>("pt_min");
  m_JET_Eta_max = cmsswSettings.getParameter<edm::ParameterSet>("jets").getParameter<double>("eta_max");
  m_JET_btag_CSVL = cmsswSettings.getParameter<edm::ParameterSet>("jets").getParameter<double>("btag_CSVL");
  m_JET_btag_CSVM = cmsswSettings.getParameter<edm::ParameterSet>("jets").getParameter<double>("btag_CSVM");
  m_JET_btag_CSVT = cmsswSettings.getParameter<edm::ParameterSet>("jets").getParameter<double>("btag_CSVT");
  m_JET_btag_TCHPT = cmsswSettings.getParameter<edm::ParameterSet>("jets").getParameter<double>("btag_TCHPT");
  m_b_tagging_efficiency = cmsswSettings.getParameter<edm::ParameterSet>("jets").getParameter<double>("b_tagging_efficiency");

  float ptBinning[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800};
  float etaBinning[] = {0., 0.8, 1.5, 2.4};

  m_gen_bjets = new TH2F("number_of_gen_bjets", "", 16, ptBinning, 3, etaBinning);
  m_gen_cjets = new TH2F("number_of_gen_cjets", "", 16, ptBinning, 3, etaBinning);
  m_gen_lightjets = new TH2F("number_of_gen_lightjets", "", 16, ptBinning, 3, etaBinning);

  m_reco_bjets = new TH2F("number_of_reco_bjets", "", 16, ptBinning, 3, etaBinning);
  m_reco_fake_bjets_among_cjets = new TH2F("number_of_reco_fake_bjets_among_cjets", "", 16, ptBinning, 3, etaBinning);
  m_reco_fake_bjets_among_lightjets = new TH2F("number_of_reco_fake_bjets_among_lightjets", "", 16, ptBinning, 3, etaBinning);

  m_selected_jets_flavor = new TH1F("selected_jets_flavor", "True flavor of selected jets", 3, 0, 3);
  m_selected_jets_flavor->GetXaxis()->SetBinLabel(1, "b jets");
  m_selected_jets_flavor->GetXaxis()->SetBinLabel(2, "c jets");
  m_selected_jets_flavor->GetXaxis()->SetBinLabel(3, "light jets");

  m_selected_jets_flavor_btagged = new TH1F("selected_jets_flavor_btagged", "True flavor of selected b-tagged jets", 3, 0, 3);
  m_selected_jets_flavor_btagged->GetXaxis()->SetBinLabel(1, "b jets");
  m_selected_jets_flavor_btagged->GetXaxis()->SetBinLabel(2, "c jets");
  m_selected_jets_flavor_btagged->GetXaxis()->SetBinLabel(3, "light jets");

  m_2b_sf = new TH1F("2b_sf", "Event weight for the 2b category", 50, 0.5, 1.5);
  m_2b_sf_flavor = new TH1F("2b_sf_flavor", "True flavor of the 2 b-tagged jets in 2b category", 3, 0, 3);
  m_2b_sf_flavor->GetXaxis()->SetBinLabel(1, "b jets");
  m_2b_sf_flavor->GetXaxis()->SetBinLabel(2, "c jets");
  m_2b_sf_flavor->GetXaxis()->SetBinLabel(3, "light jets");

  m_1b_sf = new TH1F("1b_sf", "Event weight for the 1b category", 50, 0.5, 1.5);
  m_1b_sf_flavor = new TH1F("1b_sf_flavor", "True flavor of the jets in 1b category", 3, 0, 3);
  m_1b_sf_flavor->GetXaxis()->SetBinLabel(1, "b jets");
  m_1b_sf_flavor->GetXaxis()->SetBinLabel(2, "c jets");
  m_1b_sf_flavor->GetXaxis()->SetBinLabel(3, "light jets");

  m_selected_b_jets_sf = new TH1F("selected_b_jets_sf", "Scale factor for selected b jets", 50, 0.5, 1.5);
  m_selected_c_jets_sf = new TH1F("selected_c_jets_sf", "Scale factor for selected c jets", 50, 0.5, 1.5);
  m_selected_light_jets_sf = new TH1F("selected_light_jets_sf", "Scale factor for selected light jets", 50, 0.5, 1.5);

  m_b_tagging_efficiency_provider = std::make_shared<BTaggingEfficiencyProvider>(cmsswSettings);
}

const int NO_MTT_RECO = -1;
const int FAILURE = 0;
const int SUCCESS = 1;

//
// First you start with the different physics objects selections
//

// Vertices
int mtt_analysis::VertexSel()
{
  int n_vtx = m_vertex->getSize();

  if (!n_vtx)
    return FAILURE;

  for (int i = 0; i < n_vtx; ++i)
  {
    if (m_vertex->getVtxIsFake(i)) continue;
    if (m_vertex->getVtxNdof(i) < 4) continue;

    return SUCCESS;
  }

  return FAILURE;
}

// MET
int mtt_analysis::METSel()
{
  m_mtt_MET = m_jetMet->getMETLorentzVector(0)->Pt();
  if (m_mtt_MET > m_MET_Pt_Min) {
    return SUCCESS;
  }

  return FAILURE;
}

const uint32_t NO_LEPTON = 0;
const uint32_t PASS_MUON_VETO = 1;
const uint32_t PASS_ELECTRON_VETO = 2;
const uint32_t PASS_LEPTON_VETO = PASS_MUON_VETO | PASS_ELECTRON_VETO;

const uint32_t HAS_ONE_ISOLATED_LEPTON = 8;
const uint32_t HAS_NO_ISOLATED_LEPTON = 16;
const uint32_t HAS_MANY_ISOLATED_LEPTON = 32;

int mtt_analysis::MuonSel()
{
  bool pass_cut = false;

  m_lepton_sel_flag = NO_LEPTON;

  int selected_muon_index = -1;
  int n_loose_muons = m_muon_loose->getSize();
  int n_isolated_muons = 0;
  int isolated_muon_index_in_loose_collection = -1;
  int isolated_muon_index_in_array = -1;

  std::vector<int> looseCollectionIndexes;
  for (int i = 0; i < n_loose_muons; i++) {

    // Selection from https://twiki.cern.ch/twiki/bin/view/CMS/TWikiTopRefEventSel#Muons
    if (! m_muon_loose->getMuisGlobal(i))
      continue;

    TLorentzVector *muP = m_muon_loose->getMuLorentzVector(i);

    if (fabs(muP->Pt()) <= 26)
      continue;

    if (fabs(muP->Eta()) >= 2.1)
      continue;

    if (m_muon_loose->getMunormChi2(i) >= 10.)
      continue;

    if (m_muon_loose->getTrackerLayersWithMeasurements(i) <= 5)
      continue;

    if (m_muon_loose->getGlobalTrackNumberOfValidMuonHits(i) <= 0)
      continue;

    if (m_muon_loose->getNumberOfMatchedStations(i) <= 1)
      continue;

    if (m_muon_loose->getMudB(i) >= 0.2)
      continue;

    if (m_muon_loose->getdZ(i) >= 0.5)
      continue;

    if (m_muon_loose->getMunValPixelHits(i) <= 0)
      continue;

    m_mtt_MuonPt[m_mtt_NGoodMuons]   = muP->Pt();
    m_mtt_MuonEta[m_mtt_NGoodMuons]  = muP->Eta();
    m_mtt_MuRelIso[m_mtt_NGoodMuons] = m_muon_loose->getDeltaBetaCorrectedRelativeIsolation(i);
    looseCollectionIndexes.push_back(i);

    if (m_muon_loose->getDeltaBetaCorrectedRelativeIsolation(i) < 0.12) {
      n_isolated_muons++;
      isolated_muon_index_in_loose_collection = i;
      isolated_muon_index_in_array = m_mtt_NGoodMuons;
    }

    m_mtt_NGoodMuons++;
  }

  // No muons
  if (m_mtt_NGoodMuons == 0)
    return NO_MTT_RECO;

  if (n_isolated_muons == 1) {
    // One isolated lepton = standard selection
    // Check if there are other isolated leptons in the event

    m_lepton_sel_flag |= HAS_ONE_ISOLATED_LEPTON;

    // Muon veto
    m_lepton_sel_flag |= PASS_MUON_VETO;
    int size = m_muon_loose->getSize();
    for (int i = 0; i < size; i++) {

      if (i == isolated_muon_index_in_loose_collection)
        continue;

      TLorentzVector* p4 = m_muon_loose->getMuLorentzVector(i);

      if ((m_muon_loose->getMuisGlobal(i) || m_muon_loose->getMuisTracker(i)) &&
          (p4->Pt() > 10) &&
          (fabs(p4->Eta()) < 2.5) &&
          (m_muon_loose->getDeltaBetaCorrectedRelativeIsolation(i) < 0.20)) {

        m_lepton_sel_flag = m_lepton_sel_flag & ~PASS_MUON_VETO;
        break;
      }
    }

    //electron veto for semimu channel
    m_lepton_sel_flag |= PASS_ELECTRON_VETO;

    int n_ele = m_electron_loose->getSize();
    for (int i = 0; i < n_ele; ++i) {
      TLorentzVector *eP = m_electron_loose->getEleLorentzVector(i);

      if ((eP->Pt() > 20) &&
          (fabs(eP->Eta()) < 2.5) &&
          (m_electron_loose->passVetoId(i)) &&
          (m_electron_loose->getRhoCorrectedRelativeIsolation(i) < 0.15)) {

        m_lepton_sel_flag = m_lepton_sel_flag & ~PASS_ELECTRON_VETO;
        break;
      }
    }

    selected_muon_index = isolated_muon_index_in_loose_collection;
    m_selectedLeptonIndex_in_array = isolated_muon_index_in_array;

    pass_cut = (m_lepton_sel_flag & PASS_LEPTON_VETO) == PASS_LEPTON_VETO;

  } else if (n_isolated_muons > 1) {
    // More than one isolated lepton
    // Drop the event

    return NO_MTT_RECO;

  } else {
    // No isolated muon
    // Skip lepton veto, and use the first lepton

    m_lepton_sel_flag |= HAS_NO_ISOLATED_LEPTON;
    selected_muon_index = looseCollectionIndexes[0];
    m_selectedLeptonIndex_in_array = 0;
  }

  m_refLept = m_muon_loose->getMuLorentzVector(selected_muon_index);
  m_selectedLeptonIndex_in_loose_collection = selected_muon_index;

  if (m_isMC) {
    // Get scale factor
    ScaleFactor sf = m_muon_loose->getScaleFactor(ScaleFactorService::TIGHT, ScaleFactorService::TIGHT, selected_muon_index);
    m_lepton_weight *= sf.getValue();
    m_lepton_weight_error_low += sf.getErrorLow() * sf.getErrorLow();
    m_lepton_weight_error_high += sf.getErrorHigh() * sf.getErrorHigh();
  }

  return pass_cut ? SUCCESS : FAILURE;
}

int mtt_analysis::ElectronSel()
{
  bool pass_cut = false;

  m_lepton_sel_flag = NO_LEPTON;

  int selected_electron_index = -1;
  int n_loose_electrons = m_electron_loose->getSize();
  int n_isolated_electrons = 0;
  int isolated_electron_index_in_loose_collection = -1;
  int isolated_electron_index_in_array = -1;

  std::vector<int> looseCollectionIndexes;
  for (int i = 0; i < n_loose_electrons; i++) {

    TLorentzVector *eP = m_electron_loose->getEleLorentzVector(i);

    if (fabs(eP->Pt()) <= 30)
      continue;

    if (fabs(eP->Eta()) >= 2.5)
      continue;

    if (fabs(m_electron_loose->getSuperClusterEta(i)) >= 1.4442 && fabs(m_electron_loose->getSuperClusterEta(i)) < 1.5660)
      continue;

    if (! m_electron_loose->passTightId(i))
      continue;

    m_mtt_ElectronPt[m_mtt_NGoodElectrons] = eP->Pt();
    m_mtt_ElectronEta[m_mtt_NGoodElectrons] = eP->Eta();
    m_mtt_ElRelIso[m_mtt_NGoodElectrons]   = m_electron_loose->getRhoCorrectedRelativeIsolation(i);
    looseCollectionIndexes.push_back(i);

    if (m_electron_loose->getRhoCorrectedRelativeIsolation(i) < 0.10) {
      n_isolated_electrons++;
      isolated_electron_index_in_loose_collection = i;
      isolated_electron_index_in_array = m_mtt_NGoodElectrons;
    }

    m_mtt_NGoodElectrons++;
  }

  // No electrons? bye bye
  if (m_mtt_NGoodElectrons == 0)
    return NO_MTT_RECO;

  if (n_isolated_electrons == 1) {
    // One isolated lepton = standard selection
    // Check if there are other isolated leptons in the event

    m_lepton_sel_flag |= HAS_ONE_ISOLATED_LEPTON;

    // Muon veto
    m_lepton_sel_flag |= PASS_MUON_VETO;
    int size = m_muon_loose->getSize();
    for (int i = 0; i < size; i++) {

      TLorentzVector* p4 = m_muon_loose->getMuLorentzVector(i);

      if ((m_muon_loose->getMuisGlobal(i) || m_muon_loose->getMuisTracker(i)) &&
          (p4->Pt() > 10) &&
          (fabs(p4->Eta()) < 2.5) &&
          (m_muon_loose->getDeltaBetaCorrectedRelativeIsolation(i) < 0.20)) {

        m_lepton_sel_flag = m_lepton_sel_flag & ~PASS_MUON_VETO;
        break;
      }
    }

    //electron veto for semimu channel
    m_lepton_sel_flag |= PASS_ELECTRON_VETO;

    int n_ele = m_electron_loose->getSize();
    for (int i = 0; i < n_ele; ++i) {

      if (i == isolated_electron_index_in_loose_collection)
        continue;

      TLorentzVector *eP = m_electron_loose->getEleLorentzVector(i);

      if ((eP->Pt() > 20) &&
          (fabs(eP->Eta()) < 2.5) &&
          (m_electron_loose->passVetoId(i)) &&
          (m_electron_loose->getRhoCorrectedRelativeIsolation(i) < 0.15)) {

        m_lepton_sel_flag = m_lepton_sel_flag & ~PASS_ELECTRON_VETO;
        break;
      }
    }

    selected_electron_index = isolated_electron_index_in_loose_collection;
    m_selectedLeptonIndex_in_array = isolated_electron_index_in_array;

    pass_cut = (m_lepton_sel_flag & PASS_LEPTON_VETO) == PASS_LEPTON_VETO;

  } else if (n_isolated_electrons > 1) {
    // More than one isolated lepton
    // Drop the event

    return NO_MTT_RECO;

  } else {
    // No isolated electron
    // Skip lepton veto, and use the first lepton

    m_lepton_sel_flag |= HAS_NO_ISOLATED_LEPTON;
    selected_electron_index = looseCollectionIndexes[0];
    m_selectedLeptonIndex_in_array = 0;
  }

  m_refLept = m_electron_loose->getEleLorentzVector(selected_electron_index);
  m_selectedLeptonIndex_in_loose_collection = selected_electron_index;

  if (m_isMC) {
    // Get scale factor
    ScaleFactor sf = m_electron_loose->getScaleFactor(ScaleFactorService::TIGHT, ScaleFactorService::TIGHT, selected_electron_index);
    m_lepton_weight *= sf.getValue();
    m_lepton_weight_error_low += sf.getErrorLow() * sf.getErrorLow();
    m_lepton_weight_error_high += sf.getErrorHigh() * sf.getErrorHigh();
  }

  return pass_cut ? SUCCESS : FAILURE;
}


int mtt_analysis::JetSel()
{
  AllJetsPt       = 0.;
  m_selJetsIds.clear();

  int n_jet = m_jetMet->getSize();

  if (! n_jet)
    return NO_MTT_RECO;

  std::vector<ScaleFactor> jetSF;
  std::vector<int> jetFlavor;
  std::vector<bool> jetIsBTagged;

  for (int i = 0; i < n_jet; i++)
  {
    TLorentzVector *jetP = m_jetMet->getP4(i);

    if (fabs(jetP->Pt()) < m_JET_Pt_min) continue;
    if (fabs(jetP->Eta()) >  m_JET_Eta_max) continue;

    m_mtt_JetEta[m_mtt_NJets] = jetP->Eta();
    m_mtt_JetPt[m_mtt_NJets]  = jetP->Pt();
    if (m_isMC) {
      jetSF.push_back(m_jetMet->getScaleFactor(i));
      jetFlavor.push_back(m_jetMet->getAlgoPartonFlavor(i));
    }

    ++m_mtt_NJets;

    AllJetsPt += jetP->Pt();

    if (m_mtt_NJets < 9) // Count the number of btagged jets in the selected jets
    {
      m_selJetsIds.push_back(i);
      if ((m_jetMet->getJetBTagProb_CSV(i)) > m_JET_btag_CSVL)
        ++m_mtt_NBtaggedJets_CSVL;

      if ((m_jetMet->getJetBTagProb_CSV(i)) > m_JET_btag_CSVM) {
        ++m_mtt_NBtaggedJets_CSVM;
        jetIsBTagged.push_back(true);
      } else {
        jetIsBTagged.push_back(false);
      }

      if ((m_jetMet->getJetBTagProb_CSV(i)) > m_JET_btag_CSVT)
        ++m_mtt_NBtaggedJets_CSVT;
      if ((m_jetMet->getJetBTagProb_TCHP(i)) > m_JET_btag_TCHPT)
        ++m_mtt_NBtaggedJets_TCHPT;
    }

    if (m_mtt_NJets == 1) m_mtt_1stjetpt = m_mtt_JetPt[m_mtt_NJets - 1];
    if (m_mtt_NJets == 2) m_mtt_2ndjetpt = m_mtt_JetPt[m_mtt_NJets - 1];
    if (m_mtt_NJets == 3) m_mtt_3rdjetpt = m_mtt_JetPt[m_mtt_NJets - 1];
    if (m_mtt_NJets == 4) m_mtt_4thjetpt = m_mtt_JetPt[m_mtt_NJets - 1];
  }

  // We need at least 4 good jets
  if (m_mtt_NJets < 4)
    return NO_MTT_RECO;

  auto getEfficiencyWithError = [&] (int mcFlavor, float pt, float eta) -> std::tuple<double, double, double> {
    mcFlavor = abs(mcFlavor);
    ScaleFactorService::Flavor flavor = ScaleFactorService::B;
    if (mcFlavor == 4) {
      flavor = ScaleFactorService::C;
    } else if ((mcFlavor <= 3) || (mcFlavor == 21)) {
      // If mcFlavor == 0, assume it's a light jet
      flavor = ScaleFactorService::LIGHT;
    }

    return m_b_tagging_efficiency_provider->getEfficiency(flavor, pt, eta);
  };

  if (m_isMC) {

    assert(jetIsBTagged.size() == m_selJetsIds.size());

    for (size_t i = 0; i < m_selJetsIds.size(); i++) {
      int extractorIndex = m_selJetsIds[i];

      TLorentzVector *jetP = m_jetMet->getP4(extractorIndex);

      float pt = jetP->Pt();
      float eta = fabs(jetP->Eta());

      int mcFlavor = abs(jetFlavor[i]);
      int flavorBin = 0;
      if (mcFlavor == 5) {
        // The true flavor of this jet is B
        m_gen_bjets->Fill(pt, eta);

        if (jetIsBTagged[i]) {
          // This reco jet is tagged as a B
          m_reco_bjets->Fill(pt, eta);
        }

        m_selected_b_jets_sf->Fill(jetSF[i].getValue());
      } else if (mcFlavor == 4) {
        // The true flavor of this jet is C
        m_gen_cjets->Fill(pt, eta);

        if (jetIsBTagged[i]) {
          // This reco jet is tagged as a B
          m_reco_fake_bjets_among_cjets->Fill(pt, eta);
        }

        flavorBin = 1;
        m_selected_c_jets_sf->Fill(jetSF[i].getValue());
      } else if ((mcFlavor <= 3) || (mcFlavor == 21)) {
        // The true flavor of this jet is u, d, s or g
        m_gen_lightjets->Fill(pt, eta);

        if (jetIsBTagged[i]) {
          // This reco jet is tagged as a B
          m_reco_fake_bjets_among_lightjets->Fill(pt, eta);
        }

        flavorBin = 2;
        m_selected_light_jets_sf->Fill(jetSF[i].getValue());
      }

      m_selected_jets_flavor->Fill(flavorBin);

      if (jetIsBTagged[i])
        m_selected_jets_flavor_btagged->Fill(flavorBin);
    }

    if (m_mtt_NBtaggedJets_CSVM == 1) {

      // We use method 1.a) from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods

      float error_tagged_squared_up = 0.f;
      float error_untagged_squared_up = 0.f;
      float error_tagged_squared_low = 0.f;
      float error_untagged_squared_low = 0.f;

      float eff_tagged = 1.;
      float eff_untagged = 1.;

      for (size_t i = 0; i < m_selJetsIds.size(); i++) {

        int extractorIndex = m_selJetsIds[i];
        TLorentzVector *p4 = m_jetMet->getP4(extractorIndex);

        int mcFlavor = fabs(jetFlavor[i]);
        int flavor = 0;
        if (mcFlavor == 4) {
          flavor = 1;
        } else if ((mcFlavor <= 3) || (mcFlavor == 21)) {
          // If mcFlavor == 0, assume it's a light jet
          flavor = 2;
        }
        m_1b_sf_flavor->Fill(flavor);

        float pt = p4->Pt();
        float eta = fabs(p4->Eta());

        auto eff = getEfficiencyWithError(jetFlavor[i], pt, eta);
        ScaleFactor sf = jetSF[i];

        float eff_i = std::get<0>(eff);
        float sf_i = sf.getValue();

        float error_eff_i_up = std::get<1>(eff);
        float error_sf_i_up = jetSF[i].getErrorHigh();
        float error_eff_i_low = std::get<2>(eff);
        float error_sf_i_low = jetSF[i].getErrorLow();

        // Compute b-tag weight error
        // See calculation here: https://github.com/IPNL-CMS/MttExtractorAnalysis/raw/master/doc/btagWeight.pdf
        if (jetIsBTagged[i]) {
          eff_tagged *= sf_i;
          error_tagged_squared_up += ((error_sf_i_up * error_sf_i_up) / (sf_i * sf_i));
          error_tagged_squared_low += ((error_sf_i_low * error_sf_i_low) / (sf_i * sf_i));
        } else {
          eff_untagged *= ((1. - sf_i * eff_i)) / (1. - eff_i);
          error_untagged_squared_up += std::pow(((-eff_i) / (1 - eff_i * sf_i)) * error_sf_i_up, 2) + std::pow(((1 - sf_i) / ((1 - eff_i) * (1 - sf_i * eff_i))) * error_eff_i_up, 2);
          error_untagged_squared_low += std::pow(((-eff_i) / (1 - eff_i * sf_i)) * error_sf_i_low, 2) + std::pow(((1 - sf_i) / ((1 - eff_i) * (1 - sf_i * eff_i))) * error_eff_i_low, 2);
        }
      }

      m_btag_weight = eff_tagged * eff_untagged;

      float btag_eff_squared = m_btag_weight * m_btag_weight;
      m_btag_weight_error_high = (error_tagged_squared_up + error_untagged_squared_up) * btag_eff_squared;
      m_btag_weight_error_low = (error_tagged_squared_low + error_untagged_squared_low) * btag_eff_squared;

      m_1b_sf->Fill(m_btag_weight);

    } else if (m_mtt_NBtaggedJets_CSVM > 1) {

      // In this case, the weight is simply the product
      // of the two leading B jets scale factors

      int nBTag = 0;
      m_btag_weight = 1;
      for (size_t i = 0; i < m_selJetsIds.size(); i++) {
        if (jetIsBTagged[i]) {
          nBTag++;
        } else {
          continue;
        }

        if (nBTag > 2)
          break;

        int mcFlavor = fabs(jetFlavor[i]);
        int flavor = 0;
        if (mcFlavor == 4) {
          flavor = 1;
        } else if ((mcFlavor <= 3) || (mcFlavor == 21)) {
          // If mcFlavor == 0, assume it's a light jet
          flavor = 2;
        }
        m_2b_sf_flavor->Fill(flavor);

        ScaleFactor sf = jetSF[i];

        float sf_i = sf.getValue();

        float error_sf_i_up = jetSF[i].getErrorHigh();
        float error_sf_i_low = jetSF[i].getErrorLow();

        m_btag_weight *= sf_i;
        m_btag_weight_error_high += (error_sf_i_up * error_sf_i_up) / (sf_i * sf_i);
        m_btag_weight_error_low += (error_sf_i_low * error_sf_i_low) / (sf_i * sf_i);
      }

      m_2b_sf->Fill(m_btag_weight);
      float btag_eff_squared = m_btag_weight * m_btag_weight;
      m_btag_weight_error_high *= btag_eff_squared;
      m_btag_weight_error_low *= btag_eff_squared;
    }
  }

  return SUCCESS;
}

#define PROCESS_RESULT(res, var) \
  var = res; \
  if (res == NO_MTT_RECO) { \
    m_do_mtt_reco = false; \
  }

#define   MTT_TRIGGER_NOT_FOUND   1000

void mtt_analysis::analyze(const edm::Event& event, const edm::EventSetup& iSetup, PatExtractor& extractor) {
  analyze(iSetup, extractor);
}

void mtt_analysis::analyze(const edm::EventSetup& iSetup, PatExtractor& extractor)
{
  reset();

  m_refLept  = nullptr;

  m_vertex   = std::static_pointer_cast<VertexExtractor>(extractor.getExtractor("vertex"));

  m_muon     = std::static_pointer_cast<MuonExtractor>(extractor.getExtractor("muons"));
  m_muon_loose = std::static_pointer_cast<MuonExtractor>(extractor.getExtractor("muons_loose"));

  m_electron = std::static_pointer_cast<ElectronExtractor>(extractor.getExtractor("electrons"));
  m_electron_loose = std::static_pointer_cast<ElectronExtractor>(extractor.getExtractor("electrons_loose"));

  m_jetMet   = std::static_pointer_cast<JetMETExtractor>(extractor.getExtractor("JetMET"));

  m_event    = std::static_pointer_cast<EventExtractor>(extractor.getExtractor("event"));

  std::shared_ptr<HLTExtractor> HLT = std::static_pointer_cast<HLTExtractor>(extractor.getExtractor("HLT"));
  m_trigger_passed = HLT->isTriggerFired();

  if (m_isMC)
  {
    m_MC = std::static_pointer_cast<MCExtractor>(extractor.getExtractor("MC"));
    MCidentification();
  }

  m_nPU = m_event->nPU();

  m_do_mtt_reco = true;

  int res = VertexSel();
  PROCESS_RESULT(res, m_pass_vertex_cut);

  res = METSel();
  PROCESS_RESULT(res, m_pass_met_cut);

  if (m_MAIN_doSemiMu) {
    res  = MuonSel();
  } else {
    res = ElectronSel();
  }
  PROCESS_RESULT(res, m_pass_lepton_cut);

  res = JetSel();
  PROCESS_RESULT(res, m_pass_jet_cut);

  m_mtt_isSel = m_pass_vertex_cut == SUCCESS && m_pass_met_cut == SUCCESS && m_pass_lepton_cut == SUCCESS && m_pass_jet_cut == SUCCESS;

  if (m_do_mtt_reco) {
    loopOverCombinations();
  }

  fillTree();
}

void mtt_analysis::loopOverCombinations()
{
  new ((*m_selectedLeptonP4)[0]) TLorentzVector(*m_refLept);

  std::vector<Jet> jets;
  for (size_t i = 0; i < m_selJetsIds.size(); i++) {
    Jet jet;

    int index = m_selJetsIds[i];
    TLorentzVector *jetP = m_jetMet->getP4(index);

    // LorentzVector is used in (Pt, Eta, Phi, E) coordinates
    jet.p.SetCoordinates(jetP->Pt(), jetP->Eta(), jetP->Phi(), jetP->E());
    jet.index = index;
    jet.isBTagged = isBJet(index);

    jets.push_back(jet);
  }

  LorentzVector lepton(m_refLept->Pt(), m_refLept->Eta(), m_refLept->Phi(), m_refLept->E());

  TLorentzVector* metP = m_jetMet->getMETLorentzVector(0);
  LorentzVector met(metP->Pt(), metP->Eta(), metP->Phi(), metP->E());

  for (auto& algo: m_sorting_algortihms) {
    algo->setObjects(jets, lepton, met);
    algo->work();
  }

  //FIXME
  if (m_isMC)
    checkIfSolutionIsCorrect();
}


// MC stuff

#define ID_B (5)
#define ID_T (6)

#define ID_E (11)
#define ID_NEUTRINO_E (12)
#define ID_MU (13)
#define ID_NEUTRINO_MU (14)
#define ID_TAU (15)
#define ID_NEUTRINO_TAU (16)

#define ID_W (24)

int mtt_analysis::patIndexToExtractorIndex(int patIndex) const {

  for (int i = 0; i < m_MC->getSize() ; i++) {
    if (m_MC->getPatIndex(i) == patIndex)
      return i;
  }

  return -1;
}

void mtt_analysis::MCidentification()
{
  nEle    = 0;
  nMu     = 0;
  nTau    = 0;
  nNuEle  = 0;
  nNuMu   = 0;
  nNuTau  = 0;
  nQuarkb = 0;
  nTop    = 0;
  Top.clear();
  std::vector<int> topIndexes;

  int n_MC = m_MC->getSize();

  if (!n_MC)
    return;

  for (int i = 0; i < n_MC ; ++i)
  {

    if (abs(m_MC->getType(i)) == ID_T) {

      if (std::find(topIndexes.begin(), topIndexes.end(), i) != topIndexes.end())
        continue;

      topIndexes.push_back(i);
      TLorentzVector TL_Top(m_MC->getPx(i),
          m_MC->getPy(i),
          m_MC->getPz(i),
          m_MC->getE(i));

      Top.push_back(TL_Top);
      nTop++;

      continue;
    }

    int motherIndex = patIndexToExtractorIndex(m_MC->getMom1Index(i));
    int grandMotherIndex = -1;
    if (motherIndex != -1)
      grandMotherIndex = patIndexToExtractorIndex(m_MC->getMom1Index(motherIndex));

    if (motherIndex == -1)
      continue;

    if (false) {
      std::cout << "Type: " << m_MC->getType(i) << std::endl;
      std::cout << "Mother type: " << m_MC->getType(motherIndex) << std::endl;
      if (grandMotherIndex != -1)
        std::cout << "Grandmother type: " << m_MC->getType(grandMotherIndex) << std::endl;
    }

    if (abs(m_MC->getType(motherIndex)) == ID_T || (grandMotherIndex != -1 && abs(m_MC->getType(grandMotherIndex) == ID_T)))
    {
      /// Count the number of leptons and neutrinos from Top->W
      if (fabs(m_MC->getType(i)) == 11) ++nEle;   //Electron from WTop
      if (fabs(m_MC->getType(i)) == 13) ++nMu;    //Muon	   from WTop
      if (fabs(m_MC->getType(i)) == 15) ++nTau;   //Tau	   from WTop
      if (fabs(m_MC->getType(i)) == 12) ++nNuEle; //NuEle    from WTop
      if (fabs(m_MC->getType(i)) == 14) ++nNuMu;  //NuMu	   from WTop
      if (fabs(m_MC->getType(i)) == 16) ++nNuTau; //NuTau    from WTop
    }

    /// Count the number of b quark from Top
    if (abs(m_MC->getType(i)) == ID_B && abs(m_MC->getType(motherIndex)) == ID_T)
    {
      nQuarkb++; //Quark b from Top
    }
  }

  if (nEle == 1 && nNuEle == 1 && nMu == 0 && nNuMu == 0 && nTau == 0 && nNuTau == 0 && nQuarkb > 1 && nTop == 2)
  {
    m_MC_channel = 1;
  }

  if (nEle == 0 && nNuEle == 0 && nMu == 1 && nNuMu == 1 && nTau == 0 && nNuTau == 0 && nQuarkb > 1 && nTop == 2)
  {
    m_MC_channel = 2;
  }

  if (nEle == 0 && nNuEle == 0 && nMu == 0 && nNuMu == 0 && nTau == 1 && nNuTau == 1 && nQuarkb > 1 && nTop == 2)
  {
    m_MC_channel = 3;
  }

  if (nEle == 0 && nNuEle == 0 && nMu == 0 && nNuMu == 0 && nTau == 0 && nNuTau == 0 && nQuarkb > 1 && nTop == 2)
  {
    m_MC_channel = 4;
  }

  if (nEle == 2 && nNuEle == 2 && nMu == 0 && nNuMu == 0 && nTau == 0 && nNuTau == 0 && nQuarkb > 1 && nTop == 2)
  {
    m_MC_channel = 5;
  }

  if (nEle == 0 && nNuEle == 0 && nMu == 2 && nNuMu == 2 && nTau == 0 && nNuTau == 0 && nQuarkb > 1 && nTop == 2)
  {
    m_MC_channel = 6;
  }

  if (nEle == 0 && nNuEle == 0 && nMu == 0 && nNuMu == 0 && nTau == 2 && nNuTau == 2 && nQuarkb > 1 && nTop == 2)
  {
    m_MC_channel = 7;
  }

  if (nEle == 1 && nNuEle == 1 && nMu == 1 && nNuMu == 1 && nTau == 0 && nNuTau == 0 && nQuarkb == 2 && nTop == 2)
  {
    m_MC_channel = 8 ;
  }

  if (nEle == 1 && nNuEle == 1 && nMu == 0 && nNuMu == 0 && nTau == 1 && nNuTau == 1 && nQuarkb == 2 && nTop == 2)
  {
    m_MC_channel = 9 ;
  }

  if (nEle == 0 && nNuEle == 0 && nMu == 1 && nNuMu == 1 && nTau == 1 && nNuTau == 1 && nQuarkb == 2 && nTop == 2)
  {
    m_MC_channel = 10;
  }

  if (nTop == 2) {

    TLorentzVector mc_resonance = Top[0] + Top[1];

    m_MC_mtt = mc_resonance.M();
    m_MC_pt_tt = mc_resonance.Pt();
    m_MC_eta_tt = mc_resonance.Eta();
    m_MC_beta_tt = fabs(mc_resonance.Pz() / mc_resonance.E());

    new ((*m_MC_Top1_p4)[0]) TLorentzVector(Top[0]);
    new ((*m_MC_Top2_p4)[0]) TLorentzVector(Top[1]);
  }

  if (m_MC_channel != 1 && m_MC_channel != 2) {
    return;
  }

  // Extract index of semi-leptonic event, and store them in tree. Useful if you want to know how many jets you have selected right
  if (false) {
    std::cout << "New event" << std::endl;
    for (int i = 0; i < n_MC; i++) {
      std::cout << "\t[" << i << "] Type: " << m_MC->getType(i) << std::endl;
    }
  }

  bool keepEvent = true;
  for (int i = 0; i < n_MC; i++) {

    int motherIndex = patIndexToExtractorIndex(m_MC->getMom1Index(i));
    int grandMotherIndex = -1;
    if (motherIndex != -1)
      grandMotherIndex = patIndexToExtractorIndex(m_MC->getMom1Index(motherIndex));

    if (motherIndex == -1 || grandMotherIndex == -1)
      continue;

    // Look only event coming (directly / indirectly) from a top
    if (abs(m_MC->getType(motherIndex)) == ID_T || abs(m_MC->getType(grandMotherIndex)) == ID_T)  {

      int type = abs(m_MC->getType(i));
      // W? Continue
      if (type == ID_W)
        continue;

      // Only semi-mu or semi-e events are interesting, so throw away event with a tau
      if (type == ID_TAU) {
        keepEvent = false;
        break;
      }

      switch (type) {
        case ID_E:
          if (m_leptonIndex != -1) {
            keepEvent = false;
            break;
          }
          m_leptonIndex = i;
          break;

        case ID_MU:
          if (m_leptonIndex != -1) {
            keepEvent = false;
            break;
          }
          m_leptonIndex = i;
          break;

        case ID_NEUTRINO_E:
        case ID_NEUTRINO_MU:
        case ID_NEUTRINO_TAU:
          if (m_neutrinoIndex != -1) {
            keepEvent = false;
            break;
          }
          m_neutrinoIndex = i;

          // In some case, Madgraph does not output a W in the LHE
          if (abs(m_MC->getType(grandMotherIndex)) == ID_T)
            m_leptonicTopIndex = grandMotherIndex;
          else
            m_leptonicTopIndex = motherIndex;

          break;

        case ID_B:
          if (m_leptonicBIndex == -1) {
            m_leptonicBIndex = i;
          } else {
            if (m_hadronicBIndex != -1) {
              keepEvent = false;
              break;
            }
            m_hadronicBIndex = i;
          }
          break;

        default: // Other jets
          if (m_firstJetIndex == -1) {
            m_firstJetIndex = i;
          } else {
            if (m_secondJetIndex != -1) {
              keepEvent = false;
              break;
            }
            m_secondJetIndex = i;
          }
          break;
      }

      if (! keepEvent)
        break;
    }
  }

  if (m_leptonIndex == -1 || m_neutrinoIndex == -1 || m_leptonicBIndex == -1 || m_hadronicBIndex == -1 || m_firstJetIndex == -1 || m_secondJetIndex == -1)
    keepEvent = false;

  if (! keepEvent) {
    m_leptonIndex = m_leptonicBIndex = m_hadronicBIndex = m_neutrinoIndex = m_firstJetIndex = m_secondJetIndex = m_leptonicTopIndex = -1;
    return;
  }

  // Reorder B jet indexes
  if (patIndexToExtractorIndex(m_MC->getMom1Index(m_leptonicBIndex)) != m_leptonicTopIndex) {
    // Wrong combinaison, swap
    std::swap(m_leptonicBIndex, m_hadronicBIndex);
  }

  if (false) {
    std::cout << "Lepton index: " << m_leptonIndex << std::endl;
    std::cout << "Neutrino index: " << m_neutrinoIndex << std::endl;
    std::cout << "Leptonic B index: " << m_leptonicBIndex << std::endl;
    std::cout << "Hadronic B index: " << m_hadronicBIndex << std::endl;
    std::cout << "First jet index: " << m_firstJetIndex << std::endl;
    std::cout << "Second jet index: " << m_secondJetIndex << std::endl;
  }

  // First, check if the event is associable. It is if each parton from the
  // tt system (lepton, b jets & light jets) have a associated RECO object

  // Only check jets
  m_mtt_eventIsAssociable =
    hasRecoPartner(m_leptonicBIndex) &&
    hasRecoPartner(m_hadronicBIndex) &&
    hasRecoPartner(m_firstJetIndex) &&
    hasRecoPartner(m_secondJetIndex);

  // Compute masses
  TLorentzVector mc_hadr_W = *m_MC->getP4(m_firstJetIndex) + *m_MC->getP4(m_secondJetIndex);
  m_MC_hadronicWMass = mc_hadr_W.M();

  TLorentzVector mc_hadr_top = mc_hadr_W + *m_MC->getP4(m_hadronicBIndex);
  m_MC_hadronicTopMass = mc_hadr_top.M();

  TLorentzVector mc_lept_W = *m_MC->getP4(m_neutrinoIndex) + *m_MC->getP4(m_leptonIndex);
  m_MC_leptonicWMass = mc_lept_W.M();

  TLorentzVector mc_lept_top = mc_lept_W + *m_MC->getP4(m_leptonicBIndex);
  m_MC_leptonicTopMass = mc_lept_top.M();

  // Store ref to various P4
  new ((*m_MC_lepton_p4)[0]) TLorentzVector(*m_MC->getP4(m_leptonIndex));
  new ((*m_MC_neutrino_p4)[0]) TLorentzVector(*m_MC->getP4(m_neutrinoIndex));

  new ((*m_MC_leptonic_B_p4)[0]) TLorentzVector(*m_MC->getP4(m_leptonicBIndex));
  new ((*m_MC_hadronic_B_p4)[0]) TLorentzVector(*m_MC->getP4(m_hadronicBIndex));

  new ((*m_MC_lightJet1_p4)[0]) TLorentzVector(*m_MC->getP4(m_firstJetIndex));
  new ((*m_MC_lightJet2_p4)[0]) TLorentzVector(*m_MC->getP4(m_secondJetIndex));
}

bool mtt_analysis::hasRecoPartner(int mcIndex) const {

  // Loop over all jets, and check if one has a gen particle
  // equals to mcIndex. Check also that the matched reco jet
  // pass our jet selection

  for (uint32_t i = 0; i < m_jetMet->getSize() ; i++) {
    if (m_jetMet->getJetMCIndex(i) == mcIndex) {

      if (fabs(m_jetMet->getP4(i)->Pt()) < m_JET_Pt_min)
        continue;

      if (fabs(m_jetMet->getP4(i)->Eta()) > m_JET_Eta_max)
        continue;

      return true;
    }
  }

  return false;
}

bool mtt_analysis::jetComesFromTTDecay(int mcIndex) const {
  return
    mcIndex == m_leptonicBIndex ||
    mcIndex == m_hadronicBIndex ||
    mcIndex == m_firstJetIndex ||
    mcIndex == m_secondJetIndex;
}

bool mtt_analysis::isSolutionMatched(uint32_t leptonicBIndex, uint32_t hadronicBIndex, uint32_t hadronicFirstJetIndex,
    uint32_t hadronicSecondJetIndex) {

  return (
      m_jetMet->getJetMCIndex(leptonicBIndex) == m_leptonicBIndex &&

      // Check hadronic B
      m_jetMet->getJetMCIndex(hadronicBIndex) == m_hadronicBIndex &&

      // Check light jets
      (
       (
        m_jetMet->getJetMCIndex(hadronicFirstJetIndex) == m_firstJetIndex &&
        m_jetMet->getJetMCIndex(hadronicSecondJetIndex) == m_secondJetIndex
       )
       ||
       (
        m_jetMet->getJetMCIndex(hadronicFirstJetIndex) == m_secondJetIndex &&
        m_jetMet->getJetMCIndex(hadronicSecondJetIndex) == m_firstJetIndex
       )
      )
      );
}

void mtt_analysis::checkIfSolutionIsCorrect() {

  // Check if the four select jets come from tt decay.
  // Position is not important, as long as we have the four
  //if (m_useChi2) {
    //m_mtt_recoJetsAssociatedWithChi2 = (
        //jetComesFromTTDecay(m_jetMet->getJetMCIndex(m_selectedLeptonicBIndex_AfterChi2)) &&
        //jetComesFromTTDecay(m_jetMet->getJetMCIndex(m_selectedHadronicBIndex_AfterChi2)) &&
        //jetComesFromTTDecay(m_jetMet->getJetMCIndex(m_selectedHadronicFirstJetIndex_AfterChi2)) &&
        //jetComesFromTTDecay(m_jetMet->getJetMCIndex(m_selectedHadronicSecondJetIndex_AfterChi2))
        //);

    //m_mtt_recoJetsAssociatedWellPlacedWithChi2 = isSolutionMatched(m_selectedLeptonicBIndex_AfterChi2, m_selectedHadronicBIndex_AfterChi2, m_selectedHadronicFirstJetIndex_AfterChi2, m_selectedHadronicSecondJetIndex_AfterChi2);
  //}

  //if (m_useMVA) {
    //m_mtt_recoJetsAssociatedWithMVA = (
        //jetComesFromTTDecay(m_jetMet->getJetMCIndex(m_selectedLeptonicBIndex_AfterMVA)) &&
        //jetComesFromTTDecay(m_jetMet->getJetMCIndex(m_selectedHadronicBIndex_AfterMVA)) &&
        //jetComesFromTTDecay(m_jetMet->getJetMCIndex(m_selectedHadronicFirstJetIndex_AfterMVA)) &&
        //jetComesFromTTDecay(m_jetMet->getJetMCIndex(m_selectedHadronicSecondJetIndex_AfterMVA))
        //);

    //m_mtt_recoJetsAssociatedWellPlacedWithMVA = isSolutionMatched(m_selectedLeptonicBIndex_AfterMVA, m_selectedHadronicBIndex_AfterMVA, m_selectedHadronicFirstJetIndex_AfterMVA, m_selectedHadronicSecondJetIndex_AfterMVA);
  //}
}

void mtt_analysis::fillTree()
{
  m_lepton_weight_error_low = sqrt(m_lepton_weight_error_low);
  m_lepton_weight_error_high = sqrt(m_lepton_weight_error_high);

  m_btag_weight_error_low = sqrt(m_btag_weight_error_low);
  m_btag_weight_error_high = sqrt(m_btag_weight_error_high);

  m_tree_Mtt->Fill();
}

// Here we just reset the ROOTtree parameters

void mtt_analysis::reset()
{
  m_pass_vertex_cut = -1;
  m_pass_met_cut = -1;
  m_pass_lepton_cut = -1;
  m_pass_jet_cut = -1;

  m_mtt_isSel = 0;
  m_mtt_eventIsAssociable = false;

  m_MC_lepton_p4->Clear("C");
  m_MC_neutrino_p4->Clear("C");

  m_MC_leptonic_B_p4->Clear("C");
  m_MC_hadronic_B_p4->Clear("C");

  m_MC_lightJet1_p4->Clear("C");
  m_MC_lightJet2_p4->Clear("C");

  m_MC_Top1_p4->Clear("C");
  m_MC_Top2_p4->Clear("C");

  m_selectedLeptonIndex_in_loose_collection     = -1;
  m_selectedLeptonIndex_in_array                = -1;

  m_mtt_NGoodMuons = 0;
  m_mtt_NGoodElectrons = 0;

  for (int tmp = 0; tmp < 20; ++tmp)
  {
    m_mtt_MuonPt[tmp] = 0.;
    m_mtt_MuonEta[tmp] = 1000.;
    m_mtt_MuRelIso[tmp] = -1.;

    m_mtt_ElectronPt[tmp] = 0.;
    m_mtt_ElectronEta[tmp] = 1000.;
    m_mtt_ElRelIso[tmp] = -1.;
  }


  m_mtt_NJets = 0;
  m_mtt_NBtaggedJets_CSVL = 0;
  m_mtt_NBtaggedJets_CSVM = 0;
  m_mtt_NBtaggedJets_CSVT = 0;
  m_mtt_NBtaggedJets_TCHPT = 0;

  m_mtt_1stjetpt = 0.;
  m_mtt_2ndjetpt = 0.;
  m_mtt_3rdjetpt = 0.;
  m_mtt_4thjetpt = 0.;

  for (int tmp = 0; tmp < 100; ++tmp)
  {
    m_mtt_JetEta[tmp]     = 1000.;
    m_mtt_JetPt[tmp]      = 0.;
  }

  m_mtt_MET = 0.;
  m_nPU        = 0.;
  m_MC_channel = 0;
  m_MC_mtt     = -1.;

  m_leptonIndex = -1;
  m_neutrinoIndex = -1;

  m_leptonicBIndex = -1;
  m_hadronicBIndex = -1;
  m_leptonicTopIndex = -1;

  m_firstJetIndex = -1;
  m_secondJetIndex = -1;

  m_trigger_passed = false;

  m_MC_pt_tt         = -1;
  m_MC_eta_tt        = -1;
  m_MC_beta_tt       = -1;
  m_MC_hadronicWMass = -1;
  m_MC_leptonicWMass = -1;
  m_MC_hadronicTopMass = -1;
  m_MC_leptonicTopMass = -1;

  m_selectedLeptonP4->Clear("C");

  m_lepton_weight = 1.;
  m_lepton_weight_error_low = 0.;
  m_lepton_weight_error_high = 0.;

  m_btag_weight = 1.;
  m_btag_weight_error_low = 0.;
  m_btag_weight_error_high = 0.;

  for (auto& algo: m_sorting_algortihms)
    algo->reset();
}

bool mtt_analysis::isBJet(unsigned int index) {
  // Use recommanded WP from https://indico.cern.ch/getFile.py/access?contribId=4&resId=2&materialId=slides&confId=195042
  return m_jetMet->getJetBTagProb_CSV(index) > m_JET_btag_CSVM;
}

}

DEFINE_EDM_PLUGIN(PatExtractorPluginFactory, mtt::mtt_analysis, "mtt_analysis");
