#pragma once

#include <Math/Vector4D.h> 
#include <vector>

typedef ROOT::Math::PtEtaPhiEVector LorentzVector;

class TTree;

struct Jet {
  LorentzVector p;
  int index;
  bool isBTagged;

  Jet() {
    index = -1;
    isBTagged = false;
  }
};

struct Lepton {
  LorentzVector p;
  int charge;

  Lepton() {
    charge = 0;
  }
};

const int nBins_for_gaussian_profile = 12;
const double bins_for_gaussian_profile[] = {340, 400, 450, 500, 550, 600, 650, 750, 850, 950, 1050, 1200, 1400};

class SortingAlgorithm {
  public:
    void setObjects(const std::vector<Jet>& jets, const Lepton& lepton, const LorentzVector& met);
    virtual void work() = 0;
    virtual void reset() {};
    virtual void addBranches(TTree& tree) = 0;
    virtual void end() {};

    virtual int getSelectedLeptonicBIndex() {
      return -1;
    };

    virtual int getSelectedHadronicBIndex() {
      return -1;
    };

    virtual int getSelectedHadronicFirstJetIndex() {
      return -1;
    };

    virtual int getSelectedHadronicSecondJetIndex() {
      return -1;
    };

    virtual void setRecoJetsAssociated(const bool& recoJetsAssociated) {};
    virtual void setRecoJetsAssociatedWellPlaced(const bool& recoJetsAssociatedWellPlaced) {};
    void setMttGen(float mtt_gen) {
      m_mtt_gen = mtt_gen;
    };

  protected:
    bool computeNeutrinoPz(const LorentzVector& leptonicBjet, bool* no_real_sol = nullptr);

    float m_mtt_gen;

    std::vector<Jet> m_jets;
    Lepton m_lepton;
    LorentzVector m_met;
    ROOT::Math::PxPyPzMVector m_neutrino;
};
