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

class SortingAlgorithm {
  public:
    void setObjects(const std::vector<Jet>& jets, const Lepton& lepton, const LorentzVector& met);
    virtual void work() = 0;
    virtual void reset() {};
    virtual void addBranches(TTree& tree) = 0;

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

  protected:
    bool computeNeutrinoPz(const LorentzVector& leptonicBjet, bool* no_real_sol = nullptr);

    std::vector<Jet> m_jets;
    Lepton m_lepton;
    LorentzVector m_met;
    ROOT::Math::PxPyPzMVector m_neutrino;
};
