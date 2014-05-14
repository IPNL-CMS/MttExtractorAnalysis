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

  protected:
    bool computeNeutrinoPz(const LorentzVector& leptonicBjet, bool* no_real_sol = nullptr);

    std::vector<Jet> m_jets;
    Lepton m_lepton;
    LorentzVector m_met;
    ROOT::Math::PxPyPzMVector m_neutrino;
};