#pragma once

#include "Chi2SortingAlgorithm.h"

#include <string>
#include <memory>

#include <FWCore/ParameterSet/interface/ParameterSet.h>

class SortingAlgorithmFactory {
  public:
    static std::shared_ptr<SortingAlgorithm> create(const std::string& name, const edm::ParameterSet& cfg, bool isSemiMu) {
      if (name == "chi2") {
        return std::make_shared<Chi2SortingAlgorithm>(cfg, isSemiMu);
      }

      return std::shared_ptr<SortingAlgorithm>();
    };
};
