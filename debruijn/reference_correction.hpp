#pragma once

#include "stage.hpp"

namespace debruijn_graph {

class ReferenceCorrection : public spades::AssemblyStage {
  public:
    ReferenceCorrection()
        : AssemblyStage("Reference Correction", "reference_correction") {}

    void run(conj_graph_pack &gp, const char*);
};

}
