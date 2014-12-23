#pragma once

#include "graph_pack.hpp"

#include "stage.hpp"
#include "reference_corrector.hpp"

namespace debruijn_graph {

class ReferenceCorrection : public spades::AssemblyStage {
  public:
    ReferenceCorrection()
        : AssemblyStage("Reference Correction", "reference_correction") {}

    void run(conj_graph_pack &gp, const char*);
};

}
