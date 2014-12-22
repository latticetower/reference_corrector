
#include "reference_correction.hpp"


namespace debruijn_graph {

typedef debruijn_graph::NewExtendedSequenceMapper<debruijn_graph::Graph, Index> MapperClass;


void ReferenceCorrection::run(conj_graph_pack &gp, const char*) {
    INFO("ReferenceCorrection started");
    gp.EnsureDebugInfo();
    //gp.EnsureBasicMapping();
    //TODO: to do
    ReferenceCorrector<Graph, MapperClass> corrector(gp);

    //corrector.add(new IndelChecker<Graph>(gp.g));
    int coverage_threshold = 10;
    corrector.add(new MobileElementInserionChecker<Graph>(gp, coverage_threshold));
    corrector.Process(cfg::get().ds.reference_genome);

    INFO("ReferenceCorrection ended");
}

}
