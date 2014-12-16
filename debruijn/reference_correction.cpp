
#include "reference_correction.hpp"


namespace debruijn_graph {

typedef debruijn_graph::NewExtendedSequenceMapper<debruijn_graph::Graph, Index> MapperClass;


void ReferenceCorrection::run(conj_graph_pack &gp, const char*) {
    INFO("ReferenceCorrection started");
    //TODO: to do
    ReferenceCorrector<Graph, MapperClass> corrector(gp);
    corrector.add(new IndelChecker<Graph>());
    corrector.Process(cfg::get().ds.reference_genome);
    INFO("ReferenceCorrection ended");
}

}
