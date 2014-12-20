
#include "reference_correction.hpp"


namespace debruijn_graph {

typedef debruijn_graph::NewExtendedSequenceMapper<debruijn_graph::Graph, Index> MapperClass;


void ReferenceCorrection::run(conj_graph_pack &gp, const char*) {
    INFO("ReferenceCorrection started");
    //TODO: to do
    ReferenceCorrector<Graph, MapperClass> corrector(gp);
	int coverage_threshold = 10;
    corrector.add(new IndelChecker<Graph>(gp.g, coverage_threshold));
	// probably deprecated
    //corrector.add(new MobileElementInserionChecker<Graph>(gp, gp.g));
    corrector.Process(cfg::get().ds.reference_genome);

    INFO("ReferenceCorrection ended");
}

}
