#pragma once

#include <vector>

#include "sequence_mapper.hpp"
#include "sequence/sequence.hpp"

namespace debruijn_graph {

//base class for all reference checks
class ReferenceChecker {
  public:
      ReferenceChecker() {}

      virtual void Check(MappingPath<EdgeId>& path, size_t i) {}
};

class IndelChecker : public ReferenceChecker {
  public:
      IndelChecker() {}

      void Check(MappingPath<EdgeId>& path, size_t i) {
          INFO("Check in IndelChecker called");
      }
};

template<class Graph, class Mapper>
class ReferenceCorrector {
    typedef typename Graph::EdgeId EdgeId;
    conj_graph_pack &gp;
    const Graph& g_;
    Mapper mapper_;
    EdgesPositionHandler<Graph>& edge_pos_;
    std::vector<std::unique_ptr<ReferenceChecker> > checkers_;
  public:
      ReferenceCorrector(conj_graph_pack &gp_) :
              gp(gp_), g_(gp_.g),
              mapper_(gp_.g, gp_.index, gp_.kmer_mapper),
              edge_pos_(gp_.edge_pos) {
          INFO("in ReferenceCorrector constructor");

      }

      void add(ReferenceChecker * checker) {
          checkers_.push_back(std::unique_ptr<ReferenceChecker>(checker));
          //checkers_.back()->parent_ = this;
      }

      /** method loads reference file or data and processes all */
      void LoadReference() {

          //1. Load reference to memory somehow
          //2, call Process
          //
      }

      void Process(const Sequence& sequence) const {
          INFO("in LoadReference");
          MappingPath<EdgeId> path = mapper_.MapSequence(sequence);
          int cur_pos = 0;
          TRACE("Sequence mapped on " << path.size()
              << " fragments.");
          for (size_t i = 0; i < path.size(); i++) {
              for (auto iter = checkers_.begin(); iter != checkers_.end(); ++ iter) {
                  (*iter)->Check(path, i);
              }
          }
          /*for (size_t i = 0; i < path.size(); i++) {
              EdgeId ei = path[i].first;
              MappingRange mr = path[i].second;
              int len = (int) (mr.mapped_range.end_pos - mr.mapped_range.start_pos);
              if (i > 0)
                  if (path[i - 1].first != ei)
                      if (g_.EdgeStart(ei) != g_.EdgeEnd(path[i - 1].first)) {
                          TRACE(
                              "Contig " << name
                              << " mapped on not adjacent edge. Position in contig is "
                              << path[i - 1].second.initial_range.start_pos
                              + 1
                              << "--"
                              << path[i - 1].second.initial_range.end_pos
                              << " and "
                              << mr.initial_range.start_pos + 1
                              << "--" << mr.initial_range.end_pos);
                      }
              edge_pos_.AddEdgePosition(ei, name, mr.initial_range.start_pos,
                                        mr.initial_range.end_pos,
                                        mr.mapped_range.start_pos,
                                        mr.mapped_range.end_pos);
              cur_pos += len;
          }*/
      }
};
}
