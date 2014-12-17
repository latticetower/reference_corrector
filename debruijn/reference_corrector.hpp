#pragma once

#include <vector>

#include "sequence_mapper.hpp"
#include "sequence/sequence.hpp"

namespace debruijn_graph {

//base class for all reference checks
template<class Graph>
class ReferenceChecker {
      typedef typename Graph::EdgeId EdgeId;
  public:
      ReferenceChecker() {}

      virtual void Check(MappingPath<EdgeId>& path, size_t i) {}
};

struct IndelInfo{
    size_t indel_position, start_pos, end_pos;
    IndelInfo(size_t const & i, size_t const & s, size_t const & e) : indel_position(i), start_pos(s), end_pos(e) {

    }
};

/** Very simple class, checks only indels inside edges and saves
  * them to 2 vectors (for insertions and deletions separately)
  **/
template<class Graph>
class IndelChecker : public ReferenceChecker<Graph> {
      typedef typename Graph::EdgeId EdgeId;
      Graph & g_;
      vector<IndelInfo> insertions, deletions;
      //these two vectors contain positions in reference
  public:
      IndelChecker(Graph& g): g_(g) {}

      void Check(MappingPath<EdgeId>& path, size_t i) {
          INFO("Check in IndelChecker called");
          EdgeId ei = path[i].first;
          MappingRange mr = path[i].second;

          //TODO: check and add +1 values if needed
          //check insertion in reference:
          if (path[i - 1].second.initial_range.end_pos == mr.initial_range.start_pos &&
              path[i - 1].second.mapped_range.end_pos != mr.mapped_range.start_pos) {
                //
                insertions.push_back(IndelInfo(
                        mr.initial_range.start_pos,
                        path[i - 1].second.mapped_range.end_pos,
                        mr.mapped_range.start_pos));
          }
          //check deletions in reference:
          if (path[i - 1].second.initial_range.end_pos != mr.initial_range.start_pos &&
              path[i - 1].second.mapped_range.end_pos == mr.mapped_range.start_pos) {
                //
                deletions.push_back(IndelInfo(mr.mapped_range.start_pos,
                      path[i - 1].second.initial_range.end_pos,
                      mr.initial_range.start_pos));
          }
          INFO("Check in IndelChecker called - exiting from check");
      }
};

/** mobile element insertions detector */
template<class Graph>
class MobileElementInserionChecker : public ReferenceChecker<Graph> {
      typedef typename Graph::EdgeId EdgeId;
      Graph & g_;
      //TODO: add some storage for data
  public:
      MobileElementInserionChecker(Graph& g): g_(g) {}

      void Check(MappingPath<EdgeId>& path, size_t i) {
          INFO("Check in MobileElementInserionChecker called");
          EdgeId ei = path[i].first;
          MappingRange mr = path[i].second;

          //TODO: add conditions

          INFO("Check in IndelChecker called - exiting from check");
      }
};


template<class Graph, class Mapper>
class ReferenceCorrector {
    typedef typename Graph::EdgeId EdgeId;
    conj_graph_pack &gp;
    const Graph& g_;
    Mapper mapper_;
    EdgesPositionHandler<Graph>& edge_pos_;
    std::vector<std::unique_ptr<ReferenceChecker<Graph> > > checkers_;
  public:
      ReferenceCorrector(conj_graph_pack &gp_) :
              gp(gp_), g_(gp_.g),
              mapper_(gp_.g, gp_.index, gp_.kmer_mapper),
              edge_pos_(gp_.edge_pos) {
          INFO("in ReferenceCorrector constructor");

      }

      void add(ReferenceChecker<Graph> * checker) {
          checkers_.push_back(std::unique_ptr<ReferenceChecker<Graph> >(checker));
          //checkers_.back()->parent_ = this;
      }

      void Process(const Sequence& sequence) const {
          INFO("in LoadReference");
          MappingPath<EdgeId> path = mapper_.MapSequence(sequence);
          INFO("Sequence mapped on " << path.size()
              << " fragments.");
          for (size_t i = 0; i < path.size(); i++) {
              for (auto iter = checkers_.begin(); iter != checkers_.end(); ++ iter) {
                  (*iter)->Check(path, i);
              }
          }
      }
};
}
