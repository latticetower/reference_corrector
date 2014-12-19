#pragma once

#include <vector>

#include "standard_base.hpp"
#include "sequence_mapper.hpp"
#include "sequence/sequence.hpp"

#include "omni/path_processor.hpp"
#include "omni/visualization/visualization_utils.hpp"
#include "omni/visualization/graph_colorer.hpp"

#include "stats/debruijn_stats.hpp"

namespace debruijn_graph {

//base class for all reference checks
template<class Graph>
class ReferenceChecker {
      typedef typename Graph::EdgeId EdgeId;
  public:
      ReferenceChecker() {}

      virtual void Check(MappingPath<EdgeId>& , size_t ) {}
      virtual void Write(MappingPath<EdgeId>& ) {}
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
      void Write(MappingPath<EdgeId>& ) {
      }
};

/** mobile element insertions detector */
template<class Graph>
class MobileElementInserionChecker : public ReferenceChecker<Graph> {
      typedef typename Graph::EdgeId EdgeId;
      typedef typename Graph::VertexId VertexId;

      conj_graph_pack & gp_;
      Graph & g_;
      //TODO: add some storage for data
  public:
      MobileElementInserionChecker(conj_graph_pack&gp, Graph& g) : gp_(gp), g_(g) {}

      void Check(MappingPath<EdgeId>& path, size_t i) {
          INFO("Check in MobileElementInserionChecker called");
          if (i == 0)
              return;
          EdgeId prev_edge = path[i - 1].first;
          EdgeId edge = path[i].first;
          if (edge == prev_edge)
              return; // not impl
          if (g_.EdgeStart(edge) == g_.EdgeEnd(prev_edge))
              return;
          VertexId start_v = g_.EdgeEnd(prev_edge);
          VertexId end_v = g_.EdgeStart(edge);

          //FIX: the following code is from pac_index.hpp. should modify for finding mobile elements
          PathStorageCallback<Graph> callback(g_);
          PathProcessor<Graph> path_processor(g_, 0, 4000, start_v, end_v, callback);
          //copypasted prev line from pac_index.hpp. still don't know what 0 and 4000 mean
          path_processor.Process();
          vector<vector<EdgeId> > paths = callback.paths();
          stringstream s_buf;
          for (auto p_iter = paths.begin(); p_iter != paths.end(); p_iter++) {
              size_t tlen = 0;
              for (auto path_iter = p_iter->begin();
                      path_iter != p_iter->end();
                      path_iter++) {
                  tlen += g_.length(*path_iter);
              }
              s_buf << tlen << " ";
          }
          DEBUG(s_buf.str());
          //TODO: instead of simply output to console information about path lengths, should check if edges are near or smth

          INFO("Check in IndelChecker called - exiting from check");
      }
      void Write(MappingPath<EdgeId>& path) {
          LengthIdGraphLabeler<Graph> basic_labeler(gp_.g);
          EdgePosGraphLabeler<Graph> pos_labeler(gp_.g, gp_.edge_pos);
          CompositeLabeler<Graph> labeler(basic_labeler, pos_labeler);

          auto edge_colorer = omnigraph::visualization::DefaultColorer(g_);

          WriteComponentsAlongPath(g_, path.path(), "reference_alterations/", edge_colorer, labeler);
          //  auto edge_colorer = make_shared<CompositeEdgeColorer<Graph>>("black");

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
          for (auto iter = checkers_.begin(); iter != checkers_.end(); ++ iter) {
              (*iter)->Write(path);
          }
      }
};
}
