#pragma once

#include <vector>

#include "standard_base.hpp"
#include "sequence_mapper.hpp"
#include "sequence/sequence.hpp"

#include "omni/path_processor.hpp"
#include "omni/visualization/visualization_utils.hpp"
#include "omni/visualization/graph_colorer.hpp"

#include "omni/coverage.hpp" // for CoverageIndex

#include "stats/debruijn_stats.hpp"

namespace debruijn_graph {

template<typename Graph>
using MappingElement = std::pair<const typename Graph::EdgeId, const MappingRange>;


/** helper method to draw path, smth, smth
  * omni/visualization/graph_colorer.hpp
  */
template <class Graph>
shared_ptr<omnigraph::visualization::GraphColorer<Graph>> DefaultColorer(const Graph& g,
    const vector<typename Graph::EdgeId>& graph_path,
    const vector<typename Graph::EdgeId>& reference_path) {
      using namespace omnigraph::visualization;
  shared_ptr<ElementColorer<typename Graph::EdgeId>> edge_colorer =
            make_shared<CompositeEdgeColorer<Graph>>(
                    make_shared<SetColorer<Graph>>(g, graph_path, "blue"),
                    make_shared<SetColorer<Graph>>(g, reference_path, "green"), "black");
  return DefaultColorer(g, edge_colorer);
}
/** method saves current path in graph and reference fragment to given file.
  * green edge == reference fragment, blue - current path, black - other edges in graph
  */
template <typename Graph>
void WritePathFragment(
        conj_graph_pack & gp_,
        vector<EdgeId> & graph_path,
        vector<EdgeId> & reference_path,
        string file_name = "fragment.dot") {

    INFO("Writing path fragment")
    LengthIdGraphLabeler<Graph> basic_labeler(gp_.g);
    EdgePosGraphLabeler<Graph> pos_labeler(gp_.g, gp_.edge_pos);
    CompositeLabeler<Graph> labeler(basic_labeler, pos_labeler);
    auto edge_colorer = DefaultColorer(gp_.g, graph_path, reference_path);

    set<VertexId> vertices;
    for (auto iter = graph_path.begin(); iter!= graph_path.end(); ++iter) {
        vertices.insert(gp_.g.EdgeStart(*iter));
        vertices.insert(gp_.g.EdgeEnd(*iter));
    }

    for (auto iter = reference_path.begin(); iter != reference_path.end(); ++iter) {
        vertices.insert(gp_.g.EdgeStart(*iter));
        vertices.insert(gp_.g.EdgeEnd(*iter));
    }
    GraphComponent<Graph> component(gp_.g, vertices.begin(), vertices.end());
    omnigraph::visualization::WriteComponent(component, file_name, edge_colorer, labeler);
}


//base class for all reference checks
template<class Graph>
class ReferenceChecker {
      typedef typename Graph::EdgeId EdgeId;
  public:
      ReferenceChecker() {}

      virtual void Check(MappingPath<EdgeId>& , size_t ) {}
      virtual void Write(MappingPath<EdgeId>& , string ) {}
};

struct IndelInfo{
    size_t indel_position, start_pos, end_pos;
    IndelInfo(size_t const & i, size_t const & s, size_t const & e) :
            indel_position(i), start_pos(s), end_pos(e) {

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
          if (i == 0)
              return;
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
                deletions.push_back(IndelInfo(mr.mapped_range.start_pos,
                      path[i - 1].second.initial_range.end_pos,
                      mr.initial_range.start_pos));
          }
          INFO("Check in IndelChecker called - exiting from check");
      }
      void Write(MappingPath<EdgeId>& , string ) {
      }
      DECL_LOGGER("MobileElementInserionChecker");//this class shouldn't be called now, added logger to check it
};

template<class Graph>
struct MobileElementInfo {
    typename Graph::EdgeId start_edge, end_edge;
    vector<typename Graph::EdgeId> & mobile_element_edges;

    template<class EdgeId>
    MobileElementInfo(EdgeId s, EdgeId e, vector<typename Graph::EdgeId> & edges): start_edge(s), end_edge(e), mobile_element_edges(edges) {

    }
};

/** mobile element insertions detector */
template<class Graph>
class MobileElementInserionChecker : public ReferenceChecker<Graph> {
      typedef typename Graph::EdgeId EdgeId;
      typedef typename Graph::VertexId VertexId;

      conj_graph_pack & gp_;
      Graph & g_;
      //vector<std::shared_ptr<MobileElementInfo<Graph> > > mobile_elements_;
      int coverage_threshold_;
  public:
    MobileElementInserionChecker(conj_graph_pack& gp, int coverage_threshold = 0) :
            gp_(gp), g_(gp.g), coverage_threshold_(coverage_threshold) {}

    void Check(MappingPath<EdgeId>& path, size_t i) {
        //INFO("Check in MobileElementInserionChecker called");
        if (i == 0)
            return;
            using namespace omnigraph::visualization;
        EdgeId prev_edge = path[i - 1].first;
        EdgeId edge = path[i].first;
        if (edge == prev_edge)
            return; // not impl
        if (g_.EdgeStart(edge) == g_.EdgeEnd(prev_edge))
            return;

        CollectMobileElements(prev_edge, edge);

        //TODO: instead of simply output to console information about path lengths,
        //should check if edges are near or smth

        //INFO("Check in IndelChecker called - exiting from check");
    }


    void Write(MappingPath<EdgeId>& , string output_dir ) {
        /*
        make_dir(output_dir);
        INFO("Write gets called");
        int k = 0;
        for (auto & me : mobile_elements_) {
            vector<EdgeId> reference_path;
            reference_path.push_back(me->start_edge);
            reference_path.push_back(me->end_edge);
            string file_name = output_dir + "fragment_" + ToString(k) + "/.dot";
            k++;
            WritePathFragment<Graph>(gp_, me->mobile_element_edges, me->mobile_element_edges, file_name);
        }
        */
    }

  private:


    void WriteMEToFile(vector<EdgeId>& path, int i, EdgeId start_edge, EdgeId end_edge, string output_dir = "reference_correction/") {
        make_dir(output_dir);
        INFO("WriteMEToFile gets called");
        INFO("info: " << i << gp_.g.int_id(start_edge) << " " << gp_.g.int_id(end_edge));
        //WritePathFragment(path, reference_path, file_name);
    }
    void WriteME(
                 vector<EdgeId>& path,
                 int i,
                 EdgeId start_edge,
                 EdgeId end_edge,
                 string output_dir = "") {

        INFO("WriteME gets called");
        vector<EdgeId> reference_path;
        reference_path.push_back(start_edge);
        reference_path.push_back(end_edge);
        string file_name = output_dir + "fragment_" + ToString(i) + ".dot";
        INFO(file_name);
        WritePathFragment<Graph>(gp_, path, reference_path, file_name);
    }




    /** checks all the paths between start and end.
      * if an edge on a path has coverage greater than coverage_threshold, it will be stored
      * in mobile_elements vector
     */
    void CollectMobileElements(EdgeId start_edge, EdgeId end_edge) {
        if (start_edge == end_edge)
            return;
        make_dir("reference_correction/");
        VertexId start = g_.EdgeStart(start_edge);
        VertexId end = g_.EdgeEnd(end_edge);

        PathStorageCallback<Graph> callback(g_);
        PathProcessor<Graph> path_processor(g_, 0, 4000, start, end, callback);
        //copypasted prev line from pac_index.hpp. still don't know what 0 and 4000 mean
        path_processor.Process();
        vector<vector<EdgeId>> paths = callback.paths();
        CoverageIndex<Graph> coverage_index(g_);
        int k = 0;
        bool new_mobile_element = true;

        for (vector<EdgeId>& path : paths) {
            //TODO: check logic!
            for (EdgeId edge : path) {
                INFO("coverage for current edge is " << coverage_index.coverage(edge));
                if (coverage_index.coverage(edge) >= coverage_threshold_) {
                    if (new_mobile_element) {
                        //WriteME()
                        k++;
                        WriteMEToFile(path, k, start_edge, end_edge);
                        WriteME(path, k, start_edge, end_edge);

                        /*mobile_elements_.push_back(
                            make_shared<MobileElementInfo<Graph>>(
                                start_edge, end_edge, path));
                                */
                        new_mobile_element = false;
                        //mobile_elements_.back()->mobile_element_edges.push_back(edge);

                    }
                    //mobile_elements_.back()->mobile_element_edges.push_back(edge);
                }
                //else {

                //}
            }
            new_mobile_element = true;
        }
    }


    DECL_LOGGER("MobileElementInserionChecker");
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
        using namespace omnigraph::visualization;
          INFO("in LoadReference");
          MappingPath<EdgeId> path = mapper_.MapSequence(sequence);
          INFO("Sequence mapped on " << path.size() << " fragments.");
              //FIX: nxt call to WriteComponentsAlongPath is temp
              /*LengthIdGraphLabeler<Graph> basic_labeler(gp.g);
              EdgePosGraphLabeler<Graph> pos_labeler(gp.g, gp.edge_pos);
              string file_name =  "reference_correction\\fragment_.dot";
              INFO("gp.index " << gp.index.IsAttached());
  CompositeLabeler<Graph> labeler(basic_labeler, pos_labeler);
              omnigraph::visualization::WriteComponentsAlongPath(gp.g, path.path(),
                  "reference_correction",
                  omnigraph::visualization::DefaultColorer(gp.g), labeler);
                  INFO("OK");
                  */
          for (size_t i = 0; i < path.size(); i++) {
              for (auto iter = checkers_.begin(); iter != checkers_.end(); ++ iter) {
                  (*iter)->Check(path, i);
              }
          }
          for (auto iter = checkers_.begin(); iter != checkers_.end(); ++ iter) {
              (*iter)->Write(path, "/reference_correction/");
          }
      }
};
}
