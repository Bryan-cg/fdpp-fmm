/**
 * Fast map matching.
 *
 * fmm algorithm implementation
 *
 * @author: Can Yang
 * @version: 2020.01.31
 */

#ifndef FMM_FMM_ALGORITHM_H_
#define FMM_FMM_ALGORITHM_H_

#include "network/network.hpp"
#include "network/network_graph.hpp"
#include "mm/transition_graph.hpp"
#include "mm/fmm/ubodt.hpp"
#include "python/pyfmm.hpp"
#include "config/gps_config.hpp"
#include "config/result_config.hpp"
#include "mm/stmatch/stmatch_algorithm.hpp"

#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "cxxopts/cxxopts.hpp"

namespace FMM {
namespace MM{
/**
 * Configuration class for fmm algorithm
 */
struct FastMapMatchConfig{
  /**
   * Constructor of FastMapMatch configuration
   * @param k_arg the number of candidates
   * @param r_arg the search radius, in map unit, which is the same as
   * GPS data and network data.
   * @param gps_error the gps error, in map unit
   *
   */
  FastMapMatchConfig(int k_arg = 8, double r_arg = 300, double gps_error = 50,
    double reverse_tolerance = 0.0);
  int k; /**< Number of candidates */
  double radius; /**< Search radius*/
  double gps_error; /**< GPS error */
  double reverse_tolerance;
  /**
   * Check if the configuration is valid or not
   * @return true if valid
   */
  bool validate () const;
  /**
   * Print information about this configuration
   */
  void print() const;
  /**
   * Load configuration from xml data
   * @param xml_data xml data read from an xml file
   * @return a %FastMapMatchConfig object
   */
  static FastMapMatchConfig load_from_xml(
      const boost::property_tree::ptree &xml_data);
  /**
   * Load configuration from argument data
   * @param arg_data argument data
   * @return a %FastMapMatchConfig object
   */
  static FastMapMatchConfig load_from_arg(
      const cxxopts::ParseResult &arg_data);
  /**
   * Register arguments to an option object
   */
  static void register_arg(cxxopts::Options &options);
  /**
   * Register help information to a string stream
   */
  static void register_help(std::ostringstream &oss);
};

/**
 * Fast map matching algorithm/model.
 *
 *
 */
class FastMapMatch {
 public:
  /**
   * Constructor of Fast map matching model
   * @param network road network
   * @param graph road network graph
   * @param ubodt Upperbounded origin destination table
   * @param stm spatial temporal backup matcher
   * @param stmConfig spatial temporal matcher configuration
   */
  FastMapMatch(const NETWORK::Network &network,
      const  NETWORK::NetworkGraph &graph,
      std::shared_ptr<UBODT> ubodt,
      STMATCH &stmMatch,
      STMATCHConfig &stmConfig)
      : network_(network), graph_(graph), ubodt_(ubodt), stmMatch_(stmMatch), stmConfig_(stmConfig) {
  };
  /**
   * Match a trajectory to the road network
   * @param  traj   input trajector data
   * @param  config configuration of map matching algorithm
   * @return map matching result
   */
  MatchResult match_traj(const CORE::Trajectory &traj,
                         const FastMapMatchConfig &config);
  /**
   * Match a wkt linestring to the road network.
   * @param wkt WKT representation of a trajectory
   * @param config Map matching configuration
   * @return Map matching result in POD format used in Python API
   */
  PYTHON::PyMatchResult match_wkt(
      const std::string &wkt,const FastMapMatchConfig &config);
  /**
   * Match GPS data stored in a file
   * @param  gps_config    [description]
   * @param  result_config [description]
   * @param  config        [description]
   * @return a string storing information about running time and
   * statistics.
   */
  std::string match_gps_file(
    const FMM::CONFIG::GPSConfig &gps_config,
    const FMM::CONFIG::ResultConfig &result_config,
    const FastMapMatchConfig &config,
    bool use_omp = true
  );

  /**
   * Match Trajectories stored in std::vector
   * @param  trajectories    [description]
   * @param  result_config [description]
   * @param  config        [description]
   * @return a string storing information about running time and
   * statistics.
   */
  std::string match_trajectories(
    std::vector<FMM::CORE::Trajectory> &trajectories,
    const FMM::CONFIG::ResultConfig &result_config,
    const FastMapMatchConfig &config,
    bool use_omp = true
  );
 protected:
  /**
   * Get shortest path distance between two candidates
   * @param  ca from candidate
   * @param  cb to candidate
   * @return  shortest path value
   */
  double get_sp_dist(const Candidate *ca,
                     const Candidate *cb,
                     double reverse_tolerance);
  /**
   * Update probabilities in a transition graph
   * @param tg transition graph
   * @param traj raw trajectory
   * @param config map match configuration
   */
  void update_tg(TransitionGraph *tg,
                 const CORE::Trajectory &traj,
                 double reverse_tolerance = 0);
  /**
   * Update probabilities between two layers a and b in the transition graph
   * @param level   the index of layer a
   * @param la_ptr  layer a
   * @param lb_ptr  layer b next to a
   * @param eu_dist Euclidean distance between two observed point
   * @param connected the variable is set to false if the layer is not connected
   * with the next layer
   */
  void update_layer(int level, TGLayer *la_ptr, TGLayer *lb_ptr,
                    double eu_dist, double reverse_tolerance,
                    bool *connected);
 private:
  const NETWORK::Network &network_;
  const NETWORK::NetworkGraph &graph_;
  std::shared_ptr<UBODT> ubodt_;
  STMATCHConfig &stmConfig_;
  STMATCH &stmMatch_;
};
}
}

#endif //FMM_FMM_ALGORITHM_H_
