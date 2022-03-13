#include <fmm/fmm-api.hpp>
#include "fmm_wrap.hpp"
#include "util/geometry.hpp"
#include "util/wkt_writer.hpp"
#include "util/statistics_writer.hpp"
#include "util/util.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <fmm/io/gps_reader.hpp>

using namespace FMM;
using namespace FMM::NETWORK;
using namespace FMM::CORE;
using namespace FMM::CONFIG;
namespace bg = boost::geometry;

namespace FDPP {
    namespace FORCE {
        // Edges attract or repel point trace according to cosine
        struct ElectricForce {
            ElectricForce(double length_edge_arg, double distance_arg, double c_arg = 2);

            double length_edge;
            double distance; /**< shortest distance between edge and point from trace */
            double c;        /**< Tuning parameter*/
            double calculate(double cos_theta) const;
        };

        // Adjecent vertices attract each other
        struct SpringForceAttr {
            SpringForceAttr(double distance_arg, double c1_arg = 2, double c2_arg = 1);

            double distance; /**< distance between adjecent vertices */
            double c1;       /**< Tuning paramter*/
            double c2;       /**< Tuning paramter*/
            double calculate() const;
        };

        // Non-adjecent vertices repel each other
        struct SpringForceRep {
            SpringForceRep(double distance_arg, double c3_arg = 10);

            double distance; /**< distance between non-adjecent vertices */
            double c3;       /**< Tuning paramter*/
            double calculate() const;
        };
    }
    class ForceDirectedPP {
    private:
        // Iterations FDPP
        int _iterations_fdpp;
        Network *network;
        NetworkGraph *networkGraph;
        // IO config
        ResultConfig *result_config;
        GPSConfig *gps_config;
        // FMM - STM matcher
        std::unique_ptr<fmm_wrap> fmmw;

    public:
        /**
         * @brief Construct a new Force Directed P P object
         * 
         * @param shapeFile 
         * @param tracesFile 
         * @param ubodtFile 
         * @param iterations_arg 
         */
        ForceDirectedPP(char *shapeFile, char *tracesFile, char *ubodtFile, int iterations_arg = 1);

        ~ForceDirectedPP();

        /**
         * @brief 
         * 
         * @param trajectory 
         */
        void force_directed_displacement(Trajectory &trajectory);

        /**
         * @brief 
         * 
         */
        void match();

        /**
         * @brief 
         * 
         * @param mr_pp 
         * @param mr_no_pp 
         * @param traces 
         * @param count_improved 
         * @return std::vector<MM::MatchResult> 
         */
        std::vector<MM::MatchResult> combine_fmm_fdpp_output(
                const std::vector<MM::MatchResult> &mr_pp,
                const std::vector<MM::MatchResult> &mr_no_pp,
                const std::vector<Trajectory> &traces,
                double *count_improved);

        /**
         * @brief 
         * 
         * @param point_i 
         * @param candidates 
         * @param trace_ls 
         * @return Point 
         */
        Point calculate_net_force(
                const int point_i,
                const FMM::MM::Traj_Candidates &candidates,
                const LineString &trace_ls);

        /**
         * @brief 
         * 
         * @param p1 
         * @param p2 
         * @param total_force 
         * @return std::vector<double> 
         */
        std::vector<double>
        calc_electric_force_displacement(const Point &p1, const Point &p2, const double total_force);

        /**
         * @brief 
         * 
         * @param p1 
         * @param p2 
         * @return std::vector<double> 
         */
        std::vector<double> calc_attr_spring_force_displacement(const Point &p1, const Point &p2);

        /**
         * @brief 
         * 
         * @param p1 
         * @param p2 
         * @return std::vector<double> 
         */
        std::vector<double> calc_rep_spring_force_displacement(const Point &p1, const Point &p2);
    };
}
