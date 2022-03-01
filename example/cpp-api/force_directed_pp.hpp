#include <fmm/fmm-api.hpp>
#include "fmm_wrap.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <fmm/io/gps_reader.hpp>

using namespace FMM;
using namespace FMM::NETWORK;
using namespace FMM::CORE;
using namespace FMM::CONFIG;
namespace bg = boost::geometry;

typedef bg::model::point<double, 2, bg::cs::geographic<bg::degree>> PointDeg;
typedef boost::geometry::model::linestring<PointDeg> LineStringDeg;

namespace FDPP
{
    namespace FORCE
    {
        // Edges attract or repel point trace according to cosine
        struct ElectricForce
        {
            ElectricForce(double length_edge_arg, double distance_arg, double c_arg = 1);
            double length_edge;
            double distance; /**< shortest distance between edge and point from trace */
            double c;        /**< Tuning parameter*/
            double calculate(double cos_theta) const;
            double calculate_cos_theta(const Point &p1_edge, const Point &p2_edge, const Point &p1_trace, const Point &p2_trace) const; /**< angle between edge and p1 & p2*/
        };
        // Adjecent vertices attract each other
        struct SpringForceAttr
        {
            SpringForceAttr(double distance_arg, double c1_arg = 2, double c2_arg = 1);
            double distance; /**< distance between adjecent vertices */
            double c1;       /**< Tuning paramter*/
            double c2;       /**< Tuning paramter*/
            double calculate() const;
        };
        // Non-adjecent vertices repel each other
        struct SpringForceRep
        {
            SpringForceRep(double distance_arg, double c3_arg = 10);
            double distance; /**< distance between non-adjecent vertices */
            double c3;       /**< Tuning paramter*/
            double calculate() const;
        };
    }
    class ForceDirectedPP
    {
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
        ForceDirectedPP(char *shapeFile, char *tracesFile, char *ubodtFile, int iterations_arg = 1);
        ~ForceDirectedPP();
        void force_directed_displacement(Trajectory &trajectory);
        void match();
        LineString point_to_lineString(const Point &p);
        double interpolated_distance_lineString(const LineString &ls);
        double haversine_distance_m(const Point &p1, const Point &p2);
        Point calculate_net_force(const int point_i, const FMM::MM::Traj_Candidates &candidates, const LineString &trace_ls);
        std::vector<double> calc_electric_force_displacement(const Point &p1, const Point &p2, const double total_force);
        std::vector<double> calc_attr_spring_force_displacement(const Point &p1, const Point &p2);
        std::vector<double> calc_rep_spring_force_displacement(const Point &p1, const Point &p2);
        std::vector<double> to_cart_cord(const Point &p);
        Point from_cart_cord(const std::vector<double> &coord);
        void calculate_closest_point(const LineString &ls_edge, const Point &p, double *dist, Point *c_p);
        LineStringDeg cart_to_degr_linestring(const LineString &ls);
    };
}
