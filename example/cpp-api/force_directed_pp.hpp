#include <fmm/fmm-api.hpp>

#include <iostream>
#include <cmath>
#include <algorithm> // std::min

using namespace FMM;
using namespace FMM::NETWORK;
using namespace FMM::CORE;
namespace bg = boost::geometry;

typedef bg::model::point<double, 2, bg::cs::geographic<bg::degree>> PointDeg;
typedef boost::geometry::model::linestring<PointDeg> LineStringDeg;

namespace FDPP
{
    namespace EF
    {
        struct ElectricForce
        {
            ElectricForce(double length_edge_arg, double distance_arg, double c_arg = 1);
            double length_edge;
            double distance; /**< distance between edge and point from trace perpendicular */
            double c;        /** Tuning parameter*/
            double calculate(double cos_theta) const;
            double calculate_cos_theta(const Point &p1_edge, const Point &p2_edge, const Point &p1_trace, const Point &p2_trace) const; /**< angle between edge and p1 & p2*/
        };
        struct SpringForce
        {
            SpringForce(double distance_arg, double c1_arg = 2, double c2_arg = 1);
            double distance;
            double c1; /**< Tuning paramter*/
            double c2; /**< Tuning paramter*/
            double calculate() const;
        };

    }
    class ForceDirectedPP
    {
    private:
        Network *network;
        NetworkGraph *networkGraph;

    public:
        ForceDirectedPP(char *shapeFile, char *tracesFile);
        ~ForceDirectedPP();
        void displace_linestring(const std::string &wkt);

        LineString point_to_lineString(const Point &p);
        double interpolated_distance_lineString(const LineString &ls);
        double haversine_distance_m(const Point &p1, const Point &p2);
        Point calculate_net_force(const Point &p, const Point &p_prev, const Point &p_next, const FMM::MM::Traj_Candidates &candidates);
        std::vector<double> calc_electric_force_displacement(const Point &p1, const Point &p2, const double total_force);
        std::vector<double> calc_spring_force_displacement(const Point &p1, const Point &p2);
        std::vector<double> to_cart_cord(const Point &p);
        Point from_cart_cord(const std::vector<double> &coord);
        void calculate_closest_point(const LineString &ls_edge, const Point &p, double *dist, Point *c_p);
        LineStringDeg cart_to_degr_linestring(const LineString &ls);
    };
}
