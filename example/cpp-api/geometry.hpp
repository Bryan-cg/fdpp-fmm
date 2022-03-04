#include <fmm/fmm-api.hpp>
#include <iostream>

using namespace FMM::CORE;
namespace bg = boost::geometry;

#define EARTH_R 6378137
#define d2r (M_PI / 180.0)             // degrees to radians
#define r2d (180.0 / M_PI)             // radians to degrees
#define a2m ((M_PI / 180.0) * EARTH_R) // angular units to meters
#define Haversine bg::strategy::distance::haversine<double>(EARTH_R)

typedef bg::model::point<double, 2, bg::cs::geographic<bg::degree>> PointDeg;
typedef boost::geometry::model::linestring<PointDeg> LineStringDeg;

namespace GEOM
{
    LineString point_to_lineString(const Point &p);
    std::vector<double> to_cart_cord(const Point &p);
    Point from_cart_cord(const std::vector<double> &coord);
    double interpolated_distance_lineString(const LineString &ls);
    double haversine_distance_m(const Point &p1, const Point &p2);
    LineStringDeg cart_to_degr_linestring(const LineString &ls);
    double calculate_cos_theta(const Point &p1_edge, const Point &p2_edge, const Point &p1_trace, const Point &p2_trace);
}
