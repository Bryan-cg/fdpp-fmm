#include "geometry.hpp"

LineString GEOM::point_to_lineString(const Point &p)
{
    LineString temp;
    temp.add_point(p);
    return temp;
}

// x, y, z
std::vector<double> GEOM::to_cart_cord(const Point &p)
{
    std::vector<double> result;
    result.push_back(EARTH_R * cos(p.get<1>() * d2r) * cos(p.get<0>() * d2r));
    result.push_back(EARTH_R * cos(p.get<1>() * d2r) * sin(p.get<0>() * d2r));
    result.push_back(EARTH_R * sin(p.get<1>() * d2r));
    return result;
}

Point GEOM::from_cart_cord(const std::vector<double> &coord)
{
    double lat = r2d * (asin(coord[2] / EARTH_R));
    double lon = r2d * (atan2(coord[1], coord[0]));
    Point p;
    p.set<0>(lon);
    p.set<1>(lat);
    return p;
}

double GEOM::interpolated_distance_lineString(const LineString &ls)
{
    int N = ls.get_num_points();
    Point p1 = ls.get_point(0);
    Point p2 = ls.get_point(N - 1);
    return haversine_distance_m(p1, p2);
}

double GEOM::haversine_distance_m(const Point &p1, const Point &p2)
{
    // Make degree points for distance calculations
    PointDeg _p1(p1.get<0>(), p1.get<1>());
    PointDeg _p2(p2.get<0>(), p2.get<1>());
    return bg::distance(_p1, _p2, Haversine);
}

LineStringDeg GEOM::cart_to_degr_linestring(const LineString &ls)
{
    LineStringDeg lsdg;
    for (int i = 0; i < ls.get_num_points(); i++)
    {
        Point p = ls.get_point(i);
        PointDeg pdg(p.get<0>(), p.get<1>());
        bg::append(lsdg, pdg);
    }
    return lsdg;
}

double GEOM::calculate_cos_theta(const Point &p1_edge, const Point &p2_edge, const Point &p1_trace, const Point &p2_trace)
{
    // Convert points to radians for calculations formula
    double az_edge = bg::formula::spherical_azimuth(p1_edge.get<0>() * d2r, p1_edge.get<1>() * d2r, p2_edge.get<0>() * d2r, p2_edge.get<1>() * d2r);
    double az_trace = bg::formula::spherical_azimuth(p1_trace.get<0>() * d2r, p1_trace.get<1>() * d2r, p2_trace.get<0>() * d2r, p2_trace.get<1>() * d2r);
    // std::cout << "---" << std::endl;
    // std::cout << "GEOMETRYCOLLECTION(" << p1_edge << ", " << p2_edge << ")" << std::endl;
    // std::cout << "GEOMETRYCOLLECTION(" << p1_trace << ", " << p2_trace << ")" << std::endl;
    // std::cout << az_trace * r2d << std::endl;
    // std::cout << az_edge * r2d << std::endl;
    // std::cout << (az_trace - az_edge) * r2d << std::endl;
    return cos(az_trace - az_edge);
};