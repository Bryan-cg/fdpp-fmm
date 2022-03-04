#include "util.hpp"

LineString FDPP::UTIL::add_noise(const LineString &ls)
{
    LineString result;
    int points = ls.get_num_points();
    for (int i = 0; i < points; i++)
    {
        Point p_res = ls.get_point(i);
        std::vector<double> p_res_cart = GEOM::to_cart_cord(p_res);

        // tussen de 0.1 en 2m verplaatsen
        const int dx = rand() % 10;
        const int dy = rand() % 10;
        const int dz = rand() % 10;

        p_res_cart[0] += dx;
        p_res_cart[1] += dy;
        p_res_cart[2] += dz;

        p_res = GEOM::from_cart_cord(p_res_cart);
        result.add_point(p_res);
    }
    return result;
}
