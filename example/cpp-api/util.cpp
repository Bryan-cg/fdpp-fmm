#include "util.hpp"

LineString FDPP::UTIL::add_noise(const LineString &ls) {
    LineString result;
    int N = ls.get_num_points();
    for (int i = 0; i < N; i++) {
        Point p_res = ls.get_point(i);
        std::vector<double> p_res_cart = GEOM::to_cart_cord(p_res);

        // Displacement between 0 & 10m
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

LineString FDPP::UTIL::lower_sample_freq(const LineString &ls, const int points_between) {
    LineString result;
    int N = ls.get_num_points();
    int removed = 0;
    for (int i = 0; i < N; i++) {
        if (removed == 0) {
            result.add_point(ls.get_point(i));
            removed += 1;
        } else if (removed < points_between)
            removed += 1;
        else
            removed = 0;
    }
    return result;
}
