#include "util.hpp"

constexpr int MIN = 0;
constexpr int MAX = 5;

LineString FDPP::UTIL::add_noise(const LineString &ls) {
    LineString result;
    int N = ls.get_num_points();
    for (int i = 0; i < N; i++) {
        Point p_res = ls.get_point(i);
        std::vector<double> p_res_cart = GEOM::to_cart_cord(p_res);

        // Displacement between 0 & 10m
        std::random_device rd;
        std::default_random_engine eng(rd());
        std::uniform_real_distribution<double> distribution(MIN, MAX);
        const double dx = distribution(eng);
        const double dy = distribution(eng);
        const double dz = distribution(eng);

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
