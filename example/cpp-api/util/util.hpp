#include "geometry.hpp"
#include <random>

namespace FDPP {
    namespace UTIL {
        /**
         * @brief add random noise to linestring
         * 
         * @param ls 
         * @return LineString with noise
         */
        LineString add_noise(const LineString &ls);

        /**
         * @brief delete certain number of points between to mimick lower sampling rate
         * 
         * @param ls 
         * @param points_between 
         * @return LineString 
         */
        LineString lower_sample_freq(const LineString &ls, const int points_between);
    }
}
