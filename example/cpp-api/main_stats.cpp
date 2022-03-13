//
// Created by coulier on 13/03/2022.
//
#include "util/statistics_writer.hpp"
#include <fmm/io/gps_reader.hpp>
#include <fmm/fmm-api.hpp>
#include "util/geometry.hpp"

int main(int argc, char **argv) {
    SPDLOG_INFO("Calculating accuracy");
    if (argc < 3) {
        SPDLOG_ERROR("Incorrect number of arguments: {}", argc);
    } else {
        SPDLOG_INFO("Reading traces");
        char *unmatched_traces_file = argv[1];
        char *matched_traces_file = argv[2];
        std::unique_ptr<FMM::CONFIG::GPSConfig> unmatched_gps_config(
                new FMM::CONFIG::GPSConfig(unmatched_traces_file, "id", "geom"));
        std::unique_ptr<FMM::CONFIG::GPSConfig> matched_gps_config(
                new FMM::CONFIG::GPSConfig(matched_traces_file, "id", "mgeom"));
        FMM::IO::GPSReader reader_unmatched(*unmatched_gps_config);
        FMM::IO::GPSReader reader_matched(*matched_gps_config);
        std::vector<Trajectory> traces = reader_unmatched.read_all_trajectories();
        std::vector<Trajectory> matches = reader_matched.read_all_trajectories();
        SPDLOG_INFO("Calculating accuracy & flushing calculations to output file");
        FDPP::IO::StatsWriter stats_writer("stats.csv");
        double li_total = 0.0;
        double avgE_total = 0.0;
        double frechet_total = 0.0;
        double hausdorff_total = 0.0;
        double count_unmatched = 0.0;
        int trajectories_fetched = traces.size();
        for (int i = 0; i < trajectories_fetched; i++) {
            Trajectory &trace = traces[i];
            Trajectory &match = matches[i];
            double li, avgE, frechet, hausdorff;
            GEOM::calc_accuracy(
                    GEOM::cart_to_degr_linestring(trace.geom),
                    GEOM::cart_to_degr_linestring(match.geom),
                    &li,
                    &avgE,
                    &frechet,
                    &hausdorff);
            if (li == 0.0) {
                count_unmatched += 1.0;
                continue;
            }
            stats_writer.write_result(i + 1,
                                      bg::length(GEOM::cart_to_degr_linestring(traces[i].geom), Haversine),
                                      li,
                                      avgE,
                                      hausdorff,
                                      frechet);
            li_total += li;
            avgE_total += avgE;
            frechet_total += frechet;
            hausdorff_total += hausdorff;
        }
        double total_matched = trajectories_fetched - count_unmatched;
        li_total /= total_matched;
        avgE_total /= total_matched;
        frechet_total /= total_matched;
        hausdorff_total /= total_matched;
        SPDLOG_INFO("Accuracy: li {}, avg error {}, frechet {}, hausdorff {}", li_total, avgE_total, frechet_total,
                    hausdorff_total);
        SPDLOG_INFO("unmatched traces {}", count_unmatched);
    }
    return 0;
};