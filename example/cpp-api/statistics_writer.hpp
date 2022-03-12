//
// Created by coulier on 12/03/2022.
//
#include <iostream>
#include <fstream>
#include <fmm/fmm-api.hpp>

using namespace FMM::CORE;

namespace FDPP {

    namespace IO {

        class StatsWriter {
        public:
            /**
             * Constructor
             *
             * @param result_file the filename to write result
             *
             */
            explicit StatsWriter(const std::string &result_file);

            /**
             * Write a header line for the fields exported
             */
            void write_header();

            /**
             * Write wkt linestring
             * @param ls Input trajectory wkt
             */
            void write_result(const int trace_id, const double length_trace, const double li, const double avg_e, const double hausdorff, const double frechet);

        private:
            std::ofstream m_fstream;
        };
    }
}

