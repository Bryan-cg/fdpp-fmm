
#include <iostream>
#include <fstream>
#include <omp.h>
#include <ctime>
#include <fmm/fmm-api.hpp>

using namespace FMM::CORE;

namespace FDPP {

    namespace IO {

        /**
         * An interface defined for writing the map match result
         */
        class IterationWriter {
        public:
            /**
             * Write the match result to a file
             * @param result the match result of a trajectory
             */
            virtual void write_result(const int iteration, FMM::CORE::LineString &ls) = 0;
        };

        /**
         * A writer class for writing matche result to a CSV file.
         */
        class CSVIterationWriter : public IterationWriter {
        public:
            /**
             * Constructor
             *
             * @param result_file the filename to write result
             *
             */
            explicit CSVIterationWriter(const std::string &result_file);

            /**
             * Write a header line for the fields exported
             */
            void write_header();

            /**
             * Write wkt linestring
             * @param ls Input trajectory wkt
             */
            void write_result(const int iteration, FMM::CORE::LineString &ls) override;
            /**
             * Write a header line for the fields exported
             */
        private:
            std::ofstream m_fstream;
        };

    };
}