#include "wkt_writer.hpp"
#include <omp.h>

#include <sstream>

namespace FDPP {

    namespace IO {

        CSVIterationWriter::CSVIterationWriter(const std::string &result_file) : m_fstream(result_file) {
            write_header();
        }

        void CSVIterationWriter::write_header() {
            m_fstream << "id;wkt;start_time;end_time" << '\n';
        }

        void CSVIterationWriter::write_result(const int iteration, FMM::CORE::LineString &ls) {
            std::stringstream buf;
            std::time_t start_time = std::time(nullptr) + iteration;
            const double end_time = start_time + 1;
            buf << iteration;
            buf << ";" << ls;
            buf << ";" << start_time;
            buf << ";" << end_time;
            // buf << ";" << end_time;
            buf << '\n';
            // Ensure that fstream is called corrected in OpenMP
            #pragma omp critical
            m_fstream << buf.rdbuf();
        }
    }
}
