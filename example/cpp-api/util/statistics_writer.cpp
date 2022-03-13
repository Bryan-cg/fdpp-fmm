//
// Created by coulier on 12/03/2022.
//

#include "statistics_writer.hpp"


FDPP::IO::StatsWriter::StatsWriter(const std::string &result_file) : m_fstream(result_file) {
    write_header();
}


void FDPP::IO::StatsWriter::write_header() {
    m_fstream << "id;length_trace;length_index;avg_error;hausdorff;frechet;" << '\n';
}

void
FDPP::IO::StatsWriter::write_result(const int trace_id, const double length_trace, const double li, const double avg_e,
                                    const double hausdorff, const double frechet) {
    std::stringstream buf;
    buf << trace_id;
    buf << ";" << length_trace;
    buf << ";" << li;
    buf << ";" << avg_e;
    buf << ";" << hausdorff;
    buf << ";" << frechet;
    buf << '\n';
    // Ensure that fstream is called corrected in OpenMP
#pragma omp critical
    m_fstream << buf.rdbuf();

}
