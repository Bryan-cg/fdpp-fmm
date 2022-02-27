#include "fmm_wrap.hpp"

fmm_wrap::fmm_wrap(char *shapeFile, char *tracesFile, char *ubodtFile)
{
    std::string shapeFile_str(shapeFile);
    std::string tracesFile_str(tracesFile);
    std::string ubodtFile_str(ubodtFile);

    this->network = new NETWORK::Network(shapeFile_str, "fid", "u", "v");
    this->networkGraph = new NETWORK::NetworkGraph(*(this->network));
    this->ubodt = UBODT::read_ubodt_csv(ubodtFile_str);
    this->fmm_config = new FastMapMatchConfig(8, 300, 50, 0);
    this->stm_config = new STMATCHConfig(8, 300, 50, 30.0, 1.5, 0.0);
    this->stm = new STMATCH(*network, *networkGraph);
    this->fmm = new FastMapMatch(*network, *networkGraph, ubodt, *stm, *stm_config);
    this->gpsConfig = new GPSConfig(tracesFile);
    this->resultConfig = new ResultConfig();
    resultConfig->file = "output.csv";

    SPDLOG_INFO("Network file {}", shapeFile_str);
    SPDLOG_INFO("Network node count {}", network->get_node_count());
    SPDLOG_INFO("Network edge count {}", network->get_edge_count());
    SPDLOG_INFO("networkGraph number of vertices {}", networkGraph->get_num_vertices());
    SPDLOG_INFO("ubodt row count {}", ubodt->get_num_rows());
    gpsConfig->print();
    resultConfig->print();
    SPDLOG_INFO("FMM succesfully initialized");
}

fmm_wrap::~fmm_wrap()
{
    delete network;
    delete networkGraph;
    delete fmm;
    delete fmm_config;
    delete gpsConfig;
    delete resultConfig;
    SPDLOG_INFO("Dynamic allocated memory free");
}

void fmm_wrap::match()
{
    SPDLOG_INFO("{}", fmm->match_gps_file(*gpsConfig, *resultConfig, *fmm_config, true));
}

void fmm_wrap::matchAllSTM()
{
    SPDLOG_INFO("{}", stm->match_gps_file(*gpsConfig, *resultConfig, *stm_config, true));
}
