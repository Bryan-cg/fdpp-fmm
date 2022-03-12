#include "fmm_wrap.hpp"

fmm_wrap::fmm_wrap(char *ubodtFile, const Network &network, const NetworkGraph &networkGraph) {
    std::string ubodtFile_str(ubodtFile);

    SPDLOG_INFO("Reading UBODT file");
    this->ubodt = UBODT::read_ubodt_csv(ubodtFile_str);
    this->fmm_config = new FastMapMatchConfig(8, 300, 50, 0);
    this->stm_config = new STMATCHConfig(8, 300, 50, 30.0, 1.5, 0.0);
    this->stm = new STMATCH(network, networkGraph);
    this->fmm = new FastMapMatch(network, networkGraph, ubodt, *stm, *stm_config);
    SPDLOG_INFO("ubodt row count {}", ubodt->get_num_rows());
    SPDLOG_INFO("FMM succesfully initialized");
}

fmm_wrap::~fmm_wrap() {
    delete fmm;
    delete fmm_config;
}

std::vector<MM::MatchResult> fmm_wrap::match(std::vector<FMM::CORE::Trajectory> &trajectories) {
    if (!fmm_config->validate()) {
        SPDLOG_CRITICAL("fmm_config invalid - abort matching");
        return {};
    }
    return fmm->match_trajectories(trajectories, *fmm_config, true);
}

void fmm_wrap::matchAllSTM() {
    // SPDLOG_INFO("{}", stm->match_gps_file(*gpsConfig, *resultConfig, *stm_config, true));
}
