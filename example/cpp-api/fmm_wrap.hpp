#include <fmm/fmm-api.hpp>
#include <fmm/mm/fmm/ubodt.hpp>
#include <fmm/mm/fmm/fmm_algorithm.hpp>
#include <fmm/python/pyfmm.hpp>
#include <fmm/mm/stmatch/stmatch_algorithm.hpp>
#include <fmm/io/gps_reader.hpp>

#include <iostream>
using namespace FMM;
using namespace FMM::NETWORK;
using namespace FMM::CORE;
using namespace FMM::MM;
using namespace FMM::CONFIG;

class fmm_wrap
{
private:
    //network
    Network *network;
    NetworkGraph *networkGraph;
    //ubodt precomutation for FMM
    std::shared_ptr<UBODT> ubodt;
    //FMM
    FastMapMatch *fmm;
    FastMapMatchConfig *fmm_config;
    //STM backup
    STMATCHConfig *stm_config;
    STMATCH *stm;
    //IO config
    GPSConfig *gpsConfig;
    ResultConfig *resultConfig;

public:
    fmm_wrap(char *shapeFile, char *tracesFile, char *ubodtFile);
    ~fmm_wrap();
    void matchAllSTM();
    void match();
};
