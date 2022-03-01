#include "force_directed_pp.hpp"
#include <iostream>

using namespace FDPP;
using namespace FDPP::FORCE;

int main(int argc, char **argv)
{
  SPDLOG_INFO("Fast Map Matching with force directed preprocessing");
  if (argc < 4 || argc > 5)
  {
    SPDLOG_ERROR("Incorrect number of arguments: {}", argc);
  }
  else if (argc == 4)
  {
    // shapefile, tracefile, ubodtfile
    std::unique_ptr<ForceDirectedPP> fdpp(new ForceDirectedPP(argv[1], argv[2], argv[3]));
    fdpp->match();
  }
  else
  {
    // shapefile, tracefile, ubodtfile, iterations
    std::unique_ptr<ForceDirectedPP> fdpp(new ForceDirectedPP(argv[1], argv[2], argv[3], atoi(argv[4])));
    fdpp->match();
  }
  return 0;
};
