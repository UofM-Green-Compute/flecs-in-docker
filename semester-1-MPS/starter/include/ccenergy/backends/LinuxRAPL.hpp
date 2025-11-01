#pragma once
#include "../EnergyTracker.hpp"
#include <string>
#include <vector>
#include <cstdint>

namespace ccenergy {
    class LinuxRAPLBackend:public Backend {
      public:
        void start() override;
        double stop_joules() override;
      private:
        struct Domain {
            std::string energy_path;
            std::string max_path;
            uint64_t start_uj {
            0};
            uint64_t max_uj {
            0};
        };
          std::vector < Domain > domains_;
        static std::vector < Domain > discover_domains();
        static bool read_uint64(const std::string &, uint64_t &);
    };
      std::unique_ptr < Backend > make_linux_rapl_backend();
}
