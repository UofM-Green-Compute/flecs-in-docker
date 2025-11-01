#include <ccenergy/backends/LinuxRAPL.hpp>
#include <filesystem>
#include <fstream>
namespace fs = std::filesystem;
namespace ccenergy {
    bool LinuxRAPLBackend::read_uint64(const std::string & p, uint64_t & o) {
        std::ifstream f(p);
        if (!f)
            return false;
        f >> o;
        return f.good();
    }
    std::vector < LinuxRAPLBackend::Domain > LinuxRAPLBackend::discover_domains() {
        std::vector < Domain > d;
      for (auto & e:fs::directory_iterator("/sys/class/powercap")) {
            auto n = e.path().filename().string();
            if (n.rfind("intel-rapl:", 0) == 0) {
                Domain dom;
                dom.energy_path = (e.path() / "energy_uj").string();
                dom.max_path = (e.path() / "max_energy_range_uj").string();
                d.push_back(dom);
            }
        }
        return d;
    }
    void LinuxRAPLBackend::start() {
        domains_ = discover_domains();
      for (auto & d:domains_)
            read_uint64(d.energy_path, d.start_uj);
    }
    double LinuxRAPLBackend::stop_joules() {
        double tot = 0;
      for (auto & d:domains_) {
            uint64_t e;
            if (read_uint64(d.energy_path, e))
                tot += (e - d.start_uj) / 1e6;
        }
        return tot;
    }
    std::unique_ptr < Backend > make_linux_rapl_backend() {
        return std::make_unique < LinuxRAPLBackend > ();
    }
}
