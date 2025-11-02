#pragma once
#include <string>
#include <chrono>
#include <optional>
#include <vector>
#include <memory>
#include <cstdint>
#include <functional>
#include <cstdio>

#include <fstream>
#include <filesystem>

namespace fs = std::filesystem;

namespace ccenergy {
    struct Result {
        std::string label;
        double seconds {0.0};
        double cpu_joules {0.0};
        double gpu_joules {0.0};
        double total_joules() const {
            return cpu_joules + gpu_joules;
        }
        double avg_power_watts() const {
            return seconds > 0 ? total_joules() / seconds : 0.0;
        }
    };

    struct Config {
        std::string label {"session"};
        bool measure_cpu {true};
        bool measure_gpu {false};
        bool log_to_stdout {true};
    };

    class Backend {
      public:
        virtual ~ Backend() = default;
        virtual void start() = 0;
        virtual double stop_joules() = 0;
    };

    std::unique_ptr < Backend > make_linux_rapl_backend();
    std::unique_ptr < Backend > make_nvml_backend();

    class EnergyTracker {
      public:
        explicit EnergyTracker(Config cfg = { }):cfg_(std::move(cfg)) { }
        void start();
        Result stop();
        static Result measure(const std::string & label, const std::function < void () > &fn, Config cfg = { });
      private:
        using Clock = std::chrono::steady_clock;
        Config cfg_;
        std::unique_ptr < Backend > cpu_;
        std::unique_ptr < Backend > gpu_;
        Clock::time_point start_tp_;
        static void log_result(const Result & r);
    };
    void EnergyTracker::start() {
        start_tp_ = Clock::now();
        cpu_ = make_linux_rapl_backend();
        if (cpu_)
            cpu_->start();
    }
    Result EnergyTracker::stop() {
        auto end = Clock::now();
        Result r;
        r.label = cfg_.label;
        r.seconds = std::chrono::duration < double >(end - start_tp_).count();
        if (cpu_)
            r.cpu_joules = cpu_->stop_joules();
        if (cfg_.log_to_stdout)
            log_result(r);
        return r;
    }

    Result EnergyTracker::measure(const std::string & label, const std::function < void () > &fn, Config cfg) {
        cfg.label = label;
        EnergyTracker t(cfg);
        t.start();
        fn();
        return t.stop();
    }
    void EnergyTracker::log_result(const Result & r) {
        printf("[ccenergy] %-10s time %.3fs CPU %.3fJ total %.3fJ avg %.3fW\n",
               r.label.c_str(), r.seconds, r.cpu_joules, r.total_joules(), r.avg_power_watts());
    }

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

}                               // namespace ccenergy
