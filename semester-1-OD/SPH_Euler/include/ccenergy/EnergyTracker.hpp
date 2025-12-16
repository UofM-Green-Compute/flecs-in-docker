//
// (c) 2025 Michael Sparks, University of Manchester
// You may use this under the terms of the Apache 2 License
//
//
// This file implements minimal and basic energy tracking.  Its API is
// loosely inspired by (but doesn't take code from) CodeCarbon (Hence the
// name "ccenergy") There are placeholders for some necessary, but
// un-implemented extras.
// 
// For now it accesses RAPL via /sys/class/powercap on modern linux systems.
//
// It's loosely inspired in that:
//
// * You create a tracker
// * As many times as you like:
//   * You start it
//   * You stop it
// * You report the results.
// Creation:
//
//     ccenergy::EnergyTracker energy_tracker {{ .label = "OnUpdate",
//                                               .measure_cpu = true,
//                                               .measure_gpu  = false,
//                                               .log_to_stdout = false }};
//
// Starting:
//
//             energy_tracker.start();
//
// Stopping:
//
//             auto r = energy_tracker.stop();
//
// Reporting:
//
//     std::cout << energy_tracker.mkReport() << std::endl;
//
// It doesn't implement the NVML style carbon accounting as yet.
//
// It relies on RAPL and therefore has the same limitations as RAPL in general.
// If you are using this for benchmarking, you need your software to be the
// dominant piece of code running on that machine when performing measurements.
//
// Things that are really needed:
//
// * Minimally sufficient docs - including limitations.
// * Minimally sufficient reproducible testing (manual acceptance testing is all well and
//   good, but unit tests are better).
// * Ideally a mechanism to seperate out quiesscent power usage from high level (if possible)
// * GPU power tracking. (via NVML or similar)
//

#pragma once

#include <vector>
#include <functional>
#include <fstream>
#include <filesystem>

#include <futureprint.hpp>

namespace fs = std::filesystem;

namespace ccenergy {
    // Structure for capturing and updating stats (costs) relating to a tracker.
    struct EnergyAccum {
        double seconds{0.0};
        double cpu_j{0.0};
        double gpu_j{0.0};
    };

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
    std::unique_ptr < Backend > make_nvml_backend();  // TBD

    class EnergyTracker {
      public:
        explicit EnergyTracker(Config init_config = { }) : config(std::move(init_config)) { }
        void start();
        Result stop();
        std::string mkReport();
        Result measure(const std::string & label, const std::function < void () > &fn, Config config = { });
      private:
        Config config;
        EnergyAccum energy_counters {};
        using Clock = std::chrono::steady_clock;
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
        r.label = config.label;
        r.seconds = std::chrono::duration < double >(end - start_tp_).count();
        if (cpu_)
            r.cpu_joules = cpu_->stop_joules();
        if (config.log_to_stdout)
            log_result(r);

        energy_counters.seconds += r.seconds;
        energy_counters.cpu_j   += r.cpu_joules;
        energy_counters.gpu_j   += r.gpu_joules;

        return r;
    }

    Result EnergyTracker::measure(const std::string & label, const std::function < void () > &fn, Config config) {
        config.label = label;
        EnergyTracker t(config);
        t.start();
        fn();
        return t.stop();
    }
    void EnergyTracker::log_result(const Result & r) {
        printf("[ccenergy] %-10s time %.3fs CPU %.3fJ total %.3fJ avg %.3fW\n",
               r.label.c_str(), r.seconds, r.cpu_joules, r.total_joules(), r.avg_power_watts());
    }

    std::string EnergyTracker::mkReport() {
        const double total_joules   = energy_counters.cpu_j + energy_counters.gpu_j;
        const double avg_watts = (energy_counters.seconds > 0) ? total_joules / energy_counters.seconds : 0.0;


        return fmt("[ccenergy-summary] label={} frames_seconds={:.3f} cpu_joules={:.3f} gpu_joules={:.3f} total_joules={:.3f} avg_watts={:.3f}",
                                config.label, energy_counters.seconds, energy_counters.cpu_j, energy_counters.gpu_j, total_joules, avg_watts);
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
            if (n == "intel-rapl:0") { // Rather than add all of the domains, just add the one we expect to find in our setup
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
