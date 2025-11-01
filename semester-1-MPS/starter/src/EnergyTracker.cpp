#include <ccenergy/EnergyTracker.hpp>
#include <ccenergy/backends/LinuxRAPL.hpp>
namespace ccenergy {
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
}
