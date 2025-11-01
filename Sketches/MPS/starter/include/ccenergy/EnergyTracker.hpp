#pragma once
#include <string>
#include <chrono>
#include <optional>
#include <vector>
#include <memory>
#include <cstdint>
#include <functional>
#include <cstdio>

namespace ccenergy {

    struct Result {
        std::string label;
        double seconds {
        0.0};
        double cpu_joules {
        0.0};
        double gpu_joules {
        0.0};
        double total_joules() const {
            return cpu_joules + gpu_joules;
        }
        double avg_power_watts() const {
            return seconds > 0 ? total_joules() / seconds : 0.0;
        }
    };

    struct Config {
        std::string label {
        "session"};
        bool measure_cpu {
        true};
        bool measure_gpu {
        false};
        bool log_to_stdout {
        true};
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

}                               // namespace ccenergy
