#include <ccenergy/EnergyTracker.hpp>
#include <iostream>
#include <string>

#include <chrono>

// Busy wait for the given number of milliseconds
void busy_pause(int pause_time_in_millis ){
    // Get start time
    std::chrono::time_point start = std::chrono::time_point_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now()
    );

    while(true) {
        // Get current time
        std::chrono::time_point currently = std::chrono::time_point_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now()
        );

        auto delta = currently - start;
        if (delta > std::chrono::milliseconds(pause_time_in_millis)) {
            break;
        }

    }

}

int main(int /* argc */, char ** /*argv */) {
    // Create the energy tracker
    ccenergy::EnergyTracker energy_tracker {{ .label = "OnUpdate",
                                              .measure_cpu = true,
                                              .measure_gpu  = false,
                                              .log_to_stdout = false }};

    energy_tracker.start();    // Start the energy tracker before the thing you want to measure

    busy_pause(500);           // Do the thing you want to measure

    auto r = energy_tracker.stop();      // stop the energy tracker


    std::cout << "Initial Iteration complete\n";


    for(auto iter=5; iter--;) {
        busy_pause(300);        // Do some other stuff


        energy_tracker.start();    // Start the energy tracker before the thing you want to measure

        busy_pause(500);           // Do the thing you want to measure

        auto r = energy_tracker.stop();      // stop the energy tracker

        std::cout << "Iteration complete\n";

    }

    std::cout << energy_tracker.mkReport() << std::endl;

    return 0;
}