/*
This code runs a simulation of N coupled oscillators.

Each mass has equations of motion depending on the properties
the springs it is connected to. 

First the acceleration of each particle is computed.

Than the program iterates to update the position, velocity, 
and acceleration of each entity for which that is relevant.

Entity types:
 - Particle
*/

#include <iostream>
#include <fstream> 
#include <vector>
#include <flecs.h>
#include <systems.h>


//Position of the Left Wall
double x_Lwall = 0;  // in cm

// Position of the Right Wall
double x_Rwall = 12; // in cm

struct Order {
    int i; // Which particle is it
};

struct Position { 
    double x; // In cm
};

struct Velocity { 
    double x;// In cms-1
};

struct Acceleration{
    double x; // in cms-2
};

struct Mass{
    double M; // in kg
};

struct Particle{ }; // creates a  tag

int main() {
    std::ofstream MyFile("SHM-Data.txt");

    // First Particle Data
    std::vector<double> position_1 {2}; 
    std::vector<double> velocity_1 {2};
    std::vector<double> acceleration_1 {0};

    // Second Particle Data
    std::vector<double> position_2 {7}; 
    std::vector<double> velocity_2 {-1};
    std::vector<double> acceleration_2 {0};

    //Spring Constant Data
    const std::vector<double> k_list {3, 3, 3};
    // Natural Spring Length Data
    const std::vector<double> l_list {4, 4, 4};

    flecs::world ecs;
    
    /*
        System s1 runs for any entity with the tag Particle
        Its job is to calculate the initial acceleration and
        update the vectors and entity components. 
    */
   
    flecs::system s = ecs.system<Particle>()
    .each([](Particle) {
        std::cout << "3\n";
    });

    // Create Entities
    ecs.entity("mass 1")
        // Finds and sets components
        .set<Order>({1})
        .set<Mass>({4})
        .set<Position>({position_1[0]})
        .set<Velocity>({velocity_1[0]})
        .set<Acceleration>({acceleration_1[0]})
        
        // Adds a Particle Tag
        .add<Particle>();
        
    
    ecs.entity("mass 2")
        // Finds and sets components
        .set<Order>({1})
        .set<Mass>({4})
        .set<Position>({position_2[0]})
        .set<Velocity>({velocity_2[0]})
        .set<Acceleration>({acceleration_2[0]})
        
        // Adds a Particle Tag
        .add<Particle>();
    
    s.run();
}
