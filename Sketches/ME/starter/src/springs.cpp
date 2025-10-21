/*
This code runs a simulation of N coupled oscillators.

Each mass has equations of motion depending on the properties
the springs it is connected to. 

Entity types:
Spring_Ball: /\/\/\/\/\/O
Spring_Wall: /\/\/\/\/\/|
Wall:        |

Example setup: Wall - Spring_Ball - Spring_Ball - Spring_Wall
     - |/\/\/\/\/\/O/\/\/\/\/\/O/\/\/\/\/\/|
*/
#include <iostream>
#include <fstream> 
#include <vector>
#include <flecs.h>
#include <systems.h>


struct Position { 
    double x; // In cm
};

struct Velocity { 
    double x; // In cms-1
};

struct Acceleration{
    double x; // in cms-2
};

struct Mass{
    double m; // in kg
};

struct Spring_Constant{
    double k; // in Ncm-1
};

struct Natural_Length{
    double l; // Spring's natural length in cm
};

struct Spring_Ball{ }; // creates a Spring_Ball tag

struct Spring_Wall{ }; // creates a Spring Wall tag

struct Wall{ }; // creates a Wall tag

int main() {
    std::ofstream MyFile("SHM-Data.txt");
    std::vector<double> testing_pos; 
    std::vector<double> testing_vel;
    std::vector<double> testing_acc;

    flecs::world ecs;

    flecs::system s = ecs.system<Velocity>() 
        .each([](flecs::entity e) {

            std::cout << e.name() << "\n";
            
        });

    // Create Entities
    ecs.entity("Left Wall")
        // Finds and sets components
        .set<Position>({0})

        // Adds Particle Tag
        .add<Wall>();
    
    ecs.entity("mass 1")
        // Finds and sets components
        .set<Mass>({4})
        .set<Position>({5})
        .set<Velocity>({3})
        .set<Acceleration>({0})
        .set<Spring_Constant>({3})
        .set<Natural_Length>({4})
        
        // Adds a Particle Tag
        .add<Spring_Ball>();

    ecs.entity("mass 2")
        // Finds and sets components
        .set<Mass>({4})
        .set<Position>({6})
        .set<Velocity>({0})
        .set<Acceleration>({0})
        .set<Spring_Constant>({3})
        .set<Natural_Length>({4})
        
        // Adds a Particle Tag
        .add<Spring_Ball>();
    
    ecs.entity("Right Wall")
        // Finds and sets components
        .set<Position>({12})
        .set<Spring_Constant>({3})
        .set<Natural_Length>({4})
        
        // Adds a Particle Tag
        .add<Spring_Wall>();


    s.run();

}