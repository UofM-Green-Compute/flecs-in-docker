/*
This code runs a simulation of N coupled oscillators.

Each mass has equations of motion depending on the properties
the springs it is connected to. 

First the acceleration of each particle is computed.

Than the program iterates to update the position, velocity, 
and acceleration of each entity for which that is relevant.

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

//Position of the Left Wall
const std::vector<double> x_Lwall = {0, 0};  // in cm

//Initial Position of the first mass
std::vector<double> x1 = {5, 5}; // in cm

// Initial Position of the second mass
std::vector<double> x2 = {7, 7}; // in cm   

// Position of the Right Wall
const std::vector<double> x_Rwall = {12, 12}; // in cm

//Spring Constant of the first spring
double k1 = 3; // in N cm-1

// Spring Constant of the second spring
double k2 = 3; // in N cm-1

// Spring Constant of the third spring
double k3 = 3; // in N cm-1

// Natural Length of the first spring
double l1 = 4; // in N cm-1

// Natural Length of the second spring
double l1 = 4; // in N cm-1

// Natural Length of the third spring
double l1 = 4; // in N cm-1

struct Position { 
    std::vector<double> x_before, x_after; // In cm
};

struct Position_Left { 
    std::vector<double> x_before, x_after; // In cm
};

struct Position_Right { 
    std::vector<double> x_before, x_after; // In cm
};

struct Velocity { 
    double v; // In cms-1
};

struct Acceleration{
    double a; // in cms-2
};

struct Mass{
    const double m; // in kg
};

struct Spring_Constant{
    const double k; // in Ncm-1
};

struct Spring_Constant_Right{
    const double k; // in Ncm-1
                    // Spring constant of the spring to
                    // the right of the mass
};

struct Natural_Length{
    const double l; // Spring's natural length in cm
};

struct Natural_Length_Right{
    const double l; // Spring's natural length in cm
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

    // System s runs for any entity with the tag Spring_Ball
    flecs::system s = ecs.system<Spring_Ball>() 
        .each([](flecs::entity e, const Mass& mass, Position& x1, 
            Velocity& v,Acceleration& a, Spring_Constant& k1,
            Natural_Length& l) {
               
            
            std::cout << e.name() << "\n";
            
        });

    // Create Entities
    ecs.entity("Left Wall")
        // Finds and sets components
        .set<Position>({&x_Lwall})

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

    // Create an ordered list of the entities by position
    std::vector<flecs::entity> ordered;
    ecs.each<Position>([&](flecs::entity e, const position&) {
        
        ordered.push_back(e);

    });

    s.run();

}