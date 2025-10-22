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


int N = 2;//Particle Number

//Position of the Left Wall
double p_Lwall = 0;  // in cm

// Position of the Right Wall
double p_Rwall = 12; // in cm

struct Index {
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
    std::vector<std::vector<double>> p_matrix {{2}, {7}}; 
    std::vector<std::vector<double>> v_matrix {{2}, {-1}};
    std::vector<std::vector<double>> a_matrix {{0}, {0}};

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
    flecs::system s = ecs.system<Particle, Acceleration, Position, Mass, Index>()
    .each([&](Particle, Acceleration &a, Position &p, 
        const Mass &mass, const Index &index) {
        if (index.i == 0) {
            a.x = - (k_list[index.i]/mass.M) * (p.x - p_Lwall - l_list[index.i])
                + (k_list[index.i+1]/mass.M) * (p_matrix[index.i+1][0] - p.x - l_list[index.i+1]);
        } else if (index.i == N-1) {
            a.x = - (k_list[index.i]/mass.M) * (p.x - p_matrix[index.i-1][0] - l_list[index.i])
                + (k_list[index.i+1]/mass.M) * (p_Rwall - p.x - l_list[index.i+1]);
        } else {
            a.x = - (k_list[index.i]/mass.M) * (p.x - p_matrix[index.i-1][0] - l_list[index.i])
                + (k_list[index.i+1]/mass.M) * (p_matrix[index.i+1][0] - p.x - l_list[index.i+1]);
        }
        
        a_matrix[index.i][0] = a.x;
        
    });

    // Create Entities
    ecs.entity("mass 1")
        // Finds and sets components
        .set<Index>({0})
        .set<Mass>({4})
        .set<Position>({p_matrix[0][0]})
        .set<Velocity>({v_matrix[0][0]})
        .set<Acceleration>({0})
        
        // Adds a Particle Tag
        .add<Particle>();
        
    
    ecs.entity("mass 2")
        // Finds and sets components
        .set<Index>({1})
        .set<Mass>({4})
        .set<Position>({p_matrix[1][0]})
        .set<Velocity>({v_matrix[1][0]})
        .set<Acceleration>({0})
        
        // Adds a Particle Tag
        .add<Particle>();
    
    s.run();

    std::cout << a_matrix[0][0] << "\n";
}
