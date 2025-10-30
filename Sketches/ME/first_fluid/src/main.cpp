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
#include <cmath>
#include <flecs.h>
#include <systems.h>

const int Nx = 150;  // Number of nodes in x-axis
const int Ny = 150;  // Number of nodes in y-axis

/* Creat structs to be used as components*/
struct Position { int x, y; };
struct Conserved { double phi; };

int main(int argc, char* argv[]) {
    /* Create the World */
    flecs::world world(argc, argv);

    /* Create components to be assigned to entities*/
    world.component<Position>();
    world.component<Conserved>();

    /* Create Nodes */
    std::vector<std::vector<flecs::entity>> nodes;  // place to store nodes
    nodes.reserve(Nx * Ny);  // Create the space
    for (int i = 0; i < Nx; ++i) {
        nodes.push_back(std::vector<flecs::entity>{});
        for (int j = 0; j < Ny; ++j) {
            nodes[i].push_back( 
                world.entity()
                    .set<Position>({i, j})
                    .set<Conserved>({2}));
        }
    }
    
}