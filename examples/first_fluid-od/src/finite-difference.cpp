/*
This code aims to use the finite difference method to solve simple equations in fluid dynamics

*/

/*

#include <iostream>
#include <fstream> 
#include <vector>
#include <cmath>
#include <flecs.h>
#include <systems.h>

const int Nx = 150;  // Number of nodes in x-axis
const int Ny = 150;  // Number of nodes in y-axis

// Creat structs to be used as components
struct Position { int x, y; };
struct Conserved { double phi; };

int main(int argc, char* argv[]) {
    // Create the World 
    flecs::world world(argc, argv);

    // Define grid 
    // Create components to be assigned to entities
    world.component<Position>();
    world.component<Conserved>();
    // Create Nodes 
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

// Solve matrix equation 



*/
