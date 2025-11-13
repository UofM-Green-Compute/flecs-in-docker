/*
This code runs a simulation of a simple fluid which can be 
solved analytically.
*/

#include <custom_pipeline.h>
#include <iostream>
#include <vector>
#include <fstream> 
#include <flecs.h>
#include <systems.h>
#include <cmath>

// Sets time parameters
int STEPS = 100;
int TIME = 60;

// Number of nodes in x-axis. Also equal to the length in unit of node distance
const double L = 1; // Length
const double N = 41; // Number of nodes
const int PE = 50; // Peclet Number
const double RHO = 2; // Density
const double U = 5; // Speed

// Calculates Diffusivity
double GAMMA = (RHO * U * L) / PE;

// Boundary Conditions
double Start = 0;
double End = 1;

// Pipeline structs
struct Save {};
struct RungeKutta_1 {};
struct RungeKutta_2 {};
struct RungeKutta_3 {};
struct RungeKutta_4 {};
struct Update {};

// Position Component
struct position {double x; };

// Field components
struct phi_start { double i; };
struct phihalf_first { double i; };
struct phihalf_second { double i; };
struct phinext_first { double i; };
struct phinext_final { double i; };

// Component Tags
struct Middle_Tag {};
struct End_Tag {};

double linear_function(double x, double A, double B)
{
    double m = (B-A)/L;
    double y = m * x;
    return y;
}

double fluid_function(double phi_left, double phi, double phi_right, 
    double x_left, double x, double x_right)
{
    double f;
    f = - U*((phi_right + phi_left)/(x_right-x_left)) 
    + (GAMMA/RHO) * (phi_right - 2 * phi + phi_left) / ((x_right - x) * (x - x_left));
    return f;
}

int main(int argc, char* argv[]) {
    // Prepare Save File
    std::ofstream MyFile;
    MyFile.open("starter_fluid.txt");
    if (!MyFile.is_open())
    {
        std::cout<<"Error in creating file"<<std::endl; 
        return 1; 
    }

    // Create the World
    flecs::world world(argc, argv);

    // Set up custom pipeline
    flecs::entity pipeline = world.pipeline()
        .with(flecs::System)
        .with<Save>()
        .with<RungeKutta_1>()
        .with<RungeKutta_2>()
        .with<RungeKutta_3>()
        .with<RungeKutta_4>()
        .with<Update>()
        .build();
    world.set_pipeline(pipeline);

    // Create components inside the world
    world.component<Middle_Tag>();
    world.component<End_Tag>();

    world.component<position>();

    world.component<phi_start>();
    world.component<phihalf_first>();
    world.component<phihalf_second>();
    world.component<phinext_first>();
    world.component<phinext_final>();

    std::vector<flecs::entity> nodes; // place to store nodes
    nodes.reserve(N); // Create the space
    
    nodes.push_back(
            world.entity()
                .add<End_Tag>() 
                .set<position>({0})
                .set<phi_start>({0})
                .set<phihalf_first>({0})
                .set<phihalf_second>({0})
                .set<phinext_first>({0})
            );

    for (int index = 1; index <= N-1; ++index) {
        double phi_initial;
        phi_initial = linear_function(index*L/N, Start, End);
        nodes.push_back(
            world.entity()
                .add<Middle_Tag>() 
                .set<position>({index*L/N})
                .set<phi_start>({phi_initial})
                .set<phihalf_first>({0})
                .set<phihalf_second>({0})
                .set<phinext_first>({0})
                .set<phinext_final>({0})
            );
    }

    nodes.push_back(
            world.entity()
                .add<End_Tag>() 
                .set<position>({L})
                .set<phi_start>({0})
                .set<phihalf_first>({0})
                .set<phihalf_second>({0})
                .set<phinext_first>({0})
            );
    
    world.progress();
    MyFile.close();
}
