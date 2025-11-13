/*
This code looks at the dynamical evolution of a 1D convection-conduction
system which begins out of equilibrium
*/

#include <custom_phases_no_builtin.h>
#include <iostream>
#include <fstream> 
#include <vector>

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

// Position Component
struct position {double x; };

// Field components
struct phi_start { double i; };
struct phihalf_predict { double i; };
struct phihalf_correct { double i; };
struct phi_end_predict { double i; };
struct phi_end_correct { double i; };

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

int main(int argc, char *argv[]) {
    // Create World
    flecs::world world(argc, argv);

    // Creates Phases which tell the program in which order to run the systems
    
    flecs::entity Save = world.entity()
        .add(flecs::Phase); // This Phase saves the state in its current form

    flecs::entity RungeKutta_1 = world.entity()
        .add(flecs::Phase) // This Phase calculates phihalf_predict
        .depends_on(Save);

    flecs::entity RungeKutta_2 = world.entity()
        .add(flecs::Phase) // This Phase calculates phihalf_correct
        .depends_on(RungeKutta_1);

    flecs::entity RungeKutta_3 = world.entity()
        .add(flecs::Phase) // This Phase calculates phi_end_predict
        .depends_on(RungeKutta_2);

    flecs::entity RungeKutta_4 = world.entity()
        .add(flecs::Phase) // This Phase calculates phi_end_correct
        .depends_on(RungeKutta_3);

    flecs::entity Update = world.entity()
        .add(flecs::Phase) // This phase replaces phi_start with phi_end_correct
        .depends_on(RungeKutta_4);

    // Create components inside the world
    world.component<Middle_Tag>();
    world.component<End_Tag>();

    world.component<position>();

    world.component<phi_start>();
    world.component<phihalf_predict>();
    world.component<phihalf_correct>();
    world.component<phi_end_predict>();
    world.component<phi_end_correct>();

    std::vector<flecs::entity> nodes; // place to store nodes
    nodes.reserve(N); // Create the space
    
    // Create left boundary Node
    nodes.push_back(
            world.entity()
                .add<End_Tag>() 
                .set<position>({0})
                .set<phi_start>({0})
                .set<phihalf_predict>({0})
                .set<phihalf_correct>({0})
                .set<phi_end_predict>({0})
            );
    
    // Create non-boundary Nodes
    for (int index = 1; index <= N-1; ++index) {
        double phi_initial;
        phi_initial = linear_function(index*L/N, Start, End);
        nodes.push_back(
            world.entity()
                .add<Middle_Tag>() 
                .set<position>({index*L/N})
                .set<phi_start>({phi_initial})
                .set<phihalf_predict>({0})
                .set<phihalf_correct>({0})
                .set<phi_end_predict>({0})
                .set<phi_end_correct>({0})
            );
    }

    // Create Right Boundary Nodes
    nodes.push_back(
            world.entity()
                .add<End_Tag>() 
                .set<position>({L})
                .set<phihalf_predict>({0})
                .set<phihalf_correct>({0})
                .set<phi_end_predict>({0})
            );

    // This system saves phi_start to the txt file
    world.system()
        .kind(Save);
    
    // This system finds and updates phihalf_predict
    world.system()
        .kind(RungeKutta_1);
    
    // This system finds and updates phihalf_correct
    world.system()
        .kind(RungeKutta_2);
    
    // This system finds and updates phi_end_predict
    world.system()
        .kind(RungeKutta_3);

    // This system finds and updates phi_end_correct
    world.system()
        .kind(RungeKutta_4);

    // This updates phi_start by replacing it with phi_end_correct
    world.system()
        .kind(Update);
    
    // Prepare Save File
    std::ofstream MyFile;
    MyFile.open("1D_dynamic_fluidtxt");

    // Check Save File is open/created
    if (!MyFile.is_open())
    {
        std::cout<<"Error in creating file"<<std::endl; 
        return 1; 
    }
    
    // Run through systems every time step
    for (int t_step = 0; t_step <= STEPS; ++t_step) {
        world.progress();
    }
    
    // Close Save File
    MyFile.close();
}
