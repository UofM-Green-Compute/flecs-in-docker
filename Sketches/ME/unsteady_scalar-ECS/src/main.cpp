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
int TIME = 1;

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
struct Position {double x; };

// Field components
struct Phi_start { double i; };
struct Phihalf_predict { double i; };
struct Phihalf_correct { double i; };
struct Phi_end_predict { double i; };
struct Phi_end_correct { double i; };

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
    // Prepare Save File
    std::ofstream MyFile;
    MyFile.open("1D_dynamic_fluidtxt");

    // Check Save File is open/created
    if (!MyFile.is_open())
    {
        std::cout<<"Error in creating file"<<std::endl; 
        return 1; 
    }
    // Create World
    flecs::world world(argc, argv);

    // Creates Phases which tell the program in which order to run the systems
    /*
    flecs::entity Save = world.entity()
        .add(flecs::Phase); // This Phase saves the state in its current form
    */
    flecs::entity RungeKutta_1 = world.entity()
        .add(flecs::Phase); // This Phase calculates phihalf_predict
        //.depends_on(Save);

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

    world.component<Position>();

    world.component<Phi_start>();
    world.component<Phihalf_predict>();
    world.component<Phihalf_correct>();
    world.component<Phi_end_predict>();
    world.component<Phi_end_correct>();

    // Creates Neighbour Relationships
    auto Left_Neighbour = world.entity();
    auto Right_Neighbour = world.entity();

    std::vector<flecs::entity> nodes; // place to store nodes
    nodes.reserve(N); // Create the space
    
    // Create left boundary Node
    nodes.push_back(
            world.entity()
                .add<End_Tag>() 
                .set<Position>({0})
                .set<Phi_start>({Start})
                .set<Phihalf_predict>({Start})
                .set<Phihalf_correct>({Start})
                .set<Phi_end_predict>({Start})
            );
    
    // Create non-boundary Nodes
    for (int index = 1; index < N-1; ++index) {
        double phi_initial;
        phi_initial = linear_function(index*L/N, Start, End);
        nodes.push_back(
            world.entity()
                .add<Middle_Tag>() 
                .set<Position>({index*L/N})
                .set<Phi_start>({phi_initial})
                .set<Phihalf_predict>({0})
                .set<Phihalf_correct>({0})
                .set<Phi_end_predict>({0})
                .set<Phi_end_correct>({0})
            );
    }

    // Create Right Boundary Nodes
    nodes.push_back(
            world.entity()
                .add<End_Tag>() 
                .set<Position>({L})
                .set<Phi_start>({End})
                .set<Phihalf_predict>({End})
                .set<Phihalf_correct>({End})
                .set<Phi_end_predict>({End})
            );
    
    

    // This loop adds left and right neighbours to each non-boundary node
    

    for (int index = 1; index < N-1; ++index) {
        nodes[index].add(Left_Neighbour, nodes[index-1]);
        nodes[index].add(Right_Neighbour, nodes[index+1]);
    }

    // This system saves phi_start to the txt file
    /*
    world.system<Position, Phi_start>()
        .kind(Save)
        .each([&](Position pos, Phi_start phi){
            if (pos.x == L) {
                MyFile << phi.i << "\n";
            } else {
                MyFile << phi.i << ", ";
            }

        });
    */

    // This system finds and updates phihalf_predict
    world.system<Position, Phi_start, Phihalf_predict>()
        .with<Middle_Tag>()
        .kind(RungeKutta_1)
        .each([&](Position pos, Phi_start phi0, Phihalf_predict phi_half){
            double Phi_left;
        });
    
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
    
    // Set .txt file headers
    MyFile << "t, ";
    for (int x_step = 0; x_step < N; ++x_step) {
        if (x_step == N-1) {
            MyFile << "x" << x_step << "\n";
        } else {
            MyFile << "x" << x_step << ", ";
        }
    }
    
    // Run through systems every time step
    for (int t_step = 0; t_step <= STEPS; ++t_step) {
        MyFile << t_step << ", ";
        for (int index = 0; index < N-1; ++index) {
            const Position& p = nodes[index].get<Position>();
            MyFile << p.x << ", ";
            std::cout << p.x << ", ";
        }
        const Position& p = nodes[N-1].get<Position>();
        MyFile << p.x << "\n";
        std::cout << p.x << ", ";
        world.progress();
    }
    
    // Close Save File
    MyFile.close();
}
