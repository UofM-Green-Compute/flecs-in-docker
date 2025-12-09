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
#include <ccenergy/EnergyTracker.hpp>

#include <iostream>
#include <fstream> 
#include <vector>
#include <cmath>
#include <flecs.h>
#include <systems.h>

double l = 1; // natural spring length
double k = 4*M_PI*M_PI; // spring constant in N m-2
double particle_mass = 3; // mass in kg
double omega = std::sqrt((3*k)/(particle_mass)); // angular frequency in rad s-1
int time_period = ( (2 * M_PI) / omega ); // period of oscillations in s
int period_number = 2; // number of period oscillations
float time_step = 0.00005; // Time elapsed in each time step in s
int run_time = (period_number * time_period) / time_step; // Number of loops to run 

int N = 2; // Number of particles in the system

// Initial Particle Data
std::vector<double> p_initial {1, 2}; 
std::vector<double> v_initial {-1, 1};

//Position of the Left Wall
double p_Lwall = 0;  // in m

// Position of the Right Wall
double p_Rwall = 3; // in m

struct Index {
    int i; // Which particle is it
};
struct Position { 
    double x; // In m
};
struct Velocity { 
    double x; // In ms-1
};
struct Acceleration{
    double x; // in ms-2
};
struct Mass{
    double M; // in kg
};

struct BulkTag {};

double acceleration(float position, float position_left, float position_right,
const float mass, float k_left, float k_right, float l_left, float l_right)
{
    float acc;
    acc = - (k_left/mass) * (position - position_left - l_left) 
        + (k_right/mass) * (position_right - position - l_right);
    
    return acc;
}

int main(int argc, char* argv[]) {

    std::ofstream MyFile; 
    MyFile.open("Coupled_Oscillators.txt");
    if (!MyFile.is_open())
    {
        std::cout<<"Error in creating file"<<std::endl; 
        return 1; 
    }
    MyFile << "No iterations: " << run_time << std::endl;
    MyFile << "Time (s), Position 1 (m), Velocity 1 (m s-1), Acceleration 1 (m s-2)" 
           << ", Position 2 (m), Velocity 2 (m s-1), Acceleration 2 (m s-2)" << std::endl;
    
    const std::vector<double> k_list {k, k, k}; // Spring Constants
    const std::vector<double> l_list {l, l, l}; // Natural Spring Lengths
    
    flecs::world world(argc, argv);

    // Create the energy tracker
    ccenergy::EnergyTracker energy_tracker {{ .label = "OnUpdate",
                                              .measure_cpu = true,
                                              .measure_gpu  = false,
                                              .log_to_stdout = false }};

    // Creates Phases which tell the program in which order to run the systems

    flecs::entity position_phase = world.entity()
        .add(flecs::Phase); // This Phase calculates the new position

    flecs::entity velocity_phase = world.entity()
        .add(flecs::Phase) // This Phase calculates the new velocity
        .depends_on(position_phase);
    
    flecs::entity acceleration_phase = world.entity()
        .add(flecs::Phase) // This Phase calculates the new acceleration
        .depends_on(velocity_phase);

    // Create components inside world
    world.component<Index>();
    world.component<Position>();
    world.component<Velocity>();
    world.component<Acceleration>();
    world.component<Mass>();
    world.component<BulkTag>();

    // Initialize the nodes
    std::vector<flecs::entity> nodes;
    nodes.reserve(N);

    // Create Nodes
    for (int index = 0; index <= N-1; ++index) {
        nodes.push_back(
            world.entity()
                .set<Index>({index})
                .set<Mass>({particle_mass})
                .set<Position>({p_initial[index]})
                .set<Velocity>({v_initial[index]})
                .set<Acceleration>({0})
            );
    }

    world.system<>()
        .kind(flecs::PreUpdate)
        .each([&]() {
            energy_tracker.start();
        });

    world.system<Position, Velocity>()
        .with<BulkTag>()
        .kind(position_phase)
        .each([&](Position& pos, Velocity& vel){
            pos.x += vel.x*time_step;
        });

    world.system<Velocity, Acceleration>()
        .with<BulkTag>()
        .kind(velocity_phase)
        .each([&](Velocity& vel, Acceleration& acc){
            vel.x += acc.x*time_step;
        });
    
    world.system<Index, Position, Acceleration, Mass>()
        .kind(acceleration_phase)
        .each([&](Index& ind, Position& pos, Acceleration& acc, Mass& mass){
            double p_left;
            double p_right;
            if (ind.i == 0) {
                p_left = p_Lwall;
                const Position& p_r = nodes[ind.i + 1].get<Position>();
                p_right = p_r.x;
            } else if (ind.i == N-1) {
                const Position& p_l = nodes[ind.i - 1].get<Position>();
                p_left = p_l.x;
                p_right = p_Rwall;
            } else {
                const Position& p_l = nodes[ind.i - 1].get<Position>();
                p_left = p_l.x;
                const Position& p_r = nodes[ind.i + 1].get<Position>();
                p_right = p_r.x;
            }
        acc.x = acceleration(pos.x, p_left, p_right, mass.M, k_list[ind.i], 
                                    k_list[ind.i+1], l_list[ind.i], l_list[ind.i+1]);
        });
    
    world.progress();
    const Position& p1 = nodes[0].get<Position>();
    const Position& p2 = nodes[1].get<Position>();
    const Velocity& v1 = nodes[0].get<Velocity>();
    const Velocity& v2 = nodes[1].get<Velocity>();
    const Acceleration& a1 = nodes[0].get<Acceleration>();
    const Acceleration& a2 = nodes[1].get<Acceleration>();
    MyFile << 0 << ", " << p1.x << ", " << v1.x << "," << a1.x << "," << p2.x << ", " 
    << v2.x << "," << a2.x << std::endl; 
    // Add Bulk tag to each entity in nodes
    for (flecs::entity& e : nodes) {
        	e.add<BulkTag>(); 
    }

    // Run the system
    for (int i = 1; i <= run_time; i++) {
        world.progress();
        std::cout << i << "\n";
        std::cout << "----\n";

        const Position& p1 = nodes[0].get<Position>();
        const Position& p2 = nodes[1].get<Position>();
        const Velocity& v1 = nodes[0].get<Velocity>();
        const Velocity& v2 = nodes[1].get<Velocity>();
        const Acceleration& a1 = nodes[0].get<Acceleration>();
        const Acceleration& a2 = nodes[1].get<Acceleration>();
        MyFile << i*time_step << ", " << p1.x << ", " << v1.x << "," << a1.x << "," << p2.x << ", " 
        << v2.x << "," << a2.x << std::endl; 
    }
    MyFile.close();
    std::cout << time_period << "\n";
    std::cout << energy_tracker.mkReport() << std::endl;
}
