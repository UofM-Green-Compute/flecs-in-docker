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
float time_step = 0.5; // Time elapsed in each time step in s
int run_time = (period_number * time_period) / time_step; // Number of loops to run 

int N = 2; // Number of particles in the system

// Initial Particle Data
std::vector<double> p_initial {1, 2}; 
std::vector<double> v_initial {-1, 1};

//Position of the Left Wall
double p_Lwall = 0;  // in m

// Position of the Right Wall
double p_Rwall = 3; // in m

// Index Component
struct Index {
    int i; // Which particle is it
};
// Position Components
struct PositionStart { double x;};
struct PositionHalfPredict { double x;};
struct PositionHalfCorrect { double x;};
struct PositionEndPredict { double x;};
// Velocity components
struct VelocityStart { double x;};
struct VelocityHalfPredict { double x;};
struct VelocityHalfCorrect { double x;};
struct VelocityEndPredict { double x;};
// Accerleration components
struct AccelerationStart { double x;};
struct AccelerationHalfPredict { double x;};
struct AccelerationHalfCorrect { double x;};
struct AccelerationEndPredict { double x;};
// Mass component
struct Mass{ double M;};
// Bulk tag ensure some code doesnt run on initial loop
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
    flecs::entity EnergyStart = world.entity()
        .add(flecs::Phase);
    flecs::entity RK1 = world.entity()
        .add(flecs::Phase)
        .depends_on(EnergyStart);
    flecs::entity A1 = world.entity()
        .add(flecs::Phase)
        .depends_on(RK1);
    flecs::entity RK2 = world.entity()
        .add(flecs::Phase)
        .depends_on(A1);
    flecs::entity A2 = world.entity()
        .add(flecs::Phase)
        .depends_on(RK2);
    flecs::entity RK3 = world.entity()
        .add(flecs::Phase)
        .depends_on(A2);
    flecs::entity A3 = world.entity()
        .add(flecs::Phase)
        .depends_on(RK3);
    flecs::entity RK_Update = world.entity()
        .add(flecs::Phase)
        .depends_on(A3);
    flecs::entity A_Update = world.entity()
        .add(flecs::Phase)
        .depends_on(RK_Update);
    flecs::entity EnergyEnd = world.entity()
        .add(flecs::Phase)
        .depends_on(A_Update);

    // Create components inside world
    world.component<Index>();
    world.component<PositionStart>();
    world.component<PositionHalfPredict>();
    world.component<PositionHalfCorrect>();
    world.component<PositionEndPredict>();
    world.component<VelocityStart>();
    world.component<VelocityHalfPredict>();
    world.component<VelocityHalfCorrect>();
    world.component<VelocityEndPredict>();
    world.component<AccelerationStart>();
    world.component<AccelerationHalfPredict>();
    world.component<AccelerationHalfCorrect>();
    world.component<AccelerationEndPredict>();
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
                .set<PositionStart>({p_initial[index]})
                .set<PositionHalfPredict>({0})
                .set<PositionHalfCorrect>({0})
                .set<PositionEndPredict>({0})
                .set<VelocityStart>({v_initial[index]})
                .set<VelocityHalfPredict>({0})
                .set<VelocityHalfCorrect>({0})
                .set<VelocityEndPredict>({0})
                .set<AccelerationStart>({0})
                .set<AccelerationHalfPredict>({0})
                .set<AccelerationHalfCorrect>({0})
                .set<AccelerationEndPredict>({0})
            );
    }

    world.system<>()
        .kind(EnergyStart)
        .each([&]() {
            energy_tracker.start();
        });
    
    world.system<PositionStart, PositionHalfPredict, VelocityStart, VelocityHalfPredict,
                 AccelerationStart>()
        .with<BulkTag>()
        .kind(RK1)
        .each([&](PositionStart& posStart, PositionHalfPredict& posHalf, VelocityStart& velStart, 
        VelocityHalfPredict& velHalf, AccelerationStart& accStart){
            posHalf.x = posStart.x + ((time_step/2)*velStart.x);
            velHalf.x = velStart.x + ((time_step/2)*accStart.x);
        });

    world.system<Index, PositionHalfPredict, AccelerationHalfPredict, Mass>()
        .with<BulkTag>()
        .kind(A1)
        .each([&](Index& ind, PositionHalfPredict& pos, AccelerationHalfPredict& acc, Mass& mass){
            double p_left;
            double p_right;
            if (ind.i == 0) {
                p_left = p_Lwall;
                const PositionHalfPredict& p_r = nodes[ind.i + 1].get<PositionHalfPredict>();
                p_right = p_r.x;
            } else if (ind.i == N-1) {
                const PositionHalfPredict& p_l = nodes[ind.i - 1].get<PositionHalfPredict>();
                p_left = p_l.x;
                p_right = p_Rwall;
            } else {
                const PositionHalfPredict& p_l = nodes[ind.i - 1].get<PositionHalfPredict>();
                p_left = p_l.x;
                const PositionHalfPredict& p_r = nodes[ind.i + 1].get<PositionHalfPredict>();
                p_right = p_r.x;
            }
        acc.x = acceleration(pos.x, p_left, p_right, mass.M, k_list[ind.i], 
                                    k_list[ind.i+1], l_list[ind.i], l_list[ind.i+1]);
        });

    world.system<PositionStart, PositionHalfCorrect, VelocityStart, 
                 VelocityHalfPredict, VelocityHalfCorrect, AccelerationHalfPredict>()
        .with<BulkTag>()
        .kind(RK2)
        .each([&](PositionStart& posStart, PositionHalfCorrect& posCorr,  VelocityStart& velStart, 
        VelocityHalfPredict& velPred, VelocityHalfCorrect& velCorr, 
        AccelerationHalfPredict& accPred){
            posCorr.x = posStart.x + ((time_step/2)*velPred.x);
            velCorr.x = velStart.x + ((time_step/2)*accPred.x);
        });

    world.system<Index, PositionHalfCorrect, AccelerationHalfCorrect, Mass>()
        .with<BulkTag>()
        .kind(A2)
        .each([&](Index& ind, PositionHalfCorrect& pos, AccelerationHalfCorrect& acc, Mass& mass){
            double p_left;
            double p_right;
            if (ind.i == 0) {
                p_left = p_Lwall;
                const PositionHalfCorrect& p_r = nodes[ind.i + 1].get<PositionHalfCorrect>();
                p_right = p_r.x;
            } else if (ind.i == N-1) {
                const PositionHalfCorrect& p_l = nodes[ind.i - 1].get<PositionHalfCorrect>();
                p_left = p_l.x;
                p_right = p_Rwall;
            } else {
                const PositionHalfCorrect& p_l = nodes[ind.i - 1].get<PositionHalfCorrect>();
                p_left = p_l.x;
                const PositionHalfCorrect& p_r = nodes[ind.i + 1].get<PositionHalfCorrect>();
                p_right = p_r.x;
            }
        acc.x = acceleration(pos.x, p_left, p_right, mass.M, k_list[ind.i], 
                                    k_list[ind.i+1], l_list[ind.i], l_list[ind.i+1]);
        });

    world.system<PositionStart, PositionEndPredict, VelocityStart, VelocityHalfCorrect, 
                 VelocityEndPredict, AccelerationHalfCorrect>()
        .with<BulkTag>()
        .kind(RK3)
        .each([&](PositionStart& posStart, PositionEndPredict& posEnd, VelocityStart& velStart, 
        VelocityHalfCorrect& velHalf, VelocityEndPredict& velEnd, AccelerationHalfCorrect& accHalf){
            posEnd.x = posStart.x + ((time_step/2)*velHalf.x);
            velEnd.x = velStart.x + ((time_step/2)*accHalf.x);
        });

    world.system<Index, PositionEndPredict, AccelerationEndPredict, Mass>()
        .with<BulkTag>()
        .kind(A3)
        .each([&](Index& ind, PositionEndPredict& pos, AccelerationEndPredict& acc, Mass& mass){
            double p_left;
            double p_right;
            if (ind.i == 0) {
                p_left = p_Lwall;
                const PositionEndPredict& p_r = nodes[ind.i + 1].get<PositionEndPredict>();
                p_right = p_r.x;
            } else if (ind.i == N-1) {
                const PositionEndPredict& p_l = nodes[ind.i - 1].get<PositionEndPredict>();
                p_left = p_l.x;
                p_right = p_Rwall;
            } else {
                const PositionEndPredict& p_l = nodes[ind.i - 1].get<PositionEndPredict>();
                p_left = p_l.x;
                const PositionEndPredict& p_r = nodes[ind.i + 1].get<PositionEndPredict>();
                p_right = p_r.x;
            }
        acc.x = acceleration(pos.x, p_left, p_right, mass.M, k_list[ind.i], 
                                    k_list[ind.i+1], l_list[ind.i], l_list[ind.i+1]);
        });

    world.system<PositionStart, VelocityStart, VelocityHalfPredict, VelocityHalfCorrect, 
                 VelocityEndPredict, AccelerationStart, AccelerationHalfPredict, 
                 AccelerationHalfCorrect, AccelerationEndPredict>()
        .with<BulkTag>()
        .kind(RK_Update)
        .each([&](PositionStart& pos, VelocityStart& velStart, VelocityHalfPredict& velPred, 
        VelocityHalfCorrect& velCorr, VelocityEndPredict& velEnd, AccelerationStart& accStart, 
        AccelerationHalfPredict& accPred, AccelerationHalfCorrect& accCorr, 
        AccelerationEndPredict& accEnd){
            pos.x += ((time_step/6)*(velStart.x+2*velPred.x+2*velCorr.x+velEnd.x));
            velStart.x += ((time_step/6)*(accStart.x+2*accPred.x+2*accCorr.x+accEnd.x));
        });

    world.system<Index, PositionStart, AccelerationStart, Mass>()
        .kind(A_Update)
        .each([&](Index& ind, PositionStart& pos, AccelerationStart& acc, Mass& mass){
            double p_left;
            double p_right;
            if (ind.i == 0) {
                p_left = p_Lwall;
                const PositionStart& p_r = nodes[ind.i + 1].get<PositionStart>();
                p_right = p_r.x;
            } else if (ind.i == N-1) {
                const PositionStart& p_l = nodes[ind.i - 1].get<PositionStart>();
                p_left = p_l.x;
                p_right = p_Rwall;
            } else {
                const PositionStart& p_l = nodes[ind.i - 1].get<PositionStart>();
                p_left = p_l.x;
                const PositionStart& p_r = nodes[ind.i + 1].get<PositionStart>();
                p_right = p_r.x;
            }
        acc.x = acceleration(pos.x, p_left, p_right, mass.M, k_list[ind.i], 
                                    k_list[ind.i+1], l_list[ind.i], l_list[ind.i+1]);
        });
    
     world.system<>()
        .kind(EnergyEnd)
        .each([&]() {
            auto r = energy_tracker.stop();
        });
    
    
    world.progress();
    const PositionStart& p1 = nodes[0].get<PositionStart>();
    const PositionStart& p2 = nodes[1].get<PositionStart>();
    const VelocityStart& v1 = nodes[0].get<VelocityStart>();
    const VelocityStart& v2 = nodes[1].get<VelocityStart>();
    const AccelerationStart& a1 = nodes[0].get<AccelerationStart>();
    const AccelerationStart& a2 = nodes[1].get<AccelerationStart>();
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

        const PositionStart& p1 = nodes[0].get<PositionStart>();
        const PositionStart& p2 = nodes[1].get<PositionStart>();
        const VelocityStart& v1 = nodes[0].get<VelocityStart>();
        const VelocityStart& v2 = nodes[1].get<VelocityStart>();
        const AccelerationStart& a1 = nodes[0].get<AccelerationStart>();
        const AccelerationStart& a2 = nodes[1].get<AccelerationStart>();
        MyFile << i*time_step << ", " << p1.x << ", " << v1.x << "," << a1.x << "," << p2.x << ", " 
        << v2.x << "," << a2.x << std::endl; 
    }
    MyFile.close();
    std::cout << time_period << "\n";
}

