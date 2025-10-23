/*
This code models a system of N (=2) coupled oscillators. 

Structure of the code
- 
- 
- 

Aims: 
- Produce graphs
- Animate simulations
*/
// DO NOT BUILD USING VSCODE !!!!!!!!!!

#include <iostream>
#include <fstream> 
#include <vector>
#include <cmath>
#include <flecs.h>
#include <systems.h>

double k = 2; // Spring Constant in Nm-1
double b = 0; // Damping coefficient in s-1

// Define vectors to store component data
    std::vector<double> pos; 
    std::vector<double> vel;
    std::vector<double> acc;

double particle_mass = 5; 
float w1 = std::sqrt( k / particle_mass); // Frequency of the mass-spring system 
int time_period = ( (2 * M_PI) / w1 ); // Time the spring system runs for in s 
float time_step = 0.01; // Time elapsed in each time step 
int run_time = time_period * 1/time_step; // Number of loops to run 

int N = 2; // Number of particles in the system

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
    double x; // In cms-1
};
struct Acceleration{
    double x; // in cms-2
};
struct Mass{
    double M; // in kg
};

double acceleration(float position, float position_left, float position_right,
const float mass, float k_left, float k_right, float l_left, float l_right){

    float acc;
    acc = - (k_left/mass) * (position - position_left - l_left) 
        + (k_right/mass) * (position_right - position - l_right);
    
    return acc;
}

int main() {

    // Open file for writing
    std::ofstream MyFile; 
    MyFile.open("Coupled_Oscillators.txt");
    if (!MyFile.is_open())
    {
        std::cout<<"Error in creating file"<<std::endl; 
        return 1; 
    }
    MyFile << "No iterations: " << run_time << std::endl;
    MyFile << "Time (s), Position 1 (cm), Velocity 1 (cm s-1), Acceleration 1 (cm s-2)" 
           << ", Position 2 (cm), Velocity 2 (cm s-1), Acceleration 2 (cm s-2)" << std::endl;

    // First Particle Data
    std::vector<std::vector<double>> p_matrix {{2}, {7}}; 
    std::vector<std::vector<double>> v_matrix {{2}, {-1}};
    std::vector<std::vector<double>> a_matrix {{0}, {0}};

    const std::vector<double> k_list {3, 3, 3}; // Spring Constants
    const std::vector<double> l_list {4, 4, 4}; // Natural Spring Lengths
    
    flecs::world ecs;

    flecs::system s1 = ecs.system<Acceleration, Position, Mass, Index>()
    .each([&](Acceleration &a, Position &p, const Mass &mass, const Index &index) 
    {
        if (index.i == 0) {
            a.x = acceleration(p.x, p_Lwall, p_matrix[index.i+1][0], mass.M, k_list[index.i], 
                                    k_list[index.i+1], l_list[index.i], l_list[index.i+1]);
        } else if (index.i == N-1) {
            a.x = acceleration(p.x, p_matrix[index.i-1][0], p_Rwall, mass.M, k_list[index.i], 
                                    k_list[index.i+1], l_list[index.i], l_list[index.i+1]);
        } else {
            a.x = acceleration(p.x, p_matrix[index.i-1][0], p_matrix[index.i+1][0], mass.M, k_list[index.i], 
                                    k_list[index.i+1], l_list[index.i], l_list[index.i+1]);
        } 
        a_matrix[index.i][0] = a.x;
    });

    flecs::system s2 = ecs.system<Position, Velocity, Acceleration, Mass, Index>()
    .each([&](flecs::entity e, Position &p, Velocity &v, Acceleration &a, const Mass &mass, const Index &index) 
    {
        p.x += v.x*time_step;
        v.x += a.x*time_step;
        
        if (index.i == 0) {
            a.x = acceleration(p.x, p_Lwall, p_matrix[index.i+1][0], mass.M, k_list[index.i], 
                                    k_list[index.i+1], l_list[index.i], l_list[index.i+1]);
        } else if (index.i == N-1) {
            a.x = acceleration(p.x, p_matrix[index.i-1][0], p_Rwall, mass.M, k_list[index.i], 
                                    k_list[index.i+1], l_list[index.i], l_list[index.i+1]);
        } else {
            a.x = acceleration(p.x, p_matrix[index.i-1][0], p_matrix[index.i+1][0], mass.M, k_list[index.i], 
                                    k_list[index.i+1], l_list[index.i], l_list[index.i+1]);
        } 
        
        p_matrix[index.i].push_back(p.x);
        v_matrix[index.i].push_back(v.x);
        a_matrix[index.i].push_back(a.x);
        
        std::cout << e.name() << ": {" << p.x << ", " << v.x << "," << a.x << "}\n";
    });

    // Particle prefab, a template that we can use to generate entities
    flecs::entity Particle = ecs.prefab("Particle")
        .set<Index>({})
        .set<Position>({})
        .set<Velocity>({})
        .set<Acceleration>({})
        .set<Mass>({});

    for(int i = 0; i < N; i++)
    {
        flecs::entity particle = ecs.entity().is_a(Particle); 
        particle.set(Index{i});
        particle.set(Position{p_matrix[0][0]});
        particle.set(Velocity{v_matrix[0][0]});
        particle.set(Acceleration{0}); 
        particle.set(Mass{particle_mass});
    }

    s1.run();

    // Run the system
    for(auto iter=run_time; iter--;) 
    {
        s2.run();
        std::cout << "----\n";
    }
    for(int i = 0; i < run_time; i++) 
    {
        MyFile << i*time_step << ", " << p_matrix[0][i] << ", " << v_matrix[0][i] << "," << a_matrix[0][i] << "," 
               << p_matrix[1][i] << ", " << v_matrix[1][i] << "," << a_matrix[1][i] << std::endl; 
    }
    MyFile.close();
    
}


