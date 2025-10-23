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

std::vector<double> x_init = {0,0}; // Entity 1 initial position
std::vector<double> v_init = {2,1}; // Entity 1 initial velocity
std::vector<double> m_init = {3,5}; // Entity 1 initial mass
double acceleration(int index) { return -(k*x_init[index])/m_init[index] - b * v_init[index]; } // Function to calculate acceleration

float w1 = std::sqrt( k / m_init[0]); // Entity 1 - Frequency of the mass-spring system 
int time_period = ( (2 * M_PI) / w1 ); // Time the spring system runs for in s 
int time_unit = 100; 
int run_time = time_period * time_unit; 

int no_particles = 2; 

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
    double m; // in kg
};

int main() {

    std::cout<<x_init[0]<<"AHHHHH"<<std::endl; 

    // Open file for writing
    std::ofstream MyFile; 
    MyFile.open("SHM-Data.txt");
    if (!MyFile.is_open())
    {
        std::cout<<"Error in creating file"<<std::endl; 
        return 1; 
    }
    MyFile << run_time << std::endl;

    // Define vectors to store component data
    std::vector<double> pos; 
    std::vector<double> vel;
    std::vector<double> acc;

    flecs::world ecs;
    // Create a system for Position, Velocity, Acceleration..
    flecs::system s = ecs.system<Position, Velocity, Acceleration, const Mass>()
        //.with(flecs::Prefab).optional()
        .each([&](flecs::entity e, Position& p, Velocity& v, Acceleration& a, const Mass& mass) {

            a.x = - (k*p.x)/mass.m - 2 * b * v.x ;
            v.x += a.x / time_unit;
            p.x += v.x / time_unit;

            std::cout << e.name() << ": {" << p.x << ", " << v.x << "," << a.x << "}\n";

            pos.push_back(p.x); 
            vel.push_back(v.x);
            acc.push_back(a.x);
        });

    // Particle prefab, a template that we can use to generate entities
    flecs::entity Particle = ecs.prefab("Particle")
        .set<Index>({})
        .set<Position>({})
        .set<Velocity>({})
        .set<Acceleration>({})
        .set<Mass>({});

    for(int i = 0; i < no_particles; i++)
    {
        flecs::entity particle = ecs.entity().is_a(Particle); 
        particle.set(Index{i});
        particle.set(Position{x_init[i]});
        particle.set(Velocity{v_init[i]});
        particle.set(Acceleration{acceleration(i)});
        particle.set(Mass{m_init[i]});

    }
    
    // Run the system
    for(auto iter=run_time; iter--;) 
    {
        s.run();
        std::cout << "----\n";
    }
    for(int i = 0; i < run_time; i++) 
    {
        MyFile << pos[i] << ", " << vel[i] << "," << acc[i] << "," << i << std::endl; 
    }
    MyFile.close();
    
}


