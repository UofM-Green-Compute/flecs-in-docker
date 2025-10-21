/*
This code models a one dimensional simple or damped hamonic oscillator. 
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

double x1 = 0; // Entity 1 initial position
double v1 = 2; // Entity 1 initial velocity
double m1 = 3; // Entity 1 initial mass
double a1 = - (k*x1)/m1 - b * v1; // Entity 1 initial acceleration
float w1 = std::sqrt( k / m1); // Entity 1 - Frequency of the mass-spring system 
int time_period = ( (2 * M_PI) / w1 ); // Time the spring system runs for in s 
int time_unit = 100; 
int run_time = time_period * time_unit; 


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
    std::ofstream MyFile; 
    MyFile.open("SHM-Data.txt");
    if (!MyFile.is_open())
    {
        std::cout<<"Error in creating file"<<std::endl; 
        return 1; 
    }
    else
    {
        std::cout<<"All good"<<std::endl; 
    }
    MyFile << run_time << std::endl;

    std::vector<double> pos; 
    std::vector<double> vel;
    std::vector<double> acc;

    flecs::world ecs;
    // Create a system for Position, Velocity, Acceleration..
    flecs::system s = ecs.system<Position, Velocity, Acceleration, const Mass>()
        .each([&](flecs::entity e, Position& p, Velocity& v, Acceleration& a, const Mass& mass) {
            a.x = - (k*p.x)/mass.m - 2 * b * v.x ;
            v.x += a.x / time_unit;
            p.x += v.x / time_unit;
            std::cout << e.name() << ": {" << p.x << ", " << v.x << "," << a.x << "}\n";
            pos.push_back(p.x); 
            vel.push_back(v.x);
            acc.push_back(a.x);
        });

    ecs.entity("e1")
        .set<Position>({x1})
        .set<Velocity>({v1})
        .set<Acceleration>({a1})
        .set<Mass>({m1});
    
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


