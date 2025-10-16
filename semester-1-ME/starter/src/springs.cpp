/*
This code models a one dimensional simple or damped hamonic 
oscillator. Hopefully we will also be able to produce visual
graphs and maybe even simulations.
*/
#include <iostream>
#include <fstream> 
#include <vector>
#include <flecs.h>
#include <systems.h>

// Defines Spring Properties
double k = 2; // Spring Constant in Nm-1
double b = 1; // Damping coefficient in s-1

// Defines Particle Properties
double x1 = 0; // Entity 1 initial position
double v1 = 2; // Entity 1 initial velocity
double m1 = 3; // Entity 1 initial mass
double a1 = - (k*x1)/m1 - b * v1; // Entity 1 initial acceleration

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
    std::ofstream MyFile("SHM-Data.txt");
    std::vector<double> testing_pos; 
    std::vector<double> testing_vel;
    std::vector<double> testing_acc;

    flecs::world ecs;
    // Create a system for Position, Velocity, Acceleration..
    flecs::system s = ecs.system<Position, Velocity, Acceleration, const Mass>()
        .each([&](flecs::entity e, Position& p, Velocity& v, Acceleration& a, const Mass& mass) {
            a.x = - (k*p.x)/mass.m - b * v.x ;
            v.x += a.x / 100;
            p.x += v.x / 100;
            std::cout << e.name() << ": {" << p.x << "," << v.x << "," << a.x << "}\n";
            testing_pos.push_back(p.x); 
            testing_vel.push_back(v.x);
            testing_acc.push_back(a.x);
        });

    ecs.entity("e1")
        .set<Position>({x1})
        .set<Velocity>({v1})
        .set<Acceleration>({a1})
        .set<Mass>({m1});
    
    for(auto iter=0; iter<100; iter++) {
        if (iter == 0) {
            testing_pos.push_back(x1); 
            testing_vel.push_back(v1);
            testing_acc.push_back(a1);
            std::cout << "e1: {" << x1 << ", " << v1 << "," << a1 << "}\n";
        } else {
        s.run();
        }
        std::cout << "----\n";
        MyFile << testing_pos[iter] << "," << testing_vel[iter] << "," << testing_acc[iter] << std::endl;  
    }

}