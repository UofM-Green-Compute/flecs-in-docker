/*
This code models a one dimensional simple or damped hamonic 
oscillator. Hopefully we will also be able to produce visual
graphs and maybe even simulations.
*/
#include <iostream>
#include <flecs.h>
#include <systems.h>
#include <fstream>
#include <vector>

double k = 2; // Spring Constant in N cm-1

double b = 0.01; // Damping Coefficient in s-1

double current_time = 0; // Time in s

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
            a.x = - (k*p.x)/mass.m;
            v.x += a.x/(100);
            p.x += v.x/(100);
            std::cout << e.name() << ": {" << p.x << ", " << v.x << "," << a.x << "}\n";
            testing_pos.push_back(p.x); 
            testing_vel.push_back(v.x);
            testing_acc.push_back(a.x);
        });

    ecs.entity("e1")
        .set<Position>({0})
        .set<Velocity>({1})
        .set<Acceleration>({0})
        .set<Mass>({3});
    

    MyFile << "time,position,velocity,acceleration" << std::endl;
    for(auto iter=100; iter--;) {
        current_time += 0.01;
        std::cout << "t = " << current_time << "s\n";
        s.run();
        MyFile << current_time << "," << testing_pos[iter] << "," << testing_vel[iter] << "," << testing_acc[iter] << std::endl;
        std::cout << "----\n";
    }
    MyFile.close();
}