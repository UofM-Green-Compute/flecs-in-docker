/*
This code models a one dimensional simple or damped hamonic 
oscillator. Hopefully we will also be able to produce visual
graphs and maybe even simulations.
*/
#include <iostream>
#include <flecs.h>
#include <systems.h>
#include <fstream>

double k = 2; // Spring Constant in N cm-1

double b = 10; // Damping Coefficient in s-1

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
    flecs::world ecs;
    // Create a system for Position, Velocity, Acceleration..
    flecs::system s = ecs.system<Position, Velocity, Acceleration, const Mass>()
        .each([](flecs::entity e, Position& p, Velocity& v, Acceleration& a, const Mass& mass) {
            a.x = - b*v.x - (k*p.x)/mass.m;
            v.x += a.x/(100);
            p.x += v.x/(100);
            std::cout << e.name() << ": {" << p.x << ", " << v.x << "," << a.x << "}\n";
        });

    ecs.entity("e1")
        .set<Position>({0})
        .set<Velocity>({1})
        .set<Acceleration>({0})
        .set<Mass>({3});
    
    std::ofstream MyFile("SFM-Data.txt");
    MyFile << "Hello file world" << std::endl;
    MyFile.close();

    for(auto iter=100; iter--;) {
        s.run();
        std::cout << "----\n";
    }
    
}
