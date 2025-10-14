/*
There seems to be a problem where if you uncomment this code and then try make dockerbuild,
there is confusion because there are no two main() functions (one here and one in main.cpp).

Not sure what to do about that.
*/

/* 
DO NOT BUILD USING VSCODE !!!!!!!!!!
*/

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

double k = 2; // Spring Constant in Nm-1
double gamma = 10; 

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
            a.x = - (k*p.x)/mass.m - 2 * gamma * v.x ;
            v.x += a.x / 100;
            p.x += v.x / 100;
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
    
    for(auto iter=100; iter--;) {
        s.run();
        std::cout << "----\n";
    }

    for(auto i=100; i--;) {
        MyFile << testing_pos[i] << ", " << testing_vel[i] << "," << testing_acc[i] << std::endl;
    }
    MyFile.close();
    
}


