/*
This code runs a simulation of two coupled oscillators.
The two masses are connected to walls by springs and to each other by springs.
This is a one dimensional simulation.
*/
#include <iostream>
#include <fstream> 
#include <vector>
#include <flecs.h>
#include <systems.h>


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

struct Spring_Constant{
    double k; // in Ncm-1
};

struct Damping_Coefficient{
    double b; // in s-1
};

struct Spring_Length{
    double l; // in cm
};

struct Spring_Number{
    int num; // Tells program which spring it is
}

struct Particle{ }; // creates a particle tag

struct Spring{ }; // creates a spring tag

int main() {
    std::ofstream MyFile("SHM-Data.txt");
    std::vector<double> testing_pos; 
    std::vector<double> testing_vel;
    std::vector<double> testing_acc;

    flecs::world ecs;

    flecs::query<const position, Particle> particles = 
        ecs.query<const Position, Particle>("ParticleQuery");

    flecs::query<const position, Spring> spring = 
        ecs.query<const Position, Spring>("SpringQuery");
    // Create a system to find initial acceleration
    /*
    flecs::system initial = ecs.system<>()
    */

    // Create a system for Position, Velocity, Acceleration..
    /*
    flecs::system s = ecs.system<Position, Velocity, Acceleration, const Mass>()
        .each([&](flecs::entity e, Position& p, Velocity& v, Acceleration& a, const Mass& mass) {
            a.x = - (k * p.x) / mass.m - b * v.x ;
            v.x += a.x / 100;
            p.x += v.x / 100;
            std::cout << e.name() << ": {" << p.x << "," << v.x << "," << a.x << "}\n";
            testing_pos.push_back(p.x); 
            testing_vel.push_back(v.x);
            testing_acc.push_back(a.x);
        });
    */

    // Create Entities
    ecs.entity("mass 1")
        // Finds and sets components
        .set<Position>({3})
        .set<Velocity>({-1})
        .set<Acceleration>({0})
        .set<Mass>({4})

        // Adds Particle Tag
        .add<Particle>();
    
    ecs.entity("mass 2")
        // Finds and sets components
        .set<Position>({5})
        .set<Velocity>({3})
        .set<Acceleration>({0})
        .set<Mass>({4})

        // Adds a Particle Tag
        .add<Particle>();

    ecs.entity("spring 1")
        //Finds and sets components
        .set<Spring_Constant>({2})
        .set<Damping_Coefficient>({1})
        .set<Spring_Length>({2})
        .set<Spring_Number>({1})

        // Add a Spring Tag
        .add<Spring>();
    
    ecs.entity("spring 2")
        //Finds and sets components
        .set<Spring_Constant>({2})
        .set<Damping_Coefficient>({1})
        .set<Spring_Length>({2})
        .set<Spring_Number>({2})

        // Add a Spring Tag
        .add<Spring>();

    ecs.entity("spring 3")
        //Finds and sets components
        .set<Spring_Constant>({2})
        .set<Damping_Coefficient>({1})
        .set<Spring_Length>({2})
        .set<Spring_Number>({3})

        // Add a Spring Tag
        .add<Spring>();
    
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