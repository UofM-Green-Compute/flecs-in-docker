/*
Oluwole Delano
13/11/25 
Smooth Particle Hydrodynamics using the ECS Model, rendering in ASCII
*/

#include <iostream>
#include <flecs.h>
#include <vector>
#include <random>
#include <cmath>

// Initialise Components 
struct Position { double x, y; };
struct Velocity { double dx, dy; };
struct Acceleration { double ddx, ddy; };
struct Mass { double m; };
struct ParticleTag {}; 

// Constants
int NO_PARTICLES = 5; 

// Walls of the domain 
static constexpr double HALF_WIDTH = 40.0;           // half-width  (domain x in [-W, W])
static constexpr double HALF_HEIGHT = 20.0;           // half-height (domain y in [-H, H])
static constexpr double BOX_W = 2.0 * HALF_WIDTH;    // full width
static constexpr double BOX_H = 2.0 * HALF_HEIGHT;    // full height

double CONST_H = 3; 

// Interpolant function (Gaussian)
double W(double x) { return ( 1 / CONST_H * sqrt(M_PI) ) * exp(( - pow(x,2) / pow(CONST_H,2) )); }

// Function to calculate vector distance between particles 
std::vector<double> vector_distance(flecs::entity Particle, flecs::entity Particle_i){
    std::vector<double> displacement; 

    Position p = Particle.get<Position>();
    Position p_i = Particle_i.get<Position>();

    displacement.push_back(p.x - p_i.x);
    displacement.push_back(p.y - p_i.y);

    return displacement; 
}

// Function to caluclate density rho 
double density(std::vector<flecs::entity> Particles, int particle_index){
    double density = 0; 
    for (int i = 0; i < Particles.size()-1; i++) 
    {
        // calculate distance between current particle (particle_index) and ith particle
        // add each new term to density 
    }

}

// Function to calculate force on particle a due to particle b
std::vector<double> force(flecs::entity Particle, flecs::entity Particle_i){
    std::vector<double> Force; 
    return Force; 
}

int main(int argc, char* argv[]) {

    // Create the flecs world
    flecs::world world(argc,argv);
    // Components of the world
    world.component<Position>(); 
    world.component<Velocity>(); 
    world.component<Acceleration>(); 
    world.component<Mass>(); 

    // Tools for picking random numbers
    std::mt19937 rng( std::random_device{}()  ) ; // Initialise a random number generator with random device (for actual use)

    std::uniform_real_distribution<double> UposX(-HALF_WIDTH, HALF_WIDTH); 
    std::uniform_real_distribution<double> UposY(-HALF_HEIGHT, HALF_HEIGHT); 
    std::uniform_real_distribution<double> Uvel(-0.4, 0.4); 
    std::uniform_real_distribution<double> Umass(0.5, 2.0); 

    // Initialise entities | Generate random values for initial conditions
    std::vector<flecs::entity> particles; 
    particles.reserve(NO_PARTICLES); 

    for (int i = 0; i < NO_PARTICLES; ++i) { 
        particles.push_back( 
            world.entity() 
                .add<ParticleTag>() 
                .set<Position>({UposX(rng), UposY(rng)}) 
                .set<Velocity>({Uvel(rng), Uvel(rng)}) 
                .set<Acceleration>({0.0, 0.0}) 
                .set<Mass>({Umass(rng)})); 
    } 

    world.system<Acceleration>()
        .with<ParticleTag>()
        //.kind(flecs::PreUpdate)
        .each([&](Acceleration& a){
            a.ddx = 0.0; a.ddy = 0.0; 
            std::vector<double> vec = vector_distance(particles[0],particles[1]);
            std::cout<<"Vector displacement: "<<vec[0]<<","<<vec[1]<<std::endl;
        });
             

    world.progress(); 

}

