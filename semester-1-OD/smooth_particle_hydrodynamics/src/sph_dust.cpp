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
struct ParticleIndex {int i; }; 

// Constants
int NO_PARTICLES = 5; 

// Walls of the domain 
static constexpr double HALF_WIDTH = 40.0;           // half-width  (domain x in [-W, W])
static constexpr double HALF_HEIGHT = 20.0;           // half-height (domain y in [-H, H])
static constexpr double BOX_W = 2.0 * HALF_WIDTH;    // full width
static constexpr double BOX_H = 2.0 * HALF_HEIGHT;    // full height

// Mathematical constants
double CONST_H = 3; // Parameter determining size of domain around particle

// Physical constants
double ALPHA = 0; 
double BETA = 0; 
double K_B = 1.380649 * pow(10,-23); // Units: m^2 kg s^-2 K^-1



// Interpolant function (Gaussian)
double W(double r) { return ( 1 / CONST_H * sqrt(M_PI) ) * exp(( - pow(r,2) / pow(CONST_H,2) )); }

// std::vector<double> grad_W (double x) {}

// Function to calculate vector distance between particles 
std::vector<double> vector_distance(flecs::entity Particle, flecs::entity Particle_i){
    std::vector<double> displacement; 

    Position p = Particle.get<Position>();
    Position p_i = Particle_i.get<Position>();

    displacement.push_back(p.x - p_i.x);
    displacement.push_back(p.y - p_i.y);

    return displacement; 
}

double absolute_distance(std::vector<double> displacement){
    double sum = 0; 
    for (int i = 0; i < displacement.size(); i++) { sum += displacement[i]*displacement[i]; }
    return sqrt(sum); 
}

// Function to calculate temperature ???
int T = 1; 

// Function to caluclate density rho at the position of some particle
double density(std::vector<flecs::entity> Particles, flecs::entity particle){
    double Density = 0; 
    int particle_index = particle.get<ParticleIndex>().i; 

    for (int j=0; j<Particles.size()-1; j++)
    {
        std::vector<double> distance_vec = vector_distance(Particles[particle_index],Particles[j]); 
        double R = absolute_distance(distance_vec); 
        Density += Particles[j].get<Mass>().m * W(R) ; 
    }
    return Density; 
}

// Function to calculate pressure -- rho_i can be input into some equation of state to find pressure
// Pressure calculated using the van der Waals equation of state 
double vdw_pressure(std::vector<flecs::entity> Particles, flecs::entity particle){
    double mass = particle.get<Mass>().m; 
    double Density = density(Particles, particle); 

    double alpha_bar = ALPHA / mass; 
    double beta_bar = BETA / mass; 
    double kb_bar = K_B / mass; 

    return ( (Density * kb_bar * T) / (1 - beta_bar * Density) ) - alpha_bar * (Density*Density); 
}

// Function to calculate force on particle a due to all other particles       from to particle b
std::vector<double> force(std::vector<flecs::entity> Particles, flecs::entity Particle_a){
    std::vector<double> Force; 

    // Particle position
    double x_a = Particle_a.get<Position>().x; 
    double y_a = Particle_a.get<Position>().y;

    // Particle mass
    double m_a = Particle_a.get<Mass>().m; 

    // Calculate pressure at a,b
    double p_a = vdw_pressure(Particles, Particle_a);

    // Calculate density at a
    double rho_a = density(Particles, Particle_a);

    for(int i = 0; i < Particles.size(); i++)
    {
        double x_b = Particles[i].get<Position>().x;
        double y_b = Particles[i].get<Position>().y;

        // Particle distance
        std::vector<double> disp = vector_distance(Particle_a,Particles[i]); 
        double R = absolute_distance(disp);

        double m_b = Particles[i].get<Mass>().m;

        double p_b = vdw_pressure(Particles, Particles[i]);

        double rho_b = density(Particles, Particles[i]);

        // make into one constant
        double constant = (2 * m_a * m_b / pow(CONST_H,2) ) * ( p_b / (pow(rho_b,2)) + p_a / (pow(rho_a,2)));

        Force[0] += constant * x_a - x_b * W(R); 
        Force[1] += constant * y_a - y_b * W(R);
    } 

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
                .add<ParticleIndex>(i) 
                .set<Position>({UposX(rng), UposY(rng)}) 
                .set<Velocity>({Uvel(rng), Uvel(rng)}) 
                .set<Acceleration>({0.0, 0.0}) 
                .set<Mass>({Umass(rng)})); 
    } 

    world.system<Acceleration>()
        .with<ParticleIndex>()
        //.kind(flecs::PreUpdate)
        .each([&](Acceleration& a, Mass& m, ParticleIndex& index){

        std::vector<double> Force = force(particles, particles[index.i]);

        a.ddx = Force[0] / m.m;
        a.ddy = Force[1] / m.m; 

        std::cout<<"{"<<a.ddx<<","<<a.ddy<<"}"<<std::endl; 

        });
             
    //world.progress(); 

}

