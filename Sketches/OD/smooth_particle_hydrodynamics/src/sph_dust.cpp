/*
Oluwole Delano | 13/11/25 
Simulating an ideal gas using smooth particle hydridynamics in ECS

** Questions ** :
- Should I incluide "self" particle in density sums? 
  When self particle is included in force sums it leads to weird behaviour (strange attractor)

To add:
- Sample initial velocities from Maxwell Boltzman (Metropolis Algorithm Advanced Statistics)
- Model particle interactions: Lenard Jones potential --- do i need to do this? Does EOS do this for me?

Tools:
- Measure run time: https://www.geeksforgeeks.org/cpp/how-to-measure-time-taken-by-a-program-in-c/
*/

#include <iostream>
#include <flecs.h>
#include <vector>
#include <random>
#include <cmath>
#include <thread>
#include <chrono>
#include <fstream> 
#include <time.h> 
#include <algorithm>

// Initialise Components 
struct Position { double x, y; };
struct Velocity { double dx, dy; };
struct Acceleration { double ddx, ddy; };
struct Mass { double m; };
struct ParticleIndex {int i; }; 
struct Box {int k; }; 

// Constants
int NO_PARTICLES = 20; 
const int STEPS = 10000; // Number of time steps
double DT = 0.001;    // Time step
double dt = 0; 

// Walls of the box ----
static constexpr double GX = 5;
static constexpr double GY = 5;
static constexpr double GX2 = 5;
double SLIT_WIDTH = 1; 

// Mathematical constants
double CONST_H = 0.1;                 // Parameter determining size of domain around particle
// Physical constants
int T = 100;                          // Temperature 
double ALPHA = 0;                     // Van-der-waals coefficient
double BETA = 0;                      // Van-der-waals coefficient
double K_B = 1.380649 * pow(10,-23);  // Units: m^2 kg s^-2 K^-1
double SIGMA = 10 / (7 * M_PI);       // Normalisation constant for kernel in two dimensions
int NO_DIMENSIONS = 2; 
double COURANT_NO = 1;                // Courant number


// Interpolant function (Gaussian) -- Monaghan 1992 Equation 2.6
double gaussian_W(double r) { return ( 1 / CONST_H * sqrt(M_PI) ) * exp(( - pow(r,2) / pow(CONST_H,2) )); }

// Interpolant function 2 (Spline) -- Monaghan 1992 Section 7
double spline_W(double r) { 
    double q = r / CONST_H; 
    double c = SIGMA / pow(CONST_H,NO_DIMENSIONS); 

    if( (q >= 0) && (q <= 1) )
    {
        return c * (1 - ( (3/2) * pow(q,2) ) + ( (3/4) * pow(q,3) )); 
    }
    else if( (q >= 1) && (q <= 2) )
    {
        return c * ((1/4) * pow((2-q),3)); 
    }
    else {
        return 0; 
    }
}

// Function to calculate vector distance between particles 
std::vector<double> vector_distance(flecs::entity Particle, flecs::entity Particle_i){
    std::vector<double> displacement{0.0,0.0}; 

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

// Function to caluclate density rho at the position of some particle
double density(std::vector<flecs::entity> Particles, flecs::entity particle){
    double Density = 0; 
    int particle_index = particle.get<ParticleIndex>().i; 

    for (int j=0; j<Particles.size(); j++)
    {
        std::vector<double> distance_vec = vector_distance(Particles[particle_index],Particles[j]); 
        double R = absolute_distance(distance_vec);  
        Density += Particles[j].get<Mass>().m * spline_W(R) ; 
    }
    return Density; 
}

// Gradient of Gaussian Interpolant Function
std::vector<double> grad_gaussian_W(flecs::entity p_a, flecs::entity p_b) 
{ 
    std::vector<double> vec_distance = vector_distance(p_a,p_b); 
    double R = absolute_distance(vec_distance); 
    std::vector<double> grad; 

    for(int i = 0; i < vec_distance.size(); i++) { grad.push_back( (2 / (CONST_H*CONST_H)) * vec_distance[i] * gaussian_W(R)); }

    return grad; 
} 

// Gradient of Spline Interpolant Function
std::vector<double> grad_spline_W(flecs::entity p_a, flecs::entity p_b)
{
    std::vector<double> vec_distance = vector_distance(p_a,p_b); 
    double R = absolute_distance(vec_distance); 

    double q = R / CONST_H; 
    double c = SIGMA / pow(CONST_H,NO_DIMENSIONS); 

    std::vector<double> grad;

    if( (q >= 0) && (q <= 1) )
    {
        for(int i = 0; i < vec_distance.size(); i++) 
        { grad.push_back( (3 / (CONST_H*CONST_H)) * (((3*R) / (4*CONST_H)) - 1) * vec_distance[i] ); }
    }
    else if( (q >= 1) && (q <= 2) )
    {
        for(int i = 0; i < vec_distance.size(); i++) 
        { grad.push_back( (-3 / 4) * ((2-q)*(2-q)) * vec_distance[i] ); }
    }
    else {
        for(int i = 0; i < vec_distance.size(); i++) 
        { grad.push_back(0); }
    }

    return grad; 
}

// Equations of state -- rho_i can be input into some equation of state to find pressure
double vdw_pressure(std::vector<flecs::entity> Particles, flecs::entity particle){
    // Van der Waals equation of state to find pressure
    double mass = particle.get<Mass>().m; 
    double Density = density(Particles, particle); 

    double alpha_bar = ALPHA / mass; 
    double beta_bar = BETA / mass; 
    double kb_bar = K_B / mass; 

    return ( (Density * kb_bar * T) / (1 - beta_bar * Density) ) - alpha_bar * (Density*Density); 
}

// Function to calculate force on particle a due to all other particles      
std::vector<double> force(std::vector<flecs::entity> Particles, flecs::entity Particle_a){
    std::vector<double> Force = {0.0,0.0}; 

    // Particle position
    double x_a = Particle_a.get<Position>().x; 
    double y_a = Particle_a.get<Position>().y;

    // Particle mass
    double m_a = Particle_a.get<Mass>().m; 

    // Calculate pressure at a
    double p_a = vdw_pressure(Particles, Particle_a);

    // Calculate density at a
    double rho_a = density(Particles, Particle_a);

    for(int i = 0; i < Particles.size()-1; i++)
    {
        // Particle doesn't feel force from itself
        if(!(Particle_a.get<ParticleIndex>().i == Particles[i].get<ParticleIndex>().i))
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
            double constant = (m_a * m_b) * ( p_b / (pow(rho_b,2)) + p_a / (pow(rho_a,2)));

            std::vector<double> grad_W = grad_spline_W(Particle_a, Particles[i]); 

            Force[0] += constant * grad_W[0]; 
            Force[1] += constant * grad_W[1];
        }
    } 
    return Force; 
}

int main(int argc, char* argv[]) {

    // Start measuring run time of program
    clock_t t; 
    t = clock(); 

    // Open file for writing
    std::ofstream MyFile; 
    MyFile.open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/smooth_particle_hydrodynamics/outputs/SPH_Dust.txt");
    if (!MyFile.is_open())
    {
        std::cout<<"Error in creating file"<<std::endl; 
        return 1; 
    }
    MyFile << STEPS << " | " << NO_PARTICLES <<std::endl;
    MyFile << "Particle position x_1 (cm), Particle 1 position y_1 (cm), velocity_x, velocity_y |  x_2,y_2,v_x,v_y ; ..."<< std::endl;

    // Create the flecs world
    flecs::world world(argc,argv);
    // Components of the world
    world.component<Position>(); 
    world.component<Velocity>(); 
    world.component<Acceleration>(); 
    world.component<Mass>(); 
    world.component<Box>(); 

    // Tools for picking random numbers
    std::mt19937 rng( std::random_device{}()  ) ; // Initialise a random number generator with random device (for actual use)

    std::uniform_real_distribution<double> UposX(0, GX); 
    std::uniform_real_distribution<double> UposY(0, GY); 
    std::uniform_real_distribution<double> Uvel(-5, 5);
    std::uniform_real_distribution<double> Umass(0.5, 2.0); 

    // Initialise entities | Generate random values for initial conditions
    std::vector<flecs::entity> particles; 
    particles.reserve(NO_PARTICLES); 

    for (int i = 0; i < NO_PARTICLES; ++i) { 
        particles.push_back( 
            world.entity() 
                .set<ParticleIndex>({i}) 
                .set<Position>({UposX(rng), UposY(rng)}) 
                .set<Velocity>({Uvel(rng), Uvel(rng)}) 
                .set<Acceleration>({0.0, 0.0}) 
                .set<Mass>({Umass(rng)})
                .set<Box>({0})); 
    } 

    // Write to file
    world.system<Position, Velocity>()
        .kind(flecs::PreUpdate)
        .each([&](Position& p, Velocity &v){
            MyFile << p.x << "," << p.y << "," << v.dx << "," << v.dy << "|" ;
        });

    // Check for collision
    world.system<Position, Velocity, Acceleration, Box>()
        .kind(flecs::PreUpdate)
        .each([&](Position& p, Velocity& v, Acceleration& a, Box& b){

            int cx = int(std::floor(p.x)); // Round down to nearest integer
            int upx = int(std::ceil(p.x)); 
            int cy = int(std::floor(p.y));
            double x_tolerance = v.dx * DT; 

            // Reflect particle for edge cases where position is at a wall after rounding
            if ( ((cx < 0) || (cx >= GX+GX2)) ) { v.dx = -v.dx; a.ddx = -a.ddx;}

            if ( b.k == 0 && ((p.y < (GY/2 - SLIT_WIDTH/2)) || (p.y > (GY/2 + SLIT_WIDTH/2))) )
            { if ( (cx >= GX) )  { v.dx = -v.dx; a.ddx = -a.ddx;} }
            else if (b.k == 1 && !( (p.y > (GY/2 - SLIT_WIDTH/2)) && (p.y < (GY/2 + SLIT_WIDTH/2)) ) )
            { if ( (upx <= GX) )  { v.dx = -v.dx; a.ddx = -a.ddx;} }
            else 
            { if (p.x > GX) {b.k = 1; } }

            if ( (cy < 0) || (cy >= GY) ) { v.dy = -v.dy; a.ddy = -a.ddy; } 
            
        });

    // Check that time step isn't too large
    world.system<>()
        .kind(flecs::PreUpdate)
        .each([&](){

        double time_step = 0; 

        // Particle acceleration constraint | finding min( h / sqrt(F_i) )
        double t_f = 0;
        for(int i = 0; i < particles.size(); i++){
            std::vector<double> Force = force(particles, particles[i]);
            double abs_force = absolute_distance(Force); // Function to calculate |r| can be used similarly for F
            double t1 = CONST_H / sqrt(abs_force); 
            
            if (i == 0) { t_f = t1; }
            else if(t1<t_f) { t_f = t1; }
        }

        // CFL condition
        double t_cx = 0;
        double t_cy = 0;
        for(int i = 0; i < particles.size(); i++){
            double t2 = ( COURANT_NO * CONST_H ) / abs(particles[i].get<Velocity>().dx) ; 
            if (i == 0) { t_cx = t2; }
            else if(t2<t_cx) { t_cx = t2; }
        }
        for(int i = 0; i < particles.size(); i++){
            double t2 = ( COURANT_NO * CONST_H ) / abs(particles[i].get<Velocity>().dy) ; 
            if (i == 0) { t_cy = t2; }
            else if(t2<t_cy) { t_cy = t2; }
        }
        double t_c = std::min(t_cx,t_cy); 

        dt = 0.3 * std::min(t_f,t_c); 

        if(DT > dt) { std::cout<<"ERROR: Time step too large. Decrease by "<<(DT-dt)<<std::endl; }
        });
    

    // Update Acceleration
    world.system<Position, Velocity, Acceleration, Mass, ParticleIndex>()
        .each([&](Position& p, Velocity& v, Acceleration& a, Mass& m, ParticleIndex& index){
        std::vector<double> Force = force(particles, particles[index.i]);
        a.ddx = Force[0] / m.m;
        a.ddy = Force[1] / m.m; 
        });

    // Update velocity
    world.system<Velocity, Acceleration>()
        .each([&](Velocity& v, Acceleration& a){
            v.dx += a.ddx * DT;
            v.dy += a.ddy * DT;
        });

    // Update position
    world.system<Position, Velocity>()
        .each([&](Position& p, Velocity& v){
            p.x  += v.dx * DT;
            p.y  += v.dy * DT;
        });

    for (int i = 0; i < STEPS; ++i) {

        // std::cout<<i<<std::endl; 

        world.progress();
        MyFile<<""<<std::endl; // End the line started in the write to file system

    }

    MyFile.close(); 

    t = clock() - t; 
    double time_taken = ((double)t) / CLOCKS_PER_SEC; 
    std::cout<<"RUN TIME: "<<time_taken<<"s"<<std::endl; 

}