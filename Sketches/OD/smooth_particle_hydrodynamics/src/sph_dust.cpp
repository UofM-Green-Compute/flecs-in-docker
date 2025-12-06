/*
Oluwole Delano | 13/11/25 
Simulating an ideal gas using smooth particle hydridynamics in ECS

** Questions ** :
- Should I incluide "self" particle in density sums? 
  When self particle is included in force sums it leads to weird behaviour (strange attractor)

To add:
- Classical limit of an ideal gas
- Energy checks: Total energy should be = 3/2*k_b*T*NO_PARTICLES 
- Model particle interactions: Lenard Jones potential --- do i need to do this? Does EOS do this for me?

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
int NO_PARTICLES = 50; 
const int STEPS = 1; // Number of time steps
double DT = 0.0001;    // Time step 
double PARTICLE_MASS = 1.6735575 * pow(10,-27); // Mass of Hydrogen

// Walls of the box ----
static constexpr double GX = 5;
static constexpr double GY = 5;
static constexpr double GX2 = 5;
double SLIT_WIDTH = 1; 

// Mathematical constants
double CONST_H = 0.1;                 // Parameter determining size of domain around particle
// Physical constants
int T = 10;                          // Temperature 
double ALPHA = 0;                     // Van-der-waals coefficient
double BETA = 0;                      // Van-der-waals coefficient
double K_B = 1.380649 * pow(10,-23);  // Units: m^2 kg s^-2 K^-1
double SIGMA = 10 / (7 * M_PI);       // Normalisation constant for kernel in two dimensions
int NO_DIMENSIONS = 2; 
double COURANT_NO = 1;                // Courant number

// Boltzman distribution
// Probability that a particle has velocity between v and v+dv = dv * maxwell_boltzmann()
double maxwell_boltzmann(double v)
{
    double norm = 4 * M_PI * pow((PARTICLE_MASS / (2*M_PI*K_B*T)),3.0/2.0) ;
    double exponential = std::exp( - (PARTICLE_MASS * v*v) / (2 * K_B * T) );
    return norm * v*v * exponential; 
}

// Metropolis-Hastings Algorithm
double metropolis_hastings(double x_i, double delta_i, double r)
{
    double x_trial = x_i + delta_i; 

    double w_i = maxwell_boltzmann(x_i); 
    double w_trial = maxwell_boltzmann(x_trial); 

    double alpha = std::min(1.0,(w_trial/w_i)); 

    if (r < alpha) { return x_trial; }
    else { return x_i; }
}

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
std::vector<double> vector_distance(std::vector<double> position_r, flecs::entity Particle_i){
    std::vector<double> displacement{0.0,0.0}; 

    Position p_i = Particle_i.get<Position>();

    displacement.push_back(position_r[0] - p_i.x);
    displacement.push_back(position_r[1] - p_i.y);

    return displacement; 
}

double absolute_distance(std::vector<double> displacement){
    double sum = 0; 
    for (int i = 0; i < displacement.size(); i++) { sum += displacement[i]*displacement[i]; }
    return sqrt(sum); 
}

// Function to caluclate density rho at the position of some particle
double density(std::vector<double> position_r, std::vector<flecs::entity> Particles){
    double Density = 0; 

    for (int j=0; j<Particles.size(); j++)
    {
        std::vector<double> distance_vec = vector_distance(position_r,Particles[j]); 
        double R = absolute_distance(distance_vec);  
        Density += Particles[j].get<Mass>().m * spline_W(R) ; 
    }
    return Density; 
}

// Gradient of Gaussian Interpolant Function
std::vector<double> grad_gaussian_W(flecs::entity p_a, flecs::entity p_b) 
{ 
    std::vector<double> pos_a; 
    pos_a.push_back(p_a.get<Position>().x); 
    pos_a.push_back(p_a.get<Position>().y); 

    std::vector<double> vec_distance = vector_distance(pos_a,p_b); 
    double R = absolute_distance(vec_distance); 
    std::vector<double> grad; 

    for(int i = 0; i < vec_distance.size(); i++) { grad.push_back( (2 / (CONST_H*CONST_H)) * vec_distance[i] * gaussian_W(R)); }

    return grad; 
} 

// Gradient of Spline Interpolant Function
std::vector<double> grad_spline_W(flecs::entity p_a, flecs::entity p_b)
{
    std::vector<double> pos_a; 
    pos_a.push_back(p_a.get<Position>().x); 
    pos_a.push_back(p_a.get<Position>().y); 

    std::vector<double> vec_distance = vector_distance(pos_a,p_b); 
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

    std::vector<double> particle_position; 
    particle_position.push_back(particle.get<Position>().x); 
    particle_position.push_back(particle.get<Position>().y); 

    double Density = density(particle_position,Particles); 

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
    std::vector<double> particle_a_position = {x_a,y_a};

    // Particle mass
    double m_a = Particle_a.get<Mass>().m; 

    // Calculate pressure at a
    double p_a = vdw_pressure(Particles, Particle_a);

    // Calculate density at a
    double rho_a = density(particle_a_position,Particles);

    for(int i = 0; i < Particles.size()-1; i++)
    {
        // Particle doesn't feel force from itself
        if(!(Particle_a.get<ParticleIndex>().i == Particles[i].get<ParticleIndex>().i))
        {
            double x_b = Particles[i].get<Position>().x;
            double y_b = Particles[i].get<Position>().y;
            std::vector<double> particle_b_position = {x_b,y_b}; 

            // Particle distance
            std::vector<double> disp = vector_distance(particle_a_position,Particles[i]); 
            double R = absolute_distance(disp);

            double m_b = Particles[i].get<Mass>().m;

            double p_b = vdw_pressure(Particles, Particles[i]);

            double rho_b = density(particle_b_position,Particles);

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

    // Open file for writing - Data file
    std::ofstream MyFile;
    MyFile.open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/smooth_particle_hydrodynamics/outputs/SPH_Dust.txt");
    if (!MyFile.is_open())
    {
        std::cout<<"Error in creating file"<<std::endl; 
        return 1; 
    }
    MyFile << "Particle 1 position x_1 (cm), y_1 (cm), Particle 1 velocity v_x, v_y |  x_2,y_2,v_x,v_y ; ..."<< std::endl;

    // Open file for writing - Specifications file
    std::ofstream MyFile_specs; 
    MyFile_specs.open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/smooth_particle_hydrodynamics/outputs/sph_code_specifications.txt");
    if (!MyFile_specs.is_open())
    {
        std::cout<<"Error in creating file"<<std::endl; 
        return 1; 
    }
    MyFile_specs << STEPS << " | " << NO_PARTICLES <<std::endl;
    MyFile_specs << T << " | " << PARTICLE_MASS <<std::endl;
    MyFile_specs << GX << " | " << GX2 << " | " << GY <<std::endl;
    MyFile_specs.close(); 

    // Open file for writing - Density field file
        std::ofstream MyFile_density;
        MyFile_density.open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/smooth_particle_hydrodynamics/outputs/density_field.txt");
        if (!MyFile_density.is_open())
        {
            std::cout<<"Error in opening file"<<std::endl; 
            return 1; 
        }

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
    std::uniform_real_distribution<double> Uspd(0, 5000);
    std::uniform_real_distribution<double> ZeroOne(0, 1);
    std::uniform_real_distribution<double> DELTA(-0.1, 0.1);
    // std::uniform_real_distribution<double> Umass(0.5, 2.0); 

    // Initialise entities | Generate random values for initial conditions
    std::vector<flecs::entity> particles; 
    particles.reserve(2*NO_PARTICLES); 

    int stationary_time = 1000; // Time taken for Markov Chain to reach stationary state after random initialisation
    int sample_interval = stationary_time; // 10000; // Number of time steps in bewteen taking samples of the distribution 

    std::vector<double> x_velocities; 
    x_velocities.reserve(NO_PARTICLES); 
    std::vector<double> y_velocities; 
    y_velocities.reserve(NO_PARTICLES);
    std::vector<std::vector<double>> density_matrix;
    density_matrix.reserve(GY*(GX+GX2)); 

    double v_x = Uspd(rng);
    double v_y = Uspd(rng); 
    int j = 0; 

    for(int i = 0; i < stationary_time + NO_PARTICLES*sample_interval + NO_PARTICLES; i++)
    {
        /* 
        double r_x = ZeroOne(rng); 
        double delta_x = DELTA(rng); 
        double v_x = metropolis_hastings(v_x,delta_x,r_x); 
        if(i>=stationary_time){ x_velocities.push_back(v_x); }

        double r_y = ZeroOne(rng); 
        double delta_y = DELTA(rng); 
        double v_y = metropolis_hastings(v_y,delta_y,r_y); 
        if(i>=stationary_time){ y_velocities.push_back(v_y); }
        */

        double r_x = ZeroOne(rng); 
        double delta_x = DELTA(rng); 
        double v_x = metropolis_hastings(v_x,delta_x,r_x); 

        double r_y = ZeroOne(rng); 
        double delta_y = DELTA(rng); 
        double v_y = metropolis_hastings(v_y,delta_y,r_y); 

        // std::cout<<v_x<<std::endl; 

        if (i >= stationary_time)
        {
            j+=1; 
            if(j=sample_interval){ x_velocities.push_back(v_x); y_velocities.push_back(v_y); j=0; }
                // std::cout<<v_x<<std::endl; j=0; }
        }
        
    }

    for (int i = 0; i < NO_PARTICLES; ++i) { 
        // Randomly make some of the velocities components negative
        double a = 1;  
        double rand = ZeroOne(rng); 
        if (rand > 0.5 ) { a = -1; }
        double b = 1;  
        double randb = ZeroOne(rng); 
        if (randb > 0.5 ) { b = -1; }

        particles.push_back( 
            world.entity() 
                .set<ParticleIndex>({i}) 
                .set<Position>({UposX(rng), UposY(rng)}) 
                .set<Velocity>({a * x_velocities[i], b * y_velocities[i]}) 
                .set<Acceleration>({0.0, 0.0}) 
                .set<Mass>({PARTICLE_MASS})
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
            if ( (cy < 0) || (cy >= GY) ) { v.dy = -v.dy; a.ddy = -a.ddy; } 

            if ( b.k == 0 && ((p.y < (GY/2 - SLIT_WIDTH/2)) || (p.y > (GY/2 + SLIT_WIDTH/2))) )
            { if ( (cx >= GX) )  { v.dx = -v.dx; a.ddx = -a.ddx;} }

            else if (b.k == 1 && !( (p.y > (GY/2 - SLIT_WIDTH/2)) && (p.y < (GY/2 + SLIT_WIDTH/2)) ) )
            { if ( (upx <= GX) )  { v.dx = -v.dx; a.ddx = -a.ddx;} }
            else 
            { if (p.x > GX) {b.k = 1; } else {b.k = 0; }}
        });

    // Check that time step isn't too large
    world.system<>()
        .kind(flecs::PreUpdate)
        .each([&](){

        double dt = 0;
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

        // if(DT > dt) { std::cout<<"ERROR: Time step too large. Decrease by "<<(DT-dt)<<std::endl; }
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

    // Calculate density
    world.system<>()
        .each([&](){

        for(int i = 0; i < (GX+GX2); i++)
        {
            for(int j = 0; j < GY; j++)
            {
                density_matrix[i].push_back(density({double(i),double(j)},particles)); 
                MyFile_density<<density_matrix[i][j]<<","; 
            }
            MyFile_density<<std::endl; 
        }
        });

    for (int i = 0; i < STEPS; ++i) {

        std::cout<<i<<std::endl; 

        world.progress();
        MyFile<<""<<std::endl; // End the line started in the write to file system

    }

    MyFile.close(); 

    t = clock() - t; 
    double time_taken = ((double)t) / CLOCKS_PER_SEC;
    std::cout<<"RUN TIME: "<<time_taken<<"s"<<std::endl;

}
