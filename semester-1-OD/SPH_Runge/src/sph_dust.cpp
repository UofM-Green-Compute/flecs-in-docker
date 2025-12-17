/*
Oluwole Delano | 13/11/25 
Simulating an ideal gas using smooth particle hydridynamics in ECS

** Questions / Next steps ** :
- Should I incluide "self" particle in density sums? 
  When self particle is included in force sums it leads to weird behaviour (strange attractor)
- Rescale units 
- Time step check doesnt seem to be working perfectly

** Possible extensions **
- Classical limit of an ideal gas
- Calculate total energy and check how well its conserved 
- Put in Leannard Jones potential by hand? 

** Sources **
Spline Kernel: Monaghan SPH Review article 1992
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
// Position component, plus start position at some time t_n
struct Position { double x, y; };
struct StartPosition { double x, y; };
// Velocity component, plus intermediate steps used in RK4
struct Velocity { double dx, dy; };
struct StartVelocity { double dx, dy; };
struct VelocityK1 { double dx, dy; };
struct VelocityK2 { double dx, dy; };
struct VelocityK3 { double dx, dy; };
// Acceleration component, plus intermediate steps used in RK4
struct Acceleration { double ddx, ddy; };
struct StartAcceleration { double ddx, ddy; };
struct AccelerationK1 { double ddx, ddy; };
struct AccelerationK2 { double ddx, ddy; };
struct AccelerationK3 { double ddx, ddy; };
// Other
struct Mass { double m; };
struct ParticleIndex {int i; }; 
struct Box {int k; }; 

// Walls of the box ----
static constexpr double GX = 5;
static constexpr double GY = 5;
static constexpr double GX2 = 5;
double SLIT_WIDTH = 1; 
                  
// Physical constants
int T = 10;                           // Temperature 
double ALPHA = 3.46 * pow(10,-3);     // Van-der-waals coefficient
double BETA = 23.71 * pow(10,-6);     // Van-der-waals coefficient
double K_B = 1.380649 * pow(10,-23);  // Units: m^2 kg s^-2 K^-1
double SIGMA = 10 / (7 * M_PI);       // Normalisation constant for kernel in two dimensions
double NO_DIMENSIONS = 2.0; 
double COURANT_NO = 1;                // Courant number

// Mathematical constants
double CONST_H = 0.5; // Parameter determining size of domain around particle 

// Constants
int NO_PARTICLES = 50; 
const int STEPS = 1; // Number of time steps
double DT = 0.00001;    // Time step 
double PARTICLE_MASS = 1.6735575 * pow(10,-27); // Mass of Hydrogen
double TEST_SCALE_FACTOR = 1 * pow(10,30); 

// Gambel error estimation constants
double NO_GAMBEL_DATA_POINTS = 20; // Number of data points stored for gambel error estimation graph
double NO_GAMBEL_PARTICLES = 50;  // Number of samples taken from gambel distribution  
double GAMBEL_RANGE = 50; // Range of Gambel distribution considered
double GAMBEL_MASS = (  exp(-exp(-GAMBEL_RANGE)) - exp(-(exp(GAMBEL_RANGE))) ) / NO_GAMBEL_PARTICLES ; // Integrate gambel density to find mass of each particle

// Interpolant function (Gaussian) -- Monaghan 1992 Equation 2.6
double gaussian_W(double r, double h) { return ( 1 / (h * sqrt(M_PI)) ) * exp( - (pow(r,2.0) / pow(h,2.0)) ); }

// Interpolant function 2 (Spline) -- Monaghan 1992 Section 7
double spline_W(double r,double h) { 
    double q = r / h; 
    double c = SIGMA / pow(CONST_H,NO_DIMENSIONS); 

    if( (q >= 0) && (q <= 1) )
    {
        return c * ( 1 - ((3/2) * pow(q,2)) + ((3/4) * pow(q,3)) ); 
    }
    else if( (q >= 1) && (q <= 2) )
    {
        return c * ((1/4) * pow((2-q),3)); 
    }
    else {
        return 0; 
    }
}

// Boltzman distribution
// Probability that a particle has velocity between v and v+d^3v = dv * maxwell_boltzmann()
double maxwell_boltzmann(double v)
{
    double norm = pow( (PARTICLE_MASS / (2*M_PI*K_B*T)) , (3.0/2.0) ) ;
    double exponential = std::exp( - (PARTICLE_MASS * v*v) / (2 * K_B * T) );
    return norm * exponential; 
}

// Metropolis-Hastings Algorithm -- Sampling from Maxwell-Boltzmann
double metropolis_hastings(double x_i, double delta_i, double r)
{
    double x_trial = x_i + delta_i; 

    double w_i = maxwell_boltzmann(x_i); 
    double w_trial = maxwell_boltzmann(x_trial); 

    double alpha = std::min(1.0,(w_trial/w_i)); 

    if (r < alpha) { return x_trial; }
    else { return x_i; }
}

double exact_gambel_density(double x){ return exp(-x - exp(-x)); }

// Metropolis-Hastings Algorithm -- Sampling from Gambell density
double metropolis_hastings_gambel(double x_i, double delta_i, double r)
{

    double x_trial = x_i + delta_i; 

    double w_i = exact_gambel_density(x_i); 
    double w_trial = exact_gambel_density(x_trial); 

    double alpha = std::min(1.0,(w_trial/w_i)); 

    if (r < alpha && (x_trial < GAMBEL_RANGE) && (x_trial > -GAMBEL_RANGE)) { return x_trial; }
    else { return x_i; }
}

// Function to caluclate density rho at some position
double density_gambel(double position_r, std::vector<double> configuration, double h){
    double Density = 0; 

    for (int j=0; j<configuration.size(); j++)
    {
        double R = position_r - configuration[j]; 
        Density += GAMBEL_MASS * spline_W(R,h); 
    }

    return Density; 
}

// Function to calculate vector distance between particles 
std::vector<double> vector_distance(std::vector<double> position_r, flecs::entity Particle_i){
    std::vector<double> displacement; 

    Position p_i = Particle_i.get<Position>();

    displacement.push_back(position_r[0] - p_i.x);
    displacement.push_back(position_r[1] - p_i.y);

    return displacement; 
}

// Function to calculate the absolute magnitude of a vector
double absolute_distance(std::vector<double> displacement){
    double sum = 0; 
    for (int i = 0; i < displacement.size(); i++) { sum += displacement[i]*displacement[i]; }
    return sqrt(sum); 
}

// Function to caluclate density rho at some position
double density(std::vector<double> position_r, std::vector<flecs::entity> Particles){
    double Density = 0; 

    for (int j=0; j<Particles.size(); j++)
    {
        std::vector<double> distance_vec = vector_distance(position_r,Particles[j]); 
        double R = absolute_distance(distance_vec);  
        Density += Particles[j].get<Mass>().m * gaussian_W(R,CONST_H); 
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
                                                                                                               // ! Spline XXXX 
    for(int i = 0; i < vec_distance.size(); i++) { grad.push_back( - (2 / (CONST_H*CONST_H)) * vec_distance[i] * gaussian_W(R,CONST_H)); }

    return grad; 
} 

// Gradient of Spline Interpolant Function
std::vector<double> grad_spline_W(flecs::entity p_a, flecs::entity p_b)
{
    // Get the vector position of particle a
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
        { grad.push_back( (-3 / (4 * R * CONST_H) ) * ((2-q)*(2-q)) * vec_distance[i] ); } 
    }
    else {
        for(int i = 0; i < vec_distance.size(); i++) 
        { grad.push_back(0); }
    }

    return grad; 
}

// Equations of state -- rho_i can be input into some equation of state to find pressure
// Van der Waals equation of state to find pressure -- "A review of SPH" equation 
double vdw_pressure(std::vector<flecs::entity> Particles, flecs::entity particle){
    
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

    // Particle a position
    double x_a = Particle_a.get<Position>().x; 
    double y_a = Particle_a.get<Position>().y;
    std::vector<double> particle_a_position = {x_a,y_a};

    // Particle a mass
    double m_a = Particle_a.get<Mass>().m; 

    // Pressure at particle a
    double p_a = vdw_pressure(Particles, Particle_a);

    // Density at particle a
    double rho_a = density(particle_a_position,Particles);

    for(int i = 0; i < Particles.size()-1; i++)
    {
        // Particle doesn't feel force from itself
        if(! (Particle_a.get<ParticleIndex>().i == Particles[i].get<ParticleIndex>().i) )
        {
            double x_b = Particles[i].get<Position>().x;
            double y_b = Particles[i].get<Position>().y;
            std::vector<double> particle_b_position = {x_b,y_b}; 

            // Vector distance r_a - r_i --> Absolute distance R
            std::vector<double> disp = vector_distance(particle_a_position,Particles[i]); 
            double R = absolute_distance(disp);

            // Particle i mass
            double m_b = Particles[i].get<Mass>().m;

            // Pressure at particle i
            double p_b = vdw_pressure(Particles, Particles[i]);

            // Desnity at particle i
            double rho_b = density(particle_b_position,Particles);

            // make into one constant
            double constant = - (m_a * m_b) * ( p_b / (pow(rho_b,2)) + p_a / (pow(rho_a,2)) ); 

            // Gradient of W(r_a - r_i)
            std::vector<double> grad_W = grad_gaussian_W(Particle_a, Particles[i]); 

            Force[0] += constant * grad_W[0]; 
            Force[1] += constant * grad_W[1];
            
        }
    } 
    return Force; 
}

// Function to print elements of a matrix - useful for debugging
void print_matrix(std::vector<std::vector<double>> matrix)
{
    std::cout<<"-------"<<std::endl; 
    for(int i = 0; i < matrix.size(); i++)
    {
        for(int j = 0; j < matrix[i].size(); j++) { std::cout<<matrix[i][j]<<"|"; }
        std::cout<<""<<std::endl; 
    }
    std::cout<<"-------"<<std::endl;
}

int main(int argc, char* argv[]) {

    // Start measuring run time of program
    clock_t t; 
    t = clock(); 

    // Open file for writing - Data file
    std::ofstream MyFile;
    MyFile.open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/SPH_Runge/outputs/SPH_Runge_Dust.txt");
    if (!MyFile.is_open())
    {
        std::cout<<"Error in creating file"<<std::endl; 
        return 1; 
    }
    MyFile << "Particle 1 position x_1 (cm), y_1 (cm), Particle 1 velocity v_x, v_y, Particle 1 Density |  x_2,y_2,v_x,v_y,rho_2 ; ..."<< std::endl;

    // Open file for writing - Specifications file
    std::ofstream MyFile_specs; 
    MyFile_specs.open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/SPH_Runge/outputs/SPH_Runge_specifications.txt");
    if (!MyFile_specs.is_open())
    {
        std::cout<<"Error in creating file"<<std::endl; 
        return 1; 
    }
    MyFile_specs << "Time step: " << DT << " | " << "Number of time steps: " << STEPS << " | " << "Number of particles: " << NO_PARTICLES <<std::endl;
    MyFile_specs << "System temperature: " << T << " | " << "Particle masses: " << PARTICLE_MASS <<std::endl;
    MyFile_specs << "Box 1 width: " << GX << " | " << "Box 2 width: " << GX2 << " | " << "Height of both boxes: " <<GY <<std::endl;
    MyFile_specs << "No of gambel density data points: " << NO_GAMBEL_DATA_POINTS <<std::endl;
    MyFile_specs.close(); 

    // Open file for writing - Density field file
    std::ofstream MyFile_density;
    MyFile_density.open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/SPH_Runge/outputs/density_field.txt");
    if (!MyFile_density.is_open())
    {
        std::cout<<"Error in opening density file"<<std::endl; 
        return 1; 
    }

    // Open file for writing - Gambel Density field file
    std::ofstream MyFile_gambel;
    std::string fileGambel="";
    std::string path = "/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/SPH_Runge/outputs/"; 
    if(NO_GAMBEL_PARTICLES == 20) { fileGambel=path+"Gambel_Density_ParticleNo=20.txt"; }
    else if(NO_GAMBEL_PARTICLES == 50) { fileGambel=path+"Gambel_Density_ParticleNo=50.txt"; }
    else if(NO_GAMBEL_PARTICLES == 75) { fileGambel=path+"Gambel_Density_ParticleNo=75.txt"; }
    else if(NO_GAMBEL_PARTICLES == 100) { fileGambel=path+"Gambel_Density_ParticleNo=100.txt"; }
    else if(NO_GAMBEL_PARTICLES == 200) { fileGambel=path+"Gambel_Density_ParticleNo=200.txt"; }
    else { std::cout<<"ERROR: No Gambel file matching this particle number"<<std::endl; }
    MyFile_gambel.open(fileGambel);
    if (!MyFile_gambel.is_open())
    {
        std::cout<<"Error in opening gambel density file"<<std::endl; 
        return 1; 
    }

    // Open file for writing - Particle number density file
    std::ofstream MyFile_NoDensity;
    std::string fileNoDensity="";
    if(NO_PARTICLES == 20) { fileNoDensity=path+"Particle_Number_Density_No=20.txt"; }
    else if(NO_PARTICLES == 50) { fileNoDensity=path+"Particle_Number_Density_No=50.txt"; }
    else if(NO_PARTICLES == 100) { fileNoDensity=path+"Particle_Number_Density_No=100.txt"; }
    else { std::cout<<"ERROR: No particle number density file matching this particle number"<<std::endl; }
    MyFile_NoDensity.open(fileNoDensity);
    if (!MyFile_NoDensity.is_open())
    {
        std::cout<<"Error in opening particle number density file"<<std::endl; 
        return 1; 
    }

    // Create the flecs world
    flecs::world world(argc,argv);

    // Components of the world
    world.component<Position>(); 
    world.component<StartPosition>(); 
    world.component<Velocity>(); 
    world.component<StartVelocity>();
    world.component<VelocityK1>();
    world.component<VelocityK2>(); 
    world.component<VelocityK3>();
    world.component<Acceleration>(); 
    world.component<StartAcceleration>();
    world.component<AccelerationK1>();
    world.component<AccelerationK2>();
    world.component<AccelerationK3>();
    world.component<Mass>(); 
    world.component<Box>(); 

    // Tools for picking random numbers
    std::mt19937 rng( std::random_device{}()  ) ; // Initialise a random number generator with random device (for actual use)
    std::uniform_real_distribution<double> UposX(GX, GX+GX2); 
    std::uniform_real_distribution<double> UposY(0, GY); 
    std::uniform_real_distribution<double> Uvel(-5000, 5000);
    std::uniform_real_distribution<double> Ugambel(-GAMBEL_RANGE, GAMBEL_RANGE);
    std::uniform_real_distribution<double> ZeroOne(0, 1);
    std::uniform_real_distribution<double> DELTA(-0.1, 0.1);
    std::uniform_real_distribution<double> DELTAGAMBEL(-0.1, 0.1);

    // Initialise entities | Generate random values for initial conditions
    std::vector<flecs::entity> particles; 
    particles.reserve(2*NO_PARTICLES); 

    int stationary_time = 10000000; // Give Markov Chain time to reach stationary state after random initialisation
    int sample_interval = 1000000;  // 10000; // Number of time steps in bewteen taking samples of the distribution 

    // Create vectors to store initial values of velocity components
    std::vector<double> x_velocities; 
    x_velocities.reserve(NO_PARTICLES); 
    std::vector<double> y_velocities; 
    y_velocities.reserve(NO_PARTICLES);

    // Create density matrix
    std::vector<std::vector<double>> density_matrix;
    density_matrix.reserve(GY*(GX+GX2)); 

    double v_x = Uvel(rng);
    double v_y = Uvel(rng); 
    int j = 0; 

    for(int i = 0; i < stationary_time + NO_PARTICLES*sample_interval; i++) // removed no_particles from i <
    {
        double r_x = ZeroOne(rng); 
        double delta_x = DELTA(rng); 
        double v_x = metropolis_hastings(v_x,delta_x,r_x); // remove double?? 

        double r_y = ZeroOne(rng); 
        double delta_y = DELTA(rng); 
        double v_y = metropolis_hastings(v_y,delta_y,r_y); 

        if (i >= stationary_time)
        {
            j+=1; 
            if(j==sample_interval) { x_velocities.push_back(v_x); y_velocities.push_back(v_y); j=0; }
        }
    }

    for (int i = 0; i < NO_PARTICLES; ++i) { 
        particles.push_back( 
            world.entity() 
                .set<ParticleIndex>({i}) 
                .set<Position>({UposX(rng), UposY(rng)}) 
                .set<StartPosition>({0.0, 0.0})  
                .set<Velocity>({x_velocities[i], y_velocities[i]})
                .set<StartVelocity>({0.0, 0.0})  
                .set<VelocityK1>({0.0, 0.0})
                .set<VelocityK2>({0.0, 0.0})
                .set<VelocityK3>({0.0, 0.0})
                .set<Acceleration>({0.0, 0.0})
                .set<StartAcceleration>({0.0, 0.0})  
                .set<AccelerationK1>({0.0, 0.0})
                .set<AccelerationK2>({0.0, 0.0})
                .set<AccelerationK3>({0.0, 0.0})
                .set<Mass>({PARTICLE_MASS})
                .set<Box>({1})); 
    } 

    // Note: Systems run in order they are coded in

    // Write to file
    world.system<Position, Velocity>()
        .kind(flecs::PreUpdate)
        .each([&](Position& p, Velocity &v){
            MyFile << p.x << "," << p.y << "," << v.dx << "," << v.dy << "," << density({p.x,p.y},particles) << "|" ;
        });

    // Check for collision
    world.system<Position, Velocity, Acceleration, Box>()
        .kind(flecs::PreUpdate)
        .each([&](Position& p, Velocity& v, Acceleration& a, Box& b){

            int cx = int(std::floor(p.x)); // Round down to nearest integer
            int upx = int(std::ceil(p.x)); 
            int cy = int(std::floor(p.y));
            double x_tolerance = v.dx * DT; 

            // Reflect particle when position is at a wall after rounding
            if ( ((cx < 0) || (cx >= GX+GX2)) ) { v.dx = -v.dx; }// a.ddx = -a.ddx;}
            if ( (cy < 0) || (cy >= GY) ) { v.dy = -v.dy; } // a.ddy = -a.ddy; } 

            if ( b.k == 0 && ((p.y < (GY/2 - SLIT_WIDTH/2)) || (p.y > (GY/2 + SLIT_WIDTH/2))) )
            { if ( (cx >= GX) )  { v.dx = -v.dx; } }// a.ddx = -a.ddx;} }

            else if (b.k == 1 && !( (p.y > (GY/2 - SLIT_WIDTH/2)) && (p.y < (GY/2 + SLIT_WIDTH/2)) ) )
            { if ( (upx <= GX) )  { v.dx = -v.dx; } }// a.ddx = -a.ddx;} }
            else 
            { if (p.x > GX) {b.k = 1; } else {b.k = 0; }}
        });

    // Check that time step isn't too large
    world.system<>()
        .kind(flecs::PreUpdate)
        .each([&](){

            double dt = 0;

            // Particle acceleration constraint | finding min( h / sqrt(F_i) )
            double t_f = 0;
            for(int i = 0; i < particles.size(); i++){
                std::vector<double> Force = force(particles, particles[i]);
                double abs_force = absolute_distance(Force); // Function to calculate |r| can be used similarly for F
                double t1 = CONST_H / sqrt(abs_force); 
                
                if (i == 0) { t_f = t1; }
                else if(t1<t_f) { t_f = t1; }
            }

            // CFL condition in the x and y directions
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

            // Find the minimum value of all these constraints
            double t_c = std::min(t_cx,t_cy); 

            dt = 0.3 * std::min(t_f,t_c); 

            if(DT > dt) { std::cout<<"ERROR: Time step too large. Decrease by "<<(DT-dt)<<std::endl; }
        });

    // RK1: Half predict Velocity, Position. Guess components at time t_{n+1/2}
    world.system<Position, StartPosition, Velocity, StartVelocity, Acceleration, StartAcceleration,
                VelocityK1, VelocityK2, VelocityK3, AccelerationK1, AccelerationK2, AccelerationK3>()
        .each([&](Position& p, StartPosition& pStart, Velocity& v, StartVelocity& vStart, Acceleration& a, StartAcceleration& aStart,
                VelocityK1& velk1, VelocityK2& velk2, VelocityK3& velk3, AccelerationK1& acck1, AccelerationK2& acck2, AccelerationK3& acck3)
        {
            // Components of the particle initially at time t_n
            pStart.x = p.x; 
            pStart.y = p.y;
            vStart.dx = v.dx;
            vStart.dy = v.dy; 
            aStart.ddx = a.ddx; 
            aStart.ddy = a.ddy;

            // Half predict position 
            p.x = pStart.x + (DT/2)*v.dx ;
            p.y = pStart.y + (DT/2)*v.dy ;

            // Half predict velocity
            v.dx = vStart.dx + (DT/2)*a.ddx ;
            v.dy = vStart.dy + (DT/2)*a.ddy ;

            // Store half predicted velocity
            velk1.dx = v.dx;
            velk1.dy = v.dy;
        });

    // RK1: Half predict acceleration - calculate a(t_{n+1/2}) based on the p,v just predicted
    world.system<Acceleration, Mass, ParticleIndex, AccelerationK1>()
        .each([&](Acceleration& a,Mass& m, ParticleIndex& index, AccelerationK1& acck1)
        {
            std::vector<double> Force = force(particles, particles[index.i]);

            a.ddx = Force[0] / m.m;
            a.ddy = Force[1] / m.m; 

            // Store half predicted acceleration
            acck1.ddx = a.ddx; 
            acck1.ddy = a.ddy;             
        });

    // RK2: Half correct position and velocity
    world.system<Position, StartPosition, Velocity, StartVelocity, Acceleration, StartAcceleration, VelocityK2>()
        .each([&](Position& p, StartPosition& pStart, Velocity& v, StartVelocity& vStart, Acceleration& a, StartAcceleration& aStart, VelocityK2& velk2)
        {
            // Half correct position 
            p.x = pStart.x + (DT/2)*v.dx ;
            p.y = pStart.y + (DT/2)*v.dy ;

            // Half correct velocity
            v.dx = vStart.dx + (DT/2)*a.ddx ;
            v.dy = vStart.dy + (DT/2)*a.ddy ;

            // Store half corrected velocity
            velk2.dx = v.dx;
            velk2.dy = v.dy;
        });
    
    // RK2: Half correct acceleration - correct a(t_{n+1/2}) based on the p,v just predicted
    world.system<Acceleration, Mass, ParticleIndex, AccelerationK2>()
        .each([&](Acceleration& a,Mass& m, ParticleIndex& index, AccelerationK2& acck2)
        {
            std::vector<double> Force = force(particles, particles[index.i]);

            a.ddx = Force[0] / m.m;
            a.ddy = Force[1] / m.m; 

            // Store half corrected acceleration
            acck2.ddx = a.ddx; 
            acck2.ddy = a.ddy;             
        });
    
    // RK3: Full predict position, velocity
    world.system<Position, StartPosition, Velocity, StartVelocity, Acceleration, StartAcceleration, VelocityK3>()
        .each([&](Position& p, StartPosition& pStart, Velocity& v, StartVelocity& vStart, Acceleration& a, StartAcceleration& aStart, VelocityK3& velk3)
        {
            // Full predict position 
            p.x = pStart.x + DT*v.dx ;
            p.y = pStart.y + DT*v.dy ;

            // Full predict velocity
            v.dx = vStart.dx + DT*a.ddx ;
            v.dy = vStart.dy + DT*a.ddy ;

            // Store full predicted velocity
            velk3.dx = v.dx;
            velk3.dy = v.dy;
        });

    // RK3: Full predict acceleration - predict a(t_{n+1) based on the p,v just predicted
    world.system<Acceleration, Mass, ParticleIndex, AccelerationK3>()
        .each([&](Acceleration& a,Mass& m, ParticleIndex& index, AccelerationK3& acck3)
        {
            std::vector<double> Force = force(particles, particles[index.i]);

            a.ddx = Force[0] / m.m;
            a.ddy = Force[1] / m.m; 

            // Store full predicted acceleration
            acck3.ddx = a.ddx; 
            acck3.ddy = a.ddy;             
        });
    
    // RK4: End correct position, velocity - final answer
    world.system<Position, StartPosition, Velocity, StartVelocity, Acceleration, StartAcceleration,
                VelocityK1, VelocityK2, VelocityK3, AccelerationK1, AccelerationK2, AccelerationK3>()
        .each([&](Position& p, StartPosition& pStart, Velocity& v, StartVelocity& vStart, Acceleration& a, StartAcceleration& aStart,
                    VelocityK1& velk1, VelocityK2& velk2, VelocityK3& velk3, AccelerationK1& acck1, AccelerationK2& acck2, AccelerationK3& acck3)
        {

            // Calculate position
            p.x = pStart.x + (DT/6) * ( velk1.dx + 2*velk2.dx + 2*velk3.dx + v.dx ) ;
            p.y = pStart.y + (DT/6) * ( velk1.dy + 2*velk2.dy + 2*velk3.dy + v.dy ) ;

            // Calculate velocity
            v.dx = vStart.dx + (DT/6) * ( acck1.ddx + 2*acck2.ddx + 2*acck3.ddx + a.ddx ) ;
            v.dy = vStart.dy + (DT/6) * ( acck1.ddy + 2*acck2.ddy + 2*acck3.ddy + a.ddy ) ;
        });
    
    // RK4: End correct acceleration (Final answer)
    world.system<Position, Velocity, Acceleration, Mass, ParticleIndex>()
        .each([&](Position& p, Velocity& v, Acceleration& a, Mass& m, ParticleIndex& index){
            std::vector<double> Force = force(particles, particles[index.i]);

            a.ddx = Force[0] / m.m;
            a.ddy = Force[1] / m.m; 

        });

    // *** ERROR ANALYSIS - TEST SMOTHING FUNCTIONS WITH A KNOWN DISTRIBUTION (Gambel density) *** 
    // Sample from Gambel density
    double p_x = Ugambel(rng);
    std::vector<double> particle_config; 
    particle_config.reserve(NO_GAMBEL_PARTICLES); 
    int k = 0; 

    for(int i = 0; i < stationary_time + NO_GAMBEL_PARTICLES*sample_interval; i++) 
    {
        double rx = ZeroOne(rng); 
        double deltax = DELTAGAMBEL(rng);

        p_x = metropolis_hastings_gambel(p_x,deltax,rx); 

        if (i >= stationary_time)
        {
            k+=1; 
            if(k==sample_interval) { particle_config.push_back(p_x); k=0; }
        }
    }

    // Calculate density from Gambel  
    double h = 0.1; 
    for(int i = 0; i < NO_GAMBEL_DATA_POINTS; i++){
        double sum = 0; 
        double avg = 0; 

        for(int j = 0; j < NO_GAMBEL_PARTICLES; j++)
        { sum += abs( exact_gambel_density(particle_config[j]) - density_gambel(particle_config[j],particle_config,h) ); }

        avg = sum / NO_GAMBEL_PARTICLES; 

        MyFile_gambel<<h<<","<<avg<<std::endl;
        h += 0.1; 
    }

    for (int i = 0; i < STEPS; ++i) {

        // Print current step
        std::cout<<i<<std::endl; 

        world.progress();
        MyFile<<std::endl; // End the line started in the write to file system

        // Calculate density grid at each time step
        density_matrix = {{}}; 
        for(int i = GY; i >= 0; i--)
        {
            for(int j = 0; j < (GX+GX2)+1; j++)
            {
                density_matrix[GY-i].push_back(density({double(j),double(i)},particles)); 
                if(j < (GX+GX2)) { MyFile_density<<density_matrix[GY-i][j]<<","; }
                else { MyFile_density<<density_matrix[GY-i][j]<<std::endl; } 
            } 
        }

        // Calculate the number of particles per box at each step
        int box1 = 0; 
        int box2 = 0; 
        for(int i = 0; i < particles.size(); i++)
        {
            if(particles[i].get<Box>().k == 0){box1 += 1;}
            else if(particles[i].get<Box>().k == 1){box2 += 1;}
            else{std::cout<<"What the hell is happening"<<std::endl;}
        }
        MyFile_NoDensity<<box1<<","<<box2<<std::endl; 
    }

    MyFile.close(); 

    t = clock() - t; 
    double time_taken = ((double)t) / CLOCKS_PER_SEC;
    std::cout<<"RUN TIME: "<<time_taken<<"s"<<std::endl;

}