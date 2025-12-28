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

// Walls of the box ----
static constexpr double GX = 5;
static constexpr double GY = 5;
static constexpr double GX2 = 5;

double PARTICLE_MASS = 1.6735575 * pow(10,-27); // Mass of Hydrogen

// Constants
int NO_PARTICLES = 20; 
double SIGMA = 10 / (7 * M_PI);       // Normalisation constant for kernel in two dimensions
double NO_DIMENSIONS = 2.0; 

// Gambel error estimation constants
double NO_GAMBEL_DATA_POINTS = 20; // Number of data points stored for gambel error estimation graph
double NO_GAMBEL_PARTICLES = 200;  // Number of samples taken from gambel distribution  
double GAMBEL_RANGE = 50; // Range of Gambel distribution considered
double GAMBEL_MASS = (  exp(-exp(-GAMBEL_RANGE)) - exp(-(exp(GAMBEL_RANGE))) ) / NO_GAMBEL_PARTICLES ; // Integrate gambel density to find mass of each particle
double NO_GAMBEL_AVERAGES = 0; 

// Interpolant function (Gaussian) -- Monaghan 1992 Equation 2.6
double gaussian_W(double r, double h) { return ( 1 / (h * sqrt(M_PI)) ) * exp( - (pow(r,2.0) / pow(h,2.0)) ); }

// Interpolant function 2 (Spline) -- Monaghan 1992 Section 7
double spline_W(double r,double h) { 
    double q = r / h; 
    double c = SIGMA / pow(h,NO_DIMENSIONS); 

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
        Density += GAMBEL_MASS * gaussian_W(R,h); 
    }

    return Density; 
}

// Function to calculate the absolute magnitude of a vector
double absolute_distance(std::vector<double> displacement){
    double sum = 0; 
    for (int i = 0; i < displacement.size(); i++) { sum += displacement[i]*displacement[i]; }
    return sqrt(sum); 
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

    // Open file for writing - Specifications file
    std::ofstream MyFile_specs; 
    MyFile_specs.open("/Users/oluwoledelano/ECS_Development/flecs-in-docker/Sketches/OD/SPH_Runge/outputs/SPH_Runge_specifications.txt");
    if (!MyFile_specs.is_open())
    {
        std::cout<<"Error in creating file"<<std::endl; 
        return 1; 
    }
    MyFile_specs << "Number of particles: " << NO_PARTICLES <<std::endl;
    MyFile_specs << "Particle masses: " << PARTICLE_MASS <<std::endl;
    MyFile_specs << "Box 1 width: " << GX << " | " << "Box 2 width: " << GX2 << " | " << "Height of both boxes: " << GY << " | " << std::endl;
    MyFile_specs << "No of gambel density data points: " << NO_GAMBEL_DATA_POINTS << " | " << "No gambel averages: " << NO_GAMBEL_AVERAGES << std::endl;
    MyFile_specs.close(); 

    // Open file for writing - Gambel Density file
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

    // Tools for picking random numbers
    std::mt19937 rng( std::random_device{}()  ) ; // Initialise a random number generator with random device (for actual use)
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

    // *** ERROR ANALYSIS - TEST SMOTHING FUNCTIONS WITH A KNOWN DISTRIBUTION (Gambel density) *** 
    // Sample from Gambel density
    // Do this proces NO_GAMBEL_AVERAGES times so that we can average 
    for(int l = 0; l < NO_GAMBEL_AVERAGES; l++)
    {
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

    }

    t = clock() - t; 
    double time_taken = ((double)t) / CLOCKS_PER_SEC;
    std::cout<<"RUN TIME: "<<time_taken<<"s"<<std::endl;

}