/*
This code runs a simulation of a simple fluid which can be 
solved analytically.
*/

#include <iostream>
#include <vector>
#include <fstream> 
#include <flecs.h>
#include <systems.h>
#include <cmath>


// Number of nodes in x-axis. Also equal to the length in unit of node distance
const double L = 1; // Length
const double Nx = 41; // Number of nodes
const int Pe = 50; // Peclet Number
const double Rho = 2; // Density
const double u = 5; // Speed

// Calculates Diffusivity
double Gamma = (Rho * u * L) / Pe;

// Boundary Conditions
double Start = 0;
double End = 1;

/*Create structs to be used as components for the matrix*/
struct West { std::vector<double> a; };
struct Diagonal { std::vector<double> b; };
struct East { std::vector<double> c; };
struct Qvector { std::vector<double> q; };
struct Conserved { std::vector<double> phi; };
struct MatrixTag {};

double west_function(double left, double centre, double right)
{
    double A;
    A = - (Rho * u) / (right - left) - (2 * Gamma) / ((right - left) * (centre - left));
    return A;
}

double east_function(double left, double centre, double right)
{
    double A;
    A = (Rho * u) / (right - left) - (2 * Gamma) / ((right - left) * (centre - left));
    return A;
}

int main(int argc, char* argv[]) {
    // Prepare Save File
    std::ofstream MyFile;
    MyFile.open("starter_fluid.txt");
    if (!MyFile.is_open())
    {
        std::cout<<"Error in creating file"<<std::endl; 
        return 1; 
    }

    MyFile << "Node position, Conserved Quantity" << std::endl;

    // Create the World
    flecs::world world(argc, argv);

    // Create components to be assigned to matrix entities
    world.component<West>();
    world.component<Diagonal>();
    world.component<East>();
    world.component<Qvector>();
    world.component<Conserved>();
    world.component<MatrixTag>();

    // Create the Matrix as an entity
    world.entity("Matrix")
        .add<MatrixTag>()
        .set<West>({})
        .set<Diagonal>({})
        .set<East>({})
        .set<Qvector>({})
        .set<Conserved>({});

    // System to Initialise Matrix
    world.system<West, Diagonal, East, Qvector, Conserved>()
        .with<MatrixTag>()
        .kind(flecs::PreUpdate)
        .each([](West& west, Diagonal& diag, East& east, Qvector& Qvec, Conserved& conserved){
            conserved.phi.push_back(Start);
            for(int i = 1; i <= Nx-1; ++i) {
                double A_W = west_function((i-1)/Nx, i/Nx, (i+1)/Nx);
                double A_E = east_function((i-1)/Nx, i/Nx, (i+1)/Nx);
                // Updates Matrix Diagonal
                diag.b.push_back(-A_W-A_E);
                // Initialize value of conserved Quantity
                conserved.phi.push_back(0);
                // Updates Matrix Off-Diagonals and Q vector
                if (i == 1) {   
                    east.c.push_back(A_E);
                    Qvec.q.push_back(-A_W * Start);
                } else if (i == Nx-1) {
                    west.a.push_back(A_W);
                    Qvec.q.push_back(-A_E * End);
                } else {
                    west.a.push_back(A_W);
                    east.c.push_back(A_E);
                    Qvec.q.push_back(0);    
                }

            }
            conserved.phi.push_back(End);
        });
    
    // System to Compute value of conserved quantity inside system
    world.system<West, Diagonal, East, Qvector, Conserved>()
        .with<MatrixTag>()
        .kind(flecs::OnUpdate)
        .each([&](West& west, Diagonal& diag, East& east, Qvector& Qvec, Conserved& conserved){
            for(int i = 0; i <= Nx-2; ++i) {
                if (i == 0) {
                    east.c[i] = east.c[i] / diag.b[i];
                    Qvec.q[i] = Qvec.q[i] / diag.b[i];
                } else if (i == Nx - 2) {
                    Qvec.q[i] = (Qvec.q[i] - 
                        west.a[i-1]*Qvec.q[i-1]) / (diag.b[i] - west.a[i-1]*east.c[i-1]);
                } else {
                    east.c[i] = east.c[i] / (diag.b[i] - west.a[i-1]*east.c[i-1]);
                    Qvec.q[i] = (Qvec.q[i] - 
                        west.a[i-1]*Qvec.q[i-1]) / (diag.b[i] - west.a[i-1]*east.c[i-1]);
                }
            }
            for(int i = Nx-1; i >= 1; --i) {
                if (i == Nx - 1 ) {
                    conserved.phi[i] = Qvec.q[i-1];
                } else {
                    conserved.phi[i] = Qvec.q[i-1] - east.c[i-1] * conserved.phi[i+1];
                }
            }
            
            for (int i = 0; i <= Nx; ++i) {
                    MyFile << i/Nx << ", " << conserved.phi[i] << std::endl; 
                }
            
        });

    world.progress();
    MyFile.close();
}
