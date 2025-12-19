/*
This code runs a simulation of a simple fluid which can be 
solved analytically.
*/

#include <iostream>
#include <vector>
#include <fstream> 
#include <flecs.h>
#include <systems.h>

// Number of nodes in x-axis. Also equal to the length in unit of node distance
const int Nx = 41; // Number of nodes
const int Pe = 50; // Peclet Number
const double Rho = 2; // Density
const double u = 5; // Speed

// Calculates Diffusivity
double Gamma = (Rho * u * Nx) / Pe;

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

double west_function(int left, int centre, int right)
{
    float A;
    A = - (Rho * u) / (right - left) - (2 * Gamma) / ((right - left) * (centre - left));
    return A;
}

double east_function(int left, int centre, int right)
{
    float A;
    A = (Rho * u) / (right - left) - (2 * Gamma) / ((right - left) * (right - centre));
    return A;
}

std::vector<std::vector<double>> diagonals_to_matrix(std::vector<double> main_diag, std::vector<double> upper_diag, std::vector<double> lower_diag)
{
    std::vector<std::vector<double>> matrix; 

    for(int i = 0; i < main_diag.size(); i++)
    {
        matrix.push_back({});
        for(int j = 0; j < main_diag.size(); j++) 
        { 
            if (i == j){ matrix[i].push_back(main_diag[j]); }
            else if (i == j-1 && j > 0){ matrix[i].push_back(upper_diag[j-1]); }
            else if (i == j+1){ matrix[i].push_back(lower_diag[j]); }
            else {  matrix[i].push_back(0); }
        }
    }
    return matrix; 
}

void print_vector(std::vector<double> vector)
{
    std::cout<<"---------------------------------"<<std::endl; 
    std::cout<<"{ ";
    for(int j = 0; j < vector.size(); j++) { std::cout<<vector[j]<<" "; }
    std::cout<<"}"<<std::endl;
    std::cout<<"---------------------------------"<<std::endl;
}

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
                double A_W = west_function(i-1, i, i+1);
                double A_E = east_function(i-1, i, i+1);
                
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

            std::vector<std::vector<double>> Matrix1 = diagonals_to_matrix(diag.b, east.c, west.a);
            print_matrix(Matrix1); 
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
                    MyFile << i << ", " << conserved.phi[i] << std::endl; 
                }

            std::vector<std::vector<double>> Matrix = diagonals_to_matrix(diag.b, east.c, west.a);
            print_vector(conserved.phi); 
            print_matrix(Matrix); 

        });

    world.progress();
    MyFile.close();
}
