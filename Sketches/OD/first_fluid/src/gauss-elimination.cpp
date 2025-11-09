/*
Program to trial gauss elimination
*/

#include <algorithm>
#include <vector>
#include <cmath>
#include <flecs.h>
#include <systems.h>
#include <iostream> 

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

// Matrix struct (OD used in my first version of the code)
struct Matrix { std::vector<std::vector<double>> M; }; 

// Create structs to be used as components for the matrix 
struct West { std::vector<double> a; };
struct Diagonal { std::vector<double> b; };
struct East { std::vector<double> c; };
struct Qvector { std::vector<double> q; };
struct Conserved { std::vector<double> phi; };
struct MatrixTag {};

// - Mattias 
double west_function(int left, int centre, int right)
{
    float A;
    A = - (Rho * u) / (right - left) - (2 * Gamma) / ((right - left) * (centre - left));
    return A;
}

// - Mattias 
double east_function(int left, int centre, int right)
{
    float A;
    A = (Rho * u) / (right - left) - (2 * Gamma) / ((right - left) * (right - centre));
    return A;
}

// INPUT MATRICES MUST BE PLACED IN THE FOLLOWING ORDER {main diagonal, upper diagonal, lower diagonal}
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

// Function to perform gauss elimination operations on each row of the matrix - Oluwole
std::vector<double> row_manipulation(std::vector<double> ith_row, std::vector<double> first_row, int column)
{
    for (int i = column; i < ith_row.size(); i++)
    {
        ith_row[i] = ith_row[i] - first_row[i] * ith_row[column] / first_row[column]; 
    }
    return ith_row; 
}

// Function to print elements of a matrix - Oluwole 
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

// Function to print elements of a vector - Oluwole
void print_vector(std::vector<double> vector)
{
    std::cout<<"---------------------------------"<<std::endl; 
    std::cout<<"{";
    for(int j = 0; j < vector.size(); j++) { std::cout<<vector[j]; }
    std::cout<<"}"<<std::endl;
    std::cout<<"---------------------------------"<<std::endl;
}

int main() {

    // Create world 
    flecs::world world; 

    // Create components to be assigned to matrix entities
    world.component<West>();
    world.component<Diagonal>();
    world.component<East>();
    world.component<Qvector>();
    world.component<Conserved>();
    world.component<MatrixTag>();

    // Create matrix 
    std::vector<std::vector<double>> matrix = {{1,-3,4},{2,-5,6},{-3,3,4}}; 
    //print_matrix(matrix); 
    world.entity("Matrix 1")
        .set<Matrix>({matrix});
    

    // Create the Matrix as an entity
    world.entity("Matrix")
        .add<MatrixTag>()
        .set<West>({})
        .set<Diagonal>({})
        .set<East>({})
        .set<Qvector>({})
        .set<Conserved>({});

    // System to Initialise Matrix
    flecs::system initialise_matrix = world.system<West, Diagonal, East, Qvector, Conserved>()
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
            print_vector(conserved.phi); 

            std::cout<<"Hello"<<std::endl; 
            std::vector<std::vector<double>> Matrix = diagonals_to_matrix({1,1,1}, {2,2}, {3,3});
            print_matrix(Matrix);  
        });

    ///////////////// Seems that mattias code doesnt prodice a full matrix but
    ///////////////// only the elements of it needed to use his method. I need to make a full matrix

    // System to Compute value of conserved quantity inside system
    world.system<West, Diagonal, East, Qvector, Conserved>()
        .with<MatrixTag>()
        .kind(flecs::OnUpdate)
        .each([&](West& west, Diagonal& diag, East& east, Qvector& Qvec, Conserved& conserved){

            //std::vector<std::vector<double>> Matrix = diagonals_to_matrix(diag.b, east.c, west.a); 
            //print_matrix(Matrix); 

        });
        
    flecs::system matrix_solver = world.system<Matrix>() 
    .each([&](Matrix &matrix)
    {
        int no_rows = matrix.M.size();
        int no_columns = matrix.M[0].size(); 

        for (int i = 0; i < no_columns-1; i++)
        {
            for(int j = i+1; j < no_rows; j++)
            {
            matrix.M[j] = row_manipulation(matrix.M[j], matrix.M[j], i);
            }
        }
        // print_matrix(matrix.M);
    });

    //world.progress(); 
    initialise_matrix.run(); 
    matrix_solver.run(); 
}