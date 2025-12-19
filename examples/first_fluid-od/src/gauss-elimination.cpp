/*
Program to trial gauss elimination

1. Gauss elimination doesn't know what to do when one of the elements of the matrix is 0
    There may be an issue with my algorithm... there shouldn't be any zeroes on the diagonal and hence should be
    no /0 issues
    If the element we hit is 0, swap the order of two rows
    This article looks helpful: https://www.statlect.com/matrix-algebra/Gaussian-elimination
2. Passing "Matrix[j]" into row_manipulation twice? 
3. Double check that the gauss elimination, diagonals_to_matrix functions are rigorous 

*/

/* 

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

// Matrix struct (OD: Used in my first version of the code)
struct Matrix { std::vector<std::vector<double>> M; }; 
// Create structs to be used as components for the matrix 
struct West { std::vector<double> a; };
struct Diagonal { std::vector<double> b; };
struct East { std::vector<double> c; };
struct Qvector { std::vector<double> q; };
struct Conserved { std::vector<double> phi; };
struct MatrixTag {};

// 
double west_function(int left, int centre, int right)
{
    float A;
    A = - (Rho * u) / (right - left) - (2 * Gamma) / ((right - left) * (centre - left));
    return A;
}

// 
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

// Function to print elements of a vector - Oluwole
void print_vector(std::vector<double> vector)
{
    std::cout<<"---------------------------------"<<std::endl; 
    std::cout<<"{ ";
    for(int j = 0; j < vector.size(); j++) { std::cout<<vector[j]<<" "; }
    std::cout<<"}"<<std::endl;
    std::cout<<"---------------------------------"<<std::endl;
}

// Function to print elements of a matrix 
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

std::vector<std::vector<double>> augmented_matrix(std::vector<std::vector<double>> Matrix, std::vector<double> Vector)
{
    if(Matrix.size() == Vector.size())
    {
        for (int i = 0; i < Matrix.size(); i++)
        { Matrix[i].push_back(Vector[i]); }
    }
    else
    {
        std::cout<<"*** ERROR: Incompatible matrix and vector for augmented matrix ***"<<std::endl; 
    }
    return Matrix; 
}

// Function to perform gauss elimination operations on each row of the matrix 
std::vector<double> row_manipulation(std::vector<double> ith_row, std::vector<double> row_to_manipulate, int column)
{
    // ith_row: We manipulate the matrix such that entries below ith_row[column] = 0
    // row_to_manipulate: This is the row we are manipulating. At the end of the procedure row_to_manipulate[column] = 0.
    
    for (int i = 0; i < row_to_manipulate.size(); i++)
    {
        row_to_manipulate[i] = row_to_manipulate[i] - ith_row[i] * row_to_manipulate[column] / ith_row[column]; 
    }
    int size = row_to_manipulate.size(); 
    row_to_manipulate[size-1] = row_to_manipulate[size-1] - ith_row[size-1] * row_to_manipulate[column] / ith_row[column];

    return row_to_manipulate; 
}

std::vector<double> solve_equation(std::vector<std::vector<double>> Matrix)
{
    std::vector<double> phi = {}; 
    int no_rows = Matrix.size(); 
    int no_columns = Matrix[0].size();
     
    // nth case 
    phi.push_back(Matrix[no_rows-1][no_columns-1]); 

    for(int i = no_rows-2;  i > 0; i--)
    {
        double sum = 0;
        for(int k = i+1; k < no_columns-1; k++)
        { 
            std::cout<<k<<","<<phi[(no_rows-1)-k]<<std::endl;
            sum +=  Matrix[i][k] * phi[(no_rows-1)-k]; 
        } 
        std::cout<<sum<<std::endl;
        // std::cout<<Matrix[i][no_columns-1]<<std::endl;
        phi.push_back( (Matrix[i][no_columns-1] - sum) / Matrix[i][i]); 
    }
    return phi; 
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
            std::cout<<"Initial Conserved Quantity: "<<std::endl; 
            print_vector(conserved.phi); 

        });

    // System to Compute value of conserved quantity 
    world.system<West, Diagonal, East, Qvector, Conserved>()
        .with<MatrixTag>()
        .kind(flecs::OnUpdate)
        .each([&](West& west, Diagonal& diag, East& east, Qvector& Qvec, Conserved& conserved){

            std::vector<std::vector<double>> Matrix = diagonals_to_matrix(diag.b, east.c, west.a); 
            Matrix = augmented_matrix(Matrix, Qvec.q); 

            int no_rows = Matrix.size();
            int no_columns = Matrix[0].size(); 

            for (int i = 0; i < no_columns-1; i++)
            {
                for(int j = i+1; j < no_rows; j++)
                {
                    Matrix[j] = row_manipulation(Matrix[i], Matrix[j], i);
                }
            }

            // print_matrix(Matrix); 
            conserved.phi = solve_equation(Matrix); 
            std::cout<<"Solved Conserved Quantity: "<<std::endl; 
            print_vector(conserved.phi);  

            std::cout<<"Final (Augmented) Matrix: "<<std::endl; 
            print_matrix(Matrix);

        });
        
    flecs::system matrix_solver = world.system<Matrix>() 
    .each([&](Matrix &matrix)
    {
        // print_matrix(matrix.M);
    });

    world.progress(); 
    // initialise_matrix.run(); 
    // matrix_solver.run(); 
}

*/