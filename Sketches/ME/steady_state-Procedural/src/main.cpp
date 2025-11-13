/*
This code runs a simulation of a simple fluid which can be 
solved analytically.
*/

#include <iostream>
#include <vector>
#include <fstream> 
#include <cmath>

const double L = 1;   // Length
const double N = 11;  // Number of nodes
const int PE = 50;    // Peclet Number
const double RHO = 2; // Density
const double U = 5;   // Speed

// Calculates Diffusivity
double GAMMA = (RHO * U * L) / PE;

// Boundary Conditions
double START = 0;
double END = 1;

double west_function(double left, double centre, double right){
  double A;
  A = - (RHO * U) / (right - left) - (2 * GAMMA) / ((right - left) * (centre - left));
  return A;
}

// Calculate east values
double east_function(double left, double centre, double right){
  double A;
  A = (RHO * U) / (right - left) - (2 * GAMMA) / ((right - left) * (centre - left));
  return A;
}

// Fills in the matrix diagonals
void matrix_filling(
    std::vector<double>& phi,
    std::vector<double>& west,
    std::vector<double>& diagonal,
    std::vector<double>& east,
    std::vector<double>& Q
){
    phi.push_back(START);
    for(int i = 1; i <= N-1; ++i) {
    double A_W = west_function((i-1)*L/N, i*L/N, (i+1)*L/N);
    double A_E = east_function((i-1)*L/N, i*L/N, (i+1)*L/N);
    // Updates Matrix Diagonal
    diagonal.push_back(-A_W-A_E);
    // Initialize value of conserved Quantity
    phi.push_back(0);
    // Updates Matrix Off-Diagonals and Q vector
    if (i == 1) {   
        east.push_back(A_E);
        Q.push_back(-A_W * START);
    } else if (i == N-1) {
        west.push_back(A_W);
        Q.push_back(-A_E * END);
    } else {
        west.push_back(A_W);
        east.push_back(A_E);
        Q.push_back(0);    
    }
    }
    phi.push_back(END);
}

// Solves the Matrix Equation
void matrix_solving(
    std::vector<double>& phi,
    std::vector<double>& west,
    std::vector<double>& diagonal,
    std::vector<double>& east,
    std::vector<double>& Q
){
    for(int i = 0; i <= N-2; ++i) {
        if (i == 0) {
            east[i] = east[i] / diagonal[i];
            Q[i] = Q[i] / diagonal[i];
        } else if (i == N - 2) {
            Q[i] = (Q[i] - 
                west[i-1]*Q[i-1]) / (diagonal[i] - west[i-1]*east[i-1]);
        } else {
            east[i] = east[i] / (diagonal[i] - west[i-1]*east[i-1]);
            Q[i] = (Q[i] - 
                west[i-1]*Q[i-1]) / (diagonal[i] - west[i-1]*east[i-1]);
        }
    }
    for(int i = N - 1; i >= 1; --i) {
        if (i == N - 1 ) {
            phi[i] = Q[i-1];
        } else {
            phi[i] = Q[i-1] - east[i-1] * phi[i+1];
        }
    }
}

int main() {
    // Prepare matrices and vectors
    // sparse matrix
    std::vector<double> west_steady;
    std::vector<double> diagonal_steady;
    std::vector<double> east_steady;

    // vectors
    std::vector<double> Q_steady;
    std::vector<double> phi_steady;

    // Prepare Save File
    std::ofstream MyFile;
    MyFile.open("starter_fluid.txt");
    if (!MyFile.is_open())
    {
        std::cout<<"Error in creating file"<<std::endl; 
        return 1; 
    }

    MyFile << "Node position, Conserved Quantity" << std::endl;

    matrix_filling(
        phi_steady,
        west_steady, 
        diagonal_steady, 
        east_steady, 
        Q_steady);

    matrix_solving(
        phi_steady,
        west_steady, 
        diagonal_steady, 
        east_steady, 
        Q_steady);
    
    for (int i = 0; i<=N; ++i) {
        MyFile << i*L/N << ", " << phi_steady[i] << std::endl; 
    }
    MyFile.close();
}
