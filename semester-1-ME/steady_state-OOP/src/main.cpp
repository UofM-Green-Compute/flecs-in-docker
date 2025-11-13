/*
This code runs a simulation of a simple fluid which can be 
solved analytically.
*/

#include <iostream>
#include <vector>
#include <fstream> 
#include <cmath>
#include "classes.hpp"

const double L = 1;   // Length
const double N = 41;  // Number of nodes
const int PE = 50;    // Peclet Number
const double RHO = 2; // Density
const double U = 5;   // Speed

// Calculates Diffusivity
double GAMMA = (RHO * U * L) / PE;

// Boundary Conditions
double START = 0;
double END = 1;

int main() {
    // Prepare Save File
    std::ofstream MyFile;
    MyFile.open("starter_fluid.txt");
    if (!MyFile.is_open())
    {
        std::cout<<"Error in creating file"<<std::endl; 
        return 1; 
    }

    MyFile << "Node position, Conserved Quantity" << std::endl;

    Matrix steady(L, N, PE, RHO, U, GAMMA, START, END);
    steady.matrix_filling();
    steady.matrix_solving();
    for (int i = 0; i<=N; ++i) {
        MyFile << i*L/N << ", " << steady.get_phi()[i] << std::endl; 
    }
    MyFile.close();
}
