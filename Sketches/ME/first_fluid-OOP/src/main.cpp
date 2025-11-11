/*
This code runs a simulation of a simple fluid which can be 
solved analytically.
*/

#include <iostream>
#include <vector>
#include <fstream> 
#include <cmath>
#include "classes.hpp"

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


    MyFile.close();
}
