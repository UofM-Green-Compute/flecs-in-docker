#include "classes.hpp"
#include <string>
#include <iostream>
#include <vector>

// Constructor
Matrix::Matrix(
            double new_l,
            double new_n,
            int new_peclet,
            double new_rho,
            double new_u,
            double new_gamma,
            double new_phi_start,
            double new_phi_end
          )

  : l(new_l),
    n(new_n),
    peclet(new_peclet),
    rho(new_rho),
    u(new_u),
    gamma(new_gamma),
    phi_start(new_phi_start),
    phi_end(new_phi_end),
    west(), 
    diagonal(), 
    east(), 
    qvector(), 
    phi()

{}

// Methods

// Calculate west values
double Matrix::west_function(double left, double centre, double right){
  double A;
  A = - (rho * u) / (right - left) - (2 * gamma) / ((right - left) * (centre - left));
  return A;
}

// Calculate east values
double Matrix::east_function(double left, double centre, double right){
  double A;
  A = (rho * u) / (right - left) - (2 * gamma) / ((right - left) * (centre - left));
  return A;
}

// Fills in the matrix diagonals
void Matrix::matrix_filling(){
  phi.push_back(phi_start);
  for(int i = 1; i <= n-1; ++i) {
    double A_W = west_function((i-1)*l/n, i*l/n, (i+1)*l/n);
    double A_E = east_function((i-1)*l/n, i*l/n, (i+1)*l/n);
    // Updates Matrix Diagonal
    diagonal.push_back(-A_W-A_E);
    // Initialize value of conserved Quantity
    phi.push_back(0);
    // Updates Matrix Off-Diagonals and Q vector
    if (i == 1) {   
     east.push_back(A_E);
     qvector.push_back(-A_W * phi_start);
    } else if (i == n-1) {
      west.push_back(A_W);
      qvector.push_back(-A_E * phi_end);
    } else {
      west.push_back(A_W);
      east.push_back(A_E);
      qvector.push_back(0);    
    }
  }
  phi.push_back(phi_end);
}

// Solves the Matrix Equation
void Matrix::matrix_solving(){
  for(int i = 0; i <= n-2; ++i) {
      if (i == 0) {
          east[i] = east[i] / diagonal[i];
          qvector[i] = qvector[i] / diagonal[i];
      } else if (i == n - 2) {
          qvector[i] = (qvector[i] - 
              west[i-1]*qvector[i-1]) / (diagonal[i] - west[i-1]*east[i-1]);
      } else {
          east[i] = east[i] / (diagonal[i] - west[i-1]*east[i-1]);
          qvector[i] = (qvector[i] - 
              west[i-1]*qvector[i-1]) / (diagonal[i] - west[i-1]*east[i-1]);
      }
  }
  for(int i = n-1; i >= 1; --i) {
      if (i == n - 1 ) {
          phi[i] = qvector[i-1];
      } else {
          phi[i] = qvector[i-1] - east[i-1] * phi[i+1];
      }
  }
}