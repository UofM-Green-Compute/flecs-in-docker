#include <vector>

class Matrix{

private:
    // System Attributes
    const double l;     // Length
    const double n;     // Number of nodes
    const int peclet;   // Peclet Number
    const double rho;   // Density
    const double u;     // Speed
    const double gamma; // Diffusivity
    
    // Boundary Conditions
    double phi_start;
    double phi_end;

    // sparse matrix attributes
    std::vector<double> west;
    std::vector<double> diagonal;
    std::vector<double> east;

    // vector attributes
    std::vector<double> qvector;
    std::vector<double> phi;

public:
    // Declare Constructor
    Matrix(
        double new_l,
        double new_n,
        int new_peclet,
        double new_rho,
        double new_u,
        double new_gamma,
        double new_phi_start,
        double new_phi_end
    );

    // Declare Methods
    double east_function(double left, double centre, double right);
    double west_function(double left, double centre, double right);
    void matrix_filling();
    void matrix_solving();
};