/*
This code looks at the dynamical evolution of a compressible 
2D fluid in a box with two equal size chambers connected by a small gap
*/

#include <custom_phases_no_builtin.h>
#include <iostream>
#include <fstream> 
#include <vector>

// Number of nodes in x-axis. Also equal to the length in unit of node distance
const double nodeDistance = 1; // node distance in M

const int Nx = 100; // Number of from left wall to middle wall and from middle wall to right wall
const int Ny = 40; // Number of nodes from bottom wall to hole and hole to top wall
const int Nh = 20;  // Width of hole in middle wall in units of node distance
const double TIMESTEP = 0.00001;
const int NUMBERSTEPS = 100000;
const int FILENUMBER = 11;
const int SAVESTEP = NUMBERSTEPS/(FILENUMBER-1);
const double TIME = NUMBERSTEPS*TIMESTEP;

const double R = 8.314;  // Molar Gas Constant
const double T = 300;    // Temperature of air (room temperature)
const double M = 0.0290; // Molar mass of air

const double L = 2 * Nx; // Box length in units of node distance
const double W = 2 * Ny + Nh; // Box width in units of node distance

const double RhoLeft = 1.292;
const double RhoRight = 1; // Make density in right box less than left box

// Position Components
struct Position {int x, y; }; // Node position in units of node distance and offset by + 1/2

// Velocity components
struct VelocityStart { double x, y; }; // Node velocity in m
struct VelocityHalfPredict { double x, y; };
struct VelocityHalfCorrect { double x, y; };
struct VelocityEndPredict { double x, y; };

// Velocity components
struct DensityStart { double rho; }; // Node density in kg m-3
struct DensityHalfPredict { double rho; };
struct DensityHalfCorrect { double rho; };
struct DensityEndPredict { double rho; };

// Runge Kutta functions
struct FunctionsFirst { double u, v, rho; };
struct FunctionsSecond { double u, v, rho; };
struct FunctionsThird { double u, v, rho; };
struct FunctionsFourth { double u, v, rho; };

// Wall Tags
struct LowerWall {}; // This component indicates the node has a wall below it
struct UpperWall {}; // This component indicates the node has a wall above it 
struct LeftWall {};  // This component indicates the node has a wall left of it 
struct RightWall {}; // This component indicates the node has a wall right of it 

double u_function(double u, double u_right, double u_left, double u_up, 
    double u_down, double v, double rho, double rho_right, double rho_left)
{
    /*
    This function calculates the function used in runge-kutta approximations
    to find the value of u at a new time step
    */
    double f = - (((R*T/M)*((rho_right-rho_left)/rho))
             + (u*(u_right-u_left))
             + (v*(u_up-u_down)))
             /(2*nodeDistance);
    return f;
}

double v_function(double v, double v_right, double v_left, double v_up, 
    double v_down, double u, double rho, double rho_up, double rho_down)
{
    /*
    This function calculates the function used in runge-kutta approximations
    to find the value of v at a new time step
    */
    double f = - (((R*T/M)*((rho_up-rho_down)/rho))
              + (u*(v_right-v_left))
              + (v*(v_up-v_down)))
              / (2*nodeDistance);

    return f;
}

double rho_function(double rho_right, double rho_left, double rho_up, double rho_down,
    double v_up, double v_down, double u_right, double u_left)
{
    /*
    This function calculates the function used in runge-kutta approximations
    to find the value of rho at a new time step
    */
    double f = -((rho_right * u_right)
             - (rho_left * u_left)
             + (rho_up * v_up)
             - (rho_down * v_down))
             / (2*nodeDistance);
    return f;
}

int main(int argc, char *argv[]) {
    
    // Create specifactions file
    std::ofstream Specs;
    Specs.open("specs.txt");
    // Column titles
    Specs << "Box Length, Box Width, Hole Width, time step, number of files, Density 1, Density 2" << "\n";
    // Specs
    Specs << L*nodeDistance << ", " << W*nodeDistance << ", " << Nh*nodeDistance << ", "
          << SAVESTEP*TIMESTEP << ", " << FILENUMBER + 1 << ", "<< RhoLeft << ", " <<  RhoRight << "\n";
    Specs.close();

    // Create World
    flecs::world world(argc, argv);

    // Creates Phases which tell the program in which order to run the systems
    flecs::entity RungeKutta_1 = world.entity()
        .add(flecs::Phase); // This Phase calculates phihalf_predict
    
    flecs::entity RungeKutta_2 = world.entity()
        .add(flecs::Phase) // This Phase calculates phihalf_correct
        .depends_on(RungeKutta_1);
    
    flecs::entity RungeKutta_3 = world.entity()
        .add(flecs::Phase) // This Phase calculates phi_end_predict
        .depends_on(RungeKutta_2);

    flecs::entity RungeKutta_4 = world.entity()
        .add(flecs::Phase) // This Phase calculates phi_end_correct
        .depends_on(RungeKutta_3);
    
    flecs::entity Update = world.entity()
        .add(flecs::Phase) // This phase replaces phi_start with phi_end_correct
        .depends_on(RungeKutta_4);
    
    // Create components inside the world
    world.component<Position>();

    // Velocity Components
    world.component<VelocityStart>();
    world.component<VelocityHalfPredict>();
    world.component<VelocityHalfCorrect>();
    world.component<VelocityEndPredict>();

    // Density Components
    world.component<DensityStart>();
    world.component<DensityHalfPredict>();
    world.component<DensityHalfCorrect>();
    world.component<DensityEndPredict>();

    // Wall Tags
    struct LowerWall {}; // This component indicates the node has a wall below it
    struct UpperWall {}; // This component indicates the node has a wall above it 
    struct LeftWall {};  // This component indicates the node has a wall left of it 
    struct RightWall {}; // This component indicates the node has a wall right of it 

    // Runge-Kutta Components
    world.component<FunctionsFirst>();
    world.component<FunctionsSecond>();
    world.component<FunctionsThird>();
    world.component<FunctionsFourth>();

    std::vector<std::vector<flecs::entity>> nodes; // place to store nodes
    nodes.reserve(Nx*Ny); // Create the space
    
    // Create non-boundary Nodes
    for (int n = 0; n < 2*Nx; ++n) {
        nodes.push_back({});
        for (int m = 0; m < (2*Ny + Nh); ++m) {
            nodes[n].push_back(
                world.entity()
                    .set<Position>({n, m})
                    .set<VelocityStart>({0, 0})
                    .set<VelocityHalfPredict>({0, 0})
                    .set<VelocityHalfCorrect>({0, 0})
                    .set<VelocityEndPredict>({0, 0})
                    .set<DensityHalfPredict>({0})
                    .set<DensityHalfCorrect>({0})
                    .set<DensityEndPredict>({0})
                    .set<FunctionsFirst>({0, 0, 0})
                    .set<FunctionsSecond>({0, 0, 0})
                    .set<FunctionsThird>({0, 0, 0})
                    .set<FunctionsFourth>({0, 0, 0})
            );

            if (n < Nx) {
                nodes[n][m].set<DensityStart>({RhoLeft});
                if (n == 0) {
                    nodes[n][m].add<LeftWall>();
                }
                if (m == W - 1) {
                    nodes[n][m].add<UpperWall>();
                }
                if (m == 0) {
                    nodes[n][m].add<LowerWall>();
                }
                if (n == (Nx - 1) && m < Ny) {
                    nodes[n][m].add<RightWall>();
                }
                if (n == (Nx - 1) && m >= (Ny+Nh)) {
                    nodes[n][m].add<RightWall>();
                }
            } else {
                nodes[n][m].set<DensityStart>({RhoRight});
                if (n == L-1) {
                    nodes[n][m].add<RightWall>();
                }
                if (m == W - 1) {
                    nodes[n][m].add<UpperWall>();
                }
                if (m == 0) {
                    nodes[n][m].add<LowerWall>();
                }
                if (n == Nx && m < Ny) {
                    nodes[n][m].add<LeftWall>();
                }
                if (n == Nx && m >= (Ny + Nh)) {
                    nodes[n][m].add<LeftWall>();
                }
            }  
        }
    }
    
    // This system finds and updates VelocityHalfPredict and DensityHalfPredict
    world.system<Position, VelocityStart, VelocityHalfPredict, DensityStart, 
                 DensityHalfPredict, FunctionsFirst>()
        .kind(RungeKutta_1)
        .each([&](Position& pos, VelocityStart& velocityStart, VelocityHalfPredict& velocityHalf, 
                  DensityStart& densityStart, DensityHalfPredict& densityHalf, 
                  FunctionsFirst& function){
            
            // Initialise Neighbour Node Variables
            double horizontalUp; // These are the neighbour horizontal velocities
            double horizontalDown;
            double horizontalRight;
            double horizontalLeft;

            double verticalUp; // These are the neighbour vertical velocities
            double verticalDown;
            double verticalRight;
            double verticalLeft;

            double densityUp; // These are the neighbour densities
            double densityDown;
            double densityRight;
            double densityLeft;

            // Ensure upper wall boundary conditions if necessary
            if (nodes[pos.x][pos.y].has<UpperWall>()) {
                horizontalUp = velocityStart.x;
                verticalUp = -velocityStart.y;
                densityUp = densityStart.rho;
            } else {
                const VelocityStart& velUp = nodes[pos.x][pos.y+1].get<VelocityStart>();
                const DensityStart& denUp = nodes[pos.x][pos.y+1].get<DensityStart>();
                horizontalUp = velUp.x;
                verticalUp = velUp.y;
                densityUp = denUp.rho;
            }
            // Ensure lower wall boundary conditions if necessary
            if (nodes[pos.x][pos.y].has<LowerWall>()) {
                horizontalDown = velocityStart.x;
                verticalDown = -velocityStart.y;
                densityDown = densityStart.rho;
            } else {
                const VelocityStart& velDown = nodes[pos.x][pos.y-1].get<VelocityStart>();
                const DensityStart& denDown = nodes[pos.x][pos.y-1].get<DensityStart>();
                horizontalDown = velDown.x;
                verticalDown = velDown.y;
                densityDown = denDown.rho;
            }

            // Ensure right wall boundary conditions if necessary
            if (nodes[pos.x][pos.y].has<RightWall>()) {
                horizontalRight = -velocityStart.x;
                verticalRight = velocityStart.y;
                densityRight = densityStart.rho;
            } else {
                const VelocityStart& velRight = nodes[pos.x+1][pos.y].get<VelocityStart>();
                const DensityStart& denRight = nodes[pos.x+1][pos.y].get<DensityStart>();
                horizontalRight = velRight.x;
                verticalRight = velRight.y;
                densityRight = denRight.rho;
            }

            // Ensure left wall boundary conditions if necessary
            if (nodes[pos.x][pos.y].has<LeftWall>()) {
                horizontalLeft = -velocityStart.x;
                verticalLeft = velocityStart.y;
                densityLeft = densityStart.rho;
            } else {
                const VelocityStart& velLeft = nodes[pos.x-1][pos.y].get<VelocityStart>();
                const DensityStart& denLeft = nodes[pos.x-1][pos.y].get<DensityStart>();
                horizontalLeft = velLeft.x;
                verticalLeft = velLeft.y;
                densityLeft = denLeft.rho;
            }

            // Find first Runge Kutta functions
            function.u = u_function(velocityStart.x, horizontalRight, horizontalLeft, horizontalUp, 
            horizontalDown, velocityStart.y, densityStart.rho, densityRight, densityLeft);  

            function.v = v_function(velocityStart.y, verticalRight, verticalLeft, verticalUp, 
            verticalDown, velocityStart.x, densityStart.rho, densityUp, densityDown);

            function.rho = rho_function(densityRight, densityLeft, densityUp, densityDown, 
            verticalUp, verticalDown, horizontalRight, horizontalLeft);

            // Update half predictor values
            velocityHalf.x = velocityStart.x + (TIMESTEP * function.u)/2;
            velocityHalf.y = velocityStart.y + (TIMESTEP * function.v)/2;
            densityHalf.rho = densityStart.rho + (TIMESTEP * function.rho)/2;
        });
    
    // This system finds and updates VelocityHalfCorrect and DensityHalfCorrect
    world.system<Position, VelocityStart, VelocityHalfPredict, VelocityHalfCorrect, DensityStart, 
                 DensityHalfPredict, DensityHalfCorrect, FunctionsSecond>()
        .kind(RungeKutta_2)
        .each([&](Position& pos, VelocityStart& velocityStart, VelocityHalfPredict& velocityPredict, 
                  VelocityHalfCorrect& velocityCorrect, DensityStart& densityStart, 
                  DensityHalfPredict& densityPredict, DensityHalfCorrect& densityCorrect, 
                  FunctionsSecond& function){
            
            // Initialise Neighbour Node Variables
            double horizontalUp; // These are the neighbour horizontal velocities
            double horizontalDown;
            double horizontalRight;
            double horizontalLeft;

            double verticalUp; // These are the neighbour vertical velocities
            double verticalDown;
            double verticalRight;
            double verticalLeft;

            double densityUp; // These are the neighbour densities
            double densityDown;
            double densityRight;
            double densityLeft;

            // Ensure upper wall boundary conditions if necessary
            if (nodes[pos.x][pos.y].has<UpperWall>()) {
                horizontalUp = velocityPredict.x;
                verticalUp = -velocityPredict.y;
                densityUp = densityPredict.rho;
            } else {
                const VelocityHalfPredict& velUp = nodes[pos.x][pos.y+1].get<VelocityHalfPredict>();
                const DensityHalfPredict& denUp = nodes[pos.x][pos.y+1].get<DensityHalfPredict>();
                horizontalUp = velUp.x;
                verticalUp = velUp.y;
                densityUp = denUp.rho;
            }
            // Ensure lower wall boundary conditions if necessary
            if (nodes[pos.x][pos.y].has<LowerWall>()) {
                horizontalDown = velocityPredict.x;
                verticalDown = -velocityPredict.y;
                densityDown = densityPredict.rho;
            } else {
                const VelocityHalfPredict& velDown = nodes[pos.x][pos.y-1].get<VelocityHalfPredict>();
                const DensityHalfPredict& denDown = nodes[pos.x][pos.y-1].get<DensityHalfPredict>();
                horizontalDown = velDown.x;
                verticalDown = velDown.y;
                densityDown = denDown.rho;
            }

            // Ensure right wall boundary conditions if necessary
            if (nodes[pos.x][pos.y].has<RightWall>()) {
                horizontalRight = -velocityPredict.x;
                verticalRight = velocityPredict.y;
                densityRight = densityPredict.rho;
            } else {
                const VelocityHalfPredict& velRight = nodes[pos.x+1][pos.y].get<VelocityHalfPredict>();
                const DensityHalfPredict& denRight = nodes[pos.x+1][pos.y].get<DensityHalfPredict>();
                horizontalRight = velRight.x;
                verticalRight = velRight.y;
                densityRight = denRight.rho;
            }

            // Ensure left wall boundary conditions if necessary
            if (nodes[pos.x][pos.y].has<LeftWall>()) {
                horizontalLeft = -velocityPredict.x;
                verticalLeft = velocityPredict.y;
                densityLeft = densityPredict.rho;
            } else {
                const VelocityHalfPredict& velLeft = nodes[pos.x-1][pos.y].get<VelocityHalfPredict>();
                const DensityHalfPredict& denLeft = nodes[pos.x-1][pos.y].get<DensityHalfPredict>();
                horizontalLeft = velLeft.x;
                verticalLeft = velLeft.y;
                densityLeft = denLeft.rho;
            }

            // Find second Runge Kutta functions
            function.u = u_function(velocityPredict.x, horizontalRight, horizontalLeft, horizontalUp, 
            horizontalDown, velocityPredict.y, densityPredict.rho, densityRight, densityLeft);

            function.v = v_function(velocityPredict.y, verticalRight, verticalLeft, verticalUp, 
            verticalDown, velocityPredict.x, densityPredict.rho, densityUp, densityDown);

            function.rho = rho_function(densityRight, densityLeft, densityUp, densityDown, 
                verticalUp, verticalDown, horizontalRight, horizontalLeft);

            // Update half predictor values
            velocityCorrect.x = velocityStart.x + (TIMESTEP * function.u)/2;
            velocityCorrect.y = velocityStart.y + (TIMESTEP * function.v)/2;
            densityCorrect.rho = densityStart.rho + (TIMESTEP * function.rho)/2;
        });
    
    // This system finds and updates VelocityEndPredict and DensityEndPredict
    world.system<Position, VelocityStart, VelocityHalfCorrect, VelocityEndPredict, DensityStart, 
                 DensityHalfCorrect, DensityEndPredict, FunctionsFirst>()
        .kind(RungeKutta_3)
        .each([&](Position& pos, VelocityStart& velocityStart, VelocityHalfCorrect& velocityHalf, 
                  VelocityEndPredict& velocityEnd, DensityStart& densityStart, 
                  DensityHalfCorrect& densityHalf, DensityEndPredict& densityEnd,
                  FunctionsFirst& function){
            
            // Initialise Neighbour Node Variables
            double horizontalUp; // These are the neighbour horizontal velocities
            double horizontalDown;
            double horizontalRight;
            double horizontalLeft;

            double verticalUp; // These are the neighbour vertical velocities
            double verticalDown;
            double verticalRight;
            double verticalLeft;

            double densityUp; // These are the neighbour densities
            double densityDown;
            double densityRight;
            double densityLeft;

            // Ensure upper wall boundary conditions if necessary
            if (nodes[pos.x][pos.y].has<UpperWall>()) {
                horizontalUp = velocityHalf.x;
                verticalUp = -velocityHalf.y;
                densityUp = densityHalf.rho;
            } else {
                const VelocityHalfCorrect& velUp = nodes[pos.x][pos.y+1].get<VelocityHalfCorrect>();
                const DensityHalfCorrect& denUp = nodes[pos.x][pos.y+1].get<DensityHalfCorrect>();
                horizontalUp = velUp.x;
                verticalUp = velUp.y;
                densityUp = denUp.rho;
            }
            // Ensure lower wall boundary conditions if necessary
            if (nodes[pos.x][pos.y].has<LowerWall>()) {
                horizontalDown = velocityHalf.x;
                verticalDown = -velocityHalf.y;
                densityDown = densityHalf.rho;
            } else {
                const VelocityHalfCorrect& velDown = nodes[pos.x][pos.y-1].get<VelocityHalfCorrect>();
                const DensityHalfCorrect& denDown = nodes[pos.x][pos.y-1].get<DensityHalfCorrect>();
                horizontalDown = velDown.x;
                verticalDown = velDown.y;
                densityDown = denDown.rho;
            }

            // Ensure right wall boundary conditions if necessary
            if (nodes[pos.x][pos.y].has<RightWall>()) {
                horizontalRight = -velocityHalf.x;
                verticalRight = velocityHalf.y;
                densityRight = densityHalf.rho;
            } else {
                const VelocityHalfCorrect& velRight = nodes[pos.x+1][pos.y].get<VelocityHalfCorrect>();
                const DensityHalfCorrect& denRight = nodes[pos.x+1][pos.y].get<DensityHalfCorrect>();
                horizontalRight = velRight.x;
                verticalRight = velRight.y;
                densityRight = denRight.rho;
            }

            // Ensure left wall boundary conditions if necessary
            if (nodes[pos.x][pos.y].has<LeftWall>()) {
                horizontalLeft = -velocityHalf.x;
                verticalLeft = velocityHalf.y;
                densityLeft = densityHalf.rho;
            } else {
                const VelocityHalfCorrect& velLeft = nodes[pos.x-1][pos.y].get<VelocityHalfCorrect>();
                const DensityHalfCorrect& denLeft = nodes[pos.x-1][pos.y].get<DensityHalfCorrect>();
                horizontalLeft = velLeft.x;
                verticalLeft = velLeft.y;
                densityLeft = denLeft.rho;
            }

            // Find first Runge Kutta functions
            function.u = u_function(velocityHalf.x, horizontalRight, horizontalLeft, horizontalUp, 
            horizontalDown, velocityHalf.y, densityHalf.rho, densityRight, densityLeft);

            function.v = v_function(velocityHalf.y, verticalRight, verticalLeft, verticalUp, 
            verticalDown, velocityHalf.x, densityHalf.rho, densityUp, densityDown);

            function.rho = rho_function(densityRight, densityLeft, densityUp, densityDown, 
            verticalUp, verticalDown, horizontalRight, horizontalLeft);

            // Update half predictor values
            velocityEnd.x = velocityStart.x + (TIMESTEP * function.u);
            velocityEnd.y = velocityStart.y + (TIMESTEP * function.v);
            densityEnd.rho = densityStart.rho + (TIMESTEP * function.rho);
        });

    // This system finds and updates FunctionsFourth
    world.system<Position, VelocityEndPredict, DensityEndPredict, FunctionsFourth>()
        .kind(RungeKutta_4)
        .each([&](Position& pos, VelocityEndPredict& velocity, DensityEndPredict& density, 
                  FunctionsFourth& function){
            
            // Initialise Neighbour Node Variables
            double horizontalUp; // These are the neighbour horizontal velocities
            double horizontalDown;
            double horizontalRight;
            double horizontalLeft;

            double verticalUp; // These are the neighbour vertical velocities
            double verticalDown;
            double verticalRight;
            double verticalLeft;

            double densityUp; // These are the neighbour densities
            double densityDown;
            double densityRight;
            double densityLeft;

            // Ensure upper wall boundary conditions if necessary
            if (nodes[pos.x][pos.y].has<UpperWall>()) {
                horizontalUp = velocity.x;
                verticalUp = -velocity.y;
                densityUp = density.rho;
            } else {
                const VelocityEndPredict& velUp = nodes[pos.x][pos.y+1].get<VelocityEndPredict>();
                const DensityEndPredict& denUp = nodes[pos.x][pos.y+1].get<DensityEndPredict>();
                horizontalUp = velUp.x;
                verticalUp = velUp.y;
                densityUp = denUp.rho;
            }
            // Ensure lower wall boundary conditions if necessary
            if (nodes[pos.x][pos.y].has<LowerWall>()) {
                horizontalDown = velocity.x;
                verticalDown = -velocity.y;
                densityDown = density.rho;
            } else {
                const VelocityEndPredict& velDown = nodes[pos.x][pos.y-1].get<VelocityEndPredict>();
                const DensityEndPredict& denDown = nodes[pos.x][pos.y-1].get<DensityEndPredict>();
                horizontalDown = velDown.x;
                verticalDown = velDown.y;
                densityDown = denDown.rho;
            }

            // Ensure right wall boundary conditions if necessary
            if (nodes[pos.x][pos.y].has<RightWall>()) {
                horizontalRight = -velocity.x;
                verticalRight = velocity.y;
                densityRight = density.rho;
            } else {
                const VelocityEndPredict& velRight = nodes[pos.x+1][pos.y].get<VelocityEndPredict>();
                const DensityEndPredict& denRight = nodes[pos.x+1][pos.y].get<DensityEndPredict>();
                horizontalRight = velRight.x;
                verticalRight = velRight.y;
                densityRight = denRight.rho;
            }

            // Ensure left wall boundary conditions if necessary
            if (nodes[pos.x][pos.y].has<LeftWall>()) {
                horizontalLeft = -velocity.x;
                verticalLeft = velocity.y;
                densityLeft = density.rho;
            } else {
                const VelocityEndPredict& velLeft = nodes[pos.x-1][pos.y].get<VelocityEndPredict>();
                const DensityEndPredict& denLeft = nodes[pos.x-1][pos.y].get<DensityEndPredict>();
                horizontalLeft = velLeft.x;
                verticalLeft = velLeft.y;
                densityLeft = denLeft.rho;
            }

            // Find first Runge Kutta functions
            function.u = u_function(velocity.x, horizontalRight, horizontalLeft, horizontalUp, 
            horizontalDown, velocity.y, density.rho, densityRight, densityLeft);

            function.v = v_function(velocity.y, verticalRight, verticalLeft, verticalUp, 
            verticalDown, velocity.x, density.rho, densityUp, densityDown);

            function.rho = rho_function(densityRight, densityLeft, densityUp, densityDown, 
            verticalUp, verticalDown, horizontalRight, horizontalLeft);

        });

    // This updates VelocityStart and DensityStart
    world.system<VelocityStart, DensityStart, FunctionsFirst, FunctionsSecond, FunctionsThird, 
                 FunctionsFourth>()
        .kind(Update)
        .each([](VelocityStart& velocityStart, DensityStart& densityStart, FunctionsFirst& f1, 
                 FunctionsSecond& f2, FunctionsThird& f3, FunctionsFourth& f4){
            velocityStart.x = velocityStart.x + (TIMESTEP*(f1.u + 2*f2.u + 2*f3.u + f4.u))/6;
            velocityStart.y = velocityStart.y + (TIMESTEP*(f1.v + 2*f2.v + 2*f3.v + f4.v))/6;
            densityStart.rho = densityStart.rho + (TIMESTEP*(f1.rho +2*f2.rho + 2*f3.rho + f4.rho))/6;
        });

    // Run through systems every time step
    for (int t_step = 0; t_step <= NUMBERSTEPS; t_step++) {
        if (t_step != 0) {
                world.progress();
            }
        std::cout << t_step << "\n";
        
        // Saves Data to a .txt file
        if (t_step%SAVESTEP == 0) {
            // Prepare Save Files
            std::ofstream MyFile_u;
            std::ofstream MyFile_v;
            std::ofstream MyFile_rho;
            // Create filenames
            std::string horizontalName="horizontal_t=" + std::to_string(t_step*TIMESTEP) + ".txt";
            std::string verticalName="vertical_t=" + std::to_string(t_step*TIMESTEP) + ".txt";
            std::string densityName="density_t=" + std::to_string(t_step*TIMESTEP) + ".txt";
            // Create files
            MyFile_u.open(horizontalName);
            MyFile_v.open(verticalName);
            MyFile_rho.open(densityName);


            // Check Save Files are open/created
            if (!MyFile_u.is_open())
            {
                std::cout<<"Error in creating file horizontal_velocity.txt"<<std::endl; 
                return 1; 
            }

            if (!MyFile_v.is_open())
            {
                std::cout<<"Error in creating file vertical_velocity.txt"<<std::endl; 
                return 1; 
            }

            if (!MyFile_rho.is_open())
            {
                std::cout<<"Error in creating file density.txt"<<std::endl; 
                return 1; 
            }

            for (int y = 0; y < (2*Ny + Nh); y++) {
                for (int x = 0; x < 2*Nx; x++) {
                    const VelocityStart& velocity = nodes[x][y].get<VelocityStart>();
                    const DensityStart& density = nodes[x][y].get<DensityStart>();
                    if (x == 2*Nx-1) {
                        MyFile_u << velocity.x << "\n";
                        MyFile_v << velocity.y << "\n";
                        MyFile_rho << density.rho << "\n";
                    } else {
                        MyFile_u << velocity.x << ", ";
                        MyFile_v << velocity.y << ", ";
                        MyFile_rho << density.rho << ", ";
                    }
                }
            }
            // Close Save File
            MyFile_u.close();
            MyFile_v.close();
            MyFile_rho.close();
        }
    }   
}
