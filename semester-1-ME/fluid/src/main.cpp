/*
This code looks at the dynamical evolution of a compressible 
2D fluid in a box with two equal size chambers connected by a small gap
*/

#include <custom_phases_no_builtin.h>
#include <iostream>
#include <fstream> 
#include <vector>

// Number of nodes in x-axis. Also equal to the length in unit of node distance
const double nodeDistance = 1e-5; // node distance in SI units

const double W = 1;
const int Nx = 500; // Number of from left wall to middle wall and from middle wall to right wall
const int Ny = 500; // Number of nodes from bottom wall to top wall
const int Nh = 20;  // Width of hole in middle wall in units of node distance

const double R = 8.314;  // Molar Gas Constant
const double T = 300;    // Temperature of air (room temperature)
const double M = 0.0290; // Molar mass of air

const double L = 2 * Nx; // Box length in units of node distance
const double W_b = Ny;   // Box width in units of node distance
const double W_h = Nh;   // Hole width in units of node distance

// Position Components
struct Position {int x, y; }; // Node position in units of node distance

// Velocity components
struct velocityStart { double x, y; }; // Node velocity in m
struct velocityHalfPredict { double x, y; };
struct velocityHalfCorrect { double x, y; };
struct velocityEndPredict { double x, y; };
struct velocityEndCorrect { double x, y; };

// Velocity components
struct densityStart { double rho; }; // Node density in kg m-3
struct densityHalfPredict { double rho; };
struct densityHalfCorrect { double rho; };
struct densityEndPredict { double rho; };
struct densityEnd { double rho; };

// Wall Tags
struct LowerWall {}; // This component indicates the node has a wall below it
struct UpperWall {}; // This component indicates the node has a wall above it 
struct LeftWall {};  // This component indicates the node has a wall left of it 
struct RightWall {}; // This component indicates the node has a wall right of it 

double u_function(double u, double u_right, double u_left, double u_up, 
    double u_down, double v, double rho, double rho_left, double rho_right)
{
    /*
    This function calculates the function used in runge-kutta approximations
    to find the value of u at a new time step
    */
    double f = -(1/(2*nodeDistance)) 
            * (
              ((R*T/M)*((rho_right-rho_left)/rho))
            + (u*(u_right-u_left))
            + (v*(u_up-u_down))
            );

    return f;
}

double v_function(double v, double v_right, double v_left, double v_up, 
    double v_down, double u, double rho, double rho_up, double rho_down)
{
    /*
    This function calculates the function used in runge-kutta approximations
    to find the value of v at a new time step
    */
    double f = -(1/(2*nodeDistance)) 
            * (
              ((R*T/M)*((rho_up-rho_down)/rho))
            + (u*(v_right-v_left))
            + (v*(v_up-v_down))
            );
            
    return f;
}

double rho_function(double rho_right, double rho_left, double rho_up, double rho_down,
    double v_up, double v_down, double u_right, double u_left)
{
    /*
    This function calculates the function used in runge-kutta approximations
    to find the value of rho at a new time step
    */
    double f = -(1/(2*nodeDistance)) 
            * (
              (rho_right * u_right)
            - (rho_left * u_left)
            + (rho_up * v_up)
            - (rho_down * v_down)
            );
            
    return f;
}

int main(int argc, char *argv[]) {
    // Prepare Save File
    std::ofstream MyFile;
    MyFile.open("2D_fluid.txt");

    // Check Save File is open/created
    if (!MyFile.is_open())
    {
        std::cout<<"Error in creating file"<<std::endl; 
        return 1; 
    }
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

    world.component<velocityStart>();
    world.component<velocityHalfPredict>();
    world.component<velocityHalfCorrect>();
    world.component<velocityEndPredict>();
    world.component<velocityEndCorrect>();

    // Velocity components
    world.component<densityStart>();
    world.component<densityHalfPredict>();
    world.component<densityHalfCorrect>();
    world.component<densityEndPredict>();
    world.component<densityEnd>();

    // Wall Tags
    struct LowerWall {}; // This component indicates the node has a wall below it
    struct UpperWall {}; // This component indicates the node has a wall above it 
    struct LeftWall {};  // This component indicates the node has a wall left of it 
    struct RightWall {}; // This component indicates the node has a wall right of it 

    std::vector<flecs::entity> nodes; // place to store nodes
    nodes.reserve(Nx*Ny); // Create the space
    
    // Create left boundary Node
    nodes.push_back(
            world.entity()
                .add<End_Tag>() 
                .set<Position>({0})
                .set<Index>({0})
                .set<Phi_start>({Start})
                .set<Phihalf_predict>({Start})
                .set<Phihalf_correct>({Start})
                .set<Phi_end_predict>({Start})
            );
    
    // Create non-boundary Nodes
    for (int index = 1; index < N-1; ++index) {
        double phi_initial;
        phi_initial = linear_function(index*L/(N-1), Start, End);
        nodes.push_back(
            world.entity()
                .add<Middle_Tag>() 
                .set<Position>({index*L/(N-1)})
                .set<Index>({index})
                .set<Phi_start>({phi_initial})
                .set<Phihalf_predict>({0})
                .set<Phihalf_correct>({0})
                .set<Phi_end_predict>({0})
                .set<Phi_end_correct>({0})
            );
    }

    // Create Right Boundary Nodes
    nodes.push_back(
            world.entity()
                .add<End_Tag>() 
                .set<Position>({L})
                .set<Index>({N-1})
                .set<Phi_start>({End})
                .set<Phihalf_predict>({End})
                .set<Phihalf_correct>({End})
                .set<Phi_end_predict>({End})
            );
    
    // This system finds and updates phihalf_predict
    world.system<Position, Index, Phi_start, Phihalf_predict>()
        .with<Middle_Tag>()
        .kind(RungeKutta_1)
        .each([&](Position& pos, Index& index, Phi_start& start, Phihalf_predict& half){
            const Phi_start& startMinus = nodes[index.ind - 1].get<Phi_start>();
            const Phi_start& startPlus = nodes[index.ind + 1].get<Phi_start>();

            const Position& posMinus = nodes[index.ind - 1].get<Position>();
            const Position& posPlus = nodes[index.ind + 1].get<Position>();

            double function = fluid_function(startMinus.i, start.i, startPlus.i, 
                posMinus.x, pos.x, posPlus.x);
            half.i = start.i + (TIME*function)/(2*STEPS);
        });
    
    // This system finds and updates phihalf_correct
    world.system<Position, Index, Phi_start, Phihalf_predict, Phihalf_correct>()
        .with<Middle_Tag>()
        .kind(RungeKutta_2)
        .each([&](Position& pos, Index& index, Phi_start& start, 
            Phihalf_predict& predictor, Phihalf_correct& corrector){

            const Phihalf_predict& predictorMinus = nodes[index.ind - 1].get<Phihalf_predict>();
            const Phihalf_predict& predictorPlus = nodes[index.ind + 1].get<Phihalf_predict>();

            const Position& posMinus = nodes[index.ind - 1].get<Position>();
            const Position& posPlus = nodes[index.ind + 1].get<Position>();

            double function = fluid_function(predictorMinus.i, predictor.i, predictorPlus.i, 
                posMinus.x, pos.x, posPlus.x);
            
            corrector.i = start.i + (TIME*function)/(2*STEPS);
        });
    
    // This system finds and updates phi_end_predict
    world.system<Position, Index, Phi_start, Phihalf_correct, Phi_end_predict>()
        .with<Middle_Tag>()
        .kind(RungeKutta_3)
        .each([&](Position& pos, Index& index, Phi_start& start, 
            Phihalf_correct& half, Phi_end_predict& end){

            const Phihalf_correct& halfMinus = nodes[index.ind - 1].get<Phihalf_correct>();
            const Phihalf_correct& halfPlus = nodes[index.ind + 1].get<Phihalf_correct>();

            const Position& posMinus = nodes[index.ind - 1].get<Position>();
            const Position& posPlus = nodes[index.ind + 1].get<Position>();

            double function = fluid_function(halfMinus.i, half.i, halfPlus.i, 
                posMinus.x, pos.x, posPlus.x);
            
            end.i = start.i + (TIME*function)/(2*STEPS);
        });

    // This system finds and updates phi_end_correct
    world.system<Position, Index, Phi_start, Phihalf_predict, Phihalf_correct, Phi_end_predict, Phi_end_correct>()
        .with<Middle_Tag>()
        .kind(RungeKutta_4)
        .each([&](Position& pos, Index& index, Phi_start& start, Phihalf_predict& halfPred,
            Phihalf_correct& halfCorr, Phi_end_predict& endPred, Phi_end_correct& endCorr){
            
            const Phi_start& startMinus = nodes[index.ind - 1].get<Phi_start>();
            const Phi_start& startPlus = nodes[index.ind + 1].get<Phi_start>();

            const Phihalf_predict& halfPredMinus = nodes[index.ind - 1].get<Phihalf_predict>();
            const Phihalf_predict& halfPredPlus = nodes[index.ind + 1].get<Phihalf_predict>();
            
            const Phihalf_correct& halfCorrMinus = nodes[index.ind - 1].get<Phihalf_correct>();
            const Phihalf_correct& halfCorrPlus = nodes[index.ind + 1].get<Phihalf_correct>();

            const Phi_end_predict& endPredMinus = nodes[index.ind - 1].get<Phi_end_predict>();
            const Phi_end_predict& endPredPlus = nodes[index.ind + 1].get<Phi_end_predict>();

            const Position& posMinus = nodes[index.ind - 1].get<Position>();
            const Position& posPlus = nodes[index.ind + 1].get<Position>();

            double functionStart = fluid_function(startMinus.i, start.i, startPlus.i, 
                posMinus.x, pos.x, posPlus.x);

            double functionHalfPred = fluid_function(halfPredMinus.i, halfPred.i, halfPredPlus.i, 
                posMinus.x, pos.x, posPlus.x);
            
            double functionHalfCorr = fluid_function(halfCorrMinus.i, halfCorr.i, halfCorrPlus.i, 
                posMinus.x, pos.x, posPlus.x);
            
            double functionEndPred = fluid_function(endPredMinus.i, endPred.i, endPredPlus.i, 
                posMinus.x, pos.x, posPlus.x);
            
            endCorr.i = start.i + (TIME*(functionStart + 2*functionHalfPred + 2*functionHalfCorr 
                + functionEndPred))/(6*STEPS);
        });

    // This updates phi_start by replacing it with phi_end_correct
    world.system<Phi_start, Phi_end_correct>()
        .with<Middle_Tag>()
        .kind(Update)
        .each([](Phi_start& initial, Phi_end_correct& updated){
            initial.i = updated.i;
        });

    // Set .txt file headers
    MyFile << "t, ";
    for (int x_step = 0; x_step < N; ++x_step) {
        if (x_step == N-1) {
            MyFile << "x" << x_step << "\n";
        } else {
            MyFile << "x" << x_step << ", ";
        }
    }
    
    // Run through systems every time step
    for (int t_step = 0; t_step <= STEPS; ++t_step) {
        if (t_step != 0) {
                world.progress();
            }
        std::cout << t_step << "\n";
        // Saves Data to a .txt file
        MyFile << static_cast<double>(t_step*TIME)/STEPS << ", ";
        for (int index = 0; index < N-1; ++index) {
            const Phi_start& phi = nodes[index].get<Phi_start>();
            MyFile << phi.i << ", ";
        }
        const Phi_start& phi = nodes[N-1].get<Phi_start>();
        MyFile << phi.i << "\n";
    }

    // Close Save File
    MyFile.close();
}
