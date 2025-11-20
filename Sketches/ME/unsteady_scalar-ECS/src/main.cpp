/*
This code looks at the dynamical evolution of a 1D convection-conduction
system which begins out of equilibrium
*/

#include <custom_phases_no_builtin.h>
#include <iostream>
#include <fstream> 
#include <vector>

// Sets time parameters
int STEPS = 100;
int TIME = 1;

// Number of nodes in x-axis. Also equal to the length in unit of node distance
const double L = 1; // Length
const int N = 41; // Number of nodes
const int PE = 50; // Peclet Number
const double RHO = 2; // Density
const double U = 5; // Speed

// Calculates Diffusivity
double GAMMA = (RHO * U * L) / PE;

// Boundary Conditions
double Start = 0;
double End = 1;

// Position Components
struct Position {double x; };
struct Index {int ind; };

// Field components
struct Phi_start { double i; };
struct Phihalf_predict { double i; };
struct Phihalf_correct { double i; };
struct Phi_end_predict { double i; };
struct Phi_end_correct { double i; };

// Component Tags
struct Middle_Tag {};
struct End_Tag {};

double linear_function(double x, double A, double B)
{
    double m = (B-A)/L;
    double y = m * x;
    return y;
}

double fluid_function(double phi_left, double phi, double phi_right, 
    double x_left, double x, double x_right)
{
    double f;
    f = - U*((phi_right - phi_left)/(x_right - x_left)) 
    + ((GAMMA/RHO) * (phi_right - 2 * phi + phi_left)) / ((x_right - x) * (x - x_left));
    return f;
}

int main(int argc, char *argv[]) {
    // Prepare Save File
    std::ofstream MyFile;
    MyFile.open("1D_dynamic_fluid.txt");

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
    world.component<Middle_Tag>();
    world.component<End_Tag>();

    world.component<Position>();

    world.component<Phi_start>();
    world.component<Phihalf_predict>();
    world.component<Phihalf_correct>();
    world.component<Phi_end_predict>();
    world.component<Phi_end_correct>();

    std::vector<flecs::entity> nodes; // place to store nodes
    nodes.reserve(N); // Create the space
    
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
            std::cout << half.i << "\n";

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
        // Saves Data to a .txt file
        MyFile << t_step << ", ";
        for (int index = 0; index < N-1; ++index) {
            const Phi_start& phi = nodes[index].get<Phi_start>();
            MyFile << phi.i << ", ";
        }
        const Phi_start& phi = nodes[N-1].get<Phi_start>();
        MyFile << phi.i << "\n";

        
    }
    
    world.progress();
    // Close Save File
    MyFile.close();
}
