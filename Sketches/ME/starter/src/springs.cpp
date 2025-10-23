/*
This code runs a simulation of N coupled oscillators.

Each mass has equations of motion depending on the properties
the springs it is connected to. 

First the acceleration of each particle is computed.

Than the program iterates to update the position, velocity, 
and acceleration of each entity for which that is relevant.

Entity types:
 - Particle
*/

#include <iostream>
#include <fstream> 
#include <vector>
#include <flecs.h>
#include <systems.h>

int N = 2; //Particle Number
float time_step = 0.001; //How much time passes in each iteration
int run_time = 1000; 

//Position of the Left Wall
double p_Lwall = 0;  // in cm

// Position of the Right Wall
double p_Rwall = 12; // in cm

struct Index {
    int i; // Which particle is it
};

struct Position { 
    double x; // In cm
};

struct Velocity { 
    double x;// In cms-1
};

struct Acceleration{
    double x; // in cms-2
};

struct Mass{
    double M; // in kg
};

double find_acceleration(float position, float position_left, float position_right,
const float mass, float k_left, float k_right, float l_left, float l_right){

    float acc;
    acc = - (k_left/mass) * (position - position_left - l_left) 
        + (k_right/mass) * (position_right - position - l_right);
    
    return acc;
}


int main() {
    std::ofstream MyFile("Coupled_Oscillators.txt");

    // First Particle Data
    std::vector<std::vector<double>> p_matrix {{2}, {7}}; 
    std::vector<std::vector<double>> v_matrix {{2}, {-1}};
    std::vector<std::vector<double>> a_matrix {{0}, {0}};

    //Spring Constant Data
    const std::vector<double> k_list {3, 3, 3};
    // Natural Spring Length Data
    const std::vector<double> l_list {4, 4, 4};

    flecs::world ecs;
    
    /*
        System s1 runs for any entity with the tag Particle
        Its job is to calculate the initial acceleration and
        update the vectors and entity components. 
    */
    flecs::system s1 = ecs.system<Acceleration, Position, Mass, Index>()
    .each([&](Acceleration &a, Position &p, const Mass &mass, const Index &index) {
        if (index.i == 0) {
            a.x = find_acceleration(p.x, p_Lwall, p_matrix[index.i+1][0], mass.M, k_list[index.i], 
                                    k_list[index.i+1], l_list[index.i], l_list[index.i+1]);
        } else if (index.i == N-1) {
            a.x = find_acceleration(p.x, p_matrix[index.i-1][0], p_Rwall, mass.M, k_list[index.i], 
                                    k_list[index.i+1], l_list[index.i], l_list[index.i+1]);
        } else {
            a.x = find_acceleration(p.x, p_matrix[index.i-1][0], p_matrix[index.i+1][0], mass.M, k_list[index.i], 
                                    k_list[index.i+1], l_list[index.i], l_list[index.i+1]);
        } 
        a_matrix[index.i][0] = a.x;
        
    });
    
    
    flecs::system s2 = ecs.system<Position, Velocity, Acceleration, Mass, Index>()
    .each([&](flecs::entity e, Position &p, Velocity &v, Acceleration &a, const Mass &mass, const Index &index) {
        
        p.x += v.x*time_step;
        v.x += a.x*time_step;
        
        if (index.i == 0) {
            a.x = find_acceleration(p.x, p_Lwall, p_matrix[index.i+1][0], mass.M, k_list[index.i], 
                                    k_list[index.i+1], l_list[index.i], l_list[index.i+1]);
        } else if (index.i == N-1) {
            a.x = find_acceleration(p.x, p_matrix[index.i-1][0], p_Rwall, mass.M, k_list[index.i], 
                                    k_list[index.i+1], l_list[index.i], l_list[index.i+1]);
        } else {
            a.x = find_acceleration(p.x, p_matrix[index.i-1][0], p_matrix[index.i+1][0], mass.M, k_list[index.i], 
                                    k_list[index.i+1], l_list[index.i], l_list[index.i+1]);
        } 
        
        p_matrix[index.i].push_back(p.x);
        v_matrix[index.i].push_back(v.x);
        a_matrix[index.i].push_back(a.x);
        
        std::cout << e.name() << ": {" << p.x << ", " << v.x << "," << a.x << "}\n";
        
    });
    
    // Create Entities
    ecs.entity("mass 1")
        // Finds and sets components
        .set<Index>({0})
        .set<Mass>({4})
        .set<Position>({p_matrix[0][0]})
        .set<Velocity>({v_matrix[0][0]})
        .set<Acceleration>({0});  
    
    ecs.entity("mass 2")
        // Finds and sets components
        .set<Index>({1})
        .set<Mass>({4})
        .set<Position>({p_matrix[1][0]})
        .set<Velocity>({v_matrix[1][0]})
        .set<Acceleration>({0});
    
    s1.run();

    MyFile << "Time (s), Position 1 (cm), Velocity 1 (cm s-1), Acceleration 1 (cm s-2)" 
           << "Position 2 (cm), Velocity 2 (cm s-1), Acceleration 2 (cm s-2)" << std::endl; 
    MyFile << 0 << ", " << p_matrix[0][0] << ", " << v_matrix[0][0] << "," << a_matrix[0][0] << "," 
           << p_matrix[1][0] << ", " << v_matrix[1][0] << "," << a_matrix[1][0] << std::endl; 

    for(int iter = 1; iter < run_time; iter++) {
        s2.run();
        std::cout << "----\n";
        
        MyFile << iter*time_step << ", " << p_matrix[0][iter] << ", " << v_matrix[0][iter] << "," << a_matrix[0][iter] << "," 
               << p_matrix[1][iter] << ", " << v_matrix[1][iter] << "," << a_matrix[1][iter] << std::endl; 
        
    }
    MyFile.close();
}
