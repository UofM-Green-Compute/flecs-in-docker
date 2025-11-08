/*
Program to trial gauss elimination
*/

#include <algorithm>
#include <vector>
#include <cmath>
#include <flecs.h>
#include <systems.h>
#include <iostream> 

struct Matrix { std::vector<std::vector<double>> M; }; 

std::vector<double> row_manipulation(std::vector<double> ith_row, std::vector<double> first_row, int column)
{
    for (int i = column; i < ith_row.size(); i++)
    {
        ith_row[i] = ith_row[i] - first_row[i] * ith_row[column] / first_row[column]; 
    }

    return ith_row; 
}

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

int main() {

    flecs::world ecs; 
    flecs::system matrix_solver = ecs.system<Matrix>()
    .each([&](Matrix &matrix)
    {
        int no_rows = matrix.M.size();
        int no_columns = matrix.M[0].size(); 

        for (int i = 0; i < no_columns-1; i++)
        {
            for(int j = i+1; j < no_rows; j++)
            {
            matrix.M[j] = row_manipulation(matrix.M[j], matrix.M[j], i);
            }
        }

        print_matrix(matrix.M);

    });

    // Create Matrix 
    std::vector<std::vector<double>> matrix = {{1,-3,4},{2,-5,6},{-3,3,4}}; 
    print_matrix(matrix); 

    ecs.entity("Matrix 1")
        .set<Matrix>({matrix});

    matrix_solver.run(); 
}