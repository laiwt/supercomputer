#define _USE_MATH_DEFINES

#include <chrono>
#include <iostream>
#include <cmath>

#include "include/solver.hpp"
#include "include/grid.hpp"

using namespace std;

int main(int argc, char **argv) {
    // omp_set_num_threads(8);

    double L = argc == 3 ? stod(argv[2]) : M_PI;
    int N = stoi(argv[1]);
    double T = 1.0;
    int K = 10000;
    int steps = 20;
    bool save = false;
    string filename;

    auto start_time = chrono::high_resolution_clock::now();
    Grid grid = Grid(L, L, L, N, T, K);
    Solver solver = Solver(grid);
    double error = solver.solve(steps, save);
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    double time = duration.count() / 1000000.0;

    #pragma omp parallel
    if (omp_get_thread_num() == 0) {
        int num_threads = omp_get_num_threads();
        string str_L = L == M_PI ? "pi" : to_string((int)L);
        filename = "result_" + str_L + "_" + to_string(N) + ".txt";
        ofstream f(filename, ios_base::app);
        f  << "Number of threads: " << num_threads << "   Error: " << error << "   Running time: " << time << endl;
        f.close();
    }
    return 0;
}