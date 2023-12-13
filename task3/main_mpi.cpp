#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <fstream>

#include "include/solver_mpi.hpp"
#include "include/grid.hpp"

using namespace std;

int main(int argc, char **argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double L = argc == 3 ? stod(argv[2]) : M_PI;
    int N = stoi(argv[1]);
    double T = 1.0;
    int K = 10000;
    int steps = 20;
    string filename;

    double start_time = MPI_Wtime();
    Grid grid = Grid(L, L, L, N, T, K);
    Solver_MPI solver = Solver_MPI(grid);
    double error = solver.solve(steps);
    double end_time = MPI_Wtime();
    double duration = end_time - start_time;
    double time;
    MPI_Reduce(&duration, &time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    #pragma omp parallel
    if (omp_get_thread_num() == 0) {
        int num_threads = omp_get_num_threads();
        if (rank == 0) {
            int num_processes = size;
            string str_L = L == M_PI ? "pi" : to_string((int)L);
            filename = "mpi_result_" + str_L + "_" + to_string(N) + ".txt";
            ofstream f(filename, ios_base::app);
            f << "Number of processes: " << num_processes << "   Number of threads: " << num_threads << "   Error: " << error << "   Running time: " << time << endl;
            f.close();
        }
    }

    MPI_Finalize();
    return 0;
}