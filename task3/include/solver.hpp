#pragma once

#include <omp.h>
#include <iostream>
#include <vector>

#include "grid.hpp"
#include "function.hpp"
#include "saver.hpp"

using namespace std;

class Solver {
private:
    Grid g;
    Function f;
    Saver saver;
    vector<vector<vector<vector<double>>>> u;

    double laplace(const vector<vector<vector<double>>>& u_n, int i, int j, int k) {
        return (u_n[i - 1][j][k] - 2 * u_n[i][j][k] + u_n[i + 1][j][k]) / (g.h_x * g.h_x) + (u_n[i][j - 1][k] - 2 * u_n[i][j][k] + u_n[i][j + 1][k]) / (g.h_y * g.h_y) + (u_n[i][j][k - 1] - 2 * u_n[i][j][k] + u_n[i][j][k + 1]) / (g.h_z * g.h_z);
    }

    void calculateBoundary(vector<vector<vector<double>>>& u_n, double t) {
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < g.N; i++) {
            for (int j = 0; j < g.N; j++) {
                // x - First kind
                u_n[0][i][j] = 0;
                u_n[g.N][i][j] = 0;
                // y - Periodic
                u_n[i][0][j] = f(i * g.h_x, 0, j * g.h_z, t);
                u_n[i][g.N][j] = f(i * g.h_x, g.L_y, j * g.h_z, t);
                // z - First kind
                u_n[i][j][0] = 0;
                u_n[i][j][g.N] = 0;
            }
        }
    }

    void init() {
        // Initialize boundaries
        calculateBoundary(u[0], 0);
        calculateBoundary(u[1], g.tau);

        // Initialize u_0
        #pragma omp parallel for collapse(3)
        for (int i = 1; i < g.N; i++) {
            for (int j = 1; j < g.N; j++) {
                for (int k = 1; k < g.N; k++) {
                    u[0][i][j][k] = f(i * g.h_x, j * g.h_y, k * g.h_z, 0);
                }
            }
        }

        // Initialize u_1
        #pragma omp parallel for collapse(3)
        for (int i = 1; i < g.N; i++) {
            for (int j = 1; j < g.N; j++) {
                for (int k = 1; k < g.N; k++) {
                    u[1][i][j][k] = u[0][i][j][k] + f.a_2 * g.tau * g.tau / 2.0 * laplace(u[0], i, j, k);
                }
            }
        }
    }

    double calculateError(int steps) {
        double error = 0;
        #pragma omp parallel for collapse(3) reduction(max: error)
        for (int i = 0; i <= g.N; i++) {
            for (int j = 0; j <= g.N; j++) {
                for (int k = 0; k <= g.N; k++) {
                    double analytical = f(i * g.h_x, j * g.h_y, k * g.h_z, steps * g.tau);
                    double result = u[steps % 3][i][j][k];
                    double cur_error = abs(result - analytical);
                    error = max(error, cur_error);
                }
            }
        }
        return error;
    }

    void saveData(int steps) {
        saver.openFiles();
        for (int i = 0; i <= g.N; i++) {
            for (int j = 0; j <= g.N; j++) {
                for (int k = 0; k <= g.N; k++) {
                    double analytical = f(i * g.h_x, j * g.h_y, k * g.h_z, steps * g.tau);
                    double result = u[steps % 3][i][j][k];
                    double error = abs(result - analytical);

                    saver.save("analytical", i * g.h_x, j * g.h_y, k * g.h_z, analytical);
                    saver.save("result", i * g.h_x, j * g.h_y, k * g.h_z, result);
                    saver.save("error", i * g.h_x, j * g.h_y, k * g.h_z, error);
                }
            }
        }
        saver.closeFiles();
    }
public:
    Solver(Grid g) : g(g), f(g) {
        u.resize(3);
        for (int n = 0; n < 3; n++) {
            u[n].resize(g.N + 1);
            for (int i = 0; i <= g.N; i++) {
                u[n][i].resize(g.N + 1);
                for (int j = 0; j <= g.N; j++) {
                    u[n][i][j].resize(g.N + 1);
                }
            }
        }
    }

    double solve(int steps, bool save) {
        init();
        for (int step = 2; step <= steps; step++) {
            #pragma omp parallel for collapse(3)
            for (int i = 1; i < g.N; i++) {
                for (int j = 1; j < g.N; j++) {
                    for (int k = 1; k < g.N; k++) {
                        u[step % 3][i][j][k] = 2 * u[(step + 2) % 3][i][j][k] - u[(step + 1) % 3][i][j][k] + f.a_2 * g.tau * g.tau * laplace(u[(step + 2) % 3], i, j, k);
                    }
                }
            }
            calculateBoundary(u[step % 3], step * g.tau);
        }
        double error = calculateError(steps);
        if (save) {
            saveData(steps);
        }
        return error;
    }

    ~Solver() {}
};
