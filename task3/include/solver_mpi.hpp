#pragma once

#include <omp.h>
#include <mpi.h>
#include <iostream>
#include <vector>

#include "grid.hpp"
#include "function.hpp"
#include "block.hpp"

using namespace std;

class Solver_MPI {
private:
    Grid g;
    Function f;
    Block local;
    vector<Block> blocks;
    vector<vector<double>> u;
    vector<pair<int, Block>> send_blocks;
    vector<pair<int, Block>> recv_blocks;
    vector<vector<double>> recv_data;
    int proc_rank, proc_size;

    int getIdx(int i, int j, int k, const Block& b) {
        return (i - b.x_l) * b.y_len * b.z_len + (j - b.y_l) * b.z_len + (k - b.z_l);
    }

    void split() {
        int nx, ny, nz;
        tie(nx, ny, nz) = Block::getBalancedDivision(proc_size);
        int x_len = (g.N + 1) / nx;
        int y_len = (g.N + 1) / ny;
        int z_len = (g.N + 1) / nz;

        for (int i = 0, x_l = 0; i < nx; i++) {
            int x_r = (i == nx - 1) ? g.N : x_l + x_len - 1;
            for (int j = 0, y_l = 0; j < ny; j++) {
                int y_r = (j == ny - 1) ? g.N : y_l + y_len - 1;
                for (int k = 0, z_l = 0; k < nz; k++) {
                    int z_r = (k == nz - 1) ? g.N : z_l + z_len - 1;
                    blocks.emplace_back(x_l, x_r, y_l, y_r, z_l, z_r);
                    z_l = z_r + 1;
                }
                y_l = y_r + 1;
            }
            x_l = x_r + 1;
        }
    }

    bool isBeside(const Block& b1, const Block& b2) {
        int same = 0;
        if (b1.x_l == b2.x_l && b1.x_r == b2.x_r) {
            same++;
        }
        if (b1.y_l == b2.y_l && b1.y_r == b2.y_r) {
            same++;
        }
        if (b1.z_l == b2.z_l && b1.z_r == b2.z_r) {
            same++;
        }
        return same == 2;
    }

    void getNeighbors() {
        for (int i = 0; i < proc_size; i++) {
            if (i != proc_rank) {
                Block remote = blocks[i];
                if (isBeside(local, remote)) {
                    if (local.x_l == remote.x_r + 1) {
                        send_blocks.emplace_back(i, Block(local.x_l, local.x_l, local.y_l, local.y_r, local.z_l, local.z_r));
                        recv_blocks.emplace_back(i, Block(remote.x_r, remote.x_r, local.y_l, local.y_r, local.z_l, local.z_r));
                    }
                    else if (remote.x_l == local.x_r + 1) {
                        send_blocks.emplace_back(i, Block(local.x_r, local.x_r, local.y_l, local.y_r, local.z_l, local.z_r));
                        recv_blocks.emplace_back(i, Block(remote.x_l, remote.x_l, local.y_l, local.y_r, local.z_l, local.z_r));
                    }
                    else if (local.y_l == remote.y_r + 1) {
                        send_blocks.emplace_back(i, Block(local.x_l, local.x_r, local.y_l, local.y_l, local.z_l, local.z_r));
                        recv_blocks.emplace_back(i, Block(local.x_l, local.x_r, remote.y_r, remote.y_r, local.z_l, local.z_r));
                    }
                    else if (remote.y_l == local.y_r + 1) {
                        send_blocks.emplace_back(i, Block(local.x_l, local.x_r, local.y_r, local.y_r, local.z_l, local.z_r));
                        recv_blocks.emplace_back(i, Block(local.x_l, local.x_r, remote.y_l, remote.y_l, local.z_l, local.z_r));
                    }
                    else if (local.z_l == remote.z_r + 1) {
                        send_blocks.emplace_back(i, Block(local.x_l, local.x_r, local.y_l, local.y_r, local.z_l, local.z_l));
                        recv_blocks.emplace_back(i, Block(local.x_l, local.x_r, local.y_l, local.y_r, remote.z_r, remote.z_r));
                    }
                    else if (remote.z_l == local.z_r + 1) {
                        send_blocks.emplace_back(i, Block(local.x_l, local.x_r, local.y_l, local.y_r, local.z_r, local.z_r));
                        recv_blocks.emplace_back(i, Block(local.x_l, local.x_r, local.y_l, local.y_r, remote.z_l, remote.z_l));
                    }
                }
            }
        }
    }

    vector<double> getSendData(const vector<double>& u_n, const Block& remote) {
        vector<double> send_data(remote.size);
        #pragma omp parallel for collapse(3)
        for (int i = remote.x_l; i <= remote.x_r; i++) {
            for (int j = remote.y_l; j <= remote.y_r; j++) {
                for (int k = remote.z_l; k <= remote.z_r; k++) {
                    send_data[getIdx(i, j, k, remote)] = u_n[getIdx(i, j, k, local)];
                }
            }
        }
        return send_data;
    }

    void sendRecv(const vector<double>& u_n) {
        recv_data.resize(recv_blocks.size());
        vector<MPI_Request> requests(2);
        vector<MPI_Status> statuses(2);
        for (int i = 0; i < recv_blocks.size(); i++) {
            vector<double> send_data = getSendData(u_n, send_blocks[i].second);
            recv_data[i].resize(recv_blocks[i].second.size);
            MPI_Isend(send_data.data(), send_blocks[i].second.size, MPI_DOUBLE, send_blocks[i].first, 0, MPI_COMM_WORLD, &requests[0]);
            MPI_Irecv(recv_data[i].data(), recv_blocks[i].second.size, MPI_DOUBLE, recv_blocks[i].first, 0, MPI_COMM_WORLD, &requests[1]);
            MPI_Waitall(2, requests.data(), statuses.data());
        }
    }

    double getU(const vector<double>& u_n, int i, int j, int k) {
        if (local.x_l <= i && i <= local.x_r && local.y_l <= j && j <= local.y_r && local.z_l <= k && k <= local.z_r) {
            return u_n[getIdx(i, j, k, local)];
        }
        for (int r_i = 0; r_i < recv_blocks.size(); r_i++) {
            Block remote = recv_blocks[r_i].second;
            if (remote.x_l <= i && i <= remote.x_r && remote.y_l <= j && j <= remote.y_r && remote.z_l <= k && k <= remote.z_r) {
                return recv_data[r_i][getIdx(i, j, k, remote)];
            }
        }
        throw runtime_error("Value of u not found");
    }

    double laplace(const vector<double>& u_n, int i, int j, int k) {
        return (getU(u_n, i - 1, j, k) - 2 * u_n[getIdx(i, j, k, local)] + getU(u_n, i + 1, j, k)) / (g.h_x * g.h_x) + (getU(u_n, i, j - 1, k) - 2 * u_n[getIdx(i, j, k, local)] + getU(u_n, i, j + 1, k)) / (g.h_y * g.h_y) + (getU(u_n, i, j, k - 1) - 2 * u_n[getIdx(i, j, k, local)] + getU(u_n, i, j, k + 1)) / (g.h_z * g.h_z);
    }

    void calculateBoundary(vector<double>& u_n, double t) {
        if (local.x_l == 0) {
            #pragma omp parallel for collapse(2)
            for (int i = local.y_l; i <= local.y_r; i++) {
                for (int j = local.z_l; j <= local.z_r; j++) {
                    u_n[getIdx(0, i, j, local)] = 0;
                }
            }
        }

        if (local.x_r == g.N) {
            #pragma omp parallel for collapse(2)
            for (int i = local.y_l; i <= local.y_r; i++) {
                for (int j = local.z_l; j <= local.z_r; j++) {
                    u_n[getIdx(g.N, i, j, local)] = 0;
                }
            }
        }

        if (local.y_l == 0) {
            #pragma omp parallel for collapse(2)
            for (int i = local.x_l; i <= local.x_r; i++) {
                for (int j = local.z_l; j <= local.z_r; j++) {
                    u_n[getIdx(i, 0, j, local)] = f(i * g.h_x, 0, j * g.h_z, t);
                }
            }
        }

        if (local.y_r == g.N) {
            #pragma omp parallel for collapse(2)
            for (int i = local.x_l; i <= local.x_r; i++) {
                for (int j = local.z_l; j <= local.z_r; j++) {
                    u_n[getIdx(i, g.N, j, local)] = f(i * g.h_x, g.L_y, j * g.h_z, t);
                }
            }
        }

        if (local.z_l == 0) {
            #pragma omp parallel for collapse(2)
            for (int i = local.x_l; i <= local.x_r; i++) {
                for (int j = local.y_l; j <= local.y_r; j++) {
                    u_n[getIdx(i, j, 0, local)] = 0;
                }
            }
        }

        if (local.z_r == g.N) {
            #pragma omp parallel for collapse(2)
            for (int i = local.x_l; i <= local.x_r; i++) {
                for (int j = local.y_l; j <= local.y_r; j++) {
                    u_n[getIdx(i, j, g.N, local)] = 0;
                }
            }
        }
    }

    void init() {
        // Initialize boundaries
        calculateBoundary(u[0], 0);
        calculateBoundary(u[1], g.tau);

        // 0 and N are not taken
        int x_l = max(1, local.x_l);
        int x_r = min(g.N - 1, local.x_r);
        int y_l = max(1, local.y_l);
        int y_r = min(g.N - 1, local.y_r);
        int z_l = max(1, local.z_l);
        int z_r = min(g.N - 1, local.z_r);

        // Initialize u_0
        #pragma omp parallel for collapse(3)
        for (int i = x_l; i <= x_r; i++) {
            for (int j = y_l; j <= y_r; j++) {
                for (int k = z_l; k <= z_r; k++) {
                    u[0][getIdx(i, j, k, local)] = f(i * g.h_x, j * g.h_y, k * g.h_z, 0);
                }
            }
        }

        sendRecv(u[0]);

        // Initialize u_1
        #pragma omp parallel for collapse(3)
        for (int i = x_l; i <= x_r; i++) {
            for (int j = y_l; j <= y_r; j++) {
                for (int k = z_l; k <= z_r; k++) {
                    u[1][getIdx(i, j, k, local)] = u[0][getIdx(i, j, k, local)] + f.a_2 * g.tau * g.tau / 2.0 * laplace(u[0], i, j, k);
                }
            }
        }
    }

    double calculateError(int steps) {
        double error = 0;
        #pragma omp parallel for collapse(3) reduction(max: error)
        for (int i = local.x_l; i <= local.x_r; i++) {
            for (int j = local.y_l; j <= local.y_r; j++) {
                for (int k = local.z_l; k <= local.z_r; k++) {
                    double analytical = f(i * g.h_x, j * g.h_y, k * g.h_z, steps * g.tau);
                    double result = u[steps % 3][getIdx(i, j, k, local)];
                    double cur_error = abs(result - analytical);
                    error = max(error, cur_error);
                }
            }
        }
        double error_global = 0;
        MPI_Reduce(&error, &error_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        return error_global;
    }

public:
    Solver_MPI(Grid g) : g(g), f(g) {
        MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &proc_size);
        split();
        local = blocks[proc_rank];
        u.resize(3);
        for (int i = 0; i < 3; i++) {
            u[i].resize(local.size);
        }
        getNeighbors();
    }

    double solve(int steps) {
        init();
        
        // 0 and N are not taken
        int x_l = max(1, local.x_l);
        int x_r = min(g.N - 1, local.x_r);
        int y_l = max(1, local.y_l);
        int y_r = min(g.N - 1, local.y_r);
        int z_l = max(1, local.z_l);
        int z_r = min(g.N - 1, local.z_r);

        for (int step = 2; step <= steps; step++) {
            sendRecv(u[(step + 2) % 3]);

            #pragma omp parallel for collapse(3)
            for (int i = x_l; i <= x_r; i++) {
                for (int j = y_l; j <= y_r; j++) {
                    for (int k = z_l; k <= z_r; k++) {
                        u[step % 3][getIdx(i, j, k, local)] = 2 * u[(step + 2) % 3][getIdx(i, j, k, local)] - u[(step + 1) % 3][getIdx(i, j, k, local)] + f.a_2 * g.tau * g.tau * laplace(u[(step + 2) % 3], i, j, k);
                    }
                }
            }
            calculateBoundary(u[step % 3], step * g.tau);
        }
        return calculateError(steps);
    }

    // ~Solver() {}
};
