#pragma once

struct Grid {
    double L_x, L_y, L_z;
    double h_x, h_y, h_z;
    int N;
    double T;
    int K;
    double tau;

    Grid(double L_x, double L_y, double L_z, int N, double T, int K) : L_x(L_x), L_y(L_y), L_z(L_z), N(N), T(T), K(K) {
        this -> h_x = L_x / N;
        this -> h_y = L_y / N;
        this -> h_z = L_z / N;
        this -> tau = T / K;
    }
};
