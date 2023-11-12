#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include "grid.hpp"

struct Function {
    Grid g;
    double a_t;
    double a_2;

    Function(Grid g) : g(g) {
        a_t = (M_PI / 2) * sqrt(1.0 / (g.L_x * g.L_x) + 4.0 / (g.L_y * g.L_y) + 9.0 / (g.L_z * g.L_z));
        a_2 = 1.0 / 4;
    }

    double operator()(double x, double y, double z, double t) {
        return sin(M_PI * x / g.L_x) * sin(2 * M_PI * y / g.L_y) * sin(3 * M_PI * z / g.L_z) * cos(a_t * t);
    }
};
