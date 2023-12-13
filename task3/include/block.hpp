#pragma once

#include <tuple>
#include <climits>

using namespace std;

struct Block {
    int x_l, x_r, x_len;
    int y_l, y_r, y_len;
    int z_l, z_r, z_len;
    int size;

    Block() : x_l(0), x_r(0), y_l(0), y_r(0), z_l(0), z_r(0), x_len(0), y_len(0), z_len(0), size(0) {}

    Block(int x_l, int x_r, int y_l, int y_r, int z_l, int z_r) : x_l(x_l), x_r(x_r), y_l(y_l), y_r(y_r), z_l(z_l), z_r(z_r) {
        x_len = x_r - x_l + 1;
        y_len = y_r - y_l + 1;
        z_len = z_r - z_l + 1;
        size = x_len * y_len * z_len;
    }

    static tuple<int, int, int> getBalancedDivision(int num) {
        tuple<int, int, int> res;
        int min_sum = INT_MAX;
        for (int i = 1; i * i * i <= num; i++) {
            if (num % i == 0) {
                for (int j = i; j * j <= num / i; j++) {
                    if ((num / i) % j == 0) {
                        int k = num / (i * j);
                        int sum = i + j + k;
                        if (sum < min_sum) {
                            res = make_tuple(i, j, k);
                            min_sum = sum;
                        }
                    }
                }
            }
        }
        return res;
    }
};
