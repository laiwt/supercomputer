#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

class Saver {
private:
    fstream analytical;
    fstream result;
    fstream error;
public:
    Saver() {}

    Saver(const Saver&) {}

    void openFiles() {
        analytical.open("./data/analytical.csv", ios::out | ios::trunc);
        result.open("./data/result.csv", ios::out | ios::trunc);
        error.open("./data/error.csv", ios::out | ios::trunc);
        analytical << "i,j,k,value" << endl;
        result << "i,j,k,value" << endl;
        error << "i,j,k,value" << endl;
    }

    void closeFiles() {
        analytical.close();
        result.close();
        error.close();
    }

    void save(const string& valueType, double i, double j, double k, double value) {
        if (valueType == "analytical") {
            analytical << std::fixed << std::setprecision(8) << i << "," << j << "," << k << "," << value << endl;
        }
        else if (valueType == "result") {
            result << std::fixed << std::setprecision(8) << i << "," << j << "," << k << "," << value << endl;
        }
        else if (valueType == "error") {
            error << std::fixed << std::setprecision(8) << i << "," << j << "," << k << "," << value << endl;
        }
    }

    ~Saver() {}
};
