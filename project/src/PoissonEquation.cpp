//
// Created by ighos on 16.05.2019.
//

#include "../include/PoissonEquation.h"

vector<vector<double>> mesh2D(double _L1, double _L2, double _h1, double _h2) {
    int N1 = _L1 / _h1 + 1;
    int N2 = _L2 / _h2 + 1;
    vector<vector<double>> result (N1, vector<double> (N2, 0));
    return result;
}
