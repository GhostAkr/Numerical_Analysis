//
// Created by ighos on 19.03.2019.
//

#ifndef NUMERICAL_ANALYSIS_HEATEQUATION_H
#define NUMERICAL_ANALYSIS_HEATEQUATION_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "../include/OrdinaryDE.h"

using std::cout;
using std::cin;
using std::endl;
using std::string;
using std::vector;

vector<double> mesh(size_t _N, double _L);
vector<double> Mixed(vector<double> _V0, double sigma, double h, double t);  // _V0 - previous time layer
void Fill(string path);
void explicitScheme(double _t, double _h);

// Schemes
void mixedScheme(double _sigma, double _h, double _t, string _path);  // Mixed scheme with tridiagonal function

// Special functions
double K1 (double _x, double _L);

// Other
vector<double> tridiagonalLinearSolve(const vector<vector<double>> _matrix);

#endif //NUMERICAL_ANALYSIS_HEATEQUATION_H
