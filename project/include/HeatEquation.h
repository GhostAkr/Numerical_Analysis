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
void quasilinearEquationScheme(double _h, double _t, string _path);

// Boundary conditions
vector<double> boundaryLeft(int _type, double _y0, double _y1, double _h, double _t, double _sigma, double _p, double _p1, double _L, double _layer);  // 0 --- temperature; 1 --- heat flow
vector<double> boundaryRight(int _type, double _y0, double _y1, double _h, double _t, double _sigma, double _p, double _p1, double _L, double _N, double _layer);  // 0 --- temperature; 1 --- heat flow

// Special functions
double K1 (double _x, double _L);
double K2 (double _u);

// Other
vector<double> tridiagonalLinearSolve(const vector<vector<double>> _matrix);
double exactQuasi(double _x, double _t);
double exactTest(double _x, double _t);
void exactQuasiChart(double _h, double _t, string _path);

#endif //NUMERICAL_ANALYSIS_HEATEQUATION_H
