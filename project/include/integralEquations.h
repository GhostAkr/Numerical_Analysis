//
// Created by ighos on 05.06.2019.
//

#ifndef NUMERICAL_ANALYSIS_INTEGRALEQUATIONS_H
#define NUMERICAL_ANALYSIS_INTEGRALEQUATIONS_H

#include <iostream>
#include <cmath>
#include <vector>

#include "../include/OrdinaryDE.h"

using namespace std;

// Tests
double core(double _x, double _s, int _testNum);
double rightPart(double _x, int _testNum);
vector<double> limits(int _limitType);

// Methods
void quadMethod(int _testNum, int _limitType, int _N, string _path);
void simpleMethod(int _testNum, int _limitType, int _N, string _path);
void degenMethod(int _testNum, int _limitType, int _N, string _path);

// Core decomposition
double phi(double _x, int _i);
double psi(double _s, int _i);

#endif //NUMERICAL_ANALYSIS_INTEGRALEQUATIONS_H
