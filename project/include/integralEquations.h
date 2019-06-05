//
// Created by ighos on 05.06.2019.
//

#ifndef NUMERICAL_ANALYSIS_INTEGRALEQUATIONS_H
#define NUMERICAL_ANALYSIS_INTEGRALEQUATIONS_H

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

// Tests
double core(double _x, double _s, int _testNum);
double rightPart(double _x, int _testNum);
vector<double> limits(int _limitType);

#endif //NUMERICAL_ANALYSIS_INTEGRALEQUATIONS_H
