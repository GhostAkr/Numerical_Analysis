//
// Created by ighos on 08.02.2019.
//

#ifndef NUMERICAL_ANALYSIS_ORDINARYDE_H
#define NUMERICAL_ANALYSIS_ORDINARYDE_H

#include <iostream>
#include <vector>
#include <fstream>

using std::cout;
using std::endl;
using std::vector;

// Tests

vector<double> func0(vector<double> _point);

// Methods

void EulerExplicit(vector<double> _func(vector<double>), vector<double> _startPoint);

// Other

void vectorPrint(vector<double> _sourceVector);

#endif //NUMERICAL_ANALYSIS_ORDINARYDE_H
