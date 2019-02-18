//
// Created by ighos on 08.02.2019.
//

#ifndef NUMERICAL_ANALYSIS_ORDINARYDE_H
#define NUMERICAL_ANALYSIS_ORDINARYDE_H

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;
using std::vector;

// Tests

vector<double> func0(vector<double> _point);

// Methods

void EulerExplicit(vector<double> _func(vector<double>), vector<double> _startPoint);
void EulerImplicit(vector<double> _func(vector<double>), vector<double> _startPoint);

// Linear equations

vector<double> gaussLinearSolve(vector<vector<double>> _A);
int mainElement(vector<vector<double>> _sourceColumn, int _row);

// Other

void vectorPrint(vector<double> _sourceVector);
void matrixPrint(vector<vector<double>> _sourceMatrix);
vector<double> Newtonsys(vector<double> _func(vector<double>), vector<double> _point, vector<double> _right, double _step);
vector<double> Newton(vector<double> _func(vector<double>), vector<double> _point, vector<double> _right, double _step);
vector<vector<double>> Jac(vector<double> _func(vector<double>),vector<double> _point, double _step);

#endif //NUMERICAL_ANALYSIS_ORDINARYDE_H
