//
// Created by ighos on 08.02.2019.
//

#ifndef NUMERICAL_ANALYSIS_ORDINARYDE_H
#define NUMERICAL_ANALYSIS_ORDINARYDE_H

#define STEP 0.00025
#define MESH 800

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;
using std::vector;

// TODO: CLOSE ALL FILE ITERATORS

// Tests

vector<double> func0(vector<double> _point);
vector<double> func1(vector<double> _point);
vector<double> func2(vector<double> _point);
vector<double> func3(vector<double> _point);
vector<double> funcVar1(vector<double> _point);
vector<double> funcOwn(vector<double> _point);

// Methods

void EulerExplicit(vector<double> _func(vector<double>), vector<double> _startPoint, bool _isRungeRule);
void EulerImplicit(vector<double> _func(vector<double>), vector<double> _startPoint);
void RungeKutta(vector<double> _func(vector<double>), vector<double> _startPoint, bool _isRungeRule);
void Symmetric(vector<double> _func(vector<double>), vector<double> _startPoint);
void AdamsBashfort(vector<double> _func(vector<double>), vector<double> _startPoint);
void PredCorr(vector<double> _func(vector<double>), vector<double> _startPoint);

// Methods with returning values

vector<double> EulerExplicitReturn(vector<double> _func(vector<double>), vector<double> _startPoint, int _nOfIterations, double _step);
vector<double> RungeKuttaReturn(vector<double> _func(vector<double>), vector<double> _startPoint, int _nOfIterations, double _step);

// Linear equations

vector<double> gaussLinearSolve(vector<vector<double>> _A);
int mainElement(vector<vector<double>> _sourceColumn, int _row);

// Real solutions

vector<double> real0(double _step, int _iteration);

// Other

void vectorPrint(vector<double> _sourceVector);
void matrixPrint(vector<vector<double>> _sourceMatrix);
vector<double> Newtonsys(vector<double> _func(vector<double>), vector<double> _point, vector<double> _right, double _step);
vector<double> Newton(vector<double> _func(vector<double>), vector<double> _point, vector<double> _right, double _step);
vector<vector<double>> Jac(vector<double> _func(vector<double>),vector<double> _point, double _step);
vector<double> multV(vector<double> v, double a);
vector<double> vplus(vector<double> v, vector<double> w);
vector<double> vminus(vector<double> v, vector<double> w);
double error(vector<double> _realSolution, vector<double> _apprSolution);
double normInfVect(vector<double> _vect);

#endif //NUMERICAL_ANALYSIS_ORDINARYDE_H
