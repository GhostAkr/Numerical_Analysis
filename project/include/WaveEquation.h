//
// Created by ighos on 21.04.2019.
//

#ifndef NUMERICAL_ANALYSIS_WAVEEQUATION_H
#define NUMERICAL_ANALYSIS_WAVEEQUATION_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#define Pi 3.14159265

using std::cout;
using std::endl;
using std::string;
using std::vector;

void crossScheme(double _t, double _h, string _path, int _testNum);
void crossSchemeDerivative(double _t, double _h, string _path, int _testNum);
double exactSolution(double _t, double _x, int _testNum);
void errorUpdate(vector<double> _layer, double* _error, double _t, double _h, int _testNum);
vector<double> meshDerivative(size_t _N, double _L);

// Condition functions
double gFunction(double _x, int _testNum);
double phiFunction(double _t, int _testNum);
double psiFunction(double _t, int _testNum);
double zeroFunction(double _x, int _testNum);
double zeroDerivative(double _x, int _testNum);
double muFunction(double _t, int _testNum);

#endif //NUMERICAL_ANALYSIS_WAVEEQUATION_H
