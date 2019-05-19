//
// Created by ighos on 16.05.2019.
//

#ifndef NUMERICAL_ANALYSIS_POISSONEQUATION_H
#define NUMERICAL_ANALYSIS_POISSONEQUATION_H

#include "../include/OrdinaryDE.h"
#include "../include/HeatEquation.h"

#include <iostream>
#include <vector>
#include <fstream>

using std::cout;
using std::endl;
using std::vector;
using std::string;

//Scheme
void switchDirectionsScheme(string _path, double _t, double _h1, double _h2, int _testNum);

// Task values
double fFunction(double _x1, double _x2, int _testNum);
vector<double> area(int _testNum);
vector<vector<double>> initLayer(double _L1, double _L2, double _h1, double _h2, int _testNum);

//Service functions
double lambda1(vector<vector<double>> _layer, double _h1, int _i, int _j);
double lambda2(vector<vector<double>> _layer, double _h2, int _i, int _j);
double F(vector<vector<double>> _layer, int _i, int _j, double _t, int _testNum, double _h1, double _h2);
double F1(vector<vector<double>> _layer, int _i, int _j, double _t, int _testNum, double _h1, double _h2);

// Other
vector<vector<double>> mesh2D(double _L1, double _L2, double _h1, double _h2);

#endif //NUMERICAL_ANALYSIS_POISSONEQUATION_H
