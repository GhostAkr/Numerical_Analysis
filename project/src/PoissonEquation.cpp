//
// Created by ighos on 16.05.2019.
//

#include "../include/PoissonEquation.h"

vector<vector<double>> mesh2D(double _L1, double _L2, double _h1, double _h2) {
    int N1 = _L1 / _h1 + 2;
    int N2 = _L2 / _h2 + 2;
    vector<vector<double>> result (N1, vector<double> (N2, 0));
    return result;
}

double lambda1(vector<vector<double>> _layer, double _h1, int _i, int _j) {
    double result = (_layer[_i + 1][_j] - 2.0 * _layer[_i][_j] + _layer[_i - 1][_j]) / (_h1 * _h1);
    return result;
}

double lambda2(vector<vector<double>> _layer, double _h2, int _i, int _j) {
    double result = (_layer[_i][_j + 1] - 2.0 * _layer[_i][_j] + _layer[_i][_j - 1]) / (_h2 * _h2);
    return result;
}

double fFunction(double _x1, double _x2, int _testNum) {
    switch (_testNum) {
        case 1: {
            return 0;
        }
        default: {
            cout << "Exception in fFunction" << endl;
            return -1;
        }
    }
}

vector<double> area(int _testNum) {
    vector<double> result (2, 0);
    switch (_testNum) {
        case 1: {
            result[0] = 1.0;
            result[1] = 1.0;
            return  result;
        }
        default: {
            cout << "Exception in area" << endl;
            return {};
        }
    }
}

vector<vector<double>> initLayer(double _L1, double _L2, double _h1, double _h2, int _testNum) {
    vector<vector<double>> result = mesh2D(_L1, _L2, _h1, _h2);
    switch (_testNum) {
        case 1: {
            for (int i = 0; i < result.size(); ++i) {
                for (int j = 0; j < result[0].size(); ++j) {
                    result[i][j] = 1.0;
                }
            }
            return result;
        }
        default: {
            cout << "Exception in initLayer" << endl;
            return {};
        }
    }
}

double F(vector<vector<double>> _layer, int _i, int _j, double _t, int _testNum, double _h1, double _h2) {
    double y = _layer[_i][_j];
    double phi = fFunction(_i * _h1, _j * _h2, _testNum);
    double result = 2.0 / _t * y + lambda2(_layer, _h2, _i, _j) + phi;
    return result;
}

double F1(vector<vector<double>> _layer, int _i, int _j, double _t, int _testNum, double _h1, double _h2) {
    double y = _layer[_i][_j];
    double phi = fFunction(_i * _h1, _j * _h2, _testNum);
    double result = 2.0 / _t * y + lambda1(_layer, _h1, _i, _j) + phi;
    return result;
}

vector<double> border1(double _x2, int _testNum) {
    vector<double> result (2, 0);
    switch (_testNum) {
        case 1: {
            result[0] = 1.0;
            result[1] = 1.0;
            return result;
        }
        default: {
            cout << "Exception in leftBorder" << endl;
            return {};
        }
    }
}

vector<double> border2(double _x1, int _testNum) {
    vector<double> result (2, 0);
    switch (_testNum) {
        case 1: {
            result[0] = 1.0;
            result[1] = 1.0;
            return result;
        }
        default: {
            cout << "Exception in rightBorder" << endl;
            return {};
        }
    }
}

void switchDirectionsScheme(string _path, double _t, double _h1, double _h2, int _testNum) {
    std::ofstream fOut(_path);
    if (!fOut) {
        cout << "Error while opening file" << endl;
        return;
    }
    // Constants
    vector<double> L = area(_testNum);
    int N1 = (L[0] / _h1) + 2;
    int N2 = (L[1] / _h2) + 2;
    double hh1 = _h1 * _h1;
    double hh2 = _h2 * _h2;
    double A1 = 1.0 / hh1;
    double B1 = 2.0 * (1.0 / hh1 + 1.0 / _t);
    double C1 = 1.0 / hh1;
    double A2 = 1.0 / hh2;
    double B2 = 2.0 * (1.0 / hh2 + 1.0 / _t);
    double C2 = 1.0 / hh2;
    double T = 1.0;
    // Null time layer
    vector<vector<double>> buffLayer = initLayer(L[0], L[1], _h1, _h2, _testNum);
    for (int i = 0; i < buffLayer.size(); ++i) {
        for (int j = 0; j < buffLayer[0].size(); ++j) {
            fOut << buffLayer[i][j] << " ";
        }
    }
    fOut << endl;
    // Main algorithm
    double currentTime = 0;
    while (currentTime <= T) {
        vector<vector<double>> buffHalfLayer = mesh2D(L[0], L[1], _h1, _h2);
        // First system
        for (int j = 1; j < N1 - 1; ++j) {
            vector<vector<double>> triMatrix1 (N1 - 2, vector<double> (4, 0));
            triMatrix1[0][1] = B1;
            triMatrix1[0][2] = C1;
            for(int i = 1; i <= triMatrix1.size() - 2; ++i) {
                triMatrix1[i][0] = A1;
                triMatrix1[i][1] = B1;
                triMatrix1[i][2] = C1;
            }
            triMatrix1[triMatrix1.size() - 1][0] = A1;
            triMatrix1[triMatrix1.size() - 1][1] = B1;
            for(int i = 0; i < triMatrix1.size(); ++i) {
                triMatrix1[i][3] = -F(buffLayer, i, j, _t, _testNum, _h1, _h2);
            }
            vector<double> halfSolution = tridiagonalLinearSolve(triMatrix1);
            vector<double> borders = border1(j * _h2, _testNum);
            halfSolution.insert(halfSolution.begin(), borders[0]);
            halfSolution.push_back(borders[1]);
            buffHalfLayer[j] = halfSolution;
        }
        for(int j = 0; j < buffHalfLayer[0].size(); ++j) {
            vector<double> borders = border2(j * _h1, _testNum);
            buffHalfLayer[j][0] = borders[0];
            buffHalfLayer[j][buffHalfLayer.size() - 1] = borders[1];
        }
        // Second system
        for (int j = 1; j <= N2 - 1; ++j) {
            vector<vector<double>> triMatrix2 (N2 - 1, vector<double> (4, 0));
            triMatrix2[0][1] = B2;
            triMatrix2[0][2] = C2;
            for(int i = 1; i < N2 - 2; ++i) {
                triMatrix2[i][0] = A2;
                triMatrix2[i][1] = B2;
                triMatrix2[i][2] = C2;
            }
            triMatrix2[N2 - 2][0] = A2;
            triMatrix2[N2 - 2][1] = B2;
            for(int i = 0; i < N2 - 1; ++i) {
                triMatrix2[i][3] = -F1(buffHalfLayer, j, i + 1, _t, _testNum, _h1, _h2);
            }
            vector<double> halfSolution = tridiagonalLinearSolve(triMatrix2);
            vector<double> borders = border2(j * _h1, _testNum);
            halfSolution.insert(halfSolution.begin(), borders[0]);
            halfSolution.push_back(borders[1]);
            buffHalfLayer[j] = halfSolution;
        }
        for(int j = 0; j < buffHalfLayer[0].size(); ++j) {
            vector<double> borders = border1(j * _h2, _testNum);
            buffHalfLayer[j][0] = borders[0];
            buffHalfLayer[j][buffHalfLayer.size() - 1] = borders[1];
        }
        buffLayer = buffHalfLayer;
        for (int i = 0; i < buffLayer.size(); ++i) {
            for (int j = 0; j < buffLayer[0].size(); ++j) {
                fOut << buffLayer[i][j] << " ";
            }
        }
        fOut << endl;
        currentTime += _t;
    }

    fOut.close();
}
