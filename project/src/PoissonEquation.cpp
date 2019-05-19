//
// Created by ighos on 16.05.2019.
//

#include "../include/PoissonEquation.h"

vector<vector<double>> mesh2D(double _L1, double _L2, double _h1, double _h2) {
    int N1 = _L1 / _h1 + 1;
    int N2 = _L2 / _h2 + 1;
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

void switchDirectionsScheme(string _path, double _t, double _h1, double _h2, int _testNum) {
    std::ofstream fOut(_path);
    if (!fOut) {
        cout << "Error while opening file" << endl;
        return;
    }
    // Constants
    vector<double> L = area(_testNum);
    // Null time layer
    vector<vector<double>> u0 = initLayer(L[0], L[1], _h1, _h2, _testNum);
    fOut.close();
}
