//
// Created by ighos on 16.05.2019.
//

#include "../include/PoissonEquation.h"

vector<vector<double>> mesh2D(double _L1, double _L2, double _h1, double _h2) {
    int N1 = round(_L1 / _h1 + 1);
    int N2 = round(_L2 / _h2 + 1);
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
        case 2: {
            return 0;
        }
        case 3: {
            return 4;
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
        case 2: {
            result[0] = 1.0;
            result[1] = 1.0;
            return  result;
        }
        case 3: {
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
                    result[i][j] = 3.0;
                    if (i == 0) {
                        result[i][j] = -1.0;
                    }
                    if (i == result.size() - 1) {
                        result[i][j] = -1.0;
                    }
                    if (j == 0) {
                        result[i][j] = -1.0;
                    }
                    if (j == result.size() - 1) {
                        result[i][j] = -1.0;
                    }
                }
            }
            return result;
        }
        case 2: {
            for (int i = 0; i < result.size(); ++i) {
                for (int j = 0; j < result[0].size(); ++j) {
                    result[i][j] = 1.0;
                }
            }
            return result;
        }
        case 3: {
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
    vector<double> result (3, 0);
    switch (_testNum) {
        case 1: {
            result[0] = 1.0;
            result[1] = 1.0;
            result[2] = 0;
            return result;
        }
        case 2: {
            result[0] = 1 + _x2;
            result[1] = 1 + _x2;
            result[2] = 0;
            return result;
        }
        case 3: {
            result[0] = 0;
            result[1] = -2.0;
            result[2] = 1.0;
            return result;
        }
        default: {
            cout << "Exception in leftBorder" << endl;
            return {};
        }
    }
}

vector<double> border2(double _x1, int _testNum) {
    vector<double> result (3, 0);
    switch (_testNum) {
        case 1: {
            result[0] = 1.0;
            result[1] = 1.0;
            result[2] = 0;
            return result;
        }
        case 2: {
            result[0] = 1.0;
            result[1] = 1.0;
            result[2] = 1.0;
            return result;
        }
        case 3: {
            result[0] = _x1 * _x1;
            result[1] = 1.0 + _x1 * _x1;
            result[2] = 0;
            return result;
        }
        default: {
            cout << "Exception in rightBorder" << endl;
            return {};
        }
    }
}

double exactSolutionPois(double _x1, double _x2, int _testNum) {
    switch (_testNum) {
        case 1: {
            return 1.0;
        }
        case 2: {
            return 1.0 + _x2;
        }
        case 3: {
            return _x1 * _x1 + _x2 * _x2;
        }
        default: {
            cout << "Exception in exactSolutionPois" << endl;
            return -1.0;
        }
    }
}

void switchDirectionsScheme(string _path, double _t, double _h1, double _h2, int _testNum) {
    std::ofstream fOut(_path);
    if (!fOut) {
        cout << "Error while opening file" << endl;
        return;
    }
    double nullEps = 1e-8;
    // Constants
    vector<double> L = area(_testNum);
    int N1 = round((L[0] / _h1) + 1);
    int N2 = round((L[1] / _h2) + 1);
    double hh1 = _h1 * _h1;
    double hh2 = _h2 * _h2;
    double A1 = 1.0 / hh1;
    double B1 = -2.0 * (1.0 / hh1 + 1.0 / _t);
    double C1 = 1.0 / hh1;
    double A2 = 1.0 / hh2;
    double B2 = -2.0 * (1.0 / hh2 + 1.0 / _t);
    double C2 = 1.0 / hh2;
    double T = 2;
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
    while (true) {
        vector<vector<double>> buffHalfLayer = mesh2D(L[0], L[1], _h1, _h2);
        // First system
        for (int j = 1; j < N2 - 1; ++j) {
            vector<double> borders = border1(j * _h2, _testNum);
            vector<vector<double>> triMatrix1 (N1 - 2, vector<double> (4, 0));
            triMatrix1[0][1] = B1;
            triMatrix1[0][2] = C1;
            if (fabs(borders[2] - 1.0) < nullEps) {
                triMatrix1[0][1] += A1;
            }
            for(int i = 1; i <= triMatrix1.size() - 2; ++i) {
                triMatrix1[i][0] = A1;
                triMatrix1[i][1] = B1;
                triMatrix1[i][2] = C1;
            }
            triMatrix1[triMatrix1.size() - 1][0] = A1;
            triMatrix1[triMatrix1.size() - 1][1] = B1;
            if (fabs(borders[2] - 1.0) < nullEps) {
                triMatrix1[triMatrix1.size() - 1][1] += C1;
            }
            for(int i = 0; i < triMatrix1.size(); ++i) {
                triMatrix1[i][3] = -F(buffLayer, i + 1, j, _t, _testNum, _h1, _h2);
                if (i == 0) {
                    if (fabs(borders[2] - 1.0) < nullEps) {
                        triMatrix1[i][3] += A1 * borders[0] * _h1;
                    } else {
                        triMatrix1[i][3] -= A1 * borders[0];
                    }
                }
                if (i == triMatrix1.size() - 1) {
                    if (fabs(borders[2] - 1.0) < nullEps) {
                        triMatrix1[i][3] -= C1 * borders[1] * _h1;
                    } else {
                        triMatrix1[i][3] -= C1 * borders[1];
                    }
                }
            }
            vector<double> halfSolution = tridiagonalLinearSolve(triMatrix1);
            if (fabs(borders[2] - 1.0) < nullEps) {
                halfSolution.insert(halfSolution.begin(), halfSolution[0] - borders[0] * _h1);
            } else {
                halfSolution.insert(halfSolution.begin(), borders[0]);
            }
            if (fabs(borders[2] - 1.0) < nullEps) {
                halfSolution.push_back(halfSolution[halfSolution.size() - 1] + borders[1] * _h1);
            } else {
                halfSolution.push_back(borders[1]);
            }
            for (int k = 0; k < buffHalfLayer.size(); ++k) {
                buffHalfLayer[k][j] = halfSolution[k];
            }
        }
        for(int i = 0; i < buffHalfLayer.size(); ++i) {
            vector<double> borders = border2(i * _h1, _testNum);
            if (fabs(borders[2] - 1.0) < nullEps) {
                buffHalfLayer[i][0] = buffHalfLayer[i][1] - borders[0] * _h2;
            } else {
                buffHalfLayer[i][0] = borders[0];
            }
            if (fabs(borders[2] - 1.0) < nullEps) {
                buffHalfLayer[i][buffHalfLayer.size() - 1] = buffHalfLayer[i][buffHalfLayer.size() - 2] + borders[1] * _h2;
            } else {
                buffHalfLayer[i][buffHalfLayer.size() - 1] = borders[1];
            }
        }
        vector<vector<double>> buffSecondHalfLayer = mesh2D(L[0], L[1], _h1, _h2);
        // Second system
        for (int i = 1; i < N1 - 1; ++i) {
            vector<double> borders = border2(i * _h1, _testNum);
            vector<vector<double>> triMatrix2 (N2 - 2, vector<double> (4, 0));
            triMatrix2[0][1] = B2;
            triMatrix2[0][2] = C2;
            if (fabs(borders[2] - 1.0) < nullEps) {
                triMatrix2[0][1] += A2;
            }
            for(int j = 1; j <= triMatrix2.size() - 2; ++j) {
                triMatrix2[j][0] = A2;
                triMatrix2[j][1] = B2;
                triMatrix2[j][2] = C2;
            }
            triMatrix2[triMatrix2.size() - 1][0] = A2;
            triMatrix2[triMatrix2.size() - 1][1] = B2;
            if (fabs(borders[2] - 1.0) < nullEps) {
                triMatrix2[triMatrix2.size() - 1][1] += C2;
            }
            for(int j = 0; j < triMatrix2.size(); ++j) {
                triMatrix2[j][3] = -F1(buffHalfLayer, i, j + 1, _t, _testNum, _h1, _h2);
                if (j == 0) {
                    if (fabs(borders[2] - 1.0) < nullEps) {
                        triMatrix2[j][3] += A2 * borders[0] * _h2;
                    } else {
                        triMatrix2[j][3] -= A2 * borders[0];
                    }
                }
                if (j == triMatrix2.size() - 1) {
                    if (fabs(borders[2] - 1.0) < nullEps) {
                        triMatrix2[j][3] -= C2 * borders[1] * _h2;
                    } else {
                        triMatrix2[j][3] -= C2 * borders[1];
                    }
                }
            }
            vector<double> halfSolution = tridiagonalLinearSolve(triMatrix2);
            if (fabs(borders[2] - 1.0) < nullEps) {
                halfSolution.insert(halfSolution.begin(), halfSolution[0] - borders[0] * _h1);
            } else {
                halfSolution.insert(halfSolution.begin(), borders[0]);
            }
            if (fabs(borders[2] - 1.0) < nullEps) {
                halfSolution.push_back(halfSolution[halfSolution.size() - 1] + borders[1] * _h1);
            } else {
                halfSolution.push_back(borders[1]);
            }
            buffSecondHalfLayer[i] = halfSolution;
        }
        for(int j = 0; j < buffSecondHalfLayer[0].size(); ++j) {
            vector<double> borders = border1(j * _h2, _testNum);
            if (fabs(borders[2] - 1.0) < nullEps) {
                buffSecondHalfLayer[0][j] = buffSecondHalfLayer[1][j] - borders[0] * _h1;
            } else {
                buffSecondHalfLayer[0][j] = borders[0];
            }
            if (fabs(borders[2] - 1.0) < nullEps) {
                buffSecondHalfLayer[buffSecondHalfLayer.size() - 1][j] = buffSecondHalfLayer[buffSecondHalfLayer.size() - 2][j] + borders[1] * _h1;
            } else {
                buffSecondHalfLayer[buffSecondHalfLayer.size() - 1][j] = borders[1];
            }
        }
        int isContinue = compareLayers(buffLayer, buffSecondHalfLayer, _t);
        buffLayer = buffSecondHalfLayer;
        for (int i = 0; i < buffLayer.size(); ++i) {
            for (int j = 0; j < buffLayer[0].size(); ++j) {
                fOut << buffLayer[i][j] << " ";
            }
        }
        fOut << endl;
        currentTime += _t;
        if (isContinue) {
            break;
        }
    }
    cout << "Time = " << currentTime << endl;
    double max = -1.0;
    for (int i = 0; i < buffLayer.size(); ++i) {
        for (int j = 0; j < buffLayer[0].size(); ++j) {
            double diff = fabs(buffLayer[i][j] - exactSolutionPois(i * _h1, j * _h2, _testNum));
//            cout << "Diff on [i][j] = " << "[" << i << "][" << j << "] = " << diff << endl;
            if (diff > max) {
                max = diff;
            }
        }
    }

    cout << "Error = " << max << endl;
    fOut.close();
}

int compareLayers(vector<vector<double>> _previous, vector<vector<double>> _next, double _t) {
    if (_previous.size() != _next.size() || _previous[0].size() != _next.size()) {
        cout << "Exception in compareLayers" << endl;
        return -1;
    }
    double max = 0.0;
    double eps = 1e-5;
    for (int i = 0; i < _previous.size(); ++i) {
        for (int j = 0; j < _next.size(); ++j) {
            double diff = fabs(_previous[i][j] - _next[i][j]);
            if (diff > max) {
                max = diff;
            }
        }
    }
    if (max < _t * eps) {
        return 1;
    } else {
        return 0;
    }
}
