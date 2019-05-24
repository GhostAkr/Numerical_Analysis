//
// Created by ighos on 16.05.2019.
//

#include "../include/PoissonEquation.h"

vector<vector<double>> mesh2D(double _L1, double _L2, double _h1, double _h2) {
    int N1 = _L1 / _h1 + 1;
    int N2 = _L2 / _h2 + 1;
    vector<vector<double>> result (N1 + 1, vector<double> (N2 + 1, 0));
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
                    result[i][j] = 1000.0;
                    if (i == 0) {
                        result[i][j] = 1.0;
                    }
                    if (i == result.size() - 1) {
                        result[i][j] = 1.0;
                    }
                    if (j == 0) {
                        result[i][j] = 1.0;
                    }
                    if (j == result.size() - 1) {
                        result[i][j] = 1.0;
                    }
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
//    cout << "[i][j] = " << _layer[_i][_j] << endl;
//    cout << "l2 = " << lambda2(_layer, _h2, _i, _j) << endl;
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
    int N1 = round((L[0] / _h1) + 1);
    int N2 = round((L[1] / _h2) + 1);
//    cout << "N1 = " << N1 << endl;
//    cout << "N2 = " << N2 << endl;
//    int N1 = 11;
//    int N2 = 11;
//    cout << "N2 = " << N2 << endl;
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
    while (currentTime <= T) {
        cout << "buffLayer before" << endl;
        matrixPrint(buffLayer);
        vector<vector<double>> buffHalfLayer = mesh2D(L[0], L[1], _h1, _h2);
//        cout << "buffHalfLayer" << endl;
//        matrixPrint(buffHalfLayer);
        // First system
//        cout << "N2 = " << N2 << endl;
        for (int j = 1; j < N2 - 1; ++j) {
            vector<double> borders = border1(j * _h2, _testNum);
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
//            cout << "Size 1 = " << triMatrix1.size() << endl;
//            cout << "Size 1 = " << buffLayer.size() << endl;
//            cout << "Size 2 = " << buffLayer[0].size() << endl;
            for(int i = 0; i < triMatrix1.size(); ++i) {
//                cout << "j = " << j << endl;
//                cout << "i = " << i << endl;
                triMatrix1[i][3] = -F(buffLayer, i + 1, j, _t, _testNum, _h1, _h2);
                if (i == 0) {
                    triMatrix1[i][3] -= A1 * borders[0];
                }
                if (i == triMatrix1.size() - 1) {
                    triMatrix1[i][3] -= C1 * borders[1];
                }
            }
//            cout << "triMatrix1" << endl;
//            matrixPrint(triMatrix1);
//            cout << "Test" << endl;
            vector<double> halfSolution = tridiagonalLinearSolve(triMatrix1);
//            cout << "SSSSuka = " << halfSolution.size() << endl;
            halfSolution.insert(halfSolution.begin(), borders[0]);
            halfSolution.push_back(borders[1]);
//            cout << "halfSolution" << endl;
//            vectorPrint(halfSolution);
            for (int k = 0; k < buffHalfLayer.size(); ++k) {
                buffHalfLayer[k][j] = halfSolution[k];
            }
        }
//        cout << "buffHalfLayer" << endl;
//        matrixPrint(buffHalfLayer);

        for(int i = 0; i < buffHalfLayer.size(); ++i) {
            vector<double> borders = border2(i * _h1, _testNum);
            buffHalfLayer[i][0] = borders[0];
            buffHalfLayer[i][buffHalfLayer.size() - 1] = borders[1];
        }
        // Second system
//        cout << "N1 = " << N1 << endl;
        for (int i = 1; i < N1 - 1; ++i) {
            vector<double> borders = border2(i * _h1, _testNum);
            vector<vector<double>> triMatrix2 (N2 - 2, vector<double> (4, 0));
            triMatrix2[0][1] = B2;
            triMatrix2[0][2] = C2;
            for(int j = 1; j <= triMatrix2.size() - 2; ++j) {
                triMatrix2[j][0] = A2;
                triMatrix2[j][1] = B2;
                triMatrix2[j][2] = C2;
            }
            triMatrix2[triMatrix2.size() - 1][0] = A2;
            triMatrix2[triMatrix2.size() - 1][1] = B2;
//            cout << "Size 2 = " << triMatrix2.size() << endl;
//            cout << "SSSize 1 = " << buffHalfLayer.size() << endl;
//            cout << "SSSize 2 = " << buffHalfLayer[0].size() << endl;
            for(int j = 0; j < triMatrix2.size(); ++j) {
//                cout << "j = " << j << endl;
//                cout << "i = " << i << endl;
                triMatrix2[j][3] = -F1(buffHalfLayer, i, j + 1, _t, _testNum, _h1, _h2);
                if (j == 0) {
                    triMatrix2[j][3] -= A2 * borders[0];
                }
                if (j == triMatrix2.size() - 1) {
                    triMatrix2[j][3] -= C2 * borders[1];
                }
            }
            vector<double> halfSolution = tridiagonalLinearSolve(triMatrix2);
            halfSolution.insert(halfSolution.begin(), borders[0]);
            halfSolution.push_back(borders[1]);
            buffHalfLayer[i] = halfSolution;
        }
        for(int j = 0; j < buffHalfLayer[0].size(); ++j) {
            vector<double> borders = border1(j * _h2, _testNum);
            buffHalfLayer[j][0] = borders[0];
            buffHalfLayer[j][buffHalfLayer.size() - 1] = borders[1];
        }
        buffLayer = buffHalfLayer;
//        cout << "buffLayer after" << endl;
//        matrixPrint(buffLayer);
//        cout << "buffLayer" << endl;
//        matrixPrint(buffLayer);
//        cout << "Size = " << buffLayer[1].size() << endl;
//        cout << "qqepta = " << buffLayer[1][buffLayer[0].size() - 1] << endl;
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
