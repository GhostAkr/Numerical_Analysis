//
// Created by ighos on 08.02.2019.
//

#include "../include/OrdinaryDE.h"

vector<double> func0(vector<double> _point) {
    // Constants
    double k = 20.0;
    double m = 0.3;
    if (_point.size() != 2) {  // Exception
        cout << "Incorrect point" << endl;
        return {};
    }
    vector<double> result(2);
    result[0] = _point[1];
    result[1] = -k / m * _point[0];
    return result;
}

void vectorPrint(vector<double> _sourceVector) {
    for (int i = 0; i < _sourceVector.size(); ++i) {
        cout << _sourceVector[i] << endl;
    }
    cout << endl;
}

void EulerExplicit(vector<double> _func(vector<double>), vector<double> _startPoint) {
    std::ofstream fOut("../data/solution.dat");
    if (!fOut) {  // Exception
        cout << "Error while opening file" << endl;
        return;
    }
    double step = 0.01;
    double maxMesh = 500;
    //fOut << maxMesh << endl;
    size_t nOfVars = _startPoint.size();
    //cout << nOfVars << endl;
    for (int i = 1; i < maxMesh; ++i) {
        vector<double> point(nOfVars);
        vector<double> f = _func(_startPoint);
        for (int j = 0; j < nOfVars; ++j) {
            point[j] = _startPoint[j] + step * f[j];
            fOut << point[j] << " ";
        }
        fOut << endl;
        _startPoint = point;
    }
}

//vector<vector<double>> Jac(vector<double> _func(vector<double>), vector<double> _point, vector<double> _right, double _step) {
//    vector<vector<double>> res(2, vector<double>(2));
//    double t = 0.00001;
//    vector<double> p1 = _point;
//    p1[0] += t;
//    vector<double> p2 = _point;
//    p2[1] += t;
//    vector<double> f1(2);
//    vector<double> f2(2);
//    vector<double> f3(2);
//    f1 = _func(p1);
//    f2 = _func(p2);
//    f3 = _func(_point);
//    //f3[0] -= -_right[0];
//    //f3[1] -= -_right[1];
//    res[0][0] = 1 - _step * (f1[0] - f3[0]) / t;
//    res[0][1] = -_step * (f2[0] - f3[0]) / t;
//    res[1][0] = -_step * (f1[1] - f3[1]) / t;
//    res[1][1] = 1 - _step * (f2[1] - f3[1]) / t;
//
//    return res;
//}

//vector<double> Newtonsys(vector<double> _func(vector<double>), vector<double> _point, vector<double> _right, double _step) {
//    vector<double> res(2);
//    vector<vector<double>> J;
//    vector<double> f(2);
//    double div = 0.0, norm = 0.0;
//    double eps = 0.001;
//    int iteration = 0;
//    do {
//        J = Jac(_func, _point, _right, _step);
//        div = 1 / (J[0][0] * J[1][1] - J[1][0] * J[0][1]);
//        f = _func(_point);
//        res[0] = _point[0] - div * (J[1][1] * (_point[0] -_step * f[0] - _right[0]) - J[0][1] * (_point[1] -_step * f[1] - _right[1]));
//        res[1] = _point[1] - div * (J[1][0] * (_point[0] -_step * f[0] - _right[0]) - J[0][0] * (_point[1] -_step * f[1] - _right[1]));
//        norm = sqrt((res[0] - _point[0])*(res[0] - _point[0]) + (res[1] - _point[1])*(res[1] - _point[1]));
//        iteration++;
//        _point[0] = res[0];
//        _point[1] = res[1];
//        if (iteration > 100) {
//            cout << "To many iterations " << endl;
//            vectorPrint(res);
//            break;
//        }
//    } while (norm > eps);
//
//    return res;
//}

void EulerImplicit(vector<double> _func(vector<double>), vector<double> _startPoint) {
    std::ofstream fOut("../data/solution.dat");
    if (!fOut) {  // Exception
        cout << "Error while opening file" << endl;
        return;
    }
    double step = 0.01;
    double maxMesh = 1000;
    //fOut << maxMesh << endl;
    size_t nOfVars = _startPoint.size();
    //cout << nOfVars << endl;
    vector<double> p(nOfVars, 0);
    p[0] = 0.9934;
    p[1] = -0.6623;
    for (int i = 1; i < maxMesh; ++i) {
        vector<double> point(nOfVars);
        point = Newton(_func, p, _startPoint, step);
        for (int j = 0; j < nOfVars; ++j) {
            fOut << point[j] << " ";
        }
        fOut << endl;
        _startPoint = point;
    }
}

vector<double> gaussLinearSolve(vector<vector<double>> _A) {
    double eps = 1e-14;  // For comparing with 0.0
    size_t rowsA = _A.size();
    size_t colsA = rowsA + 1;
    // First part
    for (int k = 0; k < colsA - 1; ++k) {
        int mainElem = mainElement(_A, k);
        _A[mainElem].swap(_A[k]);
        for (int i = k + 1; i < rowsA; ++i) {
            double coeff = _A[i][k] / _A[k][k];
            for (int j = k + 1; j < colsA; ++j) {
                _A[i][j] -= _A[k][j] * coeff;
                _A[i][k] = 0.0;
            }
        }
    }
//    if (!onlyDesitionCheck(_A)) {
//        cout << "Linear system has infinite number of solutions or hasn't it at all" << endl;
//        return NULL;
//    }
    // Second part
    vector<double> result(rowsA, 0);
    for (int i = rowsA - 1; i >= 0; --i) {
        double leftSum = 0.0;
        for (int j = rowsA - 1; j >= i; --j) {
            leftSum += _A[i][j] * result[j];
        }
        result[i] = (_A[i][colsA - 1] - leftSum) / _A[i][i];
        if (fabs(result[i]) < eps) {
            result[i] = 0.0;
        }
    }
    return result;
}

int mainElement(vector<vector<double>> _sourceColumn, int _row) {
    int indexOfMainElem = _row;
    int rows = _sourceColumn.size();
    double mainElem = _sourceColumn[_row][_row];
    for (int k = _row; k < rows; ++k) {
        if (_sourceColumn[k][_row] > mainElem) {
            mainElem = _sourceColumn[k][_row];
            indexOfMainElem = k;
        }
    }
    return indexOfMainElem;
}

vector<double> Newton(vector<double> _func(vector<double>), vector<double> _point, vector<double> _right, double _step) {
    size_t nOfVars = _point.size();
    vector<double> result(nOfVars);
    vector<double> fValue = _func(_point);
    vector<double> fNewton (nOfVars, 0);
    double norm = 0.0;
    double eps = 0.0;
    int iteration = 0;
    do {
        iteration++;
        // Creating of matrix for linear system
        vector<vector<double>> linearSystem = Jac(_func, _point, _step);
        for (int i = 0; i < nOfVars; ++i) {
            fNewton[i] = _point[i] - _step * fValue[i] - _right[i];
        }
        for (int i = 0; i < nOfVars; ++i) {
            linearSystem[i].push_back(-fNewton[i]);
        }
        // Solving linear equation
        vector<double> correctionVector = gaussLinearSolve(linearSystem);
        // Next point
        vector<double> nextPoint(nOfVars, 0);
        for (int i = 0; i < nOfVars; ++i) {
            nextPoint[i] = _point[i] + correctionVector[i];
        }
        // Counting norm
        double sum = 0.0;
        for (int i = 0; i < nOfVars; ++i) {
            sum += (nextPoint[i] - _point[i]) * (nextPoint[i] - _point[i]);
        }
        norm = sqrt(sum);
        // Assigning next point
        _point = nextPoint;
        if (iteration >= 1) {
            //cout << "Too many iterations" << endl;
            break;
        }
    } while (norm > eps);
    return _point;
}

vector<vector<double>> Jac(vector<double> _func(vector<double>),vector<double> _point, double _step) {
    size_t nOfVars = _point.size();
    vector<vector<double>> result(nOfVars, vector<double> (nOfVars, 0));
    double shift = 1e-7;
    for (int i = 0; i < nOfVars; ++i) {
        for (int j = 0; j < nOfVars; ++j) {
            vector<double> shiftedPoint = _point;
            shiftedPoint[j] += shift;
            if (i == j) {
                result[i][j] = 1 - _step * (_func(shiftedPoint)[i] - _func(_point)[i]) / shift;
            } else {
                result[i][j] = - _step * (_func(shiftedPoint)[i] - _func(_point)[i]) / shift;
            }
        }
    }
    return result;
}

void matrixPrint(vector<vector<double>> _sourceMatrix) {
    size_t nOfRows = _sourceMatrix.size();
    size_t nOfCols = _sourceMatrix[0].size();
    for (int i = 0; i < nOfRows; ++i) {
        for (int j = 0; j < nOfCols; ++j) {
            cout << _sourceMatrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}
