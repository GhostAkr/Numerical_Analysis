//
// Created by ighos on 05.06.2019.
//

#include "../include/integralEquations.h"
#define pi 3.14159265358979

double core(double _x, double _s, int _testNum) {
    switch (_testNum) {
        case 1: {
            return 1 - _x * cos(_x * _s);
        }
        case 2: {
            return 1 - _x * cos(_x * _s);
        }
        default: {
            cout << "Exception in core" << endl;
            return -1.0;
        }
    }
}
double rightPart(double _x, int _testNum) {
    switch (_testNum) {
        case 1: {
            return 0.5 * (1 + sin(_x));
        }
        case 2: {
            return _x * _x + sqrt(_x);
        }
        default: {
            cout << "Exception in rightPart" << endl;
            return -1.0;
        }
    }
}
vector<double> limits(int _limitType) {
    vector<double> result(2, 0);
    switch (_limitType) {
        case 1: {
            result[0] = 0.0;
            result[1] = 1.0;
            return result;
        }
        case 2: {
            result[0] = 0.1;
            result[1] = 1.0;
            return result;
        }
        default: {
            cout << "Exception in limits" << endl;
            return {};
        }
    }
}

void quadMethod(int _testNum, int _limitType, int _N, string _path) {
    vector<double> limitsInt = limits(_limitType);
    double a = limitsInt[0];
    double b = limitsInt[1];
    double h = (b - a) / _N;
    // Creating and solving of linear system
    vector<vector<double>> equation (_N + 1, vector<double> (_N + 2, 0));
    for (int i = 0; i < equation.size(); ++i) {
        double x = a + h * i;
        for (int k = 0; k < equation[0].size() - 1; ++k) {
            double s = a + h * k;
            if (k == 0) {
                equation[i][k] = -0.25 * core(x, s, _testNum) * h;
                continue;
            }
            if (k == equation[0].size() - 2) {
                equation[i][k] = -0.25 * core(x, s, _testNum) * h;
                continue;
            }
            equation[i][k] = -0.5 * core(x, s, _testNum) * h;
        }
        equation[i][i] += 1.0;
        equation[i][equation[0].size() - 1] = rightPart(x, _testNum);
    }
    vector<double> solution = gaussLinearSolve(equation);
    // Writing to file
    std::ofstream fOut(_path);
    if (!fOut) {
        cout << "Error while opening file" << endl;
        return;
    }
    for (int i = 0; i < solution.size(); ++i) {
        fOut << a + i * h << " " << solution[i] << endl;
    }
    fOut.close();
}

void simpleMethod(int _testNum, int _limitType, int _N, string _path) {
    vector<double> limitsInt = limits(_limitType);
    double a = limitsInt[0];
    double b = limitsInt[1];
    double h = (b - a) / _N;
    vector<double> solution (_N + 1, 1);
    vector<double> buffSolution (_N + 1, 0);
    int iterations = 0;
    while (true) {
        iterations++;
        for (int i = 0; i < solution.size(); ++i) {
            double x = a + i * h;
            double sum = 0.0;
            for (int k = 0; k < solution.size(); ++k) {
                double s = a + k * h;
                if (k == 0) {
                    sum += 0.25 * core(x, s, _testNum) * h * solution[k];
                    continue;
                }
                if (k == solution.size() - 1) {
                    sum += 0.25 * core(x, s, _testNum) * h * solution[k];
                    continue;
                }
                sum += 0.5 * core(x, s, _testNum) * h * solution[k];
            }
            buffSolution[i] = rightPart(x, _testNum) + sum;
        }
        solution = buffSolution;
        if (iterations == 50) {
            break;
        }
    }
    // Writing to file
    std::ofstream fOut(_path);
    if (!fOut) {
        cout << "Error while opening file" << endl;
        return;
    }
    for (int i = 0; i < solution.size(); ++i) {
        fOut << a + i * h << " " << solution[i] << endl;
    }
    fOut.close();
}

double phi(double _x, int _i) {
    if (_i == 1) {
        return 1.0 - _x;
    }
    return pow(_x, 2 * _i - 1);
}

double psi(double _s, int _i) {
    if (_i == 1) {
        return 1.0;
    }
    double fact = 1.0;
    for (int j = 1; j <= 2 * _i - 2; ++j) {
        fact *= j;
    }
    return pow(_s, 2 * _i - 2) / fact;
}

void degenMethod(int _testNum, int _limitType, int _N, string _path) {
    vector<double> limitsInt = limits(_limitType);
    double a = limitsInt[0];
    double b = limitsInt[1];
    double h = (b - a) / _N;
    int m = 5;
    double hm = (b - a) / (m - 1);
    // Coefficients calculating
    vector<vector<double>> equationC (m, vector<double> (m + 1, 0));
    for (int i = 1; i <= m; ++i) {
        // beta calculation
        double beta = 0.0;
        for (int k = 0; k < _N + 1; ++k) {
            if (k == 0) {
                beta += 0.5 * h * psi(a + h * k, i) * rightPart(a + h * k, _testNum);
                continue;
            }
            if (k == _N) {
                beta += 0.5 * h * psi(a + h * k, i) * rightPart(a + h * k, _testNum);
                continue;
            }
            beta += h * psi(a + h * k, i) * rightPart(a + h * k, _testNum);
        }
        for (int j = 1; j <= m; ++j) {
            // alpha calculation
            double alpha = 0.0;
            for (int k = 0; k < _N + 1; ++k)  {
                if (k == 0) {
                    alpha += 0.5 * h * psi(a + h * k, i) * phi(a + h * k, j);
                    continue;
                }
                if (k == _N) {
                    alpha += 0.5 * h * psi(a + h * k, i) * phi(a + h * k, j);
                    continue;
                }
                alpha += h * psi(a + h * k, i) * phi(a + h * k, j);
            }
            equationC[i - 1][j - 1] = -0.5 * alpha;
        }
        equationC[i - 1][m] = beta;
        equationC[i - 1][i - 1] += 1.0;
    }
    vector<double> solveC = gaussLinearSolve(equationC);
    vector<double> solution (_N + 1, 0);
    for (int i = 0; i < solution.size(); ++i) {
        double sum = 0.0;
        for (int k = 1; k <= m; ++k) {
            sum += solveC[i - 1] * phi(a + i * h, k);
        }
        solution[i] = rightPart(a + i * h, _testNum) + 0.5 * sum;
    }
    // Writing to file
    std::ofstream fOut(_path);
    if (!fOut) {
        cout << "Error while opening file" << endl;
        return;
    }
    for (int i = 0; i < solution.size(); ++i) {
        fOut << a + i * h << " " << solution[i] << endl;
    }
    fOut.close();
}

vector<vector<double>>kMesh(int n) {
    vector<vector<double>> result(n, vector<double>(2, 0));
    for (int i = 0; i < result.size(); ++i) {

        result[i][0] = cos(2 * pi * (i - 0.5) / n);
        result[i][1] = sin(2 * pi * (i - 0.5) / n);
    }
    return result;
}

vector<vector<double>>cMesh(int n) {
    vector<vector<double>> result(n, vector<double>(2, 0));
    for (int i = 0; i < result.size(); ++i) {
        result[i][0] = cos(2 * pi * (i - 1) / n);
        result[i][1] = sin(2 * pi * (i - 1) / n);
    }
    return result;
}


vector<double>Qker(vector<double> k, vector<double> c) {
    vector<double> result(2, 0);
    double r;
    //cout << "Test" << endl;
    r = 2 * pi * ((k[0] - c[0]) * (k[0] - c[0]) + (k[1] - c[1]) * (k[1] - c[1]));
    //cout << r << endl;
//    cout << "k[1] = " << k[1] << endl;
//    cout << "k[0] = " << k[0] << endl;
    result[0] = -(k[1]-c[1])/r;
    result[1] = (k[0] - c[0])/r;

    return result;
}

void SIE(string _path, int n) {
    std::ofstream fOut(_path);
    if (!fOut) {
        cout << "Error while opening file" << endl;
        return;
    }
    double dl = 2 * pi / n;
    vector<vector<double>> k = kMesh(n); // массив к совпадает с массивом нормалей
    vector<vector<double>> c = cMesh(n);
//    cout << "k" << endl;
//    matrixPrint(k);
//    cout << "c" << endl;
//    matrixPrint(c);
    vector<vector<double>> M(n + 1, vector<double>(n + 2, 1));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            vector<double> tmp = Qker(k[i], c[j]);
//            cout << "tmp[0] = " << tmp[0] << endl;
//            cout << "tmp[1] = " << tmp[1] << endl;
            //cout << "R" << endl;
            M[i][j] = dl * (k[i][0] * tmp[0] + k[i][1] * tmp[1]);
        }
    }
    for (int i = 0; i < n; ++i) {
        M[i][n + 1] = k[i][1];
    }
    for (int i = 0; i < M[0].size() - 1; ++i) {
        M[M.size() - 1][i] = dl;
    }
    M[M.size() - 1][M[0].size() - 1] = 0;
    matrixPrint(M);
    vector<double> solution = gaussLinearSolve(M);
    //vectorPrint(solution);
    for (int i = 0; i < solution.size() - 1; ++i) {
        fOut << c[i][0] << " " << c[i][1] << " " << solution[i] << endl;
    }
    fOut.close();
}