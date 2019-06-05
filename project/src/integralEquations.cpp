//
// Created by ighos on 05.06.2019.
//

#include "../include/integralEquations.h"

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
