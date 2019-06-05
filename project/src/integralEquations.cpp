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
    vector<double> solution (_N, 1);
    vector<double> buffSolution (_N, 0);
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
