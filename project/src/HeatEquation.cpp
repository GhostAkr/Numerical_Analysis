//
// Created by ighos on 19.03.2019.
//

#include "../include/HeatEquation.h"

vector<double> mesh(size_t _N, double _L) {
    double h = _L / _N;
    vector<double> resultMesh(_N + 1, -1);
    return resultMesh;
}

vector<double> Mixed(vector<double> _V0, double sigma, double h, double t) {
    vector<double> _V1(_V0.size());  // Next time layer
    // Coefficients for tridiagonal matrix algorithm
    vector<double> FA;
    vector<double> FB;
    double A, B, C, F, a0, a1, an;
    //int s = _V0.size() - 1;
    // Constants
    double c = 1; double ro = 1; double alpha = 2; double beta = 0.5; double gamma = 3;
    double u0 = 0.1;
    // Time layers parameters
    double T = 1;
    double N0 = T / t;
    // Coordinate parameters
    double L = 1;
    double N = L / h;
    _V1[0] = u0;  // Left border
    // Calculating A_1, B_1, C_1, F_1
    a0 = K1(h - 0.5 * h, L);  // a_1
    a1 = K1(2 * h - 0.5 * h, L);  // a_2
    A = sigma * a0 / h;
    B = sigma * a1 / h;
    C = A + B + c * ro * h / t;
    F = c * ro * _V0[1] * h / t + (1 - sigma)*(a1*(_V0[2] - _V0[1]) - a0 * (_V0[1] - _V0[0])) / h;
    FA.push_back(B / C);
    FB.push_back(F / C);
    // Calculating coefficients for tridiagonal matrix algorithm
    for (int i = 2; i <= N - 2; i++) {
        a0 = K1(i * h - 0.5 * h, L);
        a1 = K1((i + 1) * h - 0.5 * h, L);
        A = sigma * a0 / h;
        B = sigma * a1 / h;
        C = A + B + c * ro * h / t;
        F = c * ro * _V0[i] * h / t + (1 - sigma)*(a1*(_V0[i + 1] - _V0[i]) - a0 * (_V0[i] - _V0[i - 1])) / h;
        FA.push_back(-B / (A * FA[i - 1] + C));
        FB.push_back((F - A * FB[i - 1]) / (A * FA[i - 1] + C));
        if (i == N - 2) {
            i++;
            a0 = K1(i * h - 0.5 * h, L);
            a1 = K1((i + 1) * h - 0.5 * h, L);
            A = sigma * a0 / h;
            B = sigma * a1 / h;
            C = A + B + c * ro * h / t;
            F = c * ro * _V0[i] * h / t + (1 - sigma)*(a1*(_V0[i + 1] - _V0[i]) - a0 * (_V0[i] - _V0[i - 1])) / h;
            _V1[_V1.size() - 2] = (F - A * FB[FB.size() - 1]) / (C + A * FA[FA.size() - 1]);
        }
    }
    int k = 0;
    for (int i = N - 2; i >= 1; i--) {
        _V1[i] = _V1[i + 1] * FA[FA.size() - k - 1] + FB[FB.size() - k - 1];
        k++;
    }
    _V1[N] = u0;
    return _V1;
}

void explicitScheme(double _t, double _h) {
    std::ofstream fOut("../data/out.dat");
    if (!fOut) {
        cout << "Error while opening file" << endl;
        return;
    }
    // Time layers parameters
    double T = 1;
    double N0 = T / _t;
    // Coordinate parameters
    double L = 1;
    double N = L / _h;
    // Constants
    double c = 1; double ro = 1; double alpha = 2; double beta = 0.5; double gamma = 3;
    double u0 = 0.1; double k1 = 1; double k2 = 0.1; double x1 = 1.0 / 3.0; double x2 = 2.0 / 3.0;
    // Algorithm
    vector<double> vPrev = mesh(N, L);
    // Start values
    for (int i = 0; i < vPrev.size(); ++i) {
        double x = i * _h;
        vPrev[i] = u0 + x * (L - x);
        fOut << vPrev[i] << " ";
    }
    fOut << endl;
    // Next layers
    for (int j = 0; j < N0; ++j) {
        vector<double> vNext = mesh(N, L);
        vNext[0] = u0;
        fOut << vNext[0] << " ";
        cout << "Current vector is" << endl;
        vectorPrint(vPrev);
        for (int i = 1; i < vNext.size() - 1; ++i) {
            double coeff = _t / (_h * c * ro);
            double xi = i * _h;  // x_i
            double xi1 = (i - 1) * _h;  // x_{i - 1}
            double xi11 = (i + 1) * _h;  // x_{i + 1}
            double ai = 0.0;  // a_i
            double ai11 = 0.0;  // a_{i + 1}
            // Counting a_i
            double logArg = (k1 - k2) * xi - k1 * x2 + k2 * x1 / ((k1 - k2) * xi1 - k1 * x2 + k2 * x1);
            if (xi >= 0 && xi <= x1) {
                ai = _h * k1 / (xi - xi1);
            } else if (xi > x1 && xi < x2) {
                ai = _h * (k1 - k2) / ((x1 - x2) * log(fabs(logArg)));
            } else if (xi >= x2 && xi <= L) {
                ai = _h * k2 / (xi - xi1);
            } else {
                cout << "Exception while counting a_i. Index is out of range" << endl;
                return;
            }
            // Counting a_{i + 1}
            logArg = (k1 - k2) * xi11 - k1 * x2 + k2 * x1 / ((k1 - k2) * xi - k1 * x2 + k2 * x1);
            if (xi11 >= 0 && xi11 <= x1) {
                ai11 = _h * k1 / (xi11 - xi);
            } else if (xi11 > x1 && xi11 < x2) {
                ai11 = _h * (k1 - k2) / ((x1 - x2) * log(fabs(logArg)));
            } else if (xi11 >= x2 && xi11 <= L) {
                ai11 = _h * k2 / (xi11 - xi);
            } else {
                cout << "Exception while counting a_{i + 1}. Index is out of range" << endl;
                return;
            }
            vNext[i] = vPrev[i] + coeff * ((ai11 * (vPrev[i + 1] - vPrev[i]) / _h) - (ai * (vPrev[i] - vPrev[i - 1]) / _h));
            fOut << vNext[i] << " ";
        }
        vNext[vNext.size() - 1] = u0;
        fOut << vNext[vNext.size() - 1] << endl;
        vPrev = vNext;
    }
    fOut.close();
}

double K1(double _x, double _L) {
    double k1 = 1; double k2 = 1; double x1 = 1.0 / 3.0; double x2 = 2.0 / 3.0;
    if (_x >= 0 && _x <= x1) {
        return k1;
    } else if (_x > x1 && _x < x2) {
        return k1 * (_x - x2) / (x1 - x2) + k2 * (_x - x1) / (x2 - x1);
    } else if (_x >= x2 && _x <= _L) {
        return k2;
    } else {
        cout << "Exception in K1 function" << endl;
        return -1.0;
    }
}

void Fill(string path) {
    std::ofstream fOut(path);
    if (!fOut) {
        cout << "Error while opening file" << endl;
        return;
    }
    double u0 = 0.1;
    double h = 0.0001;
    double t = 0.12 ;
    double T = 1;
    double N0 = T / t;
    double L = 1;
    double N = L / h;
    vector<double> V0 = mesh(N, L);
    for (int i = 0; i < V0.size(); ++i) {
        double x = i * h;
        V0[i] = u0 + x * (L - x);
        fOut << V0[i] << " ";
    }
    fOut << endl;
    for (double i = 0; i <= T; i += t) {
        V0 = Mixed(V0, 0 , h, t);
        for (int i = 0; i < V0.size(); ++i) {
            fOut << V0[i] << " ";
        }
        fOut << endl;
    }
    fOut.close();
}

vector<double> tridiagonalLinearSolve(const vector<vector<double>> _matrix) {
    int rows = _matrix.size();
    vector<double> a (rows - 1);
    vector<double> b (rows - 1);
    a[0] = -_matrix[0][2] / _matrix[0][1];
    b[0] = _matrix[0][3] / _matrix[0][1];
    vector<double> result (rows, 0);
    for (int i = 1; i < rows - 1; ++i) {
        double A = _matrix[i][0];
        double C = _matrix[i][1];
        double B = _matrix[i][2];
        double F = _matrix[i][3];
        a[i] = -B / (A * a[i - 1] + C);
        b[i] = (F - A * b[i - 1]) / (A * a[i - 1] + C);
    }
    double F = _matrix[rows - 1][3];
    double A = _matrix[rows - 1][0];
    double C = _matrix[rows - 1][1];
    result[rows - 1] = (F - A * b[rows - 2]) / (A * a[rows - 2] + C);
    for (int i = rows - 2; i >= 0; --i) {
        result[i] = a[i] * result[i + 1] + b[i];
    }
    return result;
}

void mixedScheme(double _sigma, double _h, double _t, string _path) {
    std::ofstream fOut(_path);
    if (!fOut) {
        cout << "Error while opening file" << endl;
        return;
    }
    double max = -1.0;
    double L = 1;
    double N = L / _h;
    double T = 1;
    double N0 = T / _t;
    // Constants
    double c = 1; double ro = 1; double alpha = 2; double beta = 0.5; double gamma = 3;
    double u0 = 0.1; double k1 = 1; double k2 = 0.1; double x1 = 1.0 / 3.0; double x2 = 2.0 / 3.0;
    // First time layer
    vector<double> V0 = mesh(N, L);
    for (int i = 0; i < V0.size(); ++i) {
        double x = i * _h;
        V0[i] = u0 + x * (L - x);
        fOut << V0[i] << " ";
    }
    fOut << endl;
    vector<vector<double>> tridiagonalArgument (N - 1, vector<double> (4, 0));
    double currentTime = 0.0;
    while (currentTime <= T) {
        // Calculating A_1, B_1, C_1, F_1
        double a0 = K1(_h - 0.5 * _h, L);  // a_1
        double a1 = K1(2 * _h - 0.5 * _h, L);  // a_2
        double A = _sigma * a0 / _h;
        double B = _sigma * a1 / _h;
        double C = A + B + c * ro * _h / _t;
        double F = c * ro * V0[1] * _h / _t + (1 - _sigma)*(a1*(V0[2] - V0[1]) - a0 * (V0[1] - V0[0])) / _h;
        F += A * u0;
        // Filling of tridiagonalArgument
        tridiagonalArgument[0][1] = -C;
        tridiagonalArgument[0][2] = B;
        tridiagonalArgument[0][3] = -F;
        for (int i = 1; i < N - 1; ++i) {
            a0 = K1((i + 1) * _h - 0.5 * _h, L);
            a1 = K1((i + 2) * _h - 0.5 * _h, L);
            A = _sigma * a0 / _h;
            B = _sigma * a1 / _h;
            C = A + B + c * ro * _h / _t;
            F = c * ro * V0[i + 1] * _h / _t + (1 - _sigma)*(a1*(V0[i + 2] - V0[i + 1]) - a0 * (V0[i + 1] - V0[i])) / _h;
            if (i == 1) {
                F -= A * u0;
            }
            if (i == N - 2) {
                F += B * u0;
                tridiagonalArgument[i][0] = A;
                tridiagonalArgument[i][1] = -C;
                tridiagonalArgument[i][3] = -F;
                break;
            }
            tridiagonalArgument[i][0] = A;
            tridiagonalArgument[i][1] = -C;
            tridiagonalArgument[i][2] = B;
            tridiagonalArgument[i][3] = -F;
        }
        // Filling of next layer
        vector<double> nextLayer = tridiagonalLinearSolve(tridiagonalArgument);
        V0[0] = u0;  // Left border
        fOut << V0[0] << " ";
        for (int i = 1; i <= V0.size() - 2; ++i) {
            V0[i] = nextLayer[i - 1];
            fOut << V0[i] << " ";
        }
        V0[V0.size() - 1] = u0;  // Right border
        for (int i = 1; i < V0.size(); ++i) {
            double difference = fabs(V0[i] - V0[i - 1]);
            if (difference > max) {
                max = difference;
            }
        }
        fOut << V0[V0.size() - 1] << endl;
        currentTime += _t;
    }
    cout << "Max difference = " << max << endl;
    fOut.close();
}
