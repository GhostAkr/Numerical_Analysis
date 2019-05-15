//
// Created by ighos on 21.04.2019.
//

#include "../include/WaveEquation.h"
#include "../include/HeatEquation.h"

void crossScheme(double _t, double _h, string _path, int _testNum) {
    std::ofstream fOut(_path);
    if (!fOut) {
        cout << "Error while opening file" << endl;
        return;
    }
    // Constants and parameters
    double L = 1;
    double N = L / _h;
    double T = 10;
    double a = 1.0;
    double courantNumber = a * _t / _h;
    double courantNumberSqr = courantNumber * courantNumber;
    double error = -1.0;
    // Null time layer
    vector<double> V0 = mesh(N, L);
    // Left and right parts
    V0[0] = phiFunction(0, _testNum);
    V0[V0.size() - 1] = psiFunction(0, _testNum);
    fOut << V0[0] << " ";
    for (int i = 1; i < V0.size() - 1; ++i) {
        double x = i * _h;
        V0[i] = zeroFunction(x, _testNum);
        fOut << V0[i] << " ";
    }
    fOut << V0[V0.size() - 1] << endl;
    errorUpdate(V0, &error, 0, _h, _testNum);
    // First time layer
    vector<double> V1 = mesh (N, L);
    // Left and right parts
    V1[0] = phiFunction(_t, _testNum);
    fOut << V1[0] << " ";
    V1[V1.size() - 1] = psiFunction(_t, _testNum);
    // Middle part
    for (int i = 1; i < V1.size() - 1; ++i) {
        double x = i * _h;
        V1[i] = V0[i] + _t * gFunction(x, _testNum) + (a * a * _t * _t / 2.0) * zeroDerivative(x, _testNum);
        fOut << V1[i] << " ";
    }
    fOut << V1[V1.size() - 1] << endl;
    errorUpdate(V1, &error, _t, _h, _testNum);
    // Next layers
    double currentTime = _t;
    while (currentTime <= T) {
        currentTime += _t;
        vector<double> nextLayer = mesh(N, L);
        // Left and right parts
        nextLayer[0] = phiFunction(currentTime, _testNum);
        fOut << nextLayer[0] << " ";
        nextLayer[nextLayer.size() - 1] = psiFunction(currentTime, _testNum);
        // Middle part
        for (int n = 1; n < nextLayer.size() - 1; ++n) {
            nextLayer[n] = 2 * V1[n] - V0[n] + courantNumberSqr * (V1[n + 1] - 2 * V1[n] + V1[n - 1]);
            fOut << nextLayer[n] << " ";
        }
        fOut << nextLayer[nextLayer.size() - 1] << endl;
        errorUpdate(nextLayer, &error, currentTime, _h, _testNum);
        V0 = V1;
        V1 = nextLayer;
    }
    cout << "Error = " << error << endl;
    fOut.close();
}

void crossSchemeDerivative(double _t, double _h, string _path, int _testNum) {
    std::ofstream fOut(_path);
    if (!fOut) {
        cout << "Error while opening file" << endl;
        return;
    }
    // Constants and parameters
    double L = 1;
    double N = L / _h;
    double T = 10;
    double a = 1.0;
    double courantNumber = a * _t / _h;
    double courantNumberSqr = courantNumber * courantNumber;
    double error = -1.0;
    // Null time layer
    vector<double> V0 = mesh(N, L);
    // Left and right parts
    V0[0] = phiFunction(0, _testNum);
    fOut << V0[0] << " ";
    for (int i = 1; i < V0.size() - 1; ++i) {
//        double x = i * _h;
        double x = (i + 0.5) * _h;
        V0[i] = zeroFunction(x, _testNum);
        fOut << V0[i] << " ";
    }
//    V0[V0.size() - 1] = psiFunction(0, _testNum);
//    V0[V0.size() - 1] = 0;
    V0[V0.size() - 1] = V0[V0.size() - 2];
    fOut << V0[V0.size() - 1] << endl;
    errorUpdate(V0, &error, 0, _h, _testNum);
    // First time layer
    vector<double> V1 = mesh (N, L);
    // Left and right parts
    V1[0] = phiFunction(_t, _testNum);
    fOut << V1[0] << " ";
    // Middle part
    for (int i = 1; i < V1.size() - 1; ++i) {
//        double x = i * _h;
        double x = (i + 0.5) * _h;
        V1[i] = V0[i] + _t * gFunction(x, _testNum) + (a * a * _t * _t / 2.0) * zeroDerivative(x, _testNum);
        fOut << V1[i] << " ";
    }
//    V1[V1.size() - 1] = psiFunction(_t, _testNum);
//    V1[V1.size() - 1] = 0;
    V1[V1.size() - 1] = V1[V1.size() - 2];
    fOut << V1[V1.size() - 1] << endl;
    errorUpdate(V1, &error, _t, _h, _testNum);
    // Next layers
    double currentTime = _t;
    while (currentTime <= T) {
        currentTime += _t;
        vector<double> nextLayer = mesh(N, L);
        // Left and right parts
        nextLayer[0] = phiFunction(currentTime, _testNum);
        fOut << nextLayer[0] << " ";
        // Middle part
        for (int n = 1; n < nextLayer.size() - 1; ++n) {
            nextLayer[n] = 2 * V1[n] - V0[n] + courantNumberSqr * (V1[n + 1] - 2 * V1[n] + V1[n - 1]);
            fOut << nextLayer[n] << " ";
        }
//        nextLayer[nextLayer.size() - 1] = psiFunction(currentTime, _testNum);
//        nextLayer[nextLayer.size() - 1] = 0;
        nextLayer[nextLayer.size() - 1] = nextLayer[nextLayer.size() - 2];
        fOut << nextLayer[nextLayer.size() - 1] << endl;
        errorUpdate(nextLayer, &error, currentTime, _h, _testNum);
        V0 = V1;
        V1 = nextLayer;
    }
    cout << "Error = " << error << endl;
    fOut.close();
}

double zeroFunction(double _x, int _testNum) {
    switch (_testNum) {
        case 1: {
            return sin(Pi * _x);
        }
        case 2: {
            return _x * (1.0 - _x);
        }
        case 3: {
            return _x * (_x + 1.0);
        }
        case 4: {
            if (_x < 0.4) {
                return 0;
            } else if (_x >= 0.4 && _x <= 0.6) {
                return 1;
            } else {
                return 0;
            }
        }
        default: {
            cout << "Exception in zeroFunction" << endl;
            return -1.0;
        }
    }
}

double gFunction(double _x, int _testNum) {
    switch (_testNum) {
        case 1: {
            return 0;
        }
        case 2: {
            return 0;
        }
        case 3: {
            return cos(_x);
        }
        case 4: {
            return 0;
        }
        default: {
            cout << "Exception in gFunction" << endl;
            return -1.0;
        }
    }
}

double phiFunction(double _t, int _testNum){
    switch (_testNum) {
        case 1: {
            return 0;
        }
        case 2: {
            return 0;
        }
        case 3: {
            return 0;
        }
        case 4: {
            return 0;
        }
        default: {
            cout << "Exception in phiFunction" << endl;
            return -1.0;
        }
    }
}

double psiFunction(double _t, int _testNum) {
    switch (_testNum) {
        case 1: {
            return 0;
        }
        case 2: {
            return 0;
        }
        case 3: {
            return 2.0 * (_t + 1.0);
        }
        case 4: {
            return 0;
        }
        default: {
            cout << "Exception in psiFunction" << endl;
            return -1.0;
        }
    }
}

double muFunction(double _x, int _testNum) {
    switch (_testNum) {
        case 1: {
            return 0;
        }
        default: {
            cout << "Exception in muFunction" << endl;
            return -1.0;
        }
    }
}

double zeroDerivative(double _x, int _testNum) {
    switch (_testNum) {
        case 1: {
            return -Pi * Pi * sin(Pi * _x);
        }
        case 2: {
            return -2.0;
        }
        case 3: {
            return 2.0;
        }
        case 4: {
            return 0;
        }
        default: {
            cout << "Exception in zeroDerivative" << endl;
            return -1.0;
        }
    }
}

double exactSolution(double _t, double _x, int _testNum) {
    switch (_testNum) {
        case 1: {
            return sin(Pi * _x) * cos(Pi * _t);
        }
        case 2: {
            double eps = 1e-8;
            double k = 0.5 * (sqrt(2.0 / (Pi * Pi * eps)) - 1);
            double sum = 0.0;
            for (int n = 0; n <= k; ++n) {
                sum += 1 / pow(2 * n + 1, 3) * sin((2 * n + 1) * Pi * _x) * cos((2 * n + 1) * Pi * _t);
            }
            return 8.0 / pow(Pi, 3) * sum;
        }
        case 3: {
            return -1.0;
        }
        case 4: {
            return -1.0;
        }
        default: {
            cout << "Exception in exactSolution" << endl;
            return -1.0;
        }
    }
}

void errorUpdate(vector<double> _layer, double* _error, double _t, double _h, int _testNum) {
    for (int i = 0; i < _layer.size(); ++i) {
        double exact = exactSolution(_t, i * _h, _testNum);
        double difference = fabs(_layer[i] - exact);
        if (difference > *_error) {
            *_error = difference;
        }
    }
}

vector<double> meshDerivative(size_t _N, double _L) {
    double h = _L / _N;
    _N++;  // Adding of first point
    vector<double> resultMesh(_N + 1, -1);
    return resultMesh;
}
