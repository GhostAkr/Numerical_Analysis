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
    double L = 2;
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
//    cout << "Size of V0 = " << V0.size() << endl;
    for (int i = 1; i < V0.size() - 1; ++i) {
        double x = i * _h;
//        cout << "x0 = " << x << endl;
        V0[i] = zeroFunction(x, _testNum);
        cout << i * _h << endl;
        if (fabs(i * _h - 1) < 1e-5) {
            cout << "qq" << endl;
            V0[i] = 0.0;
        }
        cout << "V0[i] = " << V0[i] << endl;
//        if (fabs(V0[i] - V0[i - 1]) > 1e-1) {
//            V0[i] = fabs(V0[i] - V0[i - 1]) / 2.0;
//        }
//        cout << "V0[i] = " << V0[i] << endl;
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
        if (fabs(i * _h - 1) < 1e-5) {
            V1[i] = 0.0;
        }
//        if (fabs(V1[i] - V1[i - 1]) > 1e-1) {
//            V1[i] = fabs(V1[i] - V1[i - 1]) / 2.0;
//        }
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
            //cout << "n = " << n << endl;
            if (fabs(n * _h - 1) < 1e-5) {
                nextLayer[n] = 0.0;
            }
//            if (fabs(nextLayer[n] - nextLayer[n - 1]) > 1e-1) {
//                nextLayer[n] = fabs(nextLayer[n] - nextLayer[n - 1]) / 2.0;
//            }
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

double zeroFunction(double _x, int _testNum) {
    switch (_testNum) {
        case 1: {
            return sin(Pi * _x);
        }
        case 2:
            return _x * (1.0 - _x);
        case 3:
            return _x * (_x + 1.0);
        default:
            cout << "Exception in zeroFunction" << endl;
            return -1.0;
    }
}

double gFunction(double _x, int _testNum) {
    switch (_testNum) {
        case 1:
            return 0;
        case 2:
            return 0;
        case 3:
            return cos(_x);
        default:
            cout << "Exception in gFunction" << endl;
            return -1.0;
    }
}

double phiFunction(double _t, int _testNum){
    switch (_testNum) {
        case 1:
            return 0;
        case 2:
            return 0;
        case 3:
            return 0;
        default:
            cout << "Exception in phiFunction" << endl;
            return -1.0;
    }
}

double psiFunction(double _t, int _testNum) {
    switch (_testNum) {
        case 1:
            return 0;
        case 2:
            return 0;
        case 3:
            return 2.0 * (_t + 1.0);
        default:
            cout << "Exception in psiFunction" << endl;
            return -1.0;
    }
}

double zeroDerivative(double _x, int _testNum) {
    switch (_testNum) {
        case 1:
            return -Pi * Pi * sin(Pi * _x);
        case 2:
            return -2.0;
        case 3:
            return 2.0;
        default:
            cout << "Exception in zeroDerivative" << endl;
            return -1.0;
    }
}

double exactSolution(double _t, double _x, int _testNum) {
    switch (_testNum) {
        case 1:
            return sin(Pi * _x) * cos(Pi * _t);
        case 2: {
            double eps = 1e-7;
            double k = 0.5 * (sqrt(2.0 / (Pi * Pi * eps)) - 1);
            //cout << "k = " << k << endl;
            double sum = 0.0;
            for (int n = 0; n <= k; ++n) {
                sum += 1 / pow(2 * n + 1, 3) * sin((2 * n + 1) * Pi * _x) * cos((2 * n + 1) * Pi * _t);
            }
            return 8.0 / pow(Pi, 3) * sum;
        }
        case 3: {
            //cout << "Third test hasn't exact solution" << endl;
            return -1.0;
        }
        default: {
            cout << "Exception in psiFunction" << endl;
            return -1.0;
        }
    }
}

void errorUpdate(vector<double> _layer, double* _error, double _t, double _h, int _testNum) {
    for (int i = 0; i < _layer.size(); ++i) {
        double exact = exactSolution(_t, i * _h, _testNum);
        double difference = fabs(_layer[i] - exact);
        if (fabs(difference - 1.64) <= 1e-5) {
            cout << "t = " << _t << endl;
            cout << "x = " << i * _h << endl;
        }
        //cout << "difference = " << difference << endl;
        if (difference > *_error) {
            *_error = difference;
        }
    }
}
