#include "../include/OrdinaryDE.h"
#include <iostream>

int main() {
    vector<double> startPoint {1, 0};
    //cout << "Explicit Euler" << endl;
    //EulerExplicit(func0, startPoint, true);
    //cout << "Implicit Euler" << endl;
    //EulerImplicit(func0, startPoint, true);
    //cout << "Runge - Kutta" << endl;
    //RungeKutta(func0, startPoint, false);
    //cout << "Symmetric" << endl;
    //Symmetric(func0, startPoint, true);
    //cout << "Adams - Bashfort" << endl;
    //AdamsBashfort(func0, startPoint, false);
    cout << "Prediction-correction" << endl;
    PredCorr(func0, startPoint);
    return 0;
}