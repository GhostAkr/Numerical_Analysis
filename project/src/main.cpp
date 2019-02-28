#include "../include/OrdinaryDE.h"
#include <iostream>

int main() {
    vector<double> startPoint {1, 0};
    //EulerExplicit(func1, startPoint);
    //EulerImplicit(func1, startPoint);
    cout << "Runge - Kutta" << endl;
    RungeKutta(func0, startPoint);
    cout << "Symmetric" << endl;
    Symmetric(func0, startPoint);
    cout << "Adams - Bashfort" << endl;
    AdamsBashfort(func0, startPoint);
    cout << "Prediction-correction" << endl;
    PredCorr(func0, startPoint);
    return 0;
}