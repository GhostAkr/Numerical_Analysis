#include "../include/OrdinaryDE.h"
#include <iostream>

int main() {
    vector<double> startPoint {1, 0};
    //EulerExplicit(func0, startPoint);
    EulerImplicit(func0, startPoint);
    //cout << "Runge - Kutta" << endl;
    //RungeKutta(func0, startPoint);
    //cout << "Symmetric" << endl;
    //Symmetric(func0, startPoint);
    //cout << "Adams - Bashfort" << endl;
    //AdamsBashfort(func0, startPoint);
    //cout << "Prediction-correction" << endl;
    //PredCorr(func0, startPoint);
    return 0;
}