#include "../include/OrdinaryDE.h"
#include <iostream>

int main() {
    string path = "../data/solution.dat";
    string path1 = "../data/solution1.dat";
    string path2 = "../data/solution2.dat";
    string path3 = "../data/solution3.dat";
    string path4 = "../data/solution4.dat";
    string path5 = "../data/solution5.dat";
    vector<double> startPoint {1, 0};
    vector<double> startPoint1 {2, 0};
    vector<double> startPoint2 {3, 0};
    vector<double> startPoint3 {4, 0};
    vector<double> startPoint4 {5, 0};
    vector<double> startPoint5 {6, 0};
    RungeKutta(func0, startPoint, false, path);
    RungeKutta(func0, startPoint1, false, path1);
    RungeKutta(func0, startPoint2, false, path2);
    RungeKutta(func0, startPoint3, false, path3);
    RungeKutta(func0, startPoint4, false, path4);
    RungeKutta(func0, startPoint5, false, path5);
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
    //cout << "Prediction-correction" << endl;
    //PredCorr(func0, startPoint);
    return 0;
}