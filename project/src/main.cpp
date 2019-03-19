#include "../include/OrdinaryDE.h"
#include <iostream>

#define METHOD RungeKutta
#define FUNC func3
#define FLAG false

int main() {
    string path = "../data/solution.dat";
    string path1 = "../data/solution1.dat";
    string path2 = "../data/solution2.dat";
    string path3 = "../data/solution3.dat";
    string path4 = "../data/solution4.dat";
    string path5 = "../data/solution5.dat";
    //vector<double> startPoint {1, 0};


//    vector<double> startPoint {0, 0.9};
//    vector<double> startPoint1 {0, 0.5};
//    vector<double> startPoint2 {0, 1.1};
//    vector<double> startPoint3 {-0.1, -0.9};
//    vector<double> startPoint4 {0.1, -1.1};
//    vector<double> startPoint5 {0.1, -1.2};

    vector<double> startPoint {0, 0.9, 0};
    vector<double> startPoint1 {2, 0.5, 1};
    vector<double> startPoint2 {5, 1.1, 2};
    vector<double> startPoint3 {0, -0.9, -1};
    vector<double> startPoint4 {-2, -1.1, -2};
    vector<double> startPoint5 {-5, -1.2, 0};
    METHOD(FUNC, startPoint, FLAG, path);
    METHOD(FUNC, startPoint1, FLAG, path1);
    METHOD(FUNC, startPoint2, FLAG, path2);
    METHOD(FUNC, startPoint3, FLAG, path3);
    METHOD(FUNC, startPoint4, FLAG, path4);
    METHOD(FUNC, startPoint5, FLAG, path5);
    //cout << "Explicit Euler" << endl;
    //EulerExplicit(func0, startPoint, true);
    //cout << "Implicit Euler" << endl;
    //EulerImplicit(func0, startPoint, true);
    //cout << "Runge - Kutta" << endl;
    //Runy' geKutta(func0, startPoint, false);
    //cout << "Symmetric" << endl;
    //Symmetric(func0, startPoint, true);
    //cout << "Adams - Bashfort" << endl;
    //AdamsBashfort(func0, startPoint, false);
    //cout << "Prediction-correction" << endl;
    //PredCorr(func0, startPoint);
    return 0;
}