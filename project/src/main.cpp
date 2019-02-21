#include "../include/OrdinaryDE.h"
#include <iostream>

int main() {
    vector<double> startPoint {-0.1, 1.2};
    //EulerExplicit(func1, startPoint);
    //EulerImplicit(func1, startPoint);
    RungeKutta(func1, startPoint);
    return 0;
}