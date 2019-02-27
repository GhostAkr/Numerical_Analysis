#include "../include/OrdinaryDE.h"
#include <iostream>

int main() {
    vector<double> startPoint {1, 0};
    //EulerExplicit(func1, startPoint);
    //EulerImplicit(func1, startPoint);
    //RungeKutta(func0, startPoint);
    //Symmetric(func0, startPoint);
    AdamsBashfort(func0, startPoint);
    return 0;
}