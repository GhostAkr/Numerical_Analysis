#include "../include/OrdinaryDE.h"
#include <iostream>

int main() {
    vector<double> startPoint {1, 0};
    EulerExplicit(func0, startPoint);
    return 0;
}