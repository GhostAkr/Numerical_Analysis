#include "../include/integralEquations.h"

#include <iostream>

int main() {
    int testNum = 1;
    int limitType = 1;
    int N = 10;
    string path = "../data/outIntegral.dat";
    quadMethod(testNum, limitType, N, path);
    return 0;
}