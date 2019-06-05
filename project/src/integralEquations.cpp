//
// Created by ighos on 05.06.2019.
//

#include "../include/integralEquations.h"

double core(double _x, double _s, int _testNum) {
    switch (_testNum) {
        case 1: {
            return 1 - _x * cos(_x * _s);
        }
        case 2: {
            return 1 - _x * cos(_x * _s);
        }
        default: {
            cout << "Exception in core" << endl;
            return -1.0;
        }
    }
}
double rightPart(double _x, int _testNum) {
    switch (_testNum) {
        case 1: {
            return 0.5 * (1 + sin(_x));
        }
        case 2: {
            return _x * _x + sqrt(_x);
        }
        default: {
            cout << "Exception in rightPart" << endl;
            return -1.0;
        }
    }
}
vector<double> limits(int _limitType) {
    vector<double> result(2, 0);
    switch (_limitType) {
        case 1: {
            result[0] = 0.0;
            result[1] = 1.0;
            return result;
        }
        case 2: {
            result[0] = 0.1;
            result[1] = 1.0;
            return result;
        }
    }
}


