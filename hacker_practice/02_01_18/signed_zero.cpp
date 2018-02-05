#include <iostream>
#include <string>
#include <cmath>
#include <cfloat>

using namespace std;

bool isPositive(double x) {
    return signbit(x) == 0;
}

bool isNegative(double x) {
     return signbit(x) > 0;
}

void printSign(double x) {
    string sign;

    if (isPositive(x) ) {
        sign = "Positive";
    }else if (isNegative(x) ) {
        sign = "Negative";
    }else {
        sign = "Undefined";
    }
    cout << x << " is " << sign << endl;
}

void checkNaN (double x) {
    string nan_str = isnan(x) ? "NaN": "Not NaN";
    cout << x << " is " << nan_str << endl;
}

void checkInf(double x) {
    string inf_str = isinf(x) ? "": "not";
    cout << x << " is " << inf_str << " inf " << endl;
}

int main() {
    printSign(+0.0);
    printSign(-0.0);
    printSign(+1.0);
    printSign(-1.0);
    printSign(DBL_MAX);
    printSign(-1.0*DBL_MAX);
    printSign(1.0/0);
    printSign(-1.0/0);
    printSign(0.0/0);
    printSign(-0.0/0);

    checkNaN(4.0);
    checkNaN(0.0/0);
    
    checkInf(0.0/0);
    checkInf(1.0/0);
    checkInf(-1.0/0);

    return 0;
}
