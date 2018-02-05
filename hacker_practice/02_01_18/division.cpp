#include <iostream>

using namespace std;
int main() {
    double x = 1.234567890123456;
    int i = 1;

    x*= 1e-307;

    for (i = 1; i < 20; i++) {
        x /= 10.0;
        cout << x << endl;
    }
    return 0;
}
