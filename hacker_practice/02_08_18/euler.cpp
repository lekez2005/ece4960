#include <cstdio>
#include <cmath>

using namespace std;

int main() {

    double del_t[] = {0.5, 1.0, 2.0};
    double max_t = 20;

    for (int i = 0; i < sizeof(del_t)/sizeof(del_t[0]); i++) {
        double dt = del_t[i];
        int n = (int) max_t/dt + 1;


        double f_forward = 1.0;
        double f_backward = 1.0;

        printf("del_t = %g \n", dt);

        printf("%14s \t %14s \t %14s \t %14s \n", "t", "exp(-t)", "forward", "backward");

        for (int j = 0; j < n; j++) {
            double t = dt*j;
            printf("%14g \t %14g \t %14g \t %14g \n", t, exp(-t), f_forward, f_backward);
            f_forward *= (1-dt);
            f_backward /= (1 + dt); 
        }


    }

    return 0;
}
