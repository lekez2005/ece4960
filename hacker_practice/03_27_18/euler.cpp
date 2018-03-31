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
        double f_trap = 1.0;
        double f_exp = 1.0;

        printf("del_t = %g \n", dt);

        printf("%14s \t %14s \t %14s \t %14s \t %14s \t %14s \t %14s \n", "t", "forward", "backward", "trap", "e(forward)", "e(backward)", "e(trap)");

        for (int j = 0; j < n; j++) {
            double t = dt*(j+1);
            f_exp = exp(-t);
            f_forward *= (1-dt);
            f_backward /= (1 + dt); 
            f_trap = f_trap*(2 - dt)/(2 + dt);
            printf("%14g \t %14g \t %14g \t %14g \t %14g \t %14g \t %14g \n", t, f_forward, f_backward, f_trap, fabs(f_forward-f_exp), 
                    fabs(f_backward-f_exp), fabs(f_trap-f_exp));
        }


    }

    return 0;
}
