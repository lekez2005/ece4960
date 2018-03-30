#include <cstdio>
#include <cmath>

using namespace std;

double slope(double t, double x){
    return 4*std::exp(0.8*t) - 0.5*x;
}

double percent_err(double computed, double actual){
    if (fabs(actual) > 1e-12){
        return std::fabs(100/actual*(computed-actual));
    }
    return -1;
}

int main() {

    printf("\nRK4 Time steps \n\n");

    double h = 1;
    int num_steps = 20;
    int max_steps = 100;
    double t = 0;
    double tol = 1e-7;

    double x_forward = 2;
    double x_huen = 2;
    double x_huen_temp = 2;
    double x_rk4 = 2;

    printf("%12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s \n", "t", "x_true", "k1", "k2", "k3", "k4", 
            "x_rk4", "e(rk4)", "e(huen_one)", "e(f_euler)");

    for (int i = 0; i < num_steps; i++){
        x_forward += h*slope(t, x_forward);

        x_huen_temp = x_huen + h*slope(t, x_huen);
        x_huen += 0.5*h*(slope(t, x_huen) + slope(t+h, x_huen_temp));

        double k1 = slope(t, x_rk4);
        double k2 = slope(t+h/2, x_rk4+k1*h/2);
        double k3 = slope(t+h/2, x_rk4+k2*h/2);
        double k4 = slope(t+h, x_rk4+k3*h);
        x_rk4 += h/6*(k1+2*k2+2*k3+k4);

        t = t + h;
        double x_real = 4/1.3*(std::exp(0.8*t)-std::exp(-0.5*t)) + 2*std::exp(-0.5*t);
        printf("%12g \t %12g \t %12g \t %12g \t %12g \t %12g \t %12g \t %12g \t %12g \t %12g\n", t, x_real, k1, k2, k3, k4, 
                x_rk4, percent_err(x_rk4, x_real), percent_err(x_huen, x_real), percent_err(x_forward, x_real));

    }

    return 0;
}
