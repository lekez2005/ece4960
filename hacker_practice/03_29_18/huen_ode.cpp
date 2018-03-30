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

    printf("\nComparison of Forward Euler and Huen \n\n");

    double h = 1;
    int num_steps = 20;
    int max_steps = 100;
    double t = 0;
    double tol = 1e-7;

    double x_forward = 2;
    double x_huen = 2;
    double x_huen_temp = 2;
    double x_huen_tol = 2;

    printf("%12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s \n", "t", "x_true", "x_forward", "error", "x_huen_one", "error", 
            "x_huen", "iterations", "error");

    for (int i = 0; i < num_steps; i++){
        x_forward += h*slope(t, x_forward);

        x_huen_temp = x_huen + h*slope(t, x_huen);
        x_huen += 0.5*h*(slope(t, x_huen) + slope(t+h, x_huen_temp));

        double huen_old = x_huen_tol;
        int j;
        for (j = 0; j < max_steps; j++){
            x_huen_temp = x_huen_tol + h*slope(t, x_huen_tol);
            x_huen_tol = x_huen_tol + 0.5*h*(slope(t, x_huen_tol) + slope(t+h, x_huen_temp));
            if (std::fabs(x_huen_tol-huen_old) < tol){
                break;
            }
            huen_old = x_huen_tol;
        }
        t = t + h;
        double x_real = 4/1.3*(std::exp(0.8*t)-std::exp(-0.5*t)) + 2*std::exp(-0.5*t);
        printf("%12g \t %12g \t %12g \t %12g \t %12g \t %12g \t %12g \t %12d \t %12g \n", t, x_real, x_forward, percent_err(x_forward, x_real), 
                x_huen, percent_err(x_huen, x_real), x_huen_tol, j, percent_err(x_huen_tol, x_real));

    }

    return 0;
}
