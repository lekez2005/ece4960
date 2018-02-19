#include <cmath>
#include <cstdio>

#define fx_method double (*f_x)(double)
#define integ_method double (*integ)(double, double, fx_method)

double gauss_p1 = 0.5*(1-1/sqrt(3));
double gauss_p2 = 0.5*(1+1/sqrt(3));


double integrate(integ_method, fx_method, double start, double end, double h) {
    int n = round((end - start)/h);
    double result = 0;
    double x = start;
    for (int i = 0; i < n; i++){
        result += h*integ(x, h, f_x);
        x += h;
    }
    return result;
}

double rectangle(double x, double h, fx_method) {
    return f_x(x);
}

double trapezoid(double x, double h, fx_method) {
    return 0.5*(f_x(x) + f_x(x+h));
}

double midpoint(double x, double h, fx_method) {
    return f_x(x + 0.5*h);
}

double simpson(double x, double h, fx_method) {
    return ( f_x(x) + 4*f_x(x+0.5*h) + f_x(x+h))/6;
}

double two_gaussian(double x, double h, fx_method) {
    return 0.5*( f_x(x+h*gauss_p1) + f_x(x+h*gauss_p2));
}
int main() {
    double h = 0.1;
    double start = -1;
    double end = 1;
    printf("Rectangle:\t %G \n", integrate(rectangle, exp, start, end, h));
    printf("Trapezoid:\t %G \n", integrate(trapezoid, exp, start, end, h));
    printf("Mid Point:\t %G \n", integrate(midpoint, exp, start, end, h));
    printf("Simpson  :\t %G \n", integrate(simpson, exp, start, end, h));
    printf("2-Point Gauss:\t %G \n", integrate(two_gaussian, exp, start, end, h));
}

