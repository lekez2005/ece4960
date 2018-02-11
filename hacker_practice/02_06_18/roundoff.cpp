#include <cstdio>
#include <cmath>

using namespace std;

double F_PRIME = 2.0;

double rel_error(double val) {
    return abs(val-F_PRIME)/F_PRIME;
}

double f1(double x) {
    return x*x;
}

double f2(double x) {
    return x*x + 1e8;
}


int main() {

    double x_0 = 1;
    printf("h \t f'(x^2)_1 \t f'(x^2+10^10)_1 \t f'(x^2)_2 \t f'(x^2+10^10)_2 \n\n");
    for(int i = 1; i <=18; i++){
        double h = pow(10, -i);
        double x_plus = x_0 + h;
        double x_minus = x_0 - h;
        printf("%6g \t %6g \t %15g \t %12g \t %8g \n", h, rel_error((f1(x_plus) - f1(x_0))/h), 
                rel_error( (f2(x_plus) - f2(x_0))/h ),
                rel_error( (f1(x_plus) - f1(x_minus))/(2*h) ), 
                rel_error( (f2(x_plus) - f2(x_minus))/(2*h) ) );

    }
    printf("markdown output\n");
    for(int i = 1; i <=18; i++){
        double h = pow(10, -i);
        double x_plus = x_0 + h;
        double x_minus = x_0 - h;
        printf("%g | %g | %g | %g | %g |\n", h, rel_error((f1(x_plus) - f1(x_0))/h), 
                rel_error( (f2(x_plus) - f2(x_0))/h ),
                rel_error( (f1(x_plus) - f1(x_minus))/(2*h) ), 
                rel_error( (f2(x_plus) - f2(x_minus))/(2*h) ) );

    }


    return 0;
}

