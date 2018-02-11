#include <cstdio>
#include <cmath>

using namespace std;

double x_0 = 1;

double f_x(double x) {
    return pow(x, 3);
}

double f_prime_1(double h) {
    return (f_x(x_0 + h) - f_x(x_0))/h;
}

double f_prime_2(double h) {
    return f_prime_1(2*h);
}

double f_prime_3(double h) {
    return -(f_x(x_0 + 2*h)/(2*h) + 3*f_x(x_0)/(2*h)) + 2*f_x(x_0 + h)/h;
}

double rel_error(double f_prime) {
    return abs(3-f_prime)/3;
}

double eta_1(double h, double (*f_func)(double) ) {
    return rel_error(f_func(2*h))/rel_error(f_func(h));
}

double eta_2(double h, double (*f_func)(double) ) {
    return (f_func(4*h) - f_func(2*h))/(f_func(2*h) - f_func(h));
}



int main() {
    double f_prime_real = 3;

    int arr_size = 10 - 4 + 1;
    double h_array[arr_size];
    for (int i = 0; i < arr_size; i++) {
        h_array[i] = pow(10, i-4);
    }

    printf("%10s \t %10s \t %10s \t %10s \n","h", "f_prime_h", "f_prime_2h", "f_prime_second_order");
    for (int i = 0; i < arr_size; i++) {
        double h = pow(10, -i);
        printf("%10g \t %10g \t %10g \t %10g \n", h, f_prime_1(h), f_prime_2(h), f_prime_3(h));
    }

    printf("\n%10s \t %10s \t %10s \t %10s \n","h", "E(f_prime_h)", "E(f_prime_2h)", "E(f_prime_second_order)");
    for (int i = 0; i < arr_size; i++) {
        double h = pow(10, -i);
        printf("%10g \t %10g \t %10g \t %10g \n", h, rel_error(f_prime_1(h)),
                    rel_error(f_prime_2(h)), rel_error(f_prime_3(h)));
    }

    printf("\n%10s \t %10s \t %10s \t %10s \n","h", "f_1", "f_2", "f_3");
    for (int i = 0; i < arr_size; i++) {
        double h = pow(10, -i);
        printf("%10g \t %10g \t %10g \t %10g \n", h, eta_1(h, f_prime_1), eta_1(h, f_prime_2), eta_1(h, f_prime_3));
    }

    printf("\n%10s \t %10s \t %10s \t %10s \n","h", "f_1", "f_2", "f_3");
    for (int i = 0; i < arr_size; i++) {
        double h = pow(10, -i);
        printf("%10g \t %10g \t %10g \t %10g \n", h, eta_2(h, f_prime_1), eta_2(h, f_prime_2), eta_2(h, f_prime_3));
    }

    printf("markdown output \n");
    for (int i = 4; i <= 40; i++) {
        double h = pow(10, -i);
        printf("| %10g | %10g | %10g | %10g | %10g | %10g | \n", h, rel_error(f_prime_1(h)), rel_error(f_prime_2(h)), 
                rel_error(f_prime_3(h)),
                eta_1(h, f_prime_1), eta_2(h, f_prime_1));
    }

    return 0;

}
