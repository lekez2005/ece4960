#include <cstdio>
#include <cmath>

/*
double f_x(double x){
    return std::exp(50*x) - 1;
}

double f_prime(double x){
    return 50*std::exp(50*x);
}
*/
double f_x(double x){
    return x - std::cos(x);
}

double f_prime(double x){
    return 1 + std::cos(x);
}

void newton(double x_0){
    bool solved = false;
    int max_iterations = 40;
    int iterations = 1;
    double x = x_0;
    double f_k = f_x(x);
    printf("k \t\t x_k \t \t delta \t \t f_x\n"); 
    printf("%3d \t %10g \t \t -  \t %10g \n", 0, x, f_k);
    while (!solved && iterations < max_iterations){
        double f_k_1 = f_x(x);
        double del_x = -f_k_1/f_prime(x);
        x += del_x;
        printf("%3d \t %10g \t %10g \t %10g \n", iterations, x, del_x, f_k_1);
        iterations++;
        f_k = f_k_1;
        if (std::abs(del_x/x) < 1e-7)
            solved = true;
    }
}

int main(){
    printf("Newton x(0) = 1 \n");
    newton(1);
    printf("\nNewton x(0) = 10 \n");
    newton(10);
}
