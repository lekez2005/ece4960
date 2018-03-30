#include <cstdio>
#include <cmath>
#include <functional>

double f_x_exp(double x){
    return std::exp(50*x) - 1;
}

double f_prime_exp(double x){
    return 50*std::exp(50*x);
}

double f_x_cos(double x){
    return x - std::cos(x);
}

double f_prime_cos(double x){
    return 1 + std::sin(x);
}

double fixed_t(double *x, double *del_x, double *t, double curr_fx_val, 
        double (*f_x)(double), double (*f_prime)(double)){
    *del_x = -curr_fx_val/f_prime(*x);
    *t = 1;
    *x += (*t * *del_x);
}


double variable_t(double *x, double *del_x, double *t, double curr_fx_val, 
        double (*f_x)(double), double (*f_prime)(double)){
    double del_fx = -curr_fx_val/f_prime(*x);
    double temp_x = *x;
    int i = 4;
    while(i > -20){
        *t = std::pow(2, i);
        *del_x  = *t * del_fx;
        temp_x = *x + *del_x;
        double temp_fx = f_x(temp_x);
        if (std::fabs(temp_fx) <= std::fabs(curr_fx_val))
            break;
        i--;
    }
    *x = temp_x;
}

void newton(double x_0, double (*f_x)(double), double (*f_prime)(double), bool fixed){
    bool solved = false;
    int max_iterations = 50;
    int iterations = 1;
    double t = 1; // try 2^x from -5 to 1
    double del_x = 0;
    double x = x_0;
    double f_k = f_x(x);
    printf("k \t\t x_k \t \t delta_x \t t \t \t f_x\n"); 
    printf("%3d \t %10g \t \t -  \t\t - \t %10g \n", 0, x, f_k);
    while (!solved && iterations < max_iterations){
        if (fixed)
            fixed_t(&x, &del_x, &t, f_k, f_x, f_prime);
        else
            variable_t(&x, &del_x, &t, f_k, f_x, f_prime);
        f_k = f_x(x);
        printf("%3d \t %10g \t %10g \t %10g \t %10g \n", iterations, x, del_x, t, f_k);
        iterations++;
        if (std::abs(del_x/x) < 1e-7)
            solved = true;
    }
}

int main(){
    bool fixed = true;
    printf("Newton for cos, x(0) = 3 \n");
    newton(3, f_x_cos, f_prime_cos, true);

    fixed = false;
    printf("Newton for cos, x(0) = 10 \n");
    newton(10, f_x_cos, f_prime_cos, false);
    printf("\nNewton exp(50x)-1, x(0) = 10 \n");
    newton(10, f_x_exp, f_prime_exp, false);
}
