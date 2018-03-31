#include <cstdio>
#include <cmath>
#include <functional>

double f_x_exp(double x){
    return std::exp(50*x) - 1;
}

double f_prime_exp(double x){
    return 50*std::exp(50*x);
}

double f_prime_exp_quasi(double x){
    double delta_x = 0.0001*x;
    if (std::fabs(x) > 1e-20){
        return (f_x_exp(x+delta_x)-f_x_exp(x))/delta_x;
    }else{
        return 0;
    }
}

double f_x_cos(double x){
    return x - std::cos(x);
}

double f_prime_cos(double x){
    return 1 + std::sin(x);
}

int fixed_t(double *x, double *del_x, double *t, double curr_fx_val, 
        double (*f_x)(double), double (*f_prime)(double)){
    *del_x = -curr_fx_val/f_prime(*x);
    *t = 1;
    *x += (*t * *del_x);
}


int variable_t(double *x, double *del_x, double *t, double curr_fx_val, 
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
    if (i <= -20)
        return -1; // max iterations reached
    else
        return 0;
}

void newton(double x_0, double (*f_x)(double), double (*f_prime)(double), bool fixed){
    bool solved = false;
    int max_iterations = 50;
    int iterations = 1;
    double t = 1; // try 2^x from -5 to 1
    double del_x = 0;
    double x = x_0;
    double f_k = f_x(x);
    printf("%8s \t %8s \t %8s \t %8s \t %8s\n", "k", "x_k", "delta_x", "t", "f_x"); 
    printf("%8s \t %8g \t %8s \t %8s \t %8g\n", "0", x, "-", "-", f_k); 
    while (!solved && iterations < max_iterations){
        int status;
        if (fixed)
            status = fixed_t(&x, &del_x, &t, f_k, f_x, f_prime);
        else
            status = variable_t(&x, &del_x, &t, f_k, f_x, f_prime);

        if (status == -1){
            printf("Minimum step size does not result in reduction in f. Terminating ...\n");
            return;
        }
        f_k = f_x(x);
        printf("%8d \t %8g \t %8g \t %8g \t %8g \n", iterations, x, del_x, t, f_k);
        iterations++;
        if (std::abs(del_x) < 1e-15)
            solved = true;
    }
}

int main(){
    printf("\nQuasi-Newton With Line Search\n");
    bool fixed = false;
    printf("\nNewton exp(50x)-1, x(0) = 1 \n");
    newton(1, f_x_exp, f_prime_exp, false);
    printf("\nQuasi-Newton exp(50x)-1, x(0) = 1 \n");
    newton(1, f_x_exp, f_prime_exp_quasi, false);
}
