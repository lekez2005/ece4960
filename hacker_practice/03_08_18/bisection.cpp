#include <cmath>
#include <cstdio>

#define fx_method double (*f_x)(double)

bool same_sign(double x1, double x2, double tolerance){
    // not same sign if either is zero
    if (std::abs(x1) < tolerance || std::abs(x2) < tolerance)
        return false;
    return std::signbit(x1) == std::signbit(x2);
}

bool same_sign(double x1, double x2){
    return same_sign(x1, x2, 1e-12);
}


double bisect_1d(double start, double end, fx_method, int max_iterations, double tolerance){

    int k = 0;
    double start_val = f_x(start);
    
    double mid = (start + end)/2;
    printf("%s \t %6s \t %3s \t %16s \n", "k", "start", "end", "f_x(mid)");
    while(k < max_iterations){
        mid = (start + end)/2;
        double mid_val = f_x(mid);

        printf("%d \t %8g \t %8g \t %8g \n", k, start, end, mid_val);

        if (same_sign(mid_val, start_val)){
            start = mid;
            start_val = mid_val;
        }else{
            end = mid;
        }
        

        if (std::abs(mid_val) < tolerance){
            break;
        }
        k++;
    }
    return mid;
}

void bisect_2d(double *x, double *y, double (*f1_x)(double, double),
        double (*f2_x)(double, double), int max_iterations, double tolerance){
    int k = 0;


    double mid_x = (x[0] + x[1])/2;
    double mid_y = (y[0] + y[1])/2;

    printf("%s \t %6s \t %3s \t %16s %16s \n", "k", "x", "y", "f1", "f2");

    while (k < max_iterations) {

        mid_x = (x[0] + x[1])/2;
        mid_y = (y[0] + y[1])/2;
        double mid_f1_val = f1_x(mid_x, mid_y);
        double mid_f2_val = f2_x(mid_x, mid_y);
        printf("%d \t %8g \t %8g \t %8g \t %8g \n", k, mid_x, mid_y, mid_f1_val, mid_f2_val);
        if (std::abs(mid_f1_val) < tolerance && std::abs(mid_f2_val) < tolerance)
            break;
        double region_x = 0;
        double region_y = 0;
        for (int i = 0; i < 2; i++){
            for (int j = 0; j < 2; j++){
                region_x = x[i];
                region_y = y[j];
                if (! same_sign(mid_f1_val, f1_x(region_x, region_y))){
                    if (! same_sign(mid_f2_val, f2_x(region_x, region_y))){
                        goto endloop;
                    }
                }
            }
        }
        k++;
    endloop:
        x[0] = mid_x;
        x[1] = region_x;
        y[0] = mid_y;
        y[1] = region_y;
    }
}

double exp_minus(double x){
    return std::exp(x)-1;
}

double exp_diff(double x, double y){
    return std::exp(x) - std::exp(y);
}

double exp_sum(double x, double y){
    return std::exp(x) + std::exp(y) - 2;
}


int main(){
    printf("\n1-D Bisection \n\n");
    bisect_1d(-5, 10, exp_minus, 30, 1e-7);
    printf("\n\n2-D Bisection \n\n");
    double x[] = {-5, 10};
    double y[] = {-5, 10};
    bisect_2d(x, y, exp_diff, exp_sum, 30, 1e-7); 
    return 0;
}
