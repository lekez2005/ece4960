#include <cstdio>
#include <cmath>

#include <vector>

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

double rk_step(double t, double x, double h,
        int order, double *a, double *p, double *q){
    std::vector<double> k(order, 0.0);
    double res = x;

    int q_index = 0;
    
    for (int i = 0; i < order; i++){
        // calculate x evaluation point
        double x_eval = x;
        for(int j=0; j < i; j++){
            x_eval += q[q_index]*k[j]*h;
            q_index++;
        }
        // calculate t evaluation point
        double t_eval;
        if (i == 0)
            t_eval = t;
        else
            t_eval = t + p[i-1]*h;
        k[i] = slope(t_eval, x_eval);
        res += k[i]*a[i]*h;
    }
    return res;
}

// this only works for dx/dt = 4exp(0.8*t) - 0.5x

void trbdf2_step(double t, double x, double h, double *x_gamma, double *x_next){
    static double gamma = 2-std::sqrt(2);
    static double k1 = -(1-gamma)*(1-gamma)/gamma/(2-gamma);
    static double k2 = 1/gamma/(2-gamma);
    static double k3 = (1-gamma)/(2-gamma); 

    double gam_h_2 = gamma*h/2;

    *x_gamma = (x + gam_h_2*slope(t, x) + gam_h_2*4*std::exp(0.8*(t + gamma*h)))/(1+gam_h_2);

    *x_next = (k1*x + k2* *x_gamma + k3*h*4*std::exp(0.8*(t + h)))/(1 + 0.5*h*k3);

}

int main() {

    printf("\nTRBDF-2 Time steps \n\n");
    printf("%% e(trbdf2) is 100*error(trbff2)/x_real \n");
    printf("%% e(rk34) is 100*(rk4-rk3)/x_real \n");
    printf("tr_step is step size calculated using trbdf2\n\n");

    double a_rk4[] = {7.0/24, 6.0/24, 8.0/24, 3.0/24};
    double p_rk4[] = {1.0/2, 3.0/4, 1};
    double q_rk4[] = {0.5, 0, 0.75, 0, 0, 1};

    double a_rk3[] = {2.0/9, 3.0/9, 4.0/9};
    double p_rk3[] = {0.5, 0.75};
    double q_rk3[] = {0.5, 0, 0.75};

    double h = 1;
    double t = 0;
    double reltol = 1e-4;
    double abstol = 1e-6;
    int num_steps = 20;

    double gamma = 2-std::sqrt(2);

    double x_rk3 = 2;
    double x_rk4 = 2;
    double x_trbdf = 2;
    double x_trbdf_temp = 2;
    double x_gamma = 2;

    printf("%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \n", "t", "x_true", "x_trbdf2", "x_rk4", "% e(trbdf2)", "% e(rk43)", "tr_step");

    for (int i = 0; i < num_steps; i++){
        x_rk3 = rk_step(t, x_rk4, h, 3, a_rk3, p_rk3, q_rk3);
        x_rk4 = rk_step(t, x_rk4, h, 4, a_rk4, p_rk4, q_rk4);
        trbdf2_step(t, x_trbdf, h, &x_gamma, &x_trbdf_temp); 
        double err34 = std::fabs(x_rk4-x_rk3);
        double errbdf = (3*gamma*gamma-4*gamma+2)/(6*(gamma-2)) * ( slope(t, x_trbdf)/gamma -
                slope(t+gamma*h, x_gamma)/(gamma*(1-gamma)) + slope(t+h, x_trbdf_temp)/(1-gamma));
        double time_step = h*std::pow( (reltol*std::fabs(x_trbdf_temp + abstol))/std::fabs(errbdf), 1.0/3);
        t = t + h;
        double x_real = 4/1.3*(std::exp(0.8*t)-std::exp(-0.5*t)) + 2*std::exp(-0.5*t);
        printf("%10g \t %10g \t %10g \t %10g \t %10g \t %10g \t %10g \n", t, x_real, x_trbdf, x_rk4,
                std::fabs(errbdf/x_real)*100, std::fabs(err34/x_real)*100, time_step);
        x_trbdf = x_trbdf_temp;
        h = time_step;

    }

    return 0;
}
