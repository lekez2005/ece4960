#include <cstdio>
#include <cmath>
#include <Eigen/Dense>

double norm(const Eigen::Vector2d& x){
    return x[0]*x[0] + x[1]*x[1];
}

double fx(const Eigen::Vector2d &x){
    return std::pow(3*x[0]*x[0] + x[1] - 4, 2) + std::pow(x[0]*x[0]-3*x[1]+2, 2);
}


void gradient(const Eigen::Vector2d &x, Eigen::Vector2d &grad){
    double comp_1 = 3*x[0]*x[0] + x[1] - 4;
    double comp_2 = x[0]*x[0]-3*x[1]+2;
    grad[0] = x[0]*(12*comp_1 + 4*comp_2);
    grad[1] = 2*comp_1 -6*comp_2;
}

int fixed_t(Eigen::Vector2d &x, Eigen::Vector2d &delX, double *t, double curr_fx_val){
    Eigen::Vector2d grad;
    gradient(x, grad);

    *t = 1;
    delX = *t * -grad;
    x += *t * delX;
    return 0;
}


int variable_t(Eigen::Vector2d &x, Eigen::Vector2d &delX, double *t, double curr_fx_val){
    Eigen::Vector2d grad;
    gradient(x, grad);

    Eigen::Vector2d tempDelX = -grad;
    Eigen::Vector2d tempX;
    int i = 4;
    while(i > -20){
        *t = std::pow(2, i);
        delX = tempDelX * *t;
        tempX = x + delX;
        double temp_fx = fx(tempX);
        if (std::fabs(temp_fx) <= std::fabs(curr_fx_val))
            break;
        i--;
    }
    x = tempX;
    if (i <= -20)
        return -1; // max iterations reached
    else
        return 0;
}


void newton(const Eigen::MatrixXd &x0, bool fixed){
    bool solved = false;
    int max_iterations = 20;
    int iterations = 1;
    double t = 1; // try 2^x from -5 to 1
    Eigen::Vector2d x = x0;
    Eigen::Vector2d delX;

    double f_k = fx(x);
    printf("%8s \t %8s \t %8s \t %8s \t %8s \t %8s \t %8s\n", "k", "x[0]", "x[1]",
            "norm(x)", "norm(delX)", "t", "f_x"); 
    printf("%8s \t %8g \t %8g \t %8g \t %8s \t %8s \t %8g\n", "0", x[0], x[1], norm(x), 
            "-", "-", f_k); 
    while (!solved && iterations < max_iterations){
        int status;
        if (fixed)
            status = fixed_t(x, delX, &t, f_k);
        else
            status = variable_t(x, delX, &t, f_k);

        if (status == -1){
            printf("Minimum step size does not result in reduction in f. Terminating ...\n");
            return;
        }
        f_k = fx(x);
        printf("%8d \t %8g \t %8g \t %8g \t %8g \t %8g \t %8g \n", iterations, x[0], x[1],
                norm(x), norm(delX), t, f_k);
        iterations++;
        if (norm(delX) < 1e-15)
            solved = true;
    }
}



int main(){
    printf("Initial condition (0, 0) \n");
    printf("\nQuasi-Newton Without Line Search\n");
    Eigen::Vector2d x0(0, 0);
    newton(x0, true);
    printf("\nQuasi-Newton With Line Search\n");
    newton(x0, false);

    return 0;
}
    
