#include "matplotlibcpp.h"
#include <cstdio>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
namespace plt = matplotlibcpp;

typedef Eigen::Matrix<double, 1, 1> Vector1d;


double norm(const Eigen::VectorXd& x){
    double res = 0;
    for (int i = 0; i < x.size(); i++)
        res += x[i]*x[i];
    return res;
}

// params[0] = a, params[1] = m
double fx_power(const Eigen::Vector2d &params, const Vector1d &x){
    return params[0]*std::pow(x[0], params[1]);
}

void gradient_pow(const Eigen::Vector2d &params, const Eigen::MatrixXd &data, Eigen::VectorXd &grad){
    grad[0] = 0;
    grad[1] = 0;
    for (int i = 0; i < data.rows(); i++){
        double x_pow_m = std::pow(data(i, 0), params[1]);
        double x_pow_prod = 2*x_pow_m*(params[0]*x_pow_m - data(i, 1));
        grad[0] += x_pow_prod;
        grad[1] += params[0]*std::log(data(i, 0))*x_pow_prod;
    }
}

double fx(const Eigen::VectorXd params, const Eigen::MatrixXd data){
    double res = 0;
    for (int i = 0; i < data.rows(); i++)
        res += std::pow(fx_power(params, Vector1d(data(i, 0)))-data(i, 1), 2);
    return res;
}

void hessian(const Eigen::VectorXd &params, const Eigen::VectorXd &grad, const Eigen::MatrixXd &data, 
        Eigen::MatrixXd &hess){
    double step = 1e-4;
    for (int i = 0; i < params.rows(); i++){
        Eigen::VectorXd tempParams = params;
        Eigen::VectorXd tempGrad = grad;
        double delta = step*tempParams[i];
        tempParams[i] += delta;
        gradient_pow(tempParams, data, tempGrad);
        for (int j = 0; j < params.rows(); j++){
            hess(j, i) = (tempGrad[j]-grad[j])/delta;
        }
    }
}


int fixed_t(Eigen::VectorXd &params, Eigen::VectorXd &delParams, const Eigen::MatrixXd &data,
        double *t, double curr_fx_val){
    int rows = params.rows();
    Eigen::VectorXd grad(rows);
    gradient_pow(params, data, grad);
    Eigen::MatrixXd hess(rows, rows);
    hessian(params, grad, data, hess);

    delParams = -hess.inverse()*grad;
    *t = 1;
    params += *t * delParams;
    return 0;
}

int variable_t(Eigen::VectorXd &params, Eigen::VectorXd &delParams, const Eigen::MatrixXd &data,
        double *t, double curr_fx_val){
    int rows = params.rows();
    Eigen::VectorXd grad(rows);
    gradient_pow(params, data, grad);
    Eigen::MatrixXd hess(rows, rows);
    hessian(params, grad, data, hess);
    
    Eigen::Vector2d tempDelParams = -hess.inverse()*grad;
    
    Eigen::Vector2d tempParams;
    int i = -1;
    while(i > -20){
        *t = std::pow(2, i);
        delParams = tempDelParams * *t;
        tempParams = params + delParams;
        double temp_fx = fx(tempParams, data);
        if (std::fabs(temp_fx) <= std::fabs(curr_fx_val))
            break;
        i--;
    }
    params = tempParams;
    if (i <= -20)
        return -1; // max iterations reached
    else
        return 0;
}


void newton(const Eigen::VectorXd &params0, Eigen::VectorXd &soln, Eigen::MatrixXd &data, bool fixed){
    bool solved = false;
    int max_iterations = 50;
    int iterations = 1;
    double t = 1; // try 2^x from -5 to 1
    Eigen::VectorXd params = params0;
    Eigen::VectorXd delParams;

    double f_k = fx(params, data);
    printf("%8s \t %8s \t %8s \t %8s \t %8s \t %8s \t %8s\n", "k", "a", "m",
            "norm(params)", "norm(delX)", "t", "f_x"); 
    printf("%8s \t %8g \t %8g \t %8g \t %8s \t %8s \t %8g\n", "0", params[0], params[1], norm(params), 
            "-", "-", f_k); 
    int status;
    while (!solved && iterations < max_iterations){
        if (fixed)
            status = fixed_t(params, delParams, data, &t, f_k);
        else
            status = variable_t(params, delParams, data, &t, f_k);

        f_k = fx(params, data);
        printf("%8d \t %8g \t %8g \t %8g \t %8g \t %8g \t %8g \n", iterations, params[0], params[1],
                norm(params), norm(delParams), t, f_k);
        iterations++;
        if (norm(delParams) < 1e-15 || status == -1)
            solved = true;
    }
    if (status == -1){
        printf("Minimum step size does not result in reduction in f. Terminating ...\n");
    }
    soln = params;
}



int main(){

    int data_size = 6;
    Eigen::MatrixXd data(2, data_size);
    data << 1.0, 4.5, 9.0, 20, 74, 181,
         3.0, 49.4, 245, 1808, 2.2e4, 7.3e4;
    data.transposeInPlace();
    std::cout << data.transpose() << std::endl;



    Eigen::Vector2d params0(2.5, 2.25);
    Eigen::VectorXd soln0(2);
    cout << "Initial condition a = " << params0[0] << " m = " << params0[1] << endl;
    printf("\nQuasi-Newton With Line Search\n");
    newton(params0, soln0, data, false);


    Eigen::Vector2d params1(2.8, 2.8);
    Eigen::VectorXd soln1(2);
    cout << "Initial condition a = " << params1[0] << " m = " << params1[1] << endl;
    newton(params1, soln1, data, false);

    Eigen::VectorXd x_vec = data.col(0);
    std::vector<double> x(x_vec.data(), x_vec.data() + x_vec.rows());


    Eigen::VectorXd y_vec = data.col(1);
    std::vector<double> y(y_vec.data(), y_vec.data() + y_vec.rows());

    int points_plot = 100;
    vector<double> x_plot(points_plot);
    vector<double> y0_plot(points_plot);
    vector<double> y1_plot(points_plot);
    double delta_x = (x[x_vec.rows()-1]-x[0])/points_plot;
    for (int i = 0; i < points_plot; i++){
        double x_point = x[0] + delta_x*i;
        x_plot[i] = x_point;
        y0_plot[i] = soln0[0]*pow(x_point, soln0[1]);
        y1_plot[i] = soln1[0]*pow(x_point, soln1[1]);
    }

    plt::plot(x, y, "*");
    plt::named_plot("soln0", x_plot, y0_plot, "r");
    plt::named_plot("soln1", x_plot, y1_plot, "b");
    plt::legend();
    plt::save("fit_plot.png");
    //plt::show();

    return 0;
}
    
