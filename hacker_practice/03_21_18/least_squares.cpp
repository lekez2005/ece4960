#include "matplotlibcpp.h"
#include <cstdio>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
namespace plt = matplotlibcpp;

typedef Eigen::Matrix<double, 1, 1> Vector1d;

//function to be optimized
double (*opt_func)(const Eigen::VectorXd &params, const Eigen::VectorXd &x);
// how 'least squares' is defined e.g. sum_square, lagrange
double (*fx_func)(const Eigen::VectorXd &params, const Eigen::MatrixXd &data);
// how gradient is calculated
void (*grad_func)(const Eigen::VectorXd &params, const Eigen::MatrixXd &data, Eigen::VectorXd &grad);

double norm(const Eigen::VectorXd& x){
    double res = 0;
    for (int i = 0; i < x.size(); i++)
        res += x[i]*x[i];
    return res;
}

void gradient_power(const Eigen::Vector2d &params, const Eigen::MatrixXd &data,
        Eigen::VectorXd &grad, bool lagrange){
    grad[0] = 0;
    grad[1] = 0;
    for (int i = 0; i < data.rows(); i++){
        double x_pow_m = std::pow(data(i, 0), params[1]);
        double x_pow_prod = 2*x_pow_m*(params[0]*x_pow_m - data(i, 1));
        if (std::fabs(data(i, 1)) > 1e-15 && lagrange) 
            x_pow_prod /= std::pow(data(i, 1), 2);
        grad[0] += x_pow_prod;
        grad[1] += params[0]*std::log(data(i, 0))*x_pow_prod;
    }
}


void grad_power_lagrange(const Eigen::VectorXd &params, const Eigen::MatrixXd &data,
        Eigen::VectorXd &grad){
    gradient_power(params, data, grad, true);
}


void grad_power_raw(const Eigen::VectorXd &params, const Eigen::MatrixXd &data,
        Eigen::VectorXd &grad){
    gradient_power(params, data, grad, false);
}

double fx_(const Eigen::VectorXd &params, const Eigen::MatrixXd &data, bool lagrange=false){
    double res = 0;
    int rows = data.rows();
    int cols = data.cols();
    for (int i = 0; i < rows; i++){
        double square_diff = std::pow(opt_func(params, data.block(i, 0, cols-1, 1))-data(i, 1),
                    2);
        if (std::fabs(data(i, 1)) > 1e-15 && lagrange)
            square_diff /= std::pow(data(i, 1), 2);
        res += square_diff;
    }
    return res;
}

double fx_raw(const Eigen::VectorXd &params, const Eigen::MatrixXd &data){
    return fx_(params, data, false);
}


double fx_lagrange(const Eigen::VectorXd &params, const Eigen::MatrixXd &data){
    return fx_(params, data, true);
}

void hessian(const Eigen::VectorXd &params, const Eigen::VectorXd &grad, const Eigen::MatrixXd &data, 
        Eigen::MatrixXd &hess){
    double step = 1e-4;
    for (int i = 0; i < params.rows(); i++){
        Eigen::VectorXd tempParams = params;
        Eigen::VectorXd tempGrad = grad;
        double delta = step*tempParams[i];
        tempParams[i] += delta;
        grad_func(tempParams, data, tempGrad);
        for (int j = 0; j < params.rows(); j++){
            hess(j, i) = (tempGrad[j]-grad[j])/delta;
        }
    }
}


int fixed_t(Eigen::VectorXd &params, Eigen::VectorXd &delParams, const Eigen::MatrixXd &data,
        double *t, double curr_fx_val){
    int rows = params.rows();
    Eigen::VectorXd grad(rows);
    grad_func(params, data, grad);
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
    grad_func(params, data, grad);
    Eigen::MatrixXd hess(rows, rows);
    hessian(params, grad, data, hess);
    
    Eigen::Vector2d tempDelParams = -hess.inverse()*grad;
    
    Eigen::Vector2d tempParams;
    int i = -1;
    while(i > -20){
        *t = std::pow(2, i);
        delParams = tempDelParams * *t;
        tempParams = params + delParams;
        double temp_fx = fx_func(tempParams, data);
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

    double f_k = fx_func(params, data);
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

        f_k = fx_func(params, data);
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

void legend_loc(const std::string &location="best"){

    PyObject* kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "loc", PyString_FromString(location.c_str()));

    PyObject* res = PyObject_Call(plt::detail::_interpreter::get().s_python_function_legend, plt::detail::_interpreter::get().s_python_empty_tuple, kwargs);
    if(!res) throw std::runtime_error("Call to legend() failed.");

    Py_DECREF(kwargs);
    Py_DECREF(res);
}

// params[0] = a, params[1] = m
double fx_power(const Eigen::VectorXd &params, const Eigen::VectorXd &x){
    return params[0]*std::pow(x[0], params[1]);
}

int main(){

    int data_size = 6;
    Eigen::MatrixXd data(2, data_size);
    data << 1.0, 4.5, 9.0, 20, 74, 181,
         3.0, 49.4, 245, 1808, 2.2e4, 7.3e4;
    data.transposeInPlace();
    std::cout << data.transpose() << std::endl;



    Eigen::Vector2d params0(2.6, 2);
    Eigen::VectorXd soln0(2);
    cout << "Initial condition a = " << params0[0] << " m = " << params0[1] << endl;
    printf("\nQuasi-Newton With Line Search, Magnitude normalization\n");

    opt_func = fx_power;
    fx_func = fx_lagrange;
    grad_func = grad_power_lagrange;
    newton(params0, soln0, data, false);

    fx_func = fx_raw;
    grad_func = grad_power_raw;
    //Eigen::Vector2d params1(2.8, 2.8);
    Eigen::VectorXd soln1(2);
    cout << "Initial condition a = " << params0[0] << " m = " << params0[1] << endl;
    printf("\nQuasi-Newton With Line Search, No normalization\n");
    newton(params0, soln1, data, false);

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

    plt::semilogy(x, y, "*");
    plt::named_semilogy("Normalized Y", x_plot, y0_plot, "r");
    plt::named_semilogy("Raw Y", x_plot, y1_plot, "b");
    //plt::legend();
    legend_loc("right");
    plt::save("fit_plot.png");
    //plt::show();

    return 0;
}
    
