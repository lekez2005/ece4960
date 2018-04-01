#include "matplotlibcpp.h"
#include <cstdio>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
namespace plt = matplotlibcpp;

void spline(const Eigen::MatrixXd &data, Eigen::MatrixXd &params){
    int rows = data.rows();
    Eigen::VectorXd x = data.col(0);
    Eigen::VectorXd y = data.col(1);

    Eigen::MatrixXd coeffs(rows, rows);
    Eigen::VectorXd b(rows);

    // first data point
    double x10_diff = x[1]-x[0];
    coeffs(0, 0) = 2/x10_diff;
    coeffs(0, 1) = 1/x10_diff;
    b[0] = 3*(y[1] - y[0])/x10_diff;

    //last data point
    double xlast_diff = x[rows-1]-x[rows-2];
    coeffs(rows-1, rows-2) = 2/xlast_diff;
    coeffs(rows-1, rows-1) = 1/xlast_diff;
    b[rows-1] = 3*(y[rows-1] - y[rows-2])/xlast_diff;

    for (int i = 1; i < rows-1; i++){
        double x_back_diff = x[i] - x[i-1];
        double x_forw_diff = x[i+1] - x[i];
        coeffs(i, i-1) = 1/x_back_diff;
        coeffs(i, i+1) = 1/x_forw_diff;
        coeffs(i, i) = 2/x_back_diff + 2/x_forw_diff;
        b[i] = 3*(y[i] - y[i-1])/(x_back_diff * x_back_diff) + 
            3*(y[i+1] - y[i])/(x_forw_diff * x_forw_diff);
    }
    Eigen::VectorXd k = coeffs.inverse() * b; 
    params.col(0) = k;
    for (int i = 1; i < rows; i++){
        params(i, 1) = k[i-1]*(x[i]-x[i-1]) - (y[i]-y[i-1]);
        params(i, 2) = -k[i]*(x[i]-x[i-1]) + (y[i]-y[i-1]);
    }

    printf("Coefficient Matrix \n");
    cout << coeffs << endl;
    cout <<"k = " << k.transpose() << endl;
    cout << "a = " << params.col(1).transpose() << endl;
    cout << "b = " << params.col(2).transpose() << endl;
}

int get_x_index(double x, const Eigen::VectorXd x_data){
    int rows = x_data.rows();
    if (x <= x_data[1])
        return 1;
    else if(x >= x_data[rows-2])
        return rows-1;
    else{
        for (int i = 1; i < rows-2; i++){
            if ( x >= x_data[i] && x <= x_data[i+1])
                return i+1;
        }
    }
    return -1;
}

double eval_spline(double x, const Eigen::MatrixXd &data, const Eigen::MatrixXd &params){
    Eigen::VectorXd x_data = data.col(0);
    Eigen::VectorXd y = data.col(1);
    int index = get_x_index(x, x_data);
    double t = (x - x_data[index-1])/(x_data[index]-x_data[index-1]);
    double a = params(index, 1);
    double b = params(index, 2);
    return (1-t)*y[index-1] + t*y[index] + t*(1-t)*(a*(1-t) + b*t);
}


int main(){
    Eigen::MatrixXd data(2, 6);
    Eigen::MatrixXd params(6, 3);
    data << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
            0.0, 1.0, 0.0, -1.0, 0.0, 1.0;

    data.transposeInPlace();
    printf("Input data: \n");
    cout << data.transpose() << endl;

    spline(data, params);

    Eigen::VectorXd x_vec = data.col(0);
    std::vector<double> x(x_vec.data(), x_vec.data() + x_vec.rows());
    Eigen::VectorXd y_vec = data.col(1);
    std::vector<double> y(y_vec.data(), y_vec.data() + y_vec.rows());

    
    int points_plot = 100;
    vector<double> x_plot(points_plot);
    vector<double> y_plot(points_plot);
    double delta_x = (x[x_vec.rows()-1]-x[0])/points_plot;
    for (int i = 0; i < points_plot; i++){
        double x_point = x[0] + delta_x*i;
        x_plot[i] = x_point;
        y_plot[i] = eval_spline(x_point, data, params);
    }


    plt::plot(x, y, "*");
    plt::plot(x_plot, y_plot, "r");
    plt::save("spline_plot.png");
    //plt::show();


    return 0;
}



