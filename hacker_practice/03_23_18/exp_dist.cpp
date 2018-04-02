#include "matplotlibcpp.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>


namespace plt = matplotlibcpp;

double rand_exp(double lambda){
    double uniform = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
    return -std::log(1-uniform)/lambda;
}

void hist(std::vector<double> y, long bins=10,std::string color="b", double alpha=1.0, int normed=1){
	PyObject* yarray = plt::get_array(y);

    PyObject* kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "bins", PyLong_FromLong(bins));
    PyDict_SetItemString(kwargs, "color", PyString_FromString(color.c_str()));
    PyDict_SetItemString(kwargs, "alpha", PyFloat_FromDouble(alpha));
    PyDict_SetItemString(kwargs, "normed", PyLong_FromLong(normed));


    PyObject* plot_args = PyTuple_New(1);

    PyTuple_SetItem(plot_args, 0, yarray);


    PyObject* res = PyObject_Call(plt::detail::_interpreter::get().s_python_function_hist, plot_args, kwargs);


    Py_DECREF(plot_args);
    Py_DECREF(kwargs);
    if(res) Py_DECREF(res);

}


void legend_loc(const std::string &location="best"){

    PyObject* kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "loc", PyString_FromString(location.c_str()));

    PyObject* res = PyObject_Call(plt::detail::_interpreter::get().s_python_function_legend, plt::detail::_interpreter::get().s_python_empty_tuple, kwargs);
    if(!res) throw std::runtime_error("Call to legend() failed.");

    Py_DECREF(kwargs);
    Py_DECREF(res);
}

int main(){
    int N = 1000;
    double lambda = 0.2;
    double bin_size = 0.1;

    double xmin = 0;
    double xmax = 10;
    long bins = std::round((xmax-xmin)/bin_size);
    double delta_x = (xmax-xmin)/N;

    std::vector<double> x_plot(N);
    std::vector<double> data(N);
    std::vector<double> ideal_exp(N);
    std::vector<double> ideal_cdf(N);
    for(int i = 0; i < N; i++){
        double x = xmin + delta_x*i;
        x_plot[i] = x;
        data[i] = rand_exp(lambda);
        ideal_exp[i] = lambda*std::exp(-lambda*x);
        ideal_cdf[i] = 1 - std::exp(-lambda*x);
    }

    // create cdf of data
    std::vector<double> sorted_data = data;
    std::sort(sorted_data.begin(), sorted_data.end());
    std::vector<double> cdf_bins_x(bins);
    std::vector<double> cdf_bins_y(bins);
    long count = 0;
    for (int i = 0; i < bins; i++){
        cdf_bins_x[i] = xmin + (double)(i)/bins*(xmax-xmin);
        while(sorted_data[count] < cdf_bins_x[i] && count < N)
            count++;
        cdf_bins_y[i] = ((double) count)/N;
    }
    
    plt::subplot(2, 1, 1);
    hist(data, bins, "r", 0.75, 1);
    plt::plot(x_plot, ideal_exp, "k");
    plt::xlim(xmin, xmax);
    plt::title("Histogram of Sampled Data");

    plt::subplot(2, 1, 2);
    plt::named_plot("ideal", x_plot, ideal_cdf, "r");
    plt::named_plot("data", cdf_bins_x, cdf_bins_y, "b");
    legend_loc();
    plt::title("CDF of Sampled Data");

    plt::save("exponential_dist.png");
    plt::show();

    return 0;
}
