#include <cstdio>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>

int main(){
    Eigen::MatrixXd dirichlet(4, 4);
    dirichlet << -2, 1, 0, 0,
            1, -2, 1, 0,
            0, 1, -2, 1,
            0, 0, 1, -2;
    Eigen::Vector4d f(0.04, 0, 0, 0);
    Eigen::MatrixXd neumann(dirichlet);
    neumann(3, 3) = -1;

    Eigen::Vector4d diri_soln = dirichlet.inverse()*f;
    Eigen::Vector4d neumann_soln = neumann.inverse()*f;

    std::cout << "Dirichlet solution:\t" << diri_soln.transpose() << std::endl;
    std::cout << "Neumann solution:\t" << neumann_soln.transpose() << std::endl;

    return 0;

}
