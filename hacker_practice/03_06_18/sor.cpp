#include <cstdio>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/LU>

using namespace std;

int main(){

    cout << endl << "SOR Solution" << endl;

    double aArr[] = {-4, 1, 0, 0, 1,
                    4, -4, 1, 0, 0,
                    0, 1, -4, 1, 0,
                    0, 0, 1, -4, 1,
                    1, 0, 0, 1, -4};

    double b[] = {1, 0, 0, 0, 0};

    Eigen::Map<Eigen::MatrixXd> aEigen(aArr, 5, 5);
    cout << "\nEigen representation is " << endl << aEigen << endl;

    Eigen::MatrixXd Dinv = aEigen.diagonal().asDiagonal().inverse();
    cout << endl << "Diagonal inverse is " << endl << Dinv << endl;

    Eigen::Map<Eigen::VectorXd> bEigen(b, 5);
    cout << "b is " << endl << bEigen.transpose() << endl;

    Eigen::VectorXd x_k = Dinv * bEigen;
    Eigen::VectorXd x_k_1 = Dinv * bEigen;
    double w = 0.1;

    int max_iterations = 100;
    int iterations = 0;
    double tolerance = 1e-8;

    printf("k \t xk[0] \t xk_1[0] \t ||delta x|| \t ||residual|| \n");

    while (iterations < max_iterations){

        Eigen::VectorXd residual = bEigen - aEigen*x_k;
        x_k_1 = x_k + w*Dinv*residual;

        printf("%d \t %f \t %f \t %g \t %f \n", iterations+1, x_k[0], x_k_1[0], 
                (x_k_1-x_k).squaredNorm(), residual.squaredNorm());

        x_k = x_k_1;

        if (residual.squaredNorm() < tolerance)
            break;
        iterations++;
    }

    cout << "Solution after " << iterations << " iterations is " << x_k.transpose() << endl;
    //Eigen::VectorXd sol = aEigen.inverse()*bEigen;
    cout << "Solution with Direct solver " << (aEigen.inverse()*bEigen).transpose() << endl;


    return 0;
}
