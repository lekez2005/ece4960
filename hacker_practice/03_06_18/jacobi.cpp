#include <cstdio>
#include <Eigen/Core>
#include "full_matrix.h"
#include "sparse_matrix.h"
#include "gauss_solve.h"

using namespace std;

int main(){

    cout << endl << "Jacobi Solution" << endl;

    double aArr[] = {-4, 1, 0, 0, 1,
                    4, -4, 1, 0, 0,
                    0, 1, -4, 1, 0,
                    0, 0, 1, -4, 1,
                    1, 0, 0, 1, -4};

    double b[] = {1, 0, 0, 0, 0};

    SMatrix<double> A(FMatrix<double>(5, 25, aArr));
    printf("Sparse representation of A: \n");
    A.print();
    
    double soln[5] = {0};
    printf("Using direct solver, ");

    gaussSolve(A, b, soln);

    for (int i = 0; i < 5; i++){
        printf("x[%d] = %g \t", i, soln[i]);
    }
    printf("\n");

    Eigen::Map<Eigen::MatrixXd> aEigen(aArr, 5, 5);
    cout << "\nEigen representation is " << endl << aEigen << endl;

    Eigen::MatrixXd Dinv = aEigen.diagonal().asDiagonal().inverse();
    cout << endl << "Diagonal inverse is " << endl << Dinv << endl;

    Eigen::MatrixXd L = aEigen.triangularView<Eigen::StrictlyLower>();
    L = -L;
    cout << endl << "L is " << endl << L << endl;

    Eigen::MatrixXd U = aEigen.triangularView<Eigen::StrictlyUpper>();
    U = -U;
    cout << endl << "U is " << endl << U << endl;

    Eigen::MatrixXd LU = L + U;

    //Eigen::Matrix<double, 5,1> bEigen(b);
    //Eigen::Matrix<double, 5,1> bEigen(b);
    Eigen::Map<Eigen::VectorXd> bEigen(b, 5);
    cout << "b is " << endl << bEigen.transpose() << endl;

    Eigen::VectorXd x_k = Dinv * bEigen;

    int max_iterations = 50;
    int iterations = 0;
    double tolerance = 1e-8;

    while (iterations < max_iterations){
        x_k = Dinv*LU*x_k + Dinv*bEigen;
        Eigen::VectorXd residual = bEigen - aEigen*x_k;
        if (residual.squaredNorm() < tolerance)
            break;
        iterations++;
    }

    std::cout << "Solution after " << iterations << " iterations is " << x_k.transpose() << std::endl;

    Eigen::MatrixXd DLU = Dinv*LU;
    cout << "Dinv*(L+U) is " << endl << DLU << endl;

    //cout << "First norm is " << DLU.cwiseAbs().colwise().maxCoeff().max() << endl;
    Eigen::MatrixXd maxCols =  DLU.cwiseAbs().colwise().maxCoeff();
    cout << "First norm is " << maxCols.maxCoeff() << endl;

    Eigen::MatrixXd maxRows =  DLU.cwiseAbs().rowwise().maxCoeff();
    cout << "Infinity norm is " << maxRows.maxCoeff() << endl;


    return 0;
}
