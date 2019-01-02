#include <iostream>

#include <cmath>
#include <Eigen/Dense>

typedef Eigen::Matrix<double, 9, 9> Matrix9d;
typedef Eigen::Matrix<double, 9, 1> Vector9d;

int main(){
    // constants
    double q = 1.602e-19;
    double eps = 11.68*8.85e-12;
    double boltzmann = 1.38e-23;
    double ni_cm = 1e10; // intrinsic silicon in cm^-3
    double ni = ni_cm*1e6; // intrinsic silicon in m^-3
    double T = 298;
    double Vt = boltzmann*T/q;
    double Ld = std::sqrt(eps*Vt/q/ni); // remove factor of 2 to keep RHS the same

    double h = 10e-6/Ld;

    double Nd = 1e15/ni_cm;
    double Na = 1e17/ni_cm;


    int vec_length = 9;
    // homogeneous poisson reverse bias at 1V
    Matrix9d homogen_mat;
    homogen_mat << -4, 1, 0, 1, 0, 0, 0, 0, 0,
                   1, -4, 1, 0, 1, 0, 0, 0, 0,
                   0, 1, -4, 0, 0, 1, 0, 0, 0,
                   1, 0, 0, -4, 1, 0, 1, 0, 0,
                   0, 1, 0, 1, -4, 1, 0, 1, 0,
                   0, 0, 1, 0, 1, -4, 0, 0, 1,
                   0, 0, 0, 1, 0, 0, -4, 1, 0,
                   0, 0, 0, 0, 1, 0, 1, -4, 1,
                   0, 0, 0, 0, 0, 1, 0, 1, -4;
    Vector9d homogen_b;

    for (int i = 0; i < vec_length; i++){
        if (i == 0)
            homogen_b(i) = h*h*(Nd-2*1/Vt);
        else if (i == 1 || i == 3)
            homogen_b(i) = h*h*Nd;
        else
            homogen_b(i) = -h*h*Na;
    }

    std::cout << "Ld = " << Ld << std::endl;
    std::cout << "F matrix: \n"  << homogen_mat << std::endl;
    std::cout  << "B vector: \n" << homogen_b.transpose() << std::endl;
    Eigen::MatrixXd soln = homogen_mat.colPivHouseholderQr().solve(homogen_b)*Vt;
    soln.resize(3, 3);
    std::cout << "Solution (V): \n" << soln << std::endl;

    // initialize with homogeneous solution
    //Vector9d V(Eigen::Map<Vector9d>(soln.data(), 9));
    Vector9d V;
    V.fill(0.3/Vt);
    //V/=Vt;

    int max_iterations = 100;
    for (int i = 0; i < max_iterations; i++){
        Vector9d currentFv = FV(V, h, homogen_mat, homogen_b);
        Matrix9d jac = jacobian(homogen_mat, V, h); 
        Vector9d delV = -jac.colPivHouseholderQr().solve(V);
        int t = 5;
        while (t > -20){
            double scale = std::pow(2, t);
            Vector9d tempV = V + scale*delV;
            Vector9d tempFv = FV(tempV, h, homogen_mat, homogen_b);
            if (std::fabs(tempFv.sum()) < std::fabs(currentFv.sum())){
                V = tempV;
                std::cout << "t = " << t << " sum(Fv) = " << tempFv.sum();
                std::cout << " V = " << Vt*V.transpose() << std::endl;
                printf("t = %d, sum(Fv) = %g \n", t, tempFv.sum());
                break;
            }
            t--;
        }
        if (t <= -20){
            std::cout << "No t found that reduces FV" << std::endl;
            break;
        }

    }
    std::cout << V.transpose() << std::endl;



    return 0;
}
