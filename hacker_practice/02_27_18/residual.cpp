#include <cstdio>
#include <cmath>


//row major order

void product(const double *mat, const double *x, double *res){
    res[0] = mat[0]*x[0] + mat[1]*x[1];
    res[1] = mat[2]*x[0] + mat[3]*x[1];
}

void solve(const double * mat, const double *x, double *res){
    double det = (mat[0]*mat[3])-(mat[1]*mat[2]);
    double inv [] = {0, 0, 0, 0};
    inv[0] = mat[3]/det;
    inv[1] = -mat[2]/det;
    inv[2] = -mat[1]/det;
    inv[3] = mat[0]/det;
    product(inv, x, res);

}

void printArr(const double *mat, int len){
    printf("{ ");
    for (int i = 0; i < len-1; i++){
        printf("%g, ", mat[i]);
    }
    printf("%g }", mat[len-1]);
}



void residual(const double *mat, const double *x){
    double sol [] = {0, 0};
    solve(mat, x, sol);
    printf("For A = ");
    printArr(mat, 4);
    printf(", solution is ");
    printArr(sol, 2);
    double prod [] = {0, 0};
    product(mat, sol, prod);
    double resid = sqrt(pow(prod[0]-x[0], 2) + pow(prod[1]-x[1], 2));
    printf(", residual = %g \n", resid);
}



int main(){
    double mat [] = {100, 99, 99, 98.01};
    double x [] = {199, 197};

    for (int i = 2; i < 10; i++){
        mat[3] = 98.01 - pow(10, -i);
        residual(mat, x);
    }

    return 0;
}
