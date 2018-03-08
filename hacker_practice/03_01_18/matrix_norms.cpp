#include <cstdio>
#include <cmath>

#include "full_matrix.h"
#include "sparse_matrix.h"

void printArr(double *arr, int N){
    for (int i = 0; i < N; i++)
        printf("%g \t ", arr[i]);
    printf("\n");
}

double absSum(double *arr, int N){
    double res = 0;
    for (int i = 0; i < N; i++)
        res += std::abs(arr[i]);

    return res;
}

double arrMax(double *arr, int N){
    double max = arr[0];
    for (int i = 1; i < N; i++){
        if (arr[i] > max)
            max = arr[i];
    }
    return max;
}

double norm_ifty(const BMatrix<double> &inMat){
    int noRows = inMat.getRowLen();
    int colLen = inMat.getColLen();
    double maxSum = 0;
    double *tempRow = new double[colLen];
    for (int row = 0; row < noRows; row++){
        inMat.getRow(row, tempRow);
        double sum = absSum(tempRow, colLen);
        if (sum  > maxSum)
            maxSum = sum;
    }
    delete [] tempRow;
    return maxSum;
}

double norm_one(const BMatrix<double> &inMat){
    int noRows = inMat.getRowLen();
    int colLen = inMat.getColLen();
    double maxSum = 0;
    double *tempCol = new double[noRows];

    for (int col = 0; col < colLen; col++){
        inMat.getCol(col, tempCol);
        double sum = absSum(tempCol, noRows);
        if (sum  > maxSum)
            maxSum = sum;
    }
    delete [] tempCol;
    return maxSum;
}

int main(){

    double aArr [] = {1, 2, 0, 0, 3, 
                      4, 5, 6, 0, 0, 
                      0, 7, 8, 0, 9,
                      0, 0, 0, 10, 0,
                      11, 0, 0, 0, 12};

    FMatrix<double> fA(5, 25, aArr);
    fA.print();

    BMatrix<double> *mFull = new FMatrix<double>(fA);
    BMatrix<double> *mSparse = new SMatrix<double>(fA);
    printf("Max Infty norm using full matrix:\t %g\n", norm_ifty(*mFull));
    printf("Max Infty norm using sparse matrix:\t %g\n", norm_ifty(*mSparse));

    printf("Max Norm-1 using full matrix:\t %g\n", norm_one(*mFull));
    printf("Max Norm-1 using sparse matrix:\t %g\n", norm_one(*mSparse));

    return 0;
}
