#include <cstdio>
#include <iostream>
#include "full_matrix.h"
#include "sparse_matrix.h"

int main(){
    //struct Matrix matrix;
    int rowLen = 5;
    int colLen = 5;
    int value [] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}; 
    int rowPtr [] = {0, 3, 6, 9, 10, 12};
    int colInd [] = {0, 1, 4, 0, 1, 2, 1, 2, 4, 3, 0, 4};
    int numel = sizeof(value)/sizeof(value[0]);

    SMatrix<int> s(rowLen, colLen, numel, value, rowPtr, colInd);
    s.print();
    //s.productAx(NULL, NULL);

    double value_d [] = {1.6, 0.001, 3, 4, 0.02, 6, 0.00001, 8, 9, 10, 11, 12}; 
    int frowLen = 4;
    FMatrix<double> fm(frowLen, numel, value_d);
    int fcolLen = fm.getColLen();
    printf("\n");
    fm.print();

    printf("\n");
    SMatrix<double> s1(fm, 0.01);
    s1.print();

    printf("Permute full matrix \n");
    FMatrix<double> f2(s1);
    f2.rowPermute(2, 3);
    f2.print();

    printf("Permute Sparse matrix \n");
    SMatrix<double> s2(s1);
    s2.rowPermute(2, 3);
    s2.print();

    printf("\n");
    printf("Full Matrix row Scale \n");
    FMatrix<double> f3(s1);
    f3.rowScale(1, 0, 0.5);
    f3.print();

    printf("Sparse Matrix row Scale \n");
    SMatrix<double> s3(s1);
    s3.rowScale(1, 0, 0.5);
    s3.print();


    int *x = new int[rowLen];
    int *res = new int[rowLen];
    std::fill_n(x, rowLen, 2);
    FMatrix<int>(s).productAx(x, res);
    FMatrix<int> resMatrix(1, rowLen, res);
    printf("Full Matrix product: \n");
    resMatrix.print();

    printf("Sparse Matrix product \n");
    s.productAx(x, res);
    FMatrix<int> resMatrix2(1, rowLen, res);
    resMatrix2.print();

    delete [] x;


    return 0;
}
