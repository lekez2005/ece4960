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
    FMatrix<double> fm(4, numel, value_d);
    printf("\n");
    fm.print();

    printf("\n");
    SMatrix<double> s1(fm, 0.01);
    s1.print();

    return 0;
}
