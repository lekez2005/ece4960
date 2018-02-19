#include <cstdio>
#include <iostream>
#include <cmath>
#include "sparse_matrix.h"

template <class T>
int SMatrix<T>::retrieveElement(int row, int col, T *element) {
    int colOffset = rowPtr[row]; // skip to this column
    T rowSize = rowPtr[row+1] - rowPtr[row];

    *element = 0;

    for (int i = 0; i < rowSize; i++) {
        if (colInd[colOffset + i] == col) {
            *element = value[colOffset + i];
            break;
        }
    }
    return M_SUCCESS;
}


template <class T>
void SMatrix<T>::print(){
    T element = 0;
    for (int i = 0; i < rowLen; i++) {
        for (int j = 0; j < colLen; j++) {
            retrieveElement(i, j, &element);
            if (j != colLen) {
                std::cout << element << ", \t";
            }else {
                std::cout << element << std::endl;
            }
        }
        printf("\n");
    }
}

int main(){
    //struct Matrix matrix;
    int rowLen = 5;
    int colLen = 5;
    int value [] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}; 
    int rowPtr [] = {0, 3, 6, 9, 10, 12};
    int colInd [] = {0, 1, 4, 0, 1, 2, 1, 2, 4, 3, 0, 4};

    SMatrix<int> s(rowLen, colLen, value, rowPtr, colInd);
    s.print();

    double value_d [] = {1.6, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}; 
    SMatrix<double> sf(rowLen, colLen, value_d, rowPtr, colInd);
    sf.print();

    return 0;
}
