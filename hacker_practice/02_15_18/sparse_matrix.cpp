#include <cstdio>
#include <cmath>

typedef struct Matrix {
    int rowLen;
    int colLen;
    int *value;
    int *rowPtr;
    int *colInd;
} Matrix;

int retrieveElement(Matrix* matrix, int row, int col) {
    int rowLen = matrix->rowPtr[row+1] - matrix->rowPtr[row];
    int rowPtr = matrix->rowPtr[row];
    for (int i = 0; i < rowLen; i ++) {
        if (matrix->colInd[rowPtr + i] == col) {
            return matrix->value[rowPtr + i];
        }
    }
    return 0;
}

void printMatrix(Matrix* matrix) {
    for (int i = 0; i < matrix->rowLen; i++) {
        for (int j = 0; j < matrix->colLen; j++) {
            int element = retrieveElement(matrix, i, j);
            if (j != matrix->colLen) {
                printf("%d, \t", element);
            }else {
                printf("%d \n", element);
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

    Matrix matrix = {rowLen, colLen, value, rowPtr, colInd};

    printMatrix(&matrix);


    return 0;
}
