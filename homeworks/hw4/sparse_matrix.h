#ifndef SMATRIX_H
#define SMATRIX_H

#include "base_matrix.h"
#include <iostream> 

template<class T>
class SMatrix : public BMatrix<T>{
    private:
        int *rowPtr;
        int *colInd;
        int storageSize;
        T *value;
        void initialize_from_mat(const BMatrix<T> &inMat, T tolerance); 
        void initialize(int rowLen, int colLen, int valLength, T *value, int *rowPtr, int *colInd);
    public:
        int retrieveElement(int row, int col, T *element) const;
        int getRow(int row, T *res) const;
        int getCol(int col, T *res) const;

        int rowPermute(int i, int j);
        int rowScale(int i, int j, T a);
        int productAx(const T *x, T *res) const;

        bool isZero(T element, T tolerance) const;
        bool isZero(T element) const;

        ~SMatrix(){
            delete [] this->value;
            delete [] this->rowPtr;
            delete [] this->colInd;
        }
        SMatrix(int rowLen, int colLen, int valLength, T *value, int *rowPtr, int *colInd);
        SMatrix(const SMatrix<T> &inMat);
        SMatrix(const FMatrix<T> &inMat);
        SMatrix(const SMatrix<T> &inMat, T tolerance);
        SMatrix(const FMatrix<T> &inMat, T tolerance);

};

#include "sparse_matrix.cpp"

#endif
