#ifndef FULL_MATRIX_H
#define FULL_MATRIX_H

#include "base_matrix.h"
#include <stdexcept>

template<class T>
class SMatrix;

template<class T>
class FMatrix : public BMatrix<T>{
    private:
        T *value;
        void initFromMat(const BMatrix<T> &inMat);
    public:
        int retrieveElement(int row, int col, T *element) const;
        T retrieveElement(int row, int col) const;
        int getRow(int row, T *res) const;
        int getCol(int col, T *res) const;

        int getIndex(int row, int col) const;
        int productAx(const T *x, T *res) const;
        int rowPermute(int i, int j);
        int rowScale(int i, int j, T a);

        FMatrix(int rowLen, int numel, T *value);
        FMatrix(const FMatrix<T> &inMat);
        FMatrix(const SMatrix<T> &inMat);

        ~FMatrix(){
            delete [] this->value;
        }

};

#include "full_matrix.cpp"

#endif
