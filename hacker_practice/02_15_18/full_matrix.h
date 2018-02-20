#ifndef FULL_MATRIX_H
#define FULL_MATRIX_H

#include "base_matrix.h"
#include <stdexcept>

template<class T>
class FMatrix : public BMatrix<T>{
    private:
        T *value;
    public:
        int retrieveElement(int row, int col, T *element) const;
        int getRow(int row, T *res) const;
        int getCol(int col, T *res) const;
        int productAx(T *x, T *res) const;

        FMatrix(int rowLen, int numel, T *value);

};

#include "full_matrix.cpp"

#endif
