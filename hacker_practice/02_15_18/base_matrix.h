#ifndef BASE_MATRIX_H
#define BASE_MATRIX_H

#include <iostream>
#include <cstdio>

#define M_SUCCESS 0
#define M_DIM_ERROR 1

template<class T>
class BMatrix{
    protected:
        int rowLen;
        int colLen;
        int numel;
    public:
        BMatrix() {};

        virtual int retrieveElement(int row, int col, T *element) const = 0;
        virtual int getRow(int row, T *res) const = 0;
        virtual int getCol(int row, T *res) const = 0;

        int getRowLen() const { return rowLen; };
        int getColLen() const { return colLen; };
        int getNumel() const { return numel; };
        T getNullElement() const { return 0; };

        int productAx(T *x, T *res) const;

        void print() const {
            T element = getNullElement();
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
};

#endif
