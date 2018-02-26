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

        int productAx(const T *x, T *res) const;

        void print(const char* delimeter, const char* lineStart, const char* lineEnd){
            T element = getNullElement();
            for (int i = 0; i < rowLen; i++) {
                std::cout << lineStart;
                for (int j = 0; j < colLen; j++) {
                    retrieveElement(i, j, &element);
                    if (j < colLen-1) {
                        std::cout << element << delimeter;
                    }else {
                        std::cout << element;
                    }
                }
                std::cout << lineEnd;
            }
        }

        void print(){
            print(", \t", "", "\n");
        }

        void printTex(){
            print(" & ", "", " \\\\ \n");
        }

        void printMd(){
            print(" , ", "", "   \n");
        }


        template< typename A, typename B>
        static int norm(const BMatrix<A> &m1, const BMatrix<B> &m2, double *res){
            if (m1.getRowLen() != m2.getRowLen() || m1.getColLen() != m2.getColLen())
                return M_DIM_ERROR;
            *res = 0.0;
            for (int row = 0; row < m1.getRowLen(); row++){
                for (int col = 0; col < m1.getColLen(); col++) {
                    A m1Element;
                    B m2Element;
                    m1.retrieveElement(row, col, &m1Element);
                    m2.retrieveElement(row, col, &m2Element);
                    double difference = (double) ((m1Element-m2Element)*(m1Element-m2Element));
                    *res += difference;
                }
            }
            return M_SUCCESS;
        }
        int norm(const BMatrix<T> &m1, double *res){
            return this->norm<T, T>(*this, m1, res);
        }
};

#endif
