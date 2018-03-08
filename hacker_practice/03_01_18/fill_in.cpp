#include <cstdio>
#include "full_matrix.h"
#include "sparse_matrix.h"

/*
* row[j] = scale*row[i] + row[j]
*/
void eliminate(int i, int j, double scale, SMatrix<double> &mat, double *b) {
    mat.rowScale(i, j, scale);
    b[j] = scale*b[i] + b[j];
}

void rowScale(int row, double scale, SMatrix<double> &mat, double *b) {
    mat.rowScale(row, scale);
    b[row] = scale*b[row];
}

double get(int row, int col, const SMatrix<double> &mat){
    double res = 0;
    mat.retrieveElement(row, col, &res);
    return res;
}

void swapRows(int i, int j, SMatrix<double> &mat, double *b){
    double temp = b[i];
    b[i] = b[j];
    b[j] = temp;
    mat.rowPermute(i, j);
}

/* N is number of rows */
int gaussElim(int col, const SMatrix<double> &inMat,
        double *b, int N, SMatrix<double> **res){
    double bOriginal[N];
    std::copy(b, b+N, bOriginal);
    int minCount = inMat.getNumel() + 1; // +1 for when matrix is full
    int noRows = inMat.getRowLen();
    SMatrix<double> *bestPivotMat = new SMatrix<double>(inMat);

    for (int i = noRows-1; i >= col; i--){
        int pivot = i;
        double pivotVal = get(i, col, inMat);

        if (! inMat.isZero(pivotVal)) {
            SMatrix<double> pivotMat(inMat);
            double bCopy[N];
            std::copy(bOriginal, bOriginal+N, bCopy);
            rowScale(i, 1/pivotVal, pivotMat, bCopy);

            for (int j = noRows-1; j >=col; j--) {
                double element = 0;
                if ( j != i && !pivotMat.isZero(element = get(j, col, pivotMat))) {
                    eliminate(i, j, -element, pivotMat, bCopy);
                }
            }

            // find new length after pivoting
            int afterPivot = pivotMat.getCount();
            if (afterPivot < minCount) {
                delete bestPivotMat;
                bestPivotMat = new SMatrix<double>(pivotMat);
                swapRows(col, pivot, *bestPivotMat, bCopy);
                std::copy(bCopy, bCopy+N, b);
                minCount = afterPivot;
            }
            printf("Count for pivot %f is %d \n", pivotVal, afterPivot);

        }
    }
    if (*res != NULL)
        delete *res;
    *res = bestPivotMat;
    return M_SUCCESS;
}

int backSubstitute(const SMatrix<double> &inMat, const double *b, double *res){
    int N = inMat.getColLen();
    for (int row = N-1; row >= 0; row--){
        double sum = b[row];
        for (int j = N-1; j > row; j--){
            sum -= get(row, j, inMat)*res[j];
        }
        res[row] = sum/get(row, row, inMat);
    }
    return M_SUCCESS;
}

    

int main(){
    double aArr [] = {1, 2, 0, 0, 3, 
                      4, 5, 6, 0, 0, 
                      0, 7, 8, 0, 9,
                      0, 0, 0, 10, 0,
                      11, 0, 0, 0, 12};

    double b [] = {5, 4, 3, 2, 1};
    FMatrix<double> fA(5, 25, aArr);
    SMatrix<double> A(fA);
    SMatrix<double> *res = NULL;

    gaussElim(0, A, b, 5, &res);
    res->print();
    gaussElim(1, *res, b, 5, &res);
    res->print();
    gaussElim(2, *res, b, 5, &res);
    res->print();
    gaussElim(3, *res, b, 5, &res);
    res->print();
    gaussElim(4, *res, b, 5, &res);
    res->print();

    double soln[5] = {0};
    backSubstitute(*res, b, soln);
    for (int i = 0; i < 5; i++){
        printf("x[%d] = %g \n", i, soln[i]);
    }

    if (res != NULL)
        delete res;

    return 0;

}
