#include <cstdio>
#include <iostream>
#include "full_matrix.h"
#include "sparse_matrix.h"

int main(){
    
    printf("### Part I: Ground Truth \n");

    double value[] =  {1.0, 2, 0, 0, 3,
        4, 5, 6, 0, 0, 
        0, 7, 8, 0, 9,
        0, 0, 0, 10, 0,
        11, 0, 0, 0, 12};

    int rowLen = 5; 
    int colLen = 5;
    int numel = rowLen*colLen;

    FMatrix<double> aFull(rowLen, numel, value);
    SMatrix<double> aSparse(aFull);
    printf("Verify original full matrix prints correctly as a sparse matrix  \n");
    aSparse.printMd();

    double norm = 0;
    double *value2 = new double[numel];
    std::copy(value, value+numel, value2);
    value2[2] = 2.5;
    FMatrix<double> a1Full(rowLen, numel, value2);
    BMatrix<double>::norm<double, double>(aSparse, a1Full, &norm);
    printf("\nVerify norm works correctly using slightly modified test matrix  \n");
    a1Full.printMd();
    printf("Expected norm is 6.25, Actual norm is %g \n", norm);

    // permute
    FMatrix<double> aFullPerm(aFull);
    aFullPerm.rowPermute(0, 2);
    SMatrix<double> aSparsePerm(aFull);
    aSparsePerm.rowPermute(0, 2);
    printf("\nResult of permute (0, 2):  \n");
    aSparsePerm.printMd();
    aFullPerm.norm(aSparsePerm, &norm);
    printf("Norm of sparse and full matrix solutions: \t %g  \n", norm);

    printf("\nResult of permute (0, 4):   \n");
    aFullPerm.rowPermute(0, 4);
    aSparsePerm.rowPermute(0, 4);
    aSparsePerm.printMd();
    aFullPerm.norm(aSparsePerm, &norm);
    printf("Norm of sparse and full matrix solutions:   \t %g\n", norm);
   
    // test row scale
    FMatrix<double> aFullScale(aFull);
    SMatrix<double> aSparseScale(aFull);
    aFullScale.rowScale(0, 3, 3.0);
    aSparseScale.rowScale(0, 3, 3.0);
    printf("\nResult of 3.0*row[0] + row[3]  \n");
    aSparseScale.printMd();
    aFullScale.norm(aSparseScale, &norm);
    printf("Norm of sparse and full matrix solutions: \t %g  \n", norm);


    FMatrix<double> aFullScale1(aFull);
    SMatrix<double> aSparseScale1(aFull);
    aFullScale1.rowScale(4, 1, -4.4);
    aSparseScale1.rowScale(4, 1, -4.4);
    printf("\nResult of -4.4*row[4] + row[1]   \n");
    aFullScale1.printMd();
    aFullScale1.norm(aSparseScale1, &norm);
    printf("Norm of sparse and full matrix solutions: \t %g  \n", norm);

    // product Ax
    double x[] = {5, 4, 3, 2, 1};
    double *resFull = new double[rowLen];
    double *resSparse = new double[rowLen];
    FMatrix<double> aFullProduct(aFull);
    SMatrix<double> aSparseProduct(aFull);
    aFullProduct.productAx(x, resFull);
    aSparseProduct.productAx(x, resSparse);
    printf("\nFor x =   \n");
    FMatrix<double>(1, rowLen, x).print();
    printf("\nResult of A*x =   \n");
    FMatrix<double>(1, rowLen, resSparse).print();
    aFullProduct.norm(aSparseProduct, &norm);
    printf("<br /> Norm of sparse and full matrix solutions: \t %g  \n", norm);





    return 0;
}
