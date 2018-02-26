#include <cstdio>
#include <cmath>

#include <string.h>
#include <fstream>
#include <iostream>

#include <time.h>

#include "full_matrix.h"
#include "sparse_matrix.h"

int main(){


    FILE *memFile = fopen("memplus.mtx", "r");


    struct timespec start, stop;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);

    //FILE *memFile = fopen("mem2", "r"); // test smaller file
    int BUF_SIZE = 255;
    char *str_buf = new char[BUF_SIZE];

    fgets(str_buf, BUF_SIZE, memFile); // discard comment

    // read matrix dimensions
    int colLen, rowLen, numel;
    fgets(str_buf, BUF_SIZE, memFile); 
    sscanf(str_buf, "%d %d %d", &rowLen, &colLen, &numel);

    double *value = new double[numel];
    int *colInd = new int[numel];
    int *rowPtr = new int[rowLen+1];
    rowPtr[0] = 0;

    int  curRow = 0;
    int row, col;
    double val;

    int rowCount = 0;
    for (int i = 0; i < numel; i++){
        fgets(str_buf, BUF_SIZE, memFile);
        sscanf(str_buf, "%d %d %lg", &col, &row, &val);
        value[i] = val;
        colInd[i] = col-1;

        row-=1;

        if (row > curRow){
            rowPtr[curRow+1] = rowPtr[curRow] + rowCount;
            for (int j = curRow+2; j < row+1; j++) // advance rowPtr
                rowPtr[j] = rowPtr[j-1];
            rowCount = 0;
            curRow = row;
        }
        rowCount++;
    }
    for (int j = curRow; j < rowLen-1; j++) // advance rowPtr
        rowPtr[j+1] = rowPtr[j];
    rowPtr[rowLen] = rowPtr[rowLen-1] + rowCount;

    // 
    SMatrix<double> m1(rowLen, colLen, numel, value, rowPtr, colInd);
    m1.rowPermute(0, 2);
    m1.rowPermute(0, 4);
    m1.rowPermute(9, 2999);
    m1.rowPermute(4999, 9999);
    m1.rowPermute(5, 14999);

    m1.rowScale(1, 3, 3.0);
    m1.rowPermute(1, 4);
    m1.rowScale(4, 3, -3.0);

    SMatrix<double> m(rowLen, colLen, numel, value, rowPtr, colInd);
    double *x = new double[colLen];
    std::fill(x, x+colLen, 1);
    double *res = new double[colLen];
    m.productAx(x, res);

    // find sum using product
    double sumProduct = 0;
    for (int i = 0; i < colLen; i++)
        sumProduct += res[i];
    // find sum of values
    double sumValues = 0;
    for (int i = 0; i < numel; i++)
        sumValues += value[i];

    printf("### Part 2 \n");

    printf("Sum using product = %g, Sum using values = %g  \n", sumProduct, sumValues);
    printf("Residue = %g  \n", std::abs(sumProduct-sumValues));

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);

    double elapsed = (stop.tv_sec - start.tv_sec) * 1e6 + (stop.tv_nsec - start.tv_nsec) / 1e3;    // in microseconds
    printf("Time elapsed = %g ms  \n", elapsed/1000);

    printf("Peak Memory usage using /proc/self/status -> VmPeak :  \n"); 

    std::ifstream statusFile( "/proc/self/status", std::ifstream::in);
    bool found = 0;
    std::string statusLine;
    while (!found){
        if (! std::getline(statusFile, statusLine) ){
            printf("Memory usage not found  \n");
            found = true;
        }
        //std::cout << statusLine << std::endl;
        if (statusLine.find("VmPeak:") !=std::string::npos){
            std::cout << statusLine << "  \n";
            found = true;
        }
    }


    fclose(memFile);
    delete [] value;
    delete [] colInd;
    delete [] rowPtr;

    return 0;
}
