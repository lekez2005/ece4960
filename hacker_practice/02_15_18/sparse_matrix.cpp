#include <vector>

template <class T>
void SMatrix<T>::initialize(int rowLen, int colLen, int valLength, T *value, int *rowPtr, int *colInd ){
    this->rowLen = rowLen;
    this->colLen = colLen;
    this->numel = colLen*rowLen;

    this->value = new T[valLength];
    std::copy(value, value+valLength, this->value);
    this->rowPtr = new int[rowLen+1];
    std::copy(rowPtr, rowPtr+rowLen+1, this->rowPtr);
    this->colInd = new int[valLength];
    std::copy(colInd, colInd+valLength, this->colInd);
}

template <class T>
SMatrix<T>::SMatrix(int rowLen, int colLen, int valLength, T *value, int *rowPtr, int *colInd){
    this->initialize(rowLen, colLen, valLength, value, rowPtr, colInd);
}

template <class T>
SMatrix<T>::SMatrix(const BMatrix<T> &inMat){
    this->initialize_from_mat(inMat, 0);
}


template <class T>
SMatrix<T>::SMatrix(const BMatrix<T> &inMat, T tolerance){
    this->initialize_from_mat(inMat, tolerance);
}

template <class T>
void SMatrix<T>::initialize_from_mat(const BMatrix<T> &inMat, T tolerance){
    int rowLen = inMat.getRowLen();
    int colLen = inMat.getColLen();

    std::vector<T> tempValue;
    std::vector<int> tempRowPtr;
    std::vector<int> tempColInd;
    int rowIndex = 0;
    T element = this->getNullElement();

    tempRowPtr.push_back(rowIndex);
    for (int row = 0; row < rowLen; row++){
        for (int col = 0; col < colLen; col++){
            inMat.retrieveElement(row, col, &element);
            if (! isZero(element, tolerance)) {
                tempValue.push_back(element);
                tempColInd.push_back(col);
                rowIndex++;
            }
        }
        tempRowPtr.push_back(rowIndex);
    }

    this->initialize(rowLen, colLen, rowIndex, &tempValue[0], &tempRowPtr[0], &tempColInd[0]);

}

template <class T>
bool SMatrix<T>::isZero(T element, T tolerance) const{
    return element == 0;
}

template<>
bool SMatrix<double>::isZero(double element, double tolerance) const {
    return element < tolerance && element > -tolerance;
}


template <class T>
int SMatrix<T>::retrieveElement(int row, int col, T *element) const {
    int colOffset = rowPtr[row]; // skip to this column
    int rowSize = rowPtr[row+1] - rowPtr[row];

    *element = 0;

    for (int i = 0; i < rowSize; i++) {
        if (colInd[colOffset + i] == col) {
            *element = this->value[colOffset + i];
            break;
        }
    }
    return M_SUCCESS;
}

template <class T>
int SMatrix<T>::getRow(int row, T *res) const {
    if (row >= this->rowLen)
        return M_DIM_ERROR;
    int colOffset = rowPtr[row];
    int rowSize = rowPtr[row+1] - rowPtr[row];
    int col = 0;
    int colMatch = 0;
    for (; col < this->colLen; col++){
        if (col < rowSize && colInd[colOffset + colMatch] == col){
            res[col] = this->value[colOffset + colMatch];
            colMatch++;
        }else {
            res[col] = 0;
        }
    }
    return M_SUCCESS;
}


template <class T>
int SMatrix<T>::getCol(int row, T *res) const {
    return 0;
}

template <class T>
int SMatrix<T>::productAx(T *x, T *res) const {
    int rows_x = sizeof(x)/sizeof(x[0]);
    if (rows_x != this->colLen)
        return M_DIM_ERROR;
    int rows_res = sizeof(res)/sizeof(res[0]);
    if (rows_x != this->colLen) 
        return M_DIM_ERROR;

    for (int row = 0; row < this->rowLen; row++) {
        int colOffset = rowPtr[row]; 
        int rowSize = rowPtr[row+1] - rowPtr[row];
        T row_res = 0;
        for (int col = 0; col < rowSize; col++) {
            int index = colOffset + col;
            row_res += this->value[index]*colInd[index];
        }
        res[row] = row_res;
    }
    return M_SUCCESS;
}

