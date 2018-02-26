#include <vector>

template <class T>
void SMatrix<T>::initialize(int rowLen, int colLen, int valLength, T *value, int *rowPtr, int *colInd ){
    this->rowLen = rowLen;
    this->colLen = colLen;
    this->numel = colLen*rowLen;
    this->storageSize = valLength;

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
SMatrix<T>::SMatrix(const SMatrix<T> &inMat){
    this->initialize_from_mat(inMat, 0);
}


template <class T>
SMatrix<T>::SMatrix(const SMatrix<T> &inMat, T tolerance){
    this->initialize_from_mat(inMat, tolerance);
}
template <class T>
SMatrix<T>::SMatrix(const FMatrix<T> &inMat){
    this->initialize_from_mat(inMat, 0);
}


template <class T>
SMatrix<T>::SMatrix(const FMatrix<T> &inMat, T tolerance){
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
bool SMatrix<T>::isZero(T element) const{
    return this->isZero(element, 0);
}

template <class T>
bool SMatrix<T>::isZero(T element, T tolerance) const{
    return element == 0;
}

template<>
bool SMatrix<double>::isZero(double element, double tolerance) const {
    return element < tolerance && element > -tolerance;
}
template<>
bool SMatrix<double>::isZero(double element) const {
    return this->isZero(element, 1e-15);
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
int SMatrix<T>::productAx(const T *x, T *res) const {

    for (int row = 0; row < this->rowLen; row++) {
        int colOffset = rowPtr[row]; 
        int rowSize = rowPtr[row+1] - rowPtr[row];
        T row_res = 0;
        for (int colCount = 0; colCount < rowSize; colCount++) { // colCount is for non-zero items
            int col = colInd[colOffset + colCount];
            row_res += this->value[colOffset + colCount]*x[col];
        }
        res[row] = row_res;
    }
    return M_SUCCESS;
}

template<class T>
void copy_n(T *start, int n, T *destination){
    std::copy(start, start+n, destination);
}

template <class T>
int SMatrix<T>::rowPermute(int i, int j){
    if (i == j)
        return M_SUCCESS;
    if (i < 0 || j  < 0 || i >= this->rowLen || j >= this->rowLen)
        return M_DIM_ERROR;

    // get lower and higher rows
    int lowerRow = i < j ? i: j;
    int higherRow = i < j ? j: i;
    // get properties of both rows
    int colOffset = this->rowPtr[lowerRow];
    int lowLen = this->rowPtr[lowerRow+1] - this->rowPtr[lowerRow]; 
    int highLen = this->rowPtr[higherRow+1] - this->rowPtr[higherRow]; 
    int elementsBtw = this->rowPtr[higherRow+1] - this->rowPtr[lowerRow];
    // copy value
    T * tempValue = new T[elementsBtw];
    copy_n(this->value+colOffset, elementsBtw, tempValue);
    copy_n(tempValue+elementsBtw-highLen, highLen, value+colOffset);
    copy_n(tempValue+lowLen, elementsBtw-highLen, value+colOffset+highLen);
    copy_n(tempValue, lowLen, value+colOffset+elementsBtw-lowLen);
    // copy colInd
    int * tempColInd = new int[elementsBtw];
    copy_n(this->colInd+colOffset, elementsBtw, tempColInd);
    copy_n(tempColInd+elementsBtw-highLen, highLen, colInd+colOffset);
    copy_n(tempColInd+lowLen, elementsBtw-highLen, colInd+colOffset+highLen);
    copy_n(tempColInd, lowLen, colInd+colOffset+elementsBtw-lowLen);
    // adjust rowPtr
    int lenDiff = highLen - lowLen;
    for (int row = lowerRow+1; row <= higherRow; row++) {
        rowPtr[row] = rowPtr[row] + lenDiff;
    }

    delete [] tempValue;
    delete [] tempColInd;

    return M_SUCCESS;
}

template <class T>
int SMatrix<T>::rowScale(int i, int j, T a){

    T *rowI = new T[this->rowLen];
    T *rowJ = new T[this->rowLen];
    int *newCol = new int[this->rowLen]; // store col index of new row

    int colOffset = this->rowPtr[j];
    int lenBefore = rowPtr[j+1] - rowPtr[j];
    int lenAfter = 0;

    this->getRow(i, rowI);
    this->getRow(j, rowJ);
    for (int col = 0; col < this->colLen; col++){
       T res = a*rowI[col] + rowJ[col];
        if (! isZero(rowJ[col] )){
            rowJ[lenAfter] = res;
            newCol[lenAfter] = col;
            lenAfter++;
        }
    }
    if (lenBefore == lenAfter){
        copy_n(rowJ, lenAfter, this->value+colOffset);
        copy_n(newCol, lenAfter, this->colInd+colOffset);
    }else{

        int elemAfterJ = rowPtr[this->rowLen] - rowPtr[j+1];
        int lenDiff = lenAfter - lenBefore;
        // replace value array
        T *newValue = new T[this->storageSize+lenDiff];
        copy_n(this->value, colOffset, newValue);
        copy_n(rowJ, lenAfter, newValue+colOffset);
        copy_n(this->value+colOffset+lenBefore, elemAfterJ, newValue+colOffset+lenAfter);
        // replace colInd
        int *newColInd = new int[this->storageSize+lenDiff];
        copy_n(colInd, colOffset, newColInd);
        copy_n(newCol, lenAfter, newColInd+colOffset);
        copy_n(colInd+colOffset+lenBefore, elemAfterJ, newColInd+colOffset+lenAfter);
        // adjust rowPtr
        for (int row = j+1; row < this->rowLen; row++){
            rowPtr[row] += lenDiff;
        }

        delete[] this->value;
        delete[] this->colInd;
        this->storageSize += lenDiff;
        this->value = newValue;
        this->colInd = newColInd;

    }


    delete [] rowI;
    delete [] rowJ;
    delete [] newCol;

    return M_SUCCESS;

}
