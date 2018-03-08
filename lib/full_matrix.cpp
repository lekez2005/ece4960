
template<class T>
FMatrix<T>::FMatrix(int rowLen, int numel, T *value){
    if (rowLen == 0)
        throw std::invalid_argument("rowLen cannot be zero");

    this->rowLen = rowLen;
    this->numel = numel;

    this->colLen = numel/this->rowLen;
    if (this->colLen*this->rowLen != numel) {
        throw std::invalid_argument("no elements not divisible by rowLen");
    }

    this->value = new T[numel];
    std::copy( value, value+numel, this->value);

} 



template<class T>
void FMatrix<T>::initFromMat(const BMatrix<T> &inMat){
    this->numel = inMat.getNumel();
    this->rowLen = inMat.getRowLen();
    this->colLen = inMat.getColLen();

    this->value = new T[this->numel];
    int index = 0;
    T element = this->getNullElement();

    for (int row = 0; row < this->rowLen; row++){
        for(int col = 0; col < this->colLen; col++) {
            inMat.retrieveElement(row, col, &element);
            this->value[index] = element;
            index++;
        }
    }
    
}

template<class T>
FMatrix<T>::FMatrix(const FMatrix<T> &inMat){
    initFromMat(inMat);
}

template<class T>
FMatrix<T>::FMatrix(const SMatrix<T> &inMat){
    initFromMat(inMat);
}


template<class T>
T FMatrix<T>::retrieveElement(int row, int col) const{
    return this->value[getIndex(row, col)];
}


template<class T>
int FMatrix<T>::retrieveElement(int row, int col, T *element) const{
    int index = row*this->colLen + col;
    if (index >= this->numel)
        return M_DIM_ERROR;
    *element = this->value[index];
    return M_SUCCESS;
}

template<class T>
int FMatrix<T>::getIndex(int row, int col) const{
    return row*this->colLen + col;
}

template<class T>
int FMatrix<T>::getRow(int row, T *res) const{
    int offset = row*this->colLen;
    for(int i = 0; i < this->colLen; i++)
        res[i] = this->value[offset+i];
    return M_SUCCESS;
}

template<class T>
int FMatrix<T>::getCol(int col, T *res) const{
    for (int row = 0; row < this->colLen; row++)
        res[row] = this->value[row*this->colLen+col];
    return M_SUCCESS;
}


template<class T>
int FMatrix<T>::rowPermute(int i, int j){
    if (i < 0 || j  < 0 || i >= this->rowLen || j >= this->rowLen)
        return M_DIM_ERROR;
    T temp;
    for (int index = 0; index < this->colLen; index++){
        temp = this->retrieveElement(i, index);
        this->value[getIndex(i, index)] = this->retrieveElement(j, index);
        this->value[getIndex(j, index)] = temp;
    }
    return M_SUCCESS;
}

template<class T>
int FMatrix<T>::rowScale(int i, int j, T a){
    if (i < 0 || j  < 0 || i >= this->rowLen || j >= this->rowLen)
        return M_DIM_ERROR;
    for (int index = 0; index < this->colLen; index++){
        this->value[getIndex(j, index)] = a*this->retrieveElement(i, index) +
            this->retrieveElement(j, index); 
    }
    return M_SUCCESS;
}

void blasProductAx(double *value, double *x, double *res, int rowLen, int
        colLen){
    //cblas_dgemv(CblasRowMajor, CblasNoTrans, rowLen, colLen, 1, value, colLen,
     //       x, 1, 0, res, 1);
}

template<class T>
int FMatrix<T>::productAx(const T *x, T *res) const{

    for (int row = 0; row < this->rowLen; row++){
        T sum = 0;
        for (int col = 0; col < this->colLen; col++) {
            sum += x[col]*this->retrieveElement(row, col);
        }
        res[row] = sum;
    }
    return M_SUCCESS;
}

