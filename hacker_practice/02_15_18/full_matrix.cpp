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
FMatrix<T>::FMatrix(const BMatrix<T> &inMat){
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
int FMatrix<T>::retrieveElement(int row, int col, T *element) const{
    int index = row*this->colLen + col;
    if (index >= this->numel)
        return M_DIM_ERROR;
    *element = this->value[index];
    return M_SUCCESS;
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
    return M_SUCCESS;
}
