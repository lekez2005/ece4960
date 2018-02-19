#define M_SUCCESS 0

template<class T> class SMatrix{
    public:
        int rowLen;
        int colLen;
        T *value;
        int *rowPtr;
        int *colInd;

        SMatrix(int rowLen, int colLen, T *value, int *rowPtr, int *colInd):
            rowLen(rowLen),
            colLen(colLen),
            value(value),
            rowPtr(rowPtr),
            colInd(colInd)
        {}

        int retrieveElement(int row, int col, T *element);
        void print();
        int productAx(T *x, T *res);
};

