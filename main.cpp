#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>

using namespace std;
class Matrix
{
public:
    Matrix( size_t n);
    double& operator()(size_t i, size_t j);
    double operator()(size_t i, size_t j) const;
    size_t getSize();

private:
    size_t mCols;
    std::vector<double> mData;
};

Matrix::Matrix( size_t n)
    :   mCols(n),
        mData(n * n)
{
    double a[n][n];
    double f;
    //fill with random values
    for(int i=0; i<n;  i++)
    {
        for(int j=0; j<n; j++)
        {
            f= rand();
            cout<<f <<"\n";
            mData[i * mCols + j]=f;
        }
    }
}

double& Matrix::operator()(size_t i, size_t j)
{
    return mData[i * mCols + j];
}

double Matrix::operator()(size_t i, size_t j) const
{
    return mData[i * mCols + j];
}
size_t Matrix::getSize()
{
    return mCols;
}


Matrix seq_mat_mul( Matrix a, Matrix b)
{
    Matrix c(5);
    return c;
}
int main()
{
//    initialize random seed:
    srand (time(NULL));

//    initialize a Matrix of size 5
    Matrix c(5);
    cout<<"Last Element: " <<c(4,4)<< "\n";
    cout<<"Size: "<<c.getSize();
    return 0;
}
