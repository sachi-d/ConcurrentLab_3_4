#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
# include <omp.h>

using namespace std;
class Matrix
{
public:
    Matrix( size_t n);
    double& operator()(size_t i, size_t j);
    double operator()(size_t i, size_t j) const;
    size_t getSize();
    void generateRandomValues();
    void displayValues();


private:
    size_t mCols;
    std::vector<double> mData;
};

Matrix::Matrix( size_t n)
    :   mCols(n),
        mData(n * n)
{
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

void Matrix::generateRandomValues()
{
    double f;
    //fill with random values
    for(size_t i=0; i<mCols;  i++)
    {
        for(size_t j=0; j<mCols; j++)
        {
            f= rand() % 100;
            mData[i * mCols + j]= f ;
        }
    }
}

void Matrix::displayValues()
{
    for(size_t i=0; i<mCols; i++)
    {
        for(size_t j=0; j<mCols; j++)
        {
            cout<<mData[i * mCols + j]<<" , ";
        }
        cout<<"\n";
    }
    cout<<"\n";
}


Matrix seq_mat_mul( Matrix a, Matrix b)
{
    int n =a.getSize();
    Matrix c(n);
//    cout<< "n: " << n << "\n";
    clock_t begin_time = clock();
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            double local_sum=0;
            for(int k=0; k<n; k++)
            {
                local_sum+=(a(i,k)*b(k,j));
            }
            c(i,j)=local_sum;
        }
    }
    cout << "Sequential time: "<<float( clock () - begin_time ) / CLOCKS_PER_SEC  <<"\n";
    return c;
}

Matrix parallel_mat_mul(Matrix a, Matrix b)
{
    int n =a.getSize();
    Matrix c(n);
//    cout<< "n: " << n << "\n";
    clock_t begin_time = clock();
    # pragma omp parallel shared ( a, b, c, n  ) // private ( i, j, k )
    {
        # pragma omp for
        for ( int i = 0; i < n; i++ )
        {
            for (int j = 0; j < n; j++ )
            {
                double local_sum=0;
                for ( int k = 0; k < n; k++ )
                {
                    local_sum+= (a(i,k)*b(k,j));
                }
                c(i,j)=local_sum;
            }
        }

    }
    cout << "Parallel time: "<<float( clock () - begin_time ) / CLOCKS_PER_SEC  <<"\n";
    return c;
}

int main()
{
    //    initialize random seed:
    srand (time(NULL));

    int start=200, stop=1000, step=200;

    for(int n=start; n<=stop; n+=step)
    {
        cout<< "\nn: " << n << "\n";
        //    initialize a Matrix of size 5
        int my_size = n;

        Matrix a(my_size);
        a.generateRandomValues();
//    a.displayValues();

        Matrix b(my_size);
        b.generateRandomValues();
//    b.displayValues();

        Matrix c=seq_mat_mul(a,b);
        Matrix d=parallel_mat_mul(a,b);
//    c.displayValues();
    }

    return 0;
}
Matrix parallel_mat_mul(Matrix a, Matrix b, int N){
	assert((N%til))
}