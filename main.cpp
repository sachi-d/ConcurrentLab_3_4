#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using namespace std;

int mat_mul()
{
    double a[5][5];//,b[5][5],c[5][5];
    int n=10;//,i,j,k;


    /* initialize random seed: */
    srand (time(NULL));


//    cout<<f  ;
//    cin>>m>>n;
//    cout<<"Enter rows and columns of second matrix:";
//    cin>>p>>q;

    double f;
//    cout<<"\nEnter first matrix:\n";
    for(int i=0; i<n;  i++)
    {
        for(int j=0; j<n; j++)
        {
            f= (double)rand() / RAND_MAX;
            cout<<f;
            a[i][j]=f;
        }

    }


//    cout<<"\nEnter second matrix:\n";
//    for(i=0; i<n; ++i)
//        for(j=0; j<n; ++j)
////            cin>>b[i][j];
//
//            cout<<"\nThe new matrix is:\n";
//    for(i=0; i<n; ++i)
//    {
//        for(j=0; j<n; ++j)
//        {
//            c[i][j]=0;
//            for(k=0; k<n; ++k)
//                c[i][j]=c[i][j]+(a[i][k]*b[k][j]);
//            cout<<c[i][j]<<" ";
//        }
//        cout<<"\n";
//    }

    return 0;
}
int main()
{
    return mat_mul();
//    return 0;
}
