#include <vector>
#include <float.h>
#include <matplotlibcpp.h>
#include "Matrix.cpp"
using namespace std;
using namespace matplotlibcpp;


int main(){
    int n=49;
    int nn=n*n;
    int b=0;//基底
    double q=10;//電荷
    Matrix A(nn,nn) ;
    for(int i=1;i<nn+1;i++)A(i,i)=-4;
    for(int i=1;i<nn;i++) A(i+1,i)=i%n==0?0:1;
    for(int i=1;i<nn;i++)A(i,i+1)=i%n==0?0:1;
    for(int i=1;i<nn-n+1;i++)A(i+n,i)=1;
    for(int i=1;i<nn-n+1;i++)A(i,i+n)=1;
    Matrix F(nn-1,1);
    vector<vector<double>> f(nn-1,vector<double>(1,b));
    F.set(f);
    Matrix AA=A.cut(1201,1201);
    F(1200,1)=q;
    F(1201,1)=q;
    F(1153,1)=q;
    F(1248,1)=q;
    
    //F.show();
    Matrix x=AA.solveSLE_with_CG(F);

    //x.show();
    Matrix ans(n,n);
    for(int i=1;i<n+1;i++){
        for(int j=1;j<n+1;j++){
            cout<<i<<","<<j<<endl;
            if(i==j&&i==25){
                ans(j,i)=-q;
            }else if(i<25||(i==25&&j<25)){
                ans(i,j)=x((i-1)*n+j,1);
            }else if(i>25||(i==25&&j>25)){
                ans(i,j)=x((i-1)*n+j-1,1);
            }
        }
    }
    vector<vector<double>> ans1=ans.to_vector();
    vector<double> X(n*n),y(n*n),z(n*n);
    for(int i=0;i<n;i++){
        for(int k=0;k<n;k++){
            X[i*n+k]=(i+1);
            y[i*n+k]=(k+1);
            z[i*n+k]=ans1[i][k];
        }
    }
    
    scatter(X, y, z);
    legend();
    show();
}
