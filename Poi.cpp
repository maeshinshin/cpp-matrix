#include "Matrix.cpp"
#include <vector>
using namespace std;

int main(){
    int n=50;
    int nn=n*n;
    Matrix A(nn,nn) ;
    for(int i=1;i<nn+1;i++)A(i,i)=-4;
    for(int i=1;i<nn;i++) A(i+1,i)=i%n==0?0:1;
    for(int i=1;i<nn;i++)A(i,i+1)=i%n==0?0:1;
    for(int i=1;i<nn-n+1;i++)A(i+n,i)=1;
    for(int i=1;i<nn-n+1;i++)A(i,i+n)=1;
    Matrix F(nn,1);
    vector<vector<double>> f(nn,vector<double>(1,1));
    F.set(f);
    //F.show();
    Matrix x=A.solveSLE_with_pivot(F);
    //x.show();
    Matrix ans(n,n);
    for(int i=1;i<n+1;i++){
        for(int j=1;j<n+1;j++){
            ans(j,i)=-1*x((i-1)*n+j,1);
        }
    }
    ans.show();
}
