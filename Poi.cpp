#include "Matrix.cpp"
#include <vector>
using namespace std;

int main(){
    Matrix A(16,16) ;
    for(int i=1;i<17;i++)A(i,i)=-4;
    for(int i=1;i<16;i++)A(i+1,i)=1;
    for(int i=1;i<16;i++)A(i,i+1)=1;
    for(int i=1;i<13;i++)A(i+4,i)=1;
    for(int i=1;i<13;i++)A(i,i+4)=1;
    A(5,4)=0;
    A(4,5)=0;
    A(8,9)=0;
    A(9,8)=0;
    A(12,13)=0;
    A(13,12)=0;
    Matrix F(16,1);
    vector<vector<double>> f={
        {1},
        {1},
        {1},
        {1},
        {1},
        {1},
        {1},
        {1},
        {1},
        {1},
        {1},
        {1},
        {1},
        {1},
        {1},
        {1},
    };
    F.set(f);

    Matrix x=A.solveSLE_with_pivot(F);
    x.show();
    Matrix ans(4,4);
    for(int i=1;i<5;i++){
        for(int j=1;j<5;j++){
            ans(j,i)=x((i-1)*4+j,1);
        }
    }

    ans.show();
}
