#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <tuple>
#include <vector>

using namespace std;

int sign(int a,int b){
    if((a+b)%2==0){
        return 1;
    }else{
        return -1;
    }
}

class Matrix {
public:
    Matrix(int rows, int cols) : rows_num(rows), cols_num(cols), data_(rows,vector<double>(cols,0.0)) {}

    Matrix(const Matrix &other) : rows_num(other.rows_num), cols_num(other.cols_num), data_(other.rows_num,vector<double>(other.cols_num,0.0)) {}

    int getRows() const { return rows_num; }
    int getCols() const { return cols_num; }
    
    double &operator()(int i, int j) { 
        if (i<=0 ||i>rows_num ||j<=0 ||j>cols_num){
            throw out_of_range("Index out of range.");
        }
        return data_[i-1][j-1];
    }

    const double &operator()(int i, int j) const {
        if (i<=0 ||i>rows_num ||j<=0 ||j>cols_num){
            throw out_of_range("Index out of range.");
        }
        return data_[i-1][j-1]; }

    Matrix operator+(const Matrix &other) const {
        if (!((*this).getCols()==other.getCols()&&(*this).getRows()==other.getRows())){
            throw out_of_range("Index out of range.(operator+)");
        }
        Matrix result(rows_num, cols_num);
        for (int i = 1; i < rows_num+1; ++i) {
            for (int j = 1; j < cols_num+1; ++j) {
                result(i, j) = (*this)(i, j) + other(i, j);
            }
        }
        return result;
    }

    Matrix operator-(const Matrix &other) const {
        if (!((*this).getCols()==other.getCols()&&(*this).getRows()==other.getRows())){
            throw out_of_range("Index out of range.(operator-)");
        }
        Matrix result(rows_num, cols_num);
        for (int i = 1; i < rows_num+1; ++i) {
            for (int j = 1; j < cols_num+1; ++j) {
                result(i, j) = (*this)(i, j) - other(i, j);
            }
        }
        return result;
    }

    Matrix operator*(const Matrix &other) const {
        if (!((*this).getCols()==other.getRows()&&(*this).getRows()==other.getCols())){
            throw out_of_range("Index out of range.(operator*)");
        }
        Matrix result(rows_num, other.getCols());
        for (int i = 1; i < rows_num+1; ++i) {
            for (int j = 1; j < other.getCols()+1; ++j) {
                double sum = 0.0;
                for (int k = 1; k < cols_num+1; ++k) {
                    sum += (*this)(i, k) * other(k, j);
                }
                result(i, j) = sum;
            }
        }
        return result;
    }

    Matrix operator*(const vector<int> &other) const {
        Matrix result(rows_num,1);
        for (int i = 1; i < rows_num+1; ++i) {
                double sum = 0.0;
                for (int k = 1; k < cols_num+1; ++k) {
                    sum += (*this)(i, k) * other[k];
                }
                result(i, 1) = sum;
            }
        return result;
    }
    
    void set(vector<vector<double>> a){
        if (!(a.size()==rows_num&&a[0].size()==cols_num)){
            throw out_of_range("Index out of range.(set)");
        }
        for (int i =0;i<rows_num;i++){
            for (int j = 0;j<cols_num;j++){
                (*this)(i+1,j+1)=a.at(i).at(j);
            }
        }
    }

    Matrix getRow(int i) const {
        if (!(0<=i && i<=rows_num)){
            throw out_of_range("Index out of range.(getRow)");
        }
        Matrix result(1,cols_num);
        for (int j =1;j<cols_num+1;j++){
            result(1,j)=(*this)(i,j);
        }
        return result;
    }

    Matrix getCol(int i) const {
        if (!(0<=i && i<=cols_num)){
            throw out_of_range("Index out of range.(getCol)");
        }
        Matrix result(1,rows_num);
        for (int j =1;j<rows_num+1;j++){
            result(1,j)=(*this)(i,j);
        }
        return result;
    }

    void rowshow(int i) const {
        if (!(1<=i && i<=rows_num)){
            throw out_of_range("Index out of range.(rowshow)");
        }
        for (int j = 1; j < cols_num+1; j++) {
            cout << (*this)(i, j) << " ";
        }
        cout << endl;
    }

    void colshow(int i) const {
        if (!(1<=i && i<=cols_num)){
            throw out_of_range("Index out of range.(rowshow)");
        }
        for (int j = 1; j < rows_num+1; j++) {
            cout << (*this)(j, i) << " "<<endl;
        }
        cout << endl;
    }

    void show() const {
        for (int i = 1; i < rows_num+1; i++) {
            for (int j = 1; j < cols_num+1; j++) {
                cout <<setw(3)<<(*this)(i, j) << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

    Matrix change_rows(int i,int j) const{//行を入れ替える
        Matrix result(rows_num,cols_num);
        if (j<i){
            int a=i;
            i=j;
            j=a;
        }
        for (int ii =1;ii<rows_num+1;ii++){
            for (int j2=1;j2<cols_num+1;j2++){
                if(ii==i){
                    result(ii, j2)=(*this)(j,j2);
                }else if(ii==j){
                    result(ii, j2)=(*this)(i,j2);
                }else{
                    result(ii,j2)=(*this)(ii,j2);
                }
            }
        }
        return result;
    }

    Matrix preProcessing()const{
        Matrix result(*this);
        for(int k=1;k<result.rows_num+1;k++){
            if(result(k,k)==0){
                for(int i=k+1;i<rows_num+1;i++){
                    if(result(k,i)!=0){
                        if(result(i,k)!=0){
                            result.change_rows(k,i);
                        }
                    }
                }
            }
        }
        return result;
    }
    Matrix cut(int i,int j) const{
        if (!(cols_num==rows_num&&1<=i&&i<=cols_num)&&1<=j&&j<=cols_num){
            throw out_of_range("Index out of range.(cut)");
        }
        Matrix result(rows_num-1,cols_num-1);
        for (int ii =1;ii<i;ii++){
            for(int j2=1;j2<j;j2++){
                result(ii,j2)=(*this)(ii,j2); 
            }
            for(int j2=j;j2<cols_num;j2++){
                result(ii,j2)=(*this)(ii,j2+1); 
            }
        }
        for (int ii =i;ii<rows_num;ii++){
            for(int j2=1;j2<j;j2++){
                result(ii,j2)=(*this)(ii+1,j2); 
            }
            for(int j2=j;j2<cols_num;j2++){
                result(ii,j2)=(*this)(ii+1,j2+1); 
            }
        }
        return result;
    }

    double minideterminant(int i) const{
        if (!(cols_num==rows_num&&1<=i&&i<=cols_num)){
            throw out_of_range("Index out of range.(minideterminant)");
        }
        if(i==1){
            return (*this)(1,1);
        }
        double result =0;
        for (int k =1;k<i+1;k++){
            if ((1+k)%2==0){
                result+=(*this)(k,1)*(*this).cut(k,1).minideterminant(i-1);
            }else{
                result-=(*this)(k,1)*(*this).cut(k,1).minideterminant(i-1);
            }
        }
        return result;
    }


    vector<Matrix> lU1() const{ //square matrix uij=1
        Matrix l(rows_num,cols_num);  
        Matrix u(rows_num,cols_num);  
        for (int i = 1;i<cols_num+1;i++)
            u(i,i)=1; 
        for (int k=1;k<rows_num+1;k++){
            for (int i=1;i<k+1;i++){
                l(k,i)=(*this)(k,i);
                for (int i2=1;i2<i;i2++){
                    l(k,i)-=l(k,i2)*u(i2,i);
                }
            }
            for (int j=k+1;j<rows_num+1;j++){
                u(k,j)=(*this)(k,j);
                for (int j2=1;j2<k;j2++){
                    u(k,j)-=l(k,j2)*u(j2,j);
                }
                u(k,j)/=l(k,k);
            }
        }
        return {l,u} ;
    }

    Matrix lU2() const{ //square matrix uij=1
        if (!(rows_num==cols_num)){
            throw out_of_range("Index out of range.(lU2,num)");
        }
        Matrix lu(rows_num,cols_num);  
        for (int k=1;k<rows_num+1;k++){
            for (int i=1;i<k+1;i++){
                lu(k,i)=(*this)(k,i);
                for (int i2=1;i2<i;i2++){
                    lu(k,i)-=lu(k,i2)*lu(i2,i);
                }
            }
            for (int j=k+1;j<rows_num+1;j++){
                lu(k,j)=(*this)(k,j);
                for (int j2=1;j2<k;j2++){
                    lu(k,j)-=lu(k,j2)*lu(j2,j);
                }
                lu(k,j)/=lu(k,k);
            }
        }
        return lu ;
    }

    pair<Matrix,vector<int>> lU2_with_pivot() const{ //square matrix uij=1
        if (!(rows_num==cols_num)){
            throw out_of_range("Index out of range.(lU2_wp,num)");
        }
        vector<int> p(rows_num);
        for (int i=0;i<rows_num;i++){
            p[i]=i+1;
        }
        Matrix sub(*this);

        int pivot=0;
        double mx=0;
        for (int k=1;k<rows_num+1;k++){
            pivot=k;
            mx=abs(sub(pivot,k));
            for(int i=k+1;i<rows_num+1;i++){
                if(mx<abs(sub(i,k))){
                    pivot=i;
                    mx=abs(sub(pivot,k));
                }
            }
            if (pivot!=k){
                for(int i=1;i<rows_num+1;i++){
                    double tmp=sub(k,i);
                    sub(k,i)=sub(pivot,i);
                    sub(pivot,i)=tmp;
                    tmp=p[k-1];
                    p[k-1]=p[pivot-1];
                    p[pivot-1]=tmp;
                }
            }
        }

        Matrix lu=sub.lU2();
        return {lu,p} ;
    }

    Matrix solveSLE(const Matrix b) const{
        if (!(b.getRows()==cols_num)){
            throw out_of_range("Index out of range.(soloveSLE)");
        }
        Matrix y(cols_num,1);
        Matrix x(cols_num,1);
        Matrix lu=(*this).lU2();
        for (int i=1;i<cols_num+1;i++){
            y(i,1)=b(i,1);
            for(int j =1 ;j<i;j++){
                y(i,1)-=lu(i,j)*y(j,1);
            }
            y(i,1)/=lu(i,i);
        }
        for (int i = cols_num;i>0;i--){
            x(i,1)=y(i,1);
            for (int j=i+1;j<cols_num+1;j++){
                x(i,1)-=lu(i,j)*x(j,1);
            }
        }
        return x;
    }

    Matrix solveSLE_with_pivot(const Matrix b) const{
        if (!(b.getRows()==cols_num)){
            throw out_of_range("Index out of range.(soloveSLE)");
        }
        Matrix c(cols_num,1);
        Matrix y(cols_num,1);
        Matrix x(cols_num,1);
        pair<Matrix,vector<int>> lu_p=(*this).lU2_with_pivot();
        for (int i=1;i<cols_num+1;i++){
            c(i,1)=b(lu_p.second[i-1],1);
        }
        for (int i=1;i<cols_num+1;i++){
            y(i,1)=c(i,1);
            for(int j =1 ;j<i;j++){
                y(i,1)-=lu_p.first(i,j)*y(j,1);
            }
            y(i,1)/=lu_p.first(i,i);
        }
        for (int i = cols_num;i>0;i--){
            x(i,1)=y(i,1);
            for (int j=i+1;j<cols_num+1;j++){
                x(i,1)-=lu_p.first(i,j)*x(j,1);
            }
        }
        return x;
    }

    Matrix transpose() const {
        Matrix result(cols_num,rows_num);
        for (int i =1;i<rows_num+1;i++){
            for (int j =1;j<cols_num+1;j++){
                result(i,j)=(*this)(j,i);
            }
        }
        return result;
    }

    double trace() const{
        if (!(cols_num==rows_num)){
            throw out_of_range("Index out of range.(trace)");
        }
        double result=(*this)(1,1);
        for (int i = 2;i<rows_num+1;i++){
            result+=(*this)(i,i);
        }
        return result;
    }

    double ludeterminant() const{
        if (!(cols_num==rows_num)){
            throw out_of_range("Index out of range.(ludeterminant)");
        }
        return (*this).lU2().trace();
    }

    vector<vector<double>> to_vector(){
        return (*this).data_;
    }

private:
    int rows_num;
    int cols_num;
    vector<vector<double>> data_;
};

