#include <iostream>
#include <stdexcept>
#include <vector>

using namespace std;

class Matrix {
public:
    Matrix(int rows, int cols) : rows_num(rows), cols_num(cols), data_(rows,vector<double>(cols,0.0)) {}

    Matrix(const Matrix &other) : rows_num(other.rows_num), cols_num(other.cols_num), data_(other.rows_num,vector<double>(other.cols_num,0.0)) {
        for (int i = 0; i < rows_num; ++i) {
            for (int j = 0; j < cols_num; ++j) {
                (*this).data_[i][j] = other.data_[i][j];
            }
        }
    }

    int getRows() const { return rows_num; }
    int getCols() const { return cols_num; }

    
    double &operator()(int i, int j) { 
        if (i<=0 ||i>rows_num ||j<=0 ||j>cols_num){
            throw out_of_range("Index out of range.");
        }
        return data_[i-1][j-1]; }
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
                cout << (*this)(i, j) << " ";
            }
            cout << endl;
        }
        cout << endl;
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
        bool flag=true;
            if (!((flag)&&(rows_num==cols_num))){
                throw out_of_range("Index out of range.(soloveSLE)");
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

    Matrix solveSLE(const Matrix b) const{
        if (!(b.getRows()==cols_num)){
            throw out_of_range("Index out of range.(soloveSLE)");
        }
        Matrix y(cols_num,1);
        Matrix x(cols_num,1);
        Matrix lu=(*this).lU2();
        lu.show();
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

    Matrix transpose() const {
        Matrix result(cols_num,rows_num);
        for (int i =1;i<rows_num+1;i++){
            for (int j =1;j<cols_num+1;j++){
                result(i,j)=(*this)(j,i);
            }
        }
        return result;
    }

/*
    Matrix inverse1(){
        Matrix A(rows_num,cols_num);
        Matrix[] B=(*this).LU1();
        A=
        return A;
    }
*/

private:
    int rows_num;
    int cols_num;
    vector<vector<double>> data_;
};

int main() {
    Matrix A(4, 4);
    A(1, 1) = 8.0;
    A(1, 2) = 16.0;
    A(1, 3) = 24.0;
    A(1, 4) = 32.0;
    A(2, 1) = 2.0;
    A(2, 2) = 7.0;
    A(2, 3) = 12.0;
    A(2, 4) = 17.0;
    A(3, 1) = 6.0;
    A(3, 2) = 17.0;
    A(3, 3) = 32.0;
    A(3, 4) = 59.0;
    A(4, 1) = 7.0;
    A(4, 2) = 22.0;
    A(4, 3) = 46.0;
    A(4, 4) = 105.0;

    Matrix C = A +A;
    Matrix D = A -A;
    Matrix E = A * A;
    vector<Matrix> F=A.lU1();
    Matrix G = A.lU2();

    cout<<"A"<<endl;
    A.show();
    cout<<"C = A +A"<<endl;
    C.show();
    cout<<"D = A -A"<<endl;
    D.show();
    cout<<"E = A * A"<<endl;
    E.show();
    cout <<endl<<endl;

    cout<<"A=LU"<<endl;
    cout<<"L"<<endl;
    F[0].show();
    cout <<"U"<<endl;
    F[1].show();

    cout<<"LU"<<endl;
    G.show();
    Matrix H=F[0]*F[1];
    cout<<"L*U"<<endl;
    H.show();



    Matrix b(4,1);
    b.set({{160},{70},{198},{201}});
    cout<<"b"<<endl;
    b.show();
    Matrix I =A.solveSLE(b);
    cout<<"I"<<endl;
    I.show();

    
    return 0;
}
