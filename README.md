# cpp-matrix

## about this class

This class is enable you to calculate matrix.

Private variables is below,

- rows_num
  This is the number of rows in the matrix.
- cols_num
  This is the number of rows in the matrix.
- data_
  This is std::vector<std::vector<double>. as the elements of the matrix. 
  
## How to use 

Please put files `Matrix.cpp` in the same directory,
and write the following at the top of the code.

## constructor

If you create n by m matrix,
you should use constructor `M(n,m)` .
The sample is,

``` 
Matrix M(3,4);
```

## copy constructor

If you copu a matrix,
you should use copy constructor `M2(M1)` .
The sample is,

``` 
Matrix M1(3,4);
Matrix M2(M1); // copy
```

