#include "matrix.h"
#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <iterator>
#include <typeinfo>
#include <cstdlib>
#include <ctime>


//Construct 1 matrix
Matrix::Matrix(CK_BYTE rows, CK_BYTE cols, std::string name)
    : name_( name )
    , rows_( rows )
    , cols_( cols )
{
    elem_ = new bool *[rows_];
    for(CK_BYTE i = 0; i < rows_; i++)
    {
        elem_[i] = new bool [cols_];
        for(CK_BYTE j = 0; j < cols; j++)
        {
            if( i == j){
                elem_[i][j] = 1;
            }else{
                elem_[i][j] = 0;
            }
        }
    }
}

//get last N=cols columns of input Matrix
Matrix::Matrix(Matrix &A, CK_BYTE cols, std::string name)
    : name_( name )
    , rows_( 1 )
    , cols_( cols )
{
    elem_ = new bool *[rows_];
    for(CK_BYTE i = 0; i < rows_; i++)
    {
        elem_[i] = new bool [cols_];
        for(CK_BYTE j = 0; j < cols; j++){
            elem_[i][j] = 0; 
        }
    }
    for(CK_BYTE i = 0; i < A.rows_; i++){
        elem_[0][cols - A.rows_ + i] = A.elem_[i][0];
    }
}

// Конструктор копирования
Matrix::Matrix(const Matrix &other) 
    : name_(other.name_)
    , rows_(other.rows_)
    , cols_(other.cols_)
{
    // Выделяем память для нового массива
    elem_ = new bool*[rows_];
    for (short i = 0; i < rows_; ++i) {
        elem_[i] = new bool[cols_];
        // Копируем данные из другого массива
        for (short j = 0; j < cols_; ++j) {
            elem_[i][j] = other.elem_[i][j];
        }
    }
}

//Matrix = cols N to K of Matrix A
Matrix::Matrix(Matrix &A, CK_BYTE colsStart, CK_BYTE colsEnd, std::string name)
    : name_( name )
    , rows_( A.rows_ )
    , cols_( colsEnd - colsStart + 1)
{
    elem_ = new bool *[rows_];
    for(CK_BYTE i = 0; i < rows_; i++)
    {
        elem_[i] = new bool [cols_];
        for(CK_BYTE j = 0; j < cols_; j++)
        {
            elem_[i][j] = A.elem_[i][ j + colsStart ];
        }
    }
}

//Reading matrix from file (first line - rows_number space cols_number)
Matrix::Matrix(std::string s, std::string name)
    :name_( name )
{
    std::ifstream in;
    in.open(s);
    short a, b;
    in >> a >> b;
    rows_ = a;
    cols_ = b;
    elem_ = new bool *[rows_];
    for(int i = 0; i < rows_; i++)
    {
        elem_[i] = new bool [cols_];
        for(int j = 0; j < cols_; j++){
            in >> elem_[i][j]; 
        }
    }
}

Matrix::~Matrix()
{
        for(int i = 0; i < rows_; i++)
        {
            delete[] elem_[i];
        }
        delete[] elem_;
}


///////////////////////////////////////////////
///////////////////////////////////////////////
///     METHODS
///////////////////////////////////////////////
///////////////////////////////////////////////


void Matrix::SetZero()
{
    for( short i = 0; i < rows_; i++ )
    for( short j = 0; j < cols_; j++ )
        elem_[i][j] = 0;
}

//THIS = THIS * B;
void Matrix::RM(Matrix &B, Matrix &C)
{
    CK_BYTE Brows = B.rows_, Bcols = B.cols_;

    if ( cols_ == Brows )
    for( CK_BYTE i = 0; i < rows_; i++ )
    for( CK_BYTE j = 0; j < Bcols; j++ )
    {
        bool tmp = 0;
        for( CK_BYTE k = 0; k < cols_; k++ )
        {
            tmp = ( tmp != elem_[i][k] * B.elem_[k][j] );
        }
        C.elem_[i][j] = tmp;
    }

    for(int i = 0; i < rows_; i++)
    for(int j = 0; j < cols_; j++)
        elem_[i][j] = C.elem_[i][j];
}

void Matrix::RMM(Matrix &B, Matrix &C)
{
    CK_BYTE Brows = B.rows_, Bcols = B.cols_;

    if(cols_ == Brows)
    for ( CK_BYTE i = 0; i < rows_; i++ )
    for ( CK_BYTE j = 0; j < Bcols; j++ )
    {
        bool tmp = 0;
        for ( CK_BYTE k = 0; k < cols_; k++ )
        {
            tmp = ( tmp != elem_[i][k] * B.elem_[k][j] );
        }
        C.elem_[i][j] = tmp;
    }
}

//THIS = THIS * B^T;
void Matrix::RMT(Matrix &B, Matrix &C)
{
    CK_BYTE Brows = B.rows_, Bcols = B.cols_;

    for( CK_BYTE j = 0; j < Bcols; j++ )
    {
        bool tmp = 0;
        for( CK_BYTE k = 0; k < rows_; k++ )
        {
            tmp = ( tmp != elem_[k][0] * B.elem_[j][k] );
        }
        C.elem_[j][0] = tmp;
    }

    for(int i = 0; i < rows_; i++)
    for(int j = 0; j < cols_; j++)
        elem_[i][j] = C.elem_[i][j];
}



//Reduced row echelon form
bool Matrix::RREF(Matrix &A, short endParam)
{
    for (int i = 0; i < endParam; i++)
    {
        //std::cout << " Stroka " << i << std::endl;;
        // Находим ведущий элемент в текущем столбце
        int maxRow = i;
        while (maxRow < rows_ && !elem_[maxRow][i]) {
            maxRow++;
        }
        if (maxRow == rows_) return true; // Если нет ведущего элемента, переходим к следующей строке
        
        //std::cout << " Vedushi " << maxRow << std::endl;

        // Обмен строк, чтобы поставить ведущий элемент в нужную строку
        std::swap(elem_[i], elem_[maxRow]);
        std::swap(A.elem_[i], A.elem_[maxRow]);
        
        // Применяем исключение Гаусса
        for (int j = 0; j < rows_; j++)
        {
            if (elem_[j][i] && i != j) // Если элемент не нулевой, вычитаем
            {
                for (int col = 0; col < cols_; col++)
                {
                    elem_[j][col] ^= elem_[i][col]; // XOR-операция
                    if( rows_ > col )
                        A.elem_[j][col] ^= A.elem_[i][col];
                }
            }
        }
    }

    return false;
}