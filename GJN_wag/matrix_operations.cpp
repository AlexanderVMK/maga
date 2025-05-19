#include "matrix_operations.hpp"
#include "matrix.h"
#include <iostream>
#include <random>
#include <algorithm>
#include <iterator>
#include <typeinfo>
#include <cstdlib>
#include <ctime>


template <typename T>
void swap(T *a, T *b)
{
    T tmp = *a;
    *a = *b;
    *b = tmp;
}


//Random permutation
void Random_permutate( Matrix& matrix )
{
    for (int i = 0; i < matrix.rows_; i++)
    {
        std::random_device random_device;
        std::mt19937 generator(random_device());
        std::uniform_int_distribution<> distribution(0, matrix.rows_-1);
        int x = distribution(generator);
        swap(&matrix.elem_[i], &matrix.elem_[(x)]);
    }
}

void Get_syndrome( Matrix& H, Matrix& e, Matrix & s_buff )
{
    for( int i = 0; i < s_buff.rows_; i++ )
    {
        s_buff.elem_[i][0] = 0;
        for( int j = 0; j < H.cols_; j++ )
            s_buff.elem_[i][0] = ( s_buff.elem_[i][0] != H.elem_[i][j] * e.elem_[j][0] );
    }
}

void Get_reverse( Matrix& H, Matrix& e, Matrix & s_buff )
{
    for( int i = 0; i < s_buff.rows_; i++ )
    {
        s_buff.elem_[i][0] = 0;
        for( int j = 0; j < H.cols_; j++ )
            s_buff.elem_[i][0] = ( s_buff.elem_[i][0] != H.elem_[j][i] * e.elem_[j][0] );
    }
}


void Print( const Matrix& matrix )
{
    std::cout << "///////////////////////////////////\n";
    std::cout << matrix.name_ << "   rows:" << (short)matrix.rows_ << "   cols:" << (short)matrix.cols_ << std::endl;
    for ( int i = 0; i < (int)matrix.rows_ ; i++)
    {
        for ( int j = 0; j < (int)matrix.cols_; j++)
        {
            std::cout << matrix.elem_[i][j] << " "; 
        }
        std::cout << std::endl;
    }
    std::cout << "///////////////////////////////////\n";
}

void Duplicate(Matrix& to, Matrix & from)
{
    if( !(to.rows_ == from.rows_ && to.cols_ == from.cols_) )
    {
        std::cout << "\n\n\nDuplicate ERROR with " << to.name_ << " and " << from.name_ << "\n\n\n\n";
        throw "Duplicate";
    }
    
    for(int i = 0; i < to.rows_; i++)
    for(int j = 0; j < to.cols_; j++)
        to.elem_[i][j] = from.elem_[i][j];
}

void Duplicate_range(Matrix& to, Matrix & from, short from_start_rows, short count_rows, short from_start_cols, short count_cols )
{
    if( !(to.rows_ <= count_rows && to.cols_ <= count_cols) )
    {
        std::cout << "\n\n\nDuplicate ERROR with " << to.name_ << " and " << from.name_ << "\n\n\n\n";
        throw "Duplicate_range";
    }

    for(int i = from_start_rows; i < count_rows + from_start_rows; i++)
    for(int j = from_start_cols; j < count_cols + from_start_cols; j++)
        to.elem_[i - from_start_rows][j - from_start_cols] = from.elem_[i][j];
}

//for 'vectors' only
CK_BYTE Ham_weight( Matrix& matrix )
{
    CK_BYTE weight = 0;
    
    for ( CK_BYTE i = 0; i < matrix.rows_; i++ )
    for ( CK_BYTE j = 0; j < matrix.cols_; j++ )
    {
        if ( matrix.elem_[i][j] ) weight++;
    }
    
    return weight;
}

bool Is_solution( Matrix & H, Matrix & e, Matrix & s, Matrix & s_buf, short boolCount )
{
    if ( H.cols_ != e.rows_
        || s.cols_ != s_buf.cols_
        || s.rows_ != s_buf.rows_
    )
    {
        throw "Is_solution error input";
    }

    for( CK_BYTE i = 0; i < H.rows_; i++ )
    for( CK_BYTE j = 0; j < e.cols_; j++ )
    {
        bool tmp = 0;
        for( CK_BYTE k = 0; k < H.cols_; k++ )
        {
            tmp = ( tmp != H.elem_[i][k] * e.elem_[k][j] );
        }
        s_buf.elem_[i][j] = tmp;
    }

    for( int i = 0; i < boolCount; i++ )
        if( s.elem_[i][0] != s_buf.elem_[i][0] )
            return false;

    return true;
}

bool Is_solution_or_zero( Matrix & H, Matrix & e, Matrix & s, Matrix & s_buf, short boolCount )
{
    if ( H.cols_ != e.rows_
        || s.cols_ != s_buf.cols_
        || s.rows_ != s_buf.rows_
    )
    {
        std::cout << "Is_solution error input\n";
        throw "Is_solution error input";
    }

    for( CK_BYTE i = 0; i < H.rows_; i++ )
    for( CK_BYTE j = 0; j < e.cols_; j++ )
    {
        bool tmp = 0;
        for( CK_BYTE k = 0; k < H.cols_; k++ )
        {
            tmp = ( tmp != H.elem_[i][k] * e.elem_[k][j] );
        }
        s_buf.elem_[i][j] = tmp;
    }

    bool reply = true;
    bool solution = true;

    for( int i = 0; i < boolCount; i++ )
        if( s.elem_[i][0] != s_buf.elem_[i][0] )
            reply = false;

    if( !reply )
    {
        solution = false;
        reply = true;
        for( int i = 0; i < boolCount; i++ )
            if( false != s_buf.elem_[i][0] )
            reply = false;
    }

    if( reply && solution )
    {
        //std::cout << "\n solution is " << (int)solution << " boolCount " << boolCount;
        //Print( H );
        //Print( e );
        //Print( s );
        //Print( s_buf );
    }
    return reply;
}