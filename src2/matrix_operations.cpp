#include "matrix_operations.hpp"
#include "matrix.h"
#include <iostream>
#include <random>
#include <algorithm>
#include <iterator>
#include <typeinfo>
#include <cstdlib>
#include <ctime>
#include <bitset>


template <typename T>
void swap(T *a, T *b)
{
    T tmp = *a;
    *a = *b;
    *b = tmp;
}

//Random permutation
std::vector<short> Random_permutate( const Matrix & H )
{
    auto cols_num = H.cols_;

    std::vector<short> permutation( cols_num );

    for( int i = 0; i < cols_num; i++ )
        permutation[i] = i;

    for (int i = 0; i < cols_num; i++)
    {
        std::random_device random_device;
        std::mt19937 generator(random_device());
        std::uniform_int_distribution<> distribution(0, cols_num - 1);
        int x = distribution(generator);
        swap(&permutation[i], &permutation[(x)]);
    }

    uint64_t new_value;
    for( int i = 0; i < H.rows_; i++ )
    {
        new_value = 0;
        for( int j = 0; j < H.cols_; j++ )
        {
            new_value <<= 1;
            new_value += H.GetElem(i, permutation[j] );
        }

        new_value <<= 64 - H.cols_;
        H.elem_[i] = new_value;
    }

    return permutation;
}

void Get_syndrome( Matrix& H, Matrix& e, Matrix & s_buff )
{
    for( int i = 0; i < s_buff.rows_; i++ )
    {
        s_buff.SetElem(i,0,0);
        for( int j = 0; j < H.cols_; j++ )
            s_buff.SetElem(i,0, ( s_buff.GetElem(i,0) != H.GetElem(i,j) * e.GetElem(j,0) ) );
    }
}

void Get_reverse( std::vector<short> pst, Matrix& e, Matrix& e_buf )
{
    for( int i = 0; i < e.rows_; i++ )
        e_buf.elem_[pst[i]] = e.elem_[i];
}

void Duplicate(Matrix& to, Matrix & from)
{
    if( !(to.rows_ == from.rows_ && to.cols_ == from.cols_) )
    {
        std::cout << "\n\n\nDuplicate ERROR with " << to.name_ << " and " << from.name_ << "\n\n\n\n";
        throw "Duplicate";
    }
    
    for(int i = 0; i < to.rows_; i++)
        to.elem_[i] = from.elem_[i];
}

void Duplicate_range(Matrix& to, Matrix & from, short from_start_rows, short count_rows, short from_start_cols, short count_cols )
{
    if( !(to.rows_ <= count_rows && to.cols_ <= count_cols) )
    {
        std::cout << "\n\n\nDuplicate ERROR with " << to.name_ << " and " << from.name_ << "\n\n\n\n";
        throw "Duplicate_range";
    }

    for(int i = from_start_rows; i < count_rows + from_start_rows; i++)
    {
        to.elem_[i - from_start_rows] = from.elem_[i];
        to.elem_[i - from_start_rows] <<= from_start_cols;
    }
}

//for 'vectors' only
short Ham_weight( Matrix& matrix )
{
    short weight = 0;
    
    for ( short i = 0; i < matrix.rows_; i++ )
    for ( short j = 0; j < matrix.cols_; j++ )
    {
        if ( matrix.GetElem(i,j) ) weight++;
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

    for( short i = 0; i < H.rows_; i++ )
    for( short j = 0; j < e.cols_; j++ )
    {
        bool tmp = 0;
        for( short k = 0; k < H.cols_; k++ )
        {
            tmp = ( tmp != H.GetElem(i,k) * e.GetElem(k,j) );
        }
        s_buf.SetElem(i,j, tmp);
    }

    for( int i = 0; i < boolCount; i++ )
        if( s.GetElem(i,0) != s_buf.GetElem(i,0) )
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

    for( short i = 0; i < H.rows_; i++ )
    for( short j = 0; j < e.cols_; j++ )
    {
        bool tmp = 0;
        for( short k = 0; k < H.cols_; k++ )
            tmp = ( tmp != H.GetElem(i,k) * e.GetElem(k,j) );

        s_buf.SetElem(i,j, tmp);
    }

    bool
        reply = true,
        solution = true;

    for( int i = 0; i < boolCount; i++ )
    if( s.GetElem(i,0) != s_buf.GetElem(i,0) )
        reply = false;

    if( !reply )
    {
        solution = false;
        reply = true;
        for( int i = 0; i < boolCount; i++ )
        if( false != s_buf.GetElem(i,0) )
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

std::vector<Matrix> Set_Merge( std::vector<Matrix> set_one, std::vector<Matrix> set_two, Matrix & H, Matrix & s, Matrix & s_buf, short boolCount )
{
    std::vector<Matrix> set_result;
    for( auto e_1 : set_one )
    for( auto e_2 : set_two )
    {
        Matrix e( e_1 );
        for( int i = 0; i < e.rows_; i++ )
            e.elem_[i] ^= e_2.elem_[i];

        if(Is_solution_or_zero( H, e, s, s_buf, boolCount ))
            set_result.push_back( e );

    }
    return set_result;
}