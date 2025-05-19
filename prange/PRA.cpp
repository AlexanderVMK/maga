#include "matrix_operations.hpp"
#include <iostream>

int main()
{

    // Initialize step
    // buffer allocation for variables

    //some parameters
    short w = 8;
    std::string
        H_file = "H.txt",
        s_file = "s.txt";

    // Matrix H
    Matrix
        H_st(   H_file, "H_st"),
        H(      H_file, "H"),
    // Matrix s
        s_st(   s_file, "s_st"),
        s(      s_file, "s");

    //Information about problem

    std::cout << "\nStart solving ISD problem with weight: " << w << std::endl;
    H_st.Print();
    s_st.Print();

    unsigned int start_time =  clock();

    while( true )
    {
        Duplicate( H, H_st);
        Duplicate( s, s_st );

        // 1. Random Permutaion
        auto revert = Random_permutate( H );

        // 2. Gaus ( RREF = reduced row echelon form )
        if( H.RREF(s, H.rows_ ) )
            continue;

        // 3. Weight check
        if( Ham_weight( s ) != w )
            continue;

        // 4. Final steps
        unsigned int
            end_time = clock(),
            search_time = end_time - start_time;

        Matrix
            e ( H.cols_, 1, "e" ),
            e_res (H.cols_, 1, "e_result");

        for( int j = 0; j < s.rows_; j++ )
            e.elem_[j] = s.elem_[j];
            
        for( int j = s.rows_; j < H.cols_; j++ )
            e.elem_[j] = 0;

        // Revert permutation
        Get_reverse( revert, e, e_res );
        
        std::cout << "Найденный вектор ошибок: " << std::endl;
        
        e_res.Print();

        std::cout << "Время работы: " << search_time;

        return 0;
    }   
}