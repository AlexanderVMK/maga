#include "matrix_operations.hpp"
#include "GJN.hpp"
#include <iostream>

const uint64_t leftBit = 0x8000000000000000;

int main()
{

    // Initialize step
    // buffer allocation for variables

    //some parameters
    short
        l = 5,
        w = 2, // Hemming weight for error-vector
        p = 2;

    std::string 
        H_file = "H.txt",
        s_file = "s.txt";


    Matrix
    //Matrix H
        H_st(   H_file, "H_st"),
        H(      H_file, "H"),
    //Matrix s
        s_st(   s_file, "s_st"),
        s_u(    s_file, "s_u"),
    //Submatrixes H_1 and H_2
        H_2( l, H.cols_ - l, "H_2" ),
        H_1( H.rows_ - l, H.cols_ - l, "H_1" ),
    //Submatrixes s_1 and s_2
        s_1(H.rows_ - l, 1, "s_1"),
        s_2(l, 1, "s_2"),
        s_2_buf(l, 1, "s_2_buf");

    //Information about problem

    std::cout << "\nStart solving ISD problem with parameters:\nweight: " << w << "     p: " << p << "      l: " << l << std::endl;
    H_st.Print();
    s_st.Print();

    int kk = 0;

    unsigned int start_time =  clock();

    //Algorithm step
    while( true )
    {
        Duplicate( H, H_st );
        Duplicate( s_u, s_st );

        auto revert_permutation = Random_permutate( H );
        
        if( H.RREF(s_u, l) )
            continue;

        //Submatrixes H_1 and H_2
        Matrix
            H_2( l, H.cols_ - l, "H_2" ),
            H_1( H.rows_ - l, H.cols_ - l, "H_1" );

        Duplicate_range( s_1, s_u, l, H.rows_ - l, 0, 1 );
        Duplicate_range( s_2, s_u, 0, l, 0, 1 );
        Duplicate_range( H_1, H, l, H.rows_ - l, l, H.cols_ - l );
        Duplicate_range( H_2, H, 0, l, l, H.cols_ - l );

        // Заход в GJN
        std::vector<uint64_t> a = Sieving_GJN(H_1, s_1, p);

        if( !a.size() )
            continue;

        for( auto try_e : a )
        {

            for( int i = 0; i < s_2_buf.rows_; i++ )
                s_2_buf.SetElem(i,0, ( __builtin_popcountll( H_2.elem_[i] & try_e ) % 2 ) );

            for( int i = 0; i < s_2_buf.rows_; i++ )
                s_2_buf.SetElem(i, 0, (s_2_buf.GetElem(i,0) != s_2.GetElem(i,0) ) );
            
            if( Ham_weight( s_2_buf ) != w-p )
                continue;

            unsigned int
                end_time = clock(),
                search_time = end_time - start_time;

            Matrix
                e (H.cols_, 1, "e"),
                e_buf (H.cols_, 1, "e_buf");

            for( int i = 0; i < s_2_buf.rows_; i++ )
                e.elem_[i] = s_2_buf.elem_[i];
            
            try_e <<= s_2_buf.rows_;

            for( int i = s_2_buf.rows_; i < e.rows_; i++ )
            {
                e.elem_[i] = try_e & leftBit;
                try_e <<= 1;
            }

            Get_reverse( revert_permutation, e, e_buf);

            e_buf.Print();

            std::cout << "Время работы: " << search_time;

            return 0;
        }
    }
}