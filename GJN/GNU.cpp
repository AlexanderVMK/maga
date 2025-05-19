#include "matrix_operations.hpp"
#include "GJN.hpp"
#include <iostream>

int main()
{

    // Initialize step
    // buffer allocation for variables

    //some parameters
    short
        l = 10,
        w = 8,
        p = 4;


    //Matrixes H
    Matrix H_st("H.txt", "H_st");
    Matrix H("H.txt", "H");
    Matrix H_buf("H.txt", "H_buf");

    //Matrix P
    Matrix P_st(H_st.cols_, H_st.cols_, "pst");

    //Matrixes H_1 and H_2
    Matrix H_2( l, H.cols_ - l, "H_2" );
    Matrix H_1( H.rows_ - l, H.cols_ - l, "H_1" );
    
    //Matrix s
    Matrix s_st("s.txt", "s_st");
    Matrix s_u("s.txt", "s_u");
    
    //Matrixes s_1 and s_2
    Matrix s_1(H.rows_ - l, 1, "s_1");
    Matrix s_2(l, 1, "s_2");
    Matrix s_2_buf(l, 1, "s_2_buf");

    //Information about problem

    std::cout << "\nStart solving ISD problem with parameters:\nweight: " << w << "     p: " << p << "      l: " << l << std::endl;
    H_st.Print();
    s_st.Print();

    unsigned int start_time =  clock();

    int ccc = 0;
    //Algorithm step
    bool io = true;
    while( io )
    {
        Duplicate( H, H_st);
        Duplicate( s_u, s_st );

        //Random_permutate(P_st) && calculate RREF
        Random_permutate( P_st );
        H.RM(P_st, H_buf);
        if( H.RREF(s_u, l) )
            continue;

        Duplicate_range( s_1, s_u, l, H.rows_ - l, 0, 1 );
        Duplicate_range( s_2, s_u, 0, l, 0, 1 );
        Duplicate_range( H_1, H, l, H.rows_ - l, l, H.cols_ - l );
        Duplicate_range( H_2, H, 0, l, l, H.cols_ - l );

        //Vendor defined algorithme, may be replaced with another one
        std::vector<Matrix> a = Sieving_GJN(H_1, s_1, p);

        if( !a.size() )
            continue;
        
        for( int i = 0; i < a.size(); i++ )
        {
            Get_syndrome( H_2, a[i], s_2_buf );

            for( int j = 0; j < s_2_buf.rows_; j++ )
            {
                s_2_buf.SetElem(j,0, (s_2_buf.GetElem(j,0) != s_2.GetElem(j,0) ) );
            }
            
            if( Ham_weight( s_2_buf ) == w-p )
            {
                Matrix e (H.cols_, 1, "e");
                Matrix e_buf (H.cols_, 1, "e_buf");

                for( int j = 0; j < s_2_buf.rows_; j++ )
                    e.SetElem(j,0, s_2_buf.GetElem(j,0) );
                
                for( int j = 0; j < a[i].rows_; j++ )
                    e.SetElem( j+s_2_buf.rows_, 0, a[i].GetElem(j,0));

                Get_reverse(P_st, e, e_buf);

                unsigned int end_time = clock();
                unsigned int search_time = end_time - start_time;

                e_buf.Print();
                std::cout << "Время работы: " << search_time;
                return 0;
            }
        }
    }

    return 0;

}