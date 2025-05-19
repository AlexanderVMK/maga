#include "matrix_operations.hpp"
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <random>

Matrix out(1, 128, "");

bool GetRes(uint64_t op, int K)
{
    int sum = 0;
    for( int i = 0; i < K; ++i )
    if (op & (1ULL << ( 63 - i) ))
        sum += 1;

    return sum % 2;
}

std::vector<uint64_t> Create_vectors( short cols_count, short weight )
{
    /// set size is 2 ** sqrt( cols_count )
    auto sizeSet = 8;
    std::vector<uint64_t> set;

    // Генератор случайных чисел
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, cols_count - 1 );

    uint64_t rand_e;

    for ( short i = 0; i < sizeSet; ++i )
    { 
        std::vector<int> bits(cols_count, 0);

        for (int i = 0; i < weight; ++i)
        {
            int pos;
            do
            {
                pos = dis(gen);
            } while (bits[pos] == 1);
            bits[pos] = 1;
        }

        rand_e = 0;

        for (size_t i = 0; i < bits.size(); ++i)
        if (bits[i] == 1)
            rand_e |= (static_cast<uint64_t>(1) << ( i + 64 - cols_count ) );

        set.push_back( rand_e );
    }

    return set;
}

int main()
{

    // Initialize step
    // buffer allocation for variables

    //some parameters
    short w = 2;
    std::string
        H_file = "H.txt",
        s_file = "s.txt";

    // Matrix H
    Matrix
        H_st(   H_file, "H_st"),
        H(      H_file, "H"),
    //Submatrixes H_1 and H_2
        H_2( H.rows_, H.cols_ / 2, "H_2" ),
        H_1( H.rows_, H.cols_ - H_2.cols_, "H_1" ),
    // Matrix s
        s_st(   s_file, "s_st"),
        s(      s_file, "s");

    Duplicate_range( H_1, H, 0, H.rows_, 0, H_1.cols_ );
    Duplicate_range( H_2, H, 0, H.rows_, H_1.cols_, H_2.cols_ );

    //Information about problem
    std::cout << "\nStart solving ISD problem with weight: " << w << std::endl;
    H_st.Print();
    s_st.Print();

    unsigned int start_time =  clock();

    while( true )
    {

        auto set_e_1 = Create_vectors(H_1.cols_, 2 );
        auto set_e_2 = Create_vectors(H_2.cols_, 0 );

        std::vector<uint64_t>
            set_s_1,
            set_s_2;

        // Массив синдромов 1
        for( auto e : set_e_1 )
        {
            uint64_t syn = 0;
            for( int i = 0; i < H_1.rows_; i++ )
            {
                syn <<= 1;
                syn += GetRes( (uint64_t)e & H_1.elem_[i], H_1.cols_);
            }
            syn <<= 64 - H_1.rows_;
            set_s_1.push_back(syn);
        }

        // Массив синдромов 2
        for( auto e : set_e_2 )
        {
            uint64_t syn = 0;
            for( int i = 0; i < H_2.rows_; i++ )
            {
                syn <<= 1;
                syn += GetRes( (uint64_t)e & H_2.elem_[i], H_2.cols_);
            }
            syn <<= 64 - H_2.rows_;
            set_s_2.push_back(syn);
        }

        // Синдром s в uint64_t
        uint64_t syndr = 0;
        for( int i = 0; i < s.rows_; i++ )
        {
            syndr <<= 1;
            syndr += s.GetElem(i,0);
        }
        syndr <<= 64 - s.rows_;


        for( int i = 0; i < set_s_1.size(); i++ )
        for( int j = 0; j < set_s_2.size(); j++ )
        {
            if( syndr == ( set_s_1[i] ^ set_s_2[j] ) )
            {
                unsigned int
                    end_time = clock(),
                    search_time = end_time - start_time;

                set_e_2[j] >>= H_1.cols_;
                Matrix e_res(1, H.cols_, "e_res");
                e_res.elem_[0] = set_e_1[i] ^ set_e_2[j];

                std::cout << "Найденный вектор ошибок: " << std::endl;
                
                e_res.Print();

                std::cout << "Время работы: " << search_time;

                return 0;
            }
        }
    }
}