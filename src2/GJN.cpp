#include "GJN.hpp"
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <random>
#include <unordered_map>
#include <set>

uint64_t leftBit = 0x8000000000000000;

void printUint64InBits(uint64_t value) {
    for (int i = 63; i >= 0; --i) { // Проходим по всем 64 битам
        uint64_t bit = (value >> i) & 1; // Извлекаем i-й бит
        std::cout << bit; // Выводим бит
    }
    std::cout << std::endl; // Переход на новую строку после вывода всех битов
}


std::vector<uint64_t> Create_vectors( short cols_count, short weight )
{
    /// set size is 2 ** sqrt( cols_count )
    auto sizeSet = std::ceil( std::pow( 2, std::ceil( std::sqrt(cols_count) ) ) );
    std::vector<uint64_t> set;
    set.reserve( sizeSet );

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

bool CheckSolution( Matrix & H, Matrix & s, uint64_t try_e, short param, bool isLast )
{
    bool non_zero = isLast;
    bool val;

    for( int i = 0; i < param; i++ )
    {
        val = ( __builtin_popcountll( H.elem_[i] & try_e ) % 2 );

        if( val != 0 ) non_zero = true;

        if( val != __builtin_popcountll( leftBit & s.elem_[i] ) && non_zero )
            return false;
    }

    return true;
}

/// Получает из множества векторов,
/// удовлетворяющих i уравнениям системы,
/// те, что удовлетворяют i+1 уравнению.
std::vector< uint64_t > GetGoodVectors(Matrix & H, Matrix & s, std::vector<uint64_t> & old_set, short param, bool isLast )
{
    std::vector< uint64_t > good_set;
    good_set.reserve(old_set.size()/2);
    
    // TODO Неплохо бы хранить два old_set'a, чтобы проверять только 'новый' бит
    for( auto vec_e : old_set )
    if( CheckSolution(H, s, vec_e, param, isLast ) )
    {
        good_set.push_back( vec_e );
    }

    return good_set;
}

int FindFirstOneBit( uint64_t num, int x )
{
    for (int i = 63; i >= x; --i)
    {
        if ((num >> i) & 1) return i;
    }
    return -1;
}

void FindCollision( std::vector<uint64_t> & set, int w, int x, int p_prime, std::vector<uint64_t> & result )
{
    if (p_prime > 0)
    {
        // Шаг 1: Распределение чисел по коллекциям B_x, ..., B_64
        std::vector< std::vector <uint64_t> > buckets(64 - x + 1); // Индексы от x+1 до 64
        for (const uint64_t& num : set )
        {
            int first_one_bit = FindFirstOneBit(num, x);

            if (first_one_bit != -1 && first_one_bit >= x + 1)
                buckets[first_one_bit].push_back(num);
        }

        // Шаг 2: Рекурсивный вызов для каждой коллекции B_i
        for (int i = x + 1; i <= 64; ++i)
        {
            if ( !buckets[i].empty() )
            {
                FindCollision( buckets[i], w, i + 1, p_prime - 1, result );

                // Распределение векторов из B_i по оставшимся коллекциям B_i+1, ..., B_64
                for (const uint64_t & num : buckets[i])
                {
                    int next_first_one_bit = FindFirstOneBit(num, i);

                    if ( next_first_one_bit != -1 && next_first_one_bit > i )
                        buckets[next_first_one_bit].push_back(num);
                }
            }
        }
    }
    else
    {
        std::unordered_map<uint64_t, int> A; // Счетчик коллизий
        std::unordered_map<uint64_t, std::set<uint64_t>> D; // Массив множеств

        // Для каждого вектора VEC из set
        for (const uint64_t& VEC : set)
        {
            // Создаем множество Y всех комбинаций чисел с ровно w единичными битами левее x
            // Это можно сделать с помощью маски
            uint64_t
                mask = (1ULL << x) - 1,
                vec_masked = VEC & mask;

            // Генерация всех возможных комбинаций y
            for ( uint64_t y = 0; y < (1ULL << x); ++y )
            if ( __builtin_popcountll(y) == w )
            {
                // Проверяем количество единичных битов
                A[y] += 1;

                if (A[y] >= 2)
                for ( const uint64_t& a : D[y] )
                    result.push_back(VEC + a);
                D[y].insert(VEC);
            }
        }
    }
}

void Merge_vectors(Matrix & H, Matrix & s, std::vector<uint64_t> & set, std::vector<uint64_t> & new_set, short weight, short boolCounts, bool isLast )
{
    std::vector<uint64_t> result;
    FindCollision( set, weight, 0, 0, result );
    
    for( auto try_e : result )
        if( CheckSolution(H, s, try_e, boolCounts, isLast ) )
            new_set.push_back( try_e );
}

std::vector<uint64_t> Sieving_GJN(Matrix & H, Matrix & s, short weight)
{
    // Множество случайных векторов заданного веса
    auto set = Create_vectors( H.cols_, weight );

    for( short i = 0; i < s.rows_; i++ )
    {
        auto new_set = GetGoodVectors( H, s, set, i, i == s.rows_ - 1 );
        Merge_vectors( H, s, set, new_set, weight, i, i == s.rows_ - 1 );
        set = new_set;
    }

    return set;
}