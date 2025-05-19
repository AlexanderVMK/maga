#include "matrix.h"
#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <iterator>
#include <typeinfo>
#include <cstdlib>
#include <ctime>

// Функция для преобразования строки в __int128
__int128 string_to_int128(const std::string& str) {

    // Преобразование строки в __int128
    __int128 result = 0;
    for (char c : str) {
        if (c != '0' && c != '1') {
            throw std::invalid_argument("Некорректный символ в строке (допустимы только '0' и '1')");
        }
        result = (result << 1) | (c - '0'); // Сдвигаем биты и добавляем текущий бит
    }

    return result;
}

Matrix::Matrix(short rows, short cols, std::string name)
    : name_( name )
    , rows_( rows )
    , cols_( cols )
{
    elem_ = new __int128 [rows_];
    for(short i = 0; i < rows_; i++)
    {
        elem_[i] = 1;
        elem_[i] <<= ( 127 - i );
    }
}

Matrix::Matrix(const Matrix &other) 
    : name_(other.name_)
    , rows_(other.rows_)
    , cols_(other.cols_)
{
    // Выделяем память для нового массива
    elem_ = new __int128 [rows_];
    for (short i = 0; i < rows_; ++i) {
        elem_[i] = other.elem_[i];
    }
}

/// @brief Конструктор матрицы согласно матрице описанной в файле
/// @param file_name имя файла
/// @param name имя матрицы ( в финальных версиях будет упразднено, нужно для отладки )
Matrix::Matrix(std::string file_name, std::string name)
    :name_( name )
{
    std::ifstream in;
    in.open(file_name);

    short a, b;
    in >> a >> b;
    rows_ = a;
    cols_ = b;

    elem_ = new __int128 [rows_];
    std::string in_value;
    for(int i = 0; i < rows_; i++)
    {
        in >> in_value;
        elem_[i] = string_to_int128( in_value );
        elem_[i] <<= ( 128 - cols_ );
    }

    in.close();
}

Matrix::~Matrix(){ delete[] elem_; }


///////////////////////////////////////////////
///     METHODS
///////////////////////////////////////////////

/// @brief Вывод матрицы
void Matrix::Print() const
{
    std::cout << std::endl << name_ << "   rows:" << rows_ << "   cols:" << cols_ << std::endl;
    
    if ( cols_ == 1 )
    {
        for ( int i = 0; i < rows_ ; i++ ) std::cout << GetElem(i,0);
        std::cout << "\n\n";
        return;
    }

    if ( rows_ == 1 )
    {
        for ( int i = 0; i < cols_ ; i++ ) std::cout << GetElem(0,i);
        std::cout << "\n\n";
        return;
    }


    for ( int i = 0; i < rows_ ; i++ )
    {
        for ( int j = 0; j < cols_ ; j++ )
        {
            std::cout << GetElem(i,j);
        }
        std::cout << std::endl;
    }
    std::cout << "\n";
}

/// @brief Выставляет все элементы матрицы в 0
void Matrix::SetZero()
{
    for( short i = 0; i < rows_; i++ )
        elem_[i] = 0;
}

/// @brief Устанавливает значение элемента матрицы
/// @param i номер строки
/// @param j номер столбца
/// @param value значение
void Matrix::SetElem( int i, int j, bool value)
{
    if( value )
    {
        elem_[i] |= (1 << (127-j));
    }
    else
    {
        elem_[i] &= ~(1 << (127-j));
    }
}

/// @brief Инвертирует значение элемента матрицы
/// @param i номер строки
/// @param j номер столбца
void Matrix::InvertElem( int i, int j )
{
    elem_[i] ^= (1 << (127-j));
}

/// @brief Получает значение элемента матрицы
/// @param i номер строки
/// @param j номер столбца
bool Matrix::GetElem( int i, int j) const
{
    return (elem_[i] >> ( 127 - j ) ) & 1;
}

/// @brief 
/// @param B 
/// @param C 
void Matrix::RM( Matrix &B )
{
    short Brows = B.rows_, Bcols = B.cols_;
    __int128 value = 0;
    if ( cols_ == Brows )
    for( short i = 0; i < rows_; i++ )
    {
        value = 0;
        for( short j = 0; j < Bcols; j++ )
        {
            bool tmp = 0;
            for( short k = 0; k < cols_; k++ )
            {
                tmp = ( tmp != GetElem(i,k) * B.GetElem(k,j) );
            }
            value <<= 1;
            value += tmp;
        }
        value <<= 128 - cols_;
        elem_[i] = value;
    }
}

/// @brief Алгоритм Гаусса с обратным обходом
/// @param A : матрица перехода
/// @param endParam : количество итераций алгоритма ( размер единичной диагонали в матрице после выполнения алгоритма )
void Matrix::RREF(Matrix &A, short endParam)
{
    for (int i = 0; i < endParam; i++)
    {
        int maxRow = i;
        while (maxRow < rows_ && !GetElem( maxRow, i ) )
        {
            //std::cout << GetElem( maxRow, i ) << " " << i << std::endl;
            //std::cout << std::bitset<128>( elem_[i] ) << std::endl;
            maxRow++;
        }

        if (maxRow == rows_)
            return;

        // Обмен строк, чтобы поставить ведущий элемент в нужную строку
        std::swap(elem_[i], elem_[maxRow]);
        std::swap(A.elem_[i], A.elem_[maxRow]);
        
        // Применяем исключение Гаусса с обратным обходом
        for (int j = 0; j < rows_; j++)
        if ( GetElem( j, i ) && ( i != j ) )
        {
            // XOR-операция
            elem_[j] ^= elem_[i];
            A.elem_[j] ^= A.elem_[i];
        }
    }
}