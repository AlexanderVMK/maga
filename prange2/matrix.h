#pragma once
#include <string>
#include <stdint.h>

class Matrix{

public:

    ///////////////////////////////////////
    // Constructors
    ///////////////////////////////////////

    /// @brief Конструктор мастрицы, битовое представление которой имеет 1 на главной диагонали и 0 на остальных позициях
    /// @param rows: кол-во строк матрицы 
    /// @param cols: кол-во столбцов матрицы
    /// @param name: имя матрицы, в финальных версиях будет удалено, нужно для отладки
    Matrix(short rows, short cols, std::string name);

    Matrix(Matrix &A, short cols, std::string name);

    /// @brief Конструктор копирования матриц
    Matrix(const Matrix &other);

    Matrix(std::string s, std::string name);

    ~Matrix();

    ///////////////////////////////////////
    // Methods
    ///////////////////////////////////////

    void Print() const;
    
    void SetZero();

    void SetElem( int i, int j, bool value);

    void InvertElem( int i, int j );

    bool GetElem( int i, int j ) const;

    void RightMultiplication( const Matrix &B, Matrix &C );

    void RM( Matrix &B );

    //Reduced row echelon form
    void RREF(Matrix &A, short endParam);

public:
    std::string name_ = "";//TODO delete after all works
    __int128 *elem_;
    short rows_;
    short cols_;
};