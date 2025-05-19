#pragma once
#include <string>
#include <stdint.h>

class Matrix{

public:

    ///////////////////////////////////////
    // Constructors
    ///////////////////////////////////////

    Matrix(short rows, short cols, std::string name);

    Matrix(Matrix &A, short cols, std::string name);

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
    bool RREF(Matrix &A, short endParam);

public:
    std::string name_ = "";//TODO delete after all works
    uint64_t *elem_;
    short rows_;
    short cols_;
};