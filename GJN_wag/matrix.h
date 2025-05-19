#pragma once
#include <string>

typedef unsigned char CK_BYTE;

class Matrix{

public:

    ///////////////////////////////////////
    // Constructors
    ///////////////////////////////////////

    //default constructor I{rows x cols} matrix
    Matrix(CK_BYTE rows, CK_BYTE cols, std::string name);

    //get
    Matrix(Matrix &A, CK_BYTE cols, std::string name);

    //Matrix = cols N to K of Matrix A
    Matrix(Matrix &A, CK_BYTE colsStart, CK_BYTE colsEnd, std::string name);

    Matrix(const Matrix &other);

    //Reading matrix from file (first line - rows_number space cols_number)
    Matrix(std::string s, std::string name);

    ~Matrix();

    ///////////////////////////////////////
    // Methods
    ///////////////////////////////////////

    void SetZero();

    void RM(Matrix &B, Matrix &C);

    void RMT(Matrix &B, Matrix &C);

    void RMM(Matrix &B, Matrix &C);

    //Reduced row echelon form
    bool RREF(Matrix &A, short endParam);

public:
    std::string name_;//may be not needed
    bool **elem_;
    CK_BYTE rows_;
    CK_BYTE cols_;
};