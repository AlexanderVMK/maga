#pragma once
#include "matrix.h"

//Random permutation
void Random_permutate( Matrix& matrix );

void Duplicate(Matrix& to, Matrix & from);

void Duplicate_range(Matrix& to, Matrix & from, short from_start_rows, short count_rows, short from_start_cols, short count_cols );

//param[in] matrix
//return true если матрица полного ранга, false - иначе
bool Full_rank( Matrix& matrix );

short Ham_weight( Matrix& matrix );

// He => s_buff
void Get_syndrome( Matrix& H, Matrix& e, Matrix & s_buff );


// H^(-1)e => s_buff
void Get_reverse( Matrix& H, Matrix& e, Matrix & s_buff );

bool Is_solution( Matrix & H, Matrix & e, Matrix & s, Matrix & s_buf, short boolCount );
bool Is_solution_or_zero( Matrix & H, Matrix & e, Matrix & s, Matrix & s_buf, short boolCount );

std::vector<Matrix> Set_Merge( std::vector<Matrix> set_one, std::vector<Matrix> set_two, Matrix & H, Matrix & s, Matrix & s_buf, short boolCount );