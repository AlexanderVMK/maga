#pragma once
#include "matrix_operations.hpp"
#include <vector>

template <typename T>
void swap(T *a, T *b)
{
    T tmp = *a;
    *a = *b;
    *b = tmp;
}

std::vector<uint64_t> Sieving_GJN(Matrix & H, Matrix & s, short weight);
