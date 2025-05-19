WAG_ALG: matrix_all matrix WAG
	g++ -std=c++17 -o Wagner WAG.o matrix.o matrix_operations.o

WAG: WAG.cpp
	g++ -std=c++17 -c WAG.cpp

matrix_all: matrix_operations.cpp
	g++ -std=c++17 -c matrix_operations.cpp

matrix: matrix.cpp
	g++ -std=c++17 -c matrix.cpp

clear:
	rm *.o
