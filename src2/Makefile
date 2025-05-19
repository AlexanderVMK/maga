GJN_ALG: matrix_all matrix GJN_M GJN
	g++ -std=c++17 -o GJN_alg GNU.o matrix.o matrix_operations.o GJN.o

PRA_ALG: matrix_all matrix PRA
	g++ -std=c++17 -o Prange PRA.o matrix.o matrix_operations.o

PRA: PRA.cpp
	g++ -std=c++17 -c PRA.cpp

GJN_M: GNU.cpp
	g++ -std=c++17 -c GNU.cpp

GJN: GJN.cpp
	g++ -std=c++17 -c GJN.cpp

matrix_all: matrix_operations.cpp
	g++ -std=c++17 -c matrix_operations.cpp

matrix: matrix.cpp
	g++ -std=c++17 -c matrix.cpp

clear:
	rm *.o
