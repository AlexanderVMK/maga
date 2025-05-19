all: matrix_all matrix GNU GJN
	g++ -std=c++17 -o main GNU.o matrix.o matrix_operations.o GJN.o

GNU: matrix_all GNU.cpp
	g++ -std=c++17 -c GNU.cpp

GJN: matrix_all GJN.cpp
	g++ -std=c++17 -c GJN.cpp

matrix_all: matrix matrix_operations.cpp
	g++ -std=c++17 -c matrix_operations.cpp

matrix: matrix.cpp
	g++ -std=c++17 -c matrix.cpp

clear:
	rm *.o
