testmake: test.cpp 
	g++ test.cpp -o test -fopenmp -std=c++11 -O3

