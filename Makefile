#g++ (GCC) 7.2.0
#GNU Scientific Library 2.4
CC=g++
CFLAGS=-Wall -fexceptions -g -Wunreachable-code -Weffc++ -Wmain -Wfatal-errors -Wextra -Wall -std=c++14 -pg -g -g3 -O0
LIBGSL=/usr/local/lib/libgsl.a /usr/local/lib/libgsl.so /usr/local/lib/libgslcblas.a /usr/local/lib/libgslcblas.so
SRC :=./src
FILES_CPP := $(wildcard $(SRC)/*.cpp)
FILES_O := $(patsubst %.cpp,%.o,$(wildcard $(SRC)/*.cpp))

BLACKOIL : $(FILES_O)
	$(CC) -o $(FILES_O) -pg -lm $(LIBGSL)

$(SRC)/main.o : main.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(SRC)/%.o : $(SRC)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@


