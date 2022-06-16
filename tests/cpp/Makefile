ROOT := ../../
TABU_SRC := $(ROOT)/dwave/samplers/tabu/src/
TABU_INCLUDE := $(TABU_SRC)

all: catch2 test_main tests

tests: test_main.out
	./test_main 

test_main: test_main.cpp
	g++ -std=c++11 -Wall -c test_main.cpp
	g++ -std=c++11 -Wall test_main.o $(TABU_SRC)/tabu_utils.cpp tests/*.cpp -o test_main -I $(TABU_INCLUDE)

catch2:
	git submodule init
	git submodule update
