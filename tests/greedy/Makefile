ROOT := ..
SRC := $(ROOT)/greedy/src

smoke: test_descent.out
	./test_descent.out

%.out: %.cpp
	g++ -std=c++11 $< $(SRC)/*.cpp -I $(SRC) -o $@
