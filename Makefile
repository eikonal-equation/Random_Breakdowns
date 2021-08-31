ICXX = icpc
CXXFLAGS = -std=c++17 -march=native
BOOSTPATH = ../lib/boost_1_72_0/
EIGENPATH = ../lib/eigen-3.3.8/
INCLUDE = -I $(BOOSTPATH) -I $(EIGENPATH) -I ./include/
LIB = -L $(BOOSTPATH)
SRC = ./src/
TEST =

CPPFILES = $(SRC)main.cpp \
					 $(SRC)CRandomBreakdowns.cpp\
					 $(SRC)CFMM.cpp \
					 $(SRC)CGrid.cpp \
					 $(SRC)CPolicyEvaluator.cpp
main: $(CPPFILES)
	$(CXX) $(CXXFLAGS) -O3 -o $@ $^ $(INCLUDE) $(LIB)

run: main
	./main $(TEST)

debug: $(CPPFILES)
	$(CXX) $(CXXFLAGS) -O1 -g -o $@ $^ $(INCLUDE) $(LIB)

clean:
	rm main.exe

realclean: clean
	rm output/*
