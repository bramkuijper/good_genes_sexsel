EXE=good_genes.exe
EXE_DEBUG=good_genes.dbg.exe
CPP=good_genes.cpp individual.cpp main.cpp
HPP=good_genes.hpp individual.hpp parameters.hpp
CXX=g++
CXXFLAGS=-Wall -Wconversion -O3 -std=c++20
CXXFLAGS_DEBUG=-Wall -Wconversion -ggdb -std=c++20

all : $(EXE) $(EXE_DEBUG)

$(EXE) : $(CPP) $(HPP)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(CPP) 

$(EXE_DEBUG) : $(CPP) $(HPP)
	$(CXX) $(CXXFLAGS_DEBUG) -o $(EXE_DEBUG) $(CPP) 

clean :
	rm -rf $(EXE)
	rm -rf $(EXE_DEBUG)
