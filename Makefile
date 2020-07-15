main = test

# source path
SP = src
# object files path
OP = obj
# include path
INC = include
# executable path
RP = test

CXX = icpc

FLAGS = -qopenmp -I$(INC) -O3 -ipo

LIBS = -lgsl -lgslcblas -mkl -static

all : $(RP)/test.o $(OP)/dinterpl.o
	$(CXX) $(FLAGS) -o $(RP)/$(main) $(RP)/test.o $(OP)/dinterpl.o $(LIBS)
        
# test program
$(RP)/$(main).o : $(RP)/$(main).cpp $(INC)/dinterpl.h
	$(CXX) $(FLAGS) -c -o $@ $(RP)/$(main).cpp

# dinterpl
$(OP)/dinterpl.o : $(SP)/dinterpl.cpp $(INC)/dinterpl.h
	@mkdir -p $(@D)
	$(CXX) $(FLAGS) -c -o $@ $(SP)/dinterpl.cpp

# clean all object and exec files
clean :
	rm -f $(RP)/$(main) $(OP)/*.o $(RP)/test.o
