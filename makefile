.PHONY: exec clean

OBJ = main.o IntFunction.o search.o Interpolate.o

exec : $(OBJ)
	g++ -g  -o exec $^

$(OBJ) : nr.h nrutil.h Intfunction.h search.h Interpolate.hpp

clean :
	-rm -f $(OBJ)
