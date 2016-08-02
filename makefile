.PHONY: exec clean

OBJ = main.o IntFunction.o search.o Interpolate.o

exec : $(OBJ)
	g++ -g  -o exec $^

$(OBJ) : nr.h nrutil.h Intfunction.h search.h Interpolate.h

clean :
	-rm -f $(OBJ)
