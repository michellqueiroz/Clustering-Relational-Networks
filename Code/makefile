CC		= g++ -Wall -fPIC -fexceptions -DNDEBUG -DIL_STD -std=c++11 -O3

CFLAGS1_HARARY = -I/opt/ibm/ILOG/CPLEX_Studio1271/cplex/include/ -I/opt/ibm/ILOG/CPLEX_Studio1271/concert/include/
CFLAGS2_HARARY	= -L/opt/ibm/ILOG/CPLEX_Studio1271/cplex/lib/x86-64_linux/static_pic/ -lilocplex -lcplex -L/opt/ibm/ILOG/CPLEX_Studio1271/concert/lib/x86-64_linux/static_pic/ -lconcert -lm -lpthread


FILES_CC	= mip-clustering.cpp

BINARY		= a.out

all:	clean compile

clean:
	find . -name '*.o' -exec rm -f '{}' ';'
	rm -f $(BINARY);

compile:
	$(CC) $(CFLAGS1_HARARY)  $(FILES_CC) -o $(BINARY)  $(CFLAGS2_HARARY)

run:
	./$(BINARY)

