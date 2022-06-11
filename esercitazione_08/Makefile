CC = g++
CFLAGS = -Wall -O3 --std=c++11

var_MC.exe : var_MC.o random.o
	$(CC) random.o var_MC.o -o $@

var_MC.o : var_MC.cpp
	$(CC) -c var_MC.cpp -o var_MC.o $(CFLAGS)
random.o : ../generatore/random.cpp ../generatore/random.h
	$(CC) -c ../generatore/random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o var_MC.exe