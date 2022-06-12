MPI = /usr/bin/mpicxx
CC = g++
CFLAGS = -Wall -O3 --std=c++11 -Wreorder

main.exe : main.o random.o
	$(MPI) $^ -o $@

main.o : main.cpp Population.h
	$(MPI) -c $< -o $@ $(CFLAGS)
random.o : ../generatore/random.cpp ../generatore/random.h
	$(CC) -c ../generatore/random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe

move_circle :
	mv output*.dat circle

move_square :
	mv output*.dat square

start:
	/usr/bin/mpiexec -np 2 main.exe