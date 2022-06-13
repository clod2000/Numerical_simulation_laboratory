CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o
	$(CC) $^ -o $@

main.o : main.cpp Position.h Population.h
	$(CC) -c $< -o $@ $(CFLAGS)
random.o : ../generatore/random.cpp ../generatore/random.h
	$(CC) -c ../generatore/random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe

move_circle :
	mv output*.dat circle

move_square :
	mv output*.dat square