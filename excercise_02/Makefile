CC = g++
CFLAGS = -Wall -O3 --std=c++11


EXECUTABLES = main.exe main2.0.exe main3.0.exe

all: 
		make $(EXECUTABLES)



%.exe : %.o ../generatore/random.o
	$(CC) $^ -o $@

%.o : %.cpp walk.h
	$(CC) -c $< -o $@ $(CFLAGS)


random.o : ../generatore/random.cpp ../generatore/random.h
	$(CC) -c ../generatore/random.cpp -o  ../generatore/random.o $(CFLAGS)

clean :
	rm *.o main.exe main2.exe main3.exe 
