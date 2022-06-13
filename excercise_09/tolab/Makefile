CC = g++
CFLAGS = -Wall -O3 --std=c++11


EXECUTABLES = main1.1.exe main1.2.exe main1.3.exe 

all: 
		make $(EXECUTABLES)



%.exe : %.o ../generatore/random.o
	$(CC) $^ -o $@

%.o : %.cpp
	$(CC) -c $< -o $@ $(CFLAGS)


random.o : ../generatore/random.cpp ../generatore/random.h
	$(CC) -c ../generatore/random.cpp -o ../generatore/random.o $(CFLAGS)

clean :
	rm *.o main1.1.exe main1.2.exe main1.3.exe 
