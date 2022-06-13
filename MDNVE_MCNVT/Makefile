CC = g++
CFLAGS = -Wall -O3 --std=c++11

NVE_NVT.exe : MD_MC.o random.o
	$(CC) random.o MD_MC.o -o NVE_NVT.exe

temp_test.exe : temp_test.o random.o
	$(CC) random.o temp_test.o -o $@
temp_test.o: temp_test.cpp
	$(CC) -c temp_test.cpp -o  temp_test.o  $(CFLAGS)

MD_MC.o : MD_MC.cpp
	$(CC) -c MD_MC.cpp -o MD_MC.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o NVE_NVT.exe

move_4_solid:
	mv *.dat ../excercise_04/solido
move_4_liquid:
	mv *.dat ../excercise_04/liquido
move_4_gas:
	mv *.dat ../excercise_04/gas

move_7_NVT_solid:
	mv *.dat ../excercise_07/NVT/solid
move_7_NVT_liquid:
	mv *.dat ../excercise_07/NVT/liquid
move_7_NVT_gas:
	mv *.dat ../excercise_07/NVT/gas

move_7_NVE_solid:
	mv *.dat ../excercise_07/NVE/solid
move_7_NVE_liquid:
	mv *.dat ../excercise_07/NVE/liquid
move_7_NVE_gas:
	mv *.dat ../excercise_07/NVE/gas
