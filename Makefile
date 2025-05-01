


FLAGS = -O3 -ffast-math -march=native

lib:
	gcc -DTIMING -c $(FLAGS) -Iinclude ./src/fluidsim.c -o lib/fluidsim.o
	gcc -DTIMING -c $(FLAGS) -Iinclude ./src/particles.c -o lib/particles.o

all: lib
	gcc -DTIMING $(FLAGS) -fopenmp -Iinclude main.c ./lib/*.o -o bin/fluidsim

libd:
	gcc -g -O0 -c -Iinclude ./src/fluidsim.c -o lib/fluidsim.o
	gcc -g -O0 -c -Iinclude ./src/particles.c -o lib/particles.o


alld: libd
	gcc -g -fopenmp -O0 -Iinclude main.c ./lib/*.o -o bin/fluidsim

clean:
	rm ./lib/*
	rm ./bin/*