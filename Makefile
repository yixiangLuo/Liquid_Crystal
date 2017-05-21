CC=mpicc
FLAG=-lm

all: LC

LC: main.o initEnergySpace.o bingham.o calcEnergy.o para_lbfgs.o
	$(CC) $(PREFIX2) -o LC main.o initEnergySpace.o bingham.o calcEnergy.o para_lbfgs.o $(FLAG)
main.o: main.c header.h initEnergySpace.h initEnergySpace.c bingham.h bingham.c calcEnergy.h calcEnergy.c para_lbfgs.h para_lbfgs.c
	$(CC) $(PREFIX1) -c main.c
initEnergySpace.o: initEnergySpace.c initEnergySpace.h header.h gaussian_nodes.h
	$(CC) -c initEnergySpace.c
bingham.o: bingham.c bingham.h
	$(CC) -c bingham.c
calcEnergy.o: calcEnergy.c calcEnergy.h header.h bingham.h bingham.c initEnergySpace.h initEnergySpace.c
	$(CC) -c calcEnergy.c
para_lbfgs.o: para_lbfgs.c para_lbfgs.h
	$(CC) -c para_lbfgs.c

clean:
	rm -f *.o LC
