#ifndef INITENERGYSPACE_H
#define INITENERGYSPACE_H

#include "header.h"

energySpace* initEnergySpace();

void energySpaceFree(energySpace* energyDataPointer);

double* linearSpace_Energy(int size, energySpace* energyDataPointer);

void freeMatrixSpace(double** dataSpace, int row);

double** matrixSpace_Energy(int row, int col, energySpace* energyDataPointer);

#endif
