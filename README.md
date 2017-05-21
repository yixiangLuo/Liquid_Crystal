# A Parallel Routine for Liquid Crystal Simulation

## Usage

You need to install MPI before using this routine. To use it, you should modify
```
#define PointNumInRadius 32
#define PointNumInTheta 32
#define PointNumInPhi 32
```
and
```
#define RadiusSlice 2
#define ThetaSlice 2
#define PhiSlice 2
```
in `header.h`.

The first three parameters define how many discrete points you have on each axis. These integers should be little equal than 64 or equal to 100. The last three parameters define how many MPI tasks you have on each axis. The total MPI tasks number should be their product.

Then type `make` in command line. To run the routine, type
```
mpirun -np INTEGER LC
```
where `INTEGER` is the total tasks number and must be equal to `RadiusSlice`\*`ThetaSlice`\*`PhiSlice`.

## Algorithm, Results and Timings

Please see these details in `A_Parallel_Routine_for_Liquid_Crystal_Simulation.pdf`
