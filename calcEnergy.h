#ifndef CALCENERGY_H
#define CALCENERGY_H

#include <mpi.h>

#define Radius_Low(rank) ((rank)/(ThetaSlice*PhiSlice)*(PointNumInRadius/RadiusSlice))
#define Radius_Up(rank) (((rank)/(ThetaSlice*PhiSlice)+1)*(PointNumInRadius/RadiusSlice))
#define Theta_Low(rank) ((rank)%(ThetaSlice*PhiSlice)/PhiSlice*(PointNumInTheta/ThetaSlice))
#define Theta_Up(rank) (((rank)%(ThetaSlice*PhiSlice)/PhiSlice+1)*(PointNumInTheta/ThetaSlice))
#define Phi_Low(rank) ((rank)%PhiSlice*(PointNumInPhi/PhiSlice))
#define Phi_Up(rank) (((rank)%PhiSlice+1)*(PointNumInPhi/PhiSlice))


typedef struct mpi_comm_data{
	MPI_Comm* mpi_comm_groups;
	double* mpi_msg_send;
	double* mpi_msg_rec;
	int local_total_num;
} mpi_comm_data;

mpi_comm_data* mpi_init_space();

void mpi_free_space(mpi_comm_data* mpi_comm_space);

double* linearSpace_MPI(int size, mpi_comm_data* mpi_comm_space);

void exchange_vals(int rank, double** src, double** des, mpi_comm_data* mpi_comm_space, int same);


//void makeB();

//void makeFreeVars();

//void transVarsBack();

double eval_energy();

void preCalculate(int rank);

void calc_Q_dQdB(int i);

void calc_Fbulk(int rank);

void calc_Felas(int rank);

void calc_Fpena(int rank);

void setInitValue();

#endif
