#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "header.h"
#include "initEnergySpace.h"
#include "bingham.h"
#include "calcEnergy.h"

#include "para_lbfgs.h"

energySpace* energyDataPointer;
binghamDataSpace* binghamDataPointer;
mpi_comm_data* mpi_comm_space;
clock_t comm_time = 0;
clock_t comm_start, comm_end;

void outPut();

double para_norm(double* vector, int dim) {
	int i;
	double norm, norm_reduce;
	norm = 0;
	for(i=0;i<dim;i++){
		norm += vector[i]*vector[i];
	}
  MPI_Allreduce(&norm, &norm_reduce, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	norm = sqrt(norm_reduce);
	return norm;
}

static lbfgsfloatval_t evaluate(
    void *instance,
    const double* x,
    double* grad,
    const int n,
    const lbfgsfloatval_t step
    )
{
	int i;
	double energy;

	energy = eval_energy();

	for(i=0;i<5*PointTotalNum;i++){
		grad[i]=energyDataPointer->dEnergy_dB[i];
	}

  return energy;
}

static int progress(
    void *instance,
    const double* freeVars,
    const double* dEnergy_dFreeVars,
    const double Energy,
    const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
    )
{
  int i, rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(k%10 == 0 && rank==0){
		printf("Iteration %d:  ",k);
		printf("Energy = %-16.15f  ",Energy);
		printf("normX = %-10.8f  normdF = %-10.9f  step = %-10.8f\n", xnorm, gnorm, step);
	}
  if (k%500 == 0){
		//outPut();
	}
  return 0;
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

	int rank, stop;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	mpi_comm_space = mpi_init_space();

  energyDataPointer=initEnergySpace();
	binghamDataPointer=initiateBingham();

  setInitValue();

	int ret = 0;
	double Energy;
  lbfgs_parameter_t param;

  lbfgs_parameter_init(&param);
	param.m = 10;
	param.delta = 0;
	param.epsilon = 1e-5;
	param.max_iterations = 10000;
	param.max_linesearch = 50;
	param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
	param.ftol = 1e-4;
	param.wolfe = 0.9;

	clock_t start, finish;
	start = clock();

  ret = lbfgs(5*PointTotalNum, energyDataPointer->Bdiag_Angle_Elements, &Energy, evaluate, progress, NULL, &param);

	finish = clock();
  double norm = para_norm(energyDataPointer->Bdiag_Angle_Elements,5*PointTotalNum);

	FILE* fp;
	char filename[100];
	char number[10];

  if(rank==0){
		strcpy(filename, "LC_Timing_");
		sprintf(number, "%d_" , PointNumInRadius);
		strcat(filename, number);
		sprintf(number, "%d_" , PointNumInTheta);
		strcat(filename, number);
		sprintf(number, "%d__" , PointNumInPhi);
		strcat(filename, number);
		sprintf(number, "%d_" , RadiusSlice);
		strcat(filename, number);
		sprintf(number, "%d_" , ThetaSlice);
		strcat(filename, number);
		sprintf(number, "%d.txt" , PhiSlice);
		strcat(filename, number);

		if((fp=fopen(filename,"wt"))==NULL){
			printf("Cannot open file strike any key exit!");
			exit(1);
		}

  	fprintf(fp, "\nL-BFGS optimization terminated with status code = %d\n", ret);
    fprintf(fp, "Energy = %-20.19f normX = %-20.19f\n", Energy, norm);

		fprintf(fp, "total time: %.6f seconds\n", (double)(finish - start) / CLOCKS_PER_SEC);
		fprintf(fp, "communication time: %.6f seconds. %.2f%% \n", (double)(comm_time)/CLOCKS_PER_SEC, (double)comm_time / (finish - start)*100);

		fclose(fp);

		//outPut();
  }

	mpi_free_space(mpi_comm_space);

	binghamDataFree(binghamDataPointer);
	energySpaceFree(energyDataPointer);

  MPI_Finalize();

	return 0;
}


void outPut(){
	int i, j;
	double* p;
	FILE* fp;

	/*p=energyDataPointer->freeVars;
	if((fp=fopen("result/freeVars.txt","wt"))==NULL){
		printf("Cannot open file strike any key exit!");
		exit(1);
	}
	for(i=0; i<FreeVarsTotalNum; i++){
		fprintf(fp, "%.16f ", *p);
		p++;
	}
	fclose(fp); */

	if((fp=fopen("result/energy.txt","wt"))==NULL){
		printf("Cannot open file strike any key exit!");
		exit(1);
	}
	fprintf(fp, "%.16f ", energyDataPointer->energy);
	fclose(fp);

	p=energyDataPointer->Bdiag_Angle_Elements;
	if((fp=fopen("result/B.txt","wt"))==NULL){
		printf("Cannot open file strike any key exit!");
		exit(1);
	}
	for(i=0; i<5*PointTotalNum; i++){
		fprintf(fp, "%.16f ", *p);
		if((i+1)%PointTotalNum==0) fprintf(fp, "\n");
		p++;
	}
	fclose(fp);

	p=energyDataPointer->Qdiag_elements;
	if((fp=fopen("result/Q_eigVal.txt","wt"))==NULL){
		printf("Cannot open file strike any key exit!");
		exit(1);
	}
	for(i=0; i<2*PointTotalNum; i++){
		fprintf(fp, "%.16f ", *p);
		if((i+1)%PointTotalNum==0) fprintf(fp, "\n");
		p++;
	}
	fclose(fp);

	p=energyDataPointer->Q_elements;
	if((fp=fopen("result/Q.txt","wt"))==NULL){
		printf("Cannot open file strike any key exit!");
		exit(1);
	}
	for(i=0; i<5*PointTotalNum; i++){
		fprintf(fp, "%.16f ", *p);
		if((i+1)%PointTotalNum==0) fprintf(fp, "\n");
		p++;
	}
	fclose(fp);

	p=energyDataPointer->Z_elements;
	if((fp=fopen("result/Z.txt","wt"))==NULL){
		printf("Cannot open file strike any key exit!");
		exit(1);
	}
	for(i=0; i<PointTotalNum; i++){
		fprintf(fp, "%.16f ", *p);
		if((i+1)%PointTotalNum==0) fprintf(fp, "\n");
		p++;
	}
	fclose(fp);

	p=energyDataPointer->dEnergy_dB;
	if((fp=fopen("result/dEnergy_dB.txt","wt"))==NULL){
		printf("Cannot open file strike any key exit!");
		exit(1);
	}
	for(i=0; i<5*PointTotalNum; i++){
		fprintf(fp, "%.16f ", *p);
		if((i+1)%PointTotalNum==0) fprintf(fp, "\n");
		p++;
	}
	fclose(fp);

	if((fp=fopen("result/dZ_dBdiag.txt","a+"))==NULL){
		printf("Cannot open file strike any key exit!");
		exit(1);
	}
	for(i=0;i<2;i++){
		p=energyDataPointer->dZ_dBdiag[i];
		for(j=0; j<PointTotalNum; j++){
			fprintf(fp, "%.16f ", *p);
			p++;
		}
	}
	fclose(fp);

	if((fp=fopen("result/dQdiag_dBdiag.txt","a+"))==NULL){
		printf("Cannot open file strike any key exit!");
		exit(1);
	}
	for(i=0;i<4;i++){
		p=energyDataPointer->dQdiag_dBdiag[i];
		for(j=0; j<PointTotalNum; j++){
			fprintf(fp, "%.16f ", *p);
			p++;
		}
	}
	fclose(fp);

	if((fp=fopen("result/dQ_dBdiagAngle.txt","a+"))==NULL){
		printf("Cannot open file strike any key exit!");
		exit(1);
	}
	for(i=0;i<PointTotalNum;i++){
		p=energyDataPointer->dQ_dBdiagAngle[i];
		for(j=0; j<25; j++){
			fprintf(fp, "%.16f ", *p);
			p++;
		}
	}
	fclose(fp);

}
