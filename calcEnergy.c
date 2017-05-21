#include "calcEnergy.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "header.h"
#include "bingham.h"
#include "initEnergySpace.h"

extern energySpace* energyDataPointer;
extern binghamDataSpace* binghamDataPointer;
extern mpi_comm_data* mpi_comm_space;
extern clock_t comm_start, comm_end, comm_time;



mpi_comm_data* mpi_init_space(){
	mpi_comm_data* mpi_comm_space = (mpi_comm_data*) malloc(sizeof(mpi_comm_data));
	mpi_comm_space->mpi_comm_groups = (MPI_Comm*) malloc(3*sizeof(MPI_Comm));
	mpi_comm_space->local_total_num = PointTotalNum/(RadiusSlice*ThetaSlice*PhiSlice);
	mpi_comm_space->mpi_msg_send = linearSpace_MPI(5*mpi_comm_space->local_total_num, mpi_comm_space);
	mpi_comm_space->mpi_msg_rec = linearSpace_MPI(max3(RadiusSlice,ThetaSlice,PhiSlice)*5*mpi_comm_space->local_total_num, mpi_comm_space);

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int color_radius = rank%(ThetaSlice*PhiSlice);
	int color_theta = rank/(ThetaSlice*PhiSlice)*size+rank%PhiSlice;
	int color_phi = rank/(ThetaSlice*PhiSlice)*size+rank%(ThetaSlice*PhiSlice)/PhiSlice;

	MPI_Comm_split(MPI_COMM_WORLD, color_radius, rank, &mpi_comm_space->mpi_comm_groups[0]);
	MPI_Comm_split(MPI_COMM_WORLD, color_theta, rank, &mpi_comm_space->mpi_comm_groups[1]);
	MPI_Comm_split(MPI_COMM_WORLD, color_phi, rank, &mpi_comm_space->mpi_comm_groups[2]);

	return mpi_comm_space;
}

void mpi_free_space(mpi_comm_data* mpi_comm_space){
	if(!mpi_comm_space) return;
	MPI_Comm_free(&mpi_comm_space->mpi_comm_groups[0]);
	MPI_Comm_free(&mpi_comm_space->mpi_comm_groups[1]);
	MPI_Comm_free(&mpi_comm_space->mpi_comm_groups[2]);
	free((void*) mpi_comm_space->mpi_comm_groups);
	free((void*) mpi_comm_space->mpi_msg_send);
	free((void*) mpi_comm_space->mpi_msg_rec);
	free((void *)mpi_comm_space);
}

double* linearSpace_MPI(int size, mpi_comm_data* mpi_comm_space){
	double* dataSpace = (double*) malloc(size*sizeof(double));
	if(dataSpace==NULL){
		mpi_free_space(mpi_comm_space);
		printf("fail to allocate temp space: menmory not enough!\n");
		exit(1);
	}
	return dataSpace;
}

void exchange_vals(int rank, double** src, double** des, mpi_comm_data* mpi_comm_space, int same){
	int i, n, m, l, q_index, n_index, m_index, l_index, block_index;
	int radius_l = Radius_Low(rank);
	int radius_u = Radius_Up(rank);
	int theta_l = Theta_Low(rank);
	int theta_u = Theta_Up(rank);
	int phi_l = Phi_Low(rank);
	int phi_u = Phi_Up(rank);
	int local_total_num = mpi_comm_space->local_total_num;

	double* dEnergy_dQdr_msg = (double*) malloc(5*local_total_num*sizeof(double));
	double* dEnergy_dQdr_rev = (double*) malloc(RadiusSlice*5*local_total_num*sizeof(double));

	for(i=0;i<3;i++){
		if(same==0 || (same==1 && i==0))
			for(n=radius_l; n<radius_u; n++){
				n_index = (n-radius_l)*(theta_u-theta_l)*(phi_u-phi_l);
				for(m=theta_l; m<theta_u; m++){
					m_index = (m-theta_l)*(phi_u-phi_l);
					for(l=phi_l; l<phi_u; l++){
						l_index = l-phi_l;
						for(q_index=0; q_index<5; q_index++){
							mpi_comm_space->mpi_msg_send[n_index+m_index+l_index+q_index*local_total_num] = src[i][Point_index(n,m,l)+q_index*PointTotalNum];
						}
					}
				}
			}

		comm_start = clock();
		MPI_Allgather(mpi_comm_space->mpi_msg_send, 5*local_total_num, MPI_DOUBLE, mpi_comm_space->mpi_msg_rec, 5*local_total_num, MPI_DOUBLE, mpi_comm_space->mpi_comm_groups[i]);
		comm_end = clock();
		comm_time += comm_end-comm_start;

		for(n=(i==0?0:radius_l); n<(i==0?PointNumInRadius:radius_u); n++){
			n_index = n%(radius_u-radius_l)*(theta_u-theta_l)*(phi_u-phi_l);
			if(i==0) block_index = n/(radius_u-radius_l);
			for(m=(i==1?0:theta_l); m<(i==1?PointNumInTheta:theta_u); m++){
				m_index = m%(theta_u-theta_l)*(phi_u-phi_l);
				if(i==1) block_index = m/(theta_u-theta_l);
				for(l=(i==2?0:phi_l); l<(i==2?PointNumInPhi:phi_u); l++){
					l_index = l%(phi_u-phi_l);
					if(i==2) block_index = l/(phi_u-phi_l);
					for(q_index=0; q_index<5; q_index++){
						des[i][Point_index(n,m,l)+q_index*PointTotalNum] = mpi_comm_space->mpi_msg_rec[block_index*5*local_total_num+n_index+m_index+l_index+q_index*local_total_num];
					}
				}
			}
		}
	}
}

/*void makeB(){
	int i, j, k, q_index;
	int indexBasis_B, indexBasis_Var;

	//径向对称
	/ *for(i=0; i<PointNumInRadius; i++){
		for(q_index=0; q_index<2; q_index++){
			for(j=0; j<PointNumInTheta; j++){
				indexBasis_B=i*PointNumInTheta*PointNumInPhi+j*PointNumInPhi+q_index*PointTotalNum;
				for(k=0; k<PointNumInPhi; k++){
					energyDataPointer->Bdiag_Angle_Elements[indexBasis_B+k]=energyDataPointer->freeVars[i];
				}
			}
		}
	}* /

	//旋转对称
	for(i=0; i<PointNumInRadius; i++){
		for(q_index=0; q_index<4; q_index++){
			indexBasis_Var=i*PointNumInTheta+q_index*FreeVarsUnitNum;
			for(j=0; j<PointNumInTheta; j++){
				//if(fabs(energyDataPointer->pointInRadius[i]*cos(energyDataPointer->pointInTheta[j]))<SplitCorePoint) continue;	//split core 初次求解，限定中部横向，配合初值使用
				indexBasis_B=i*PointNumInTheta*PointNumInPhi+j*PointNumInPhi+q_index*PointTotalNum;
				for(k=0; k<PointNumInPhi; k++){
					energyDataPointer->Bdiag_Angle_Elements[indexBasis_B+k]=energyDataPointer->freeVars[indexBasis_Var+j];
				}
			}
		}
		//q_index=4
		indexBasis_Var=i*PointNumInTheta+4*FreeVarsUnitNum;
		for(j=0; j<PointNumInTheta; j++){
			//if(fabs(energyDataPointer->pointInRadius[i]*cos(energyDataPointer->pointInTheta[j]))<SplitCorePoint) continue;	//split core 初次求解，限定中部横向，配合初值使用
			indexBasis_B=i*PointNumInTheta*PointNumInPhi+j*PointNumInPhi+4*PointTotalNum;
			for(k=0; k<PointNumInPhi; k++){
				energyDataPointer->Bdiag_Angle_Elements[indexBasis_B+k]=energyDataPointer->freeVars[indexBasis_Var+j]-energyDataPointer->pointInPhi[k];
			}
		}
	}
}*/

/*void makeFreeVars(){
	int i, j, k, q_index;
	int indexBasis_B, indexBasis_Var;

	//径向对称
	/ *for(i=0; i<PointNumInRadius; i++){
		energyDataPointer->freeVars[i]=energyDataPointer->Bdiag_Angle_Elements[i*PointNumInTheta*PointNumInPhi];
	}* /

	//旋转对称
	for(i=0; i<PointNumInRadius; i++){
		for(q_index=0; q_index<4; q_index++){
			indexBasis_Var=i*PointNumInTheta+q_index*FreeVarsUnitNum;
			for(j=0; j<PointNumInTheta; j++){
				indexBasis_B=i*PointNumInTheta*PointNumInPhi+j*PointNumInPhi+q_index*PointTotalNum;
				energyDataPointer->freeVars[indexBasis_Var+j]=energyDataPointer->Bdiag_Angle_Elements[indexBasis_B];
			}
		}
		//q_index=4
		indexBasis_Var=i*PointNumInTheta+4*FreeVarsUnitNum;
		for(j=0; j<PointNumInTheta; j++){
			indexBasis_B=i*PointNumInTheta*PointNumInPhi+j*PointNumInPhi+4*PointTotalNum;
			energyDataPointer->freeVars[indexBasis_Var+j]=energyDataPointer->Bdiag_Angle_Elements[indexBasis_B]+energyDataPointer->pointInPhi[0];
		}
	}
}*/

/*void transVarsBack(){
	int i, j, k, q_index;
	int indexBasis_B, indexBasis_Var;

	for(i=0; i<FreeVarsTotalNum; i++){
		energyDataPointer->dEnergy_freeVars[i]=0;
	}

	//径向对称
	/ *for(i=0; i<PointNumInRadius; i++){
		for(q_index=0; q_index<2; q_index++){
			for(j=0; j<PointNumInTheta; j++){
				indexBasis_B=i*PointNumInTheta*PointNumInPhi+j*PointNumInPhi+q_index*PointTotalNum;
				for(k=0; k<PointNumInPhi; k++){
					energyDataPointer->dEnergy_freeVars[i]+=energyDataPointer->dEnergy_dB[indexBasis_B+k];
				}
			}
		}
	}* /

	//旋转对称
	for(i=0; i<PointNumInRadius; i++){
		for(q_index=0; q_index<5; q_index++){
			indexBasis_Var=i*PointNumInTheta+q_index*FreeVarsUnitNum;
			for(j=0; j<PointNumInTheta; j++){
				//if(fabs(energyDataPointer->pointInRadius[i]*cos(energyDataPointer->pointInTheta[j]))<SplitCorePoint) continue;	//split core 初次求解，限定中部横向，配合初值使用
				indexBasis_B=i*PointNumInTheta*PointNumInPhi+j*PointNumInPhi+q_index*PointTotalNum;
				for(k=0; k<PointNumInPhi; k++){
					energyDataPointer->dEnergy_freeVars[indexBasis_Var+j]+=energyDataPointer->dEnergy_dB[indexBasis_B+k];
				}
			}
		}
	}
}*/

double eval_energy(){
	int rank;
  double energy=0;
	double *Q_elem[3];
	Q_elem[0]=energyDataPointer->Q_elements;
	Q_elem[1]=energyDataPointer->Q_elements;
	Q_elem[2]=energyDataPointer->Q_elements;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	preCalculate(rank);

	exchange_vals(rank, Q_elem, Q_elem, mpi_comm_space, 1);

	calc_Fbulk(rank);
	calc_Felas(rank);
	calc_Fpena(rank);
	comm_start = clock();
  MPI_Allreduce(&energyDataPointer->energy, &energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	comm_end = clock();
	comm_time += comm_end-comm_start;

  return energy;
}

void preCalculate(int rank){
	int i;
	int n, m, l, q_index;
	int radius_l = Radius_Low(rank);
	int radius_u = Radius_Up(rank);
	int theta_l = Theta_Low(rank);
	int theta_u = Theta_Up(rank);
	int phi_l = Phi_Low(rank);
	int phi_u = Phi_Up(rank);

	Moments M;
	double bingham_B[3]={0};

	//计算Q的特征值、Z以及它们关于B的特征值的偏导数
	for(n=radius_l; n<radius_u; n++){
		for(m=theta_l; m<theta_u; m++){
			for(l=phi_l; l<phi_u; l++){
				i=Point_index(n,m,l);

				bingham_B[0]=energyDataPointer->Bdiag_Angle_Elements[i];
				bingham_B[1]=energyDataPointer->Bdiag_Angle_Elements[PointTotalNum+i];
				bingham_B[2]=-energyDataPointer->Bdiag_Angle_Elements[i]-energyDataPointer->Bdiag_Angle_Elements[PointTotalNum+i];

		    M = bingham(bingham_B, binghamDataPointer);

				energyDataPointer->Z_elements[i]=M.Z;
				energyDataPointer->dZ_dBdiag[0][i]=2*M.Z_x12+M.Z_x22-M.Z;
				energyDataPointer->dZ_dBdiag[1][i]=M.Z_x12+2*M.Z_x22-M.Z;
				energyDataPointer->Qdiag_elements[i]=M.Q_x12-1.0/3;
				energyDataPointer->Qdiag_elements[PointTotalNum+i]=M.Q_x22-1.0/3;
				energyDataPointer->dQdiag_dBdiag[0][i]=2*M.Q_x14+M.Q_x12x22-(2*M.Q_x12+M.Q_x22)*M.Q_x12;					//dQ1_dB1
				energyDataPointer->dQdiag_dBdiag[1][i]=M.Q_x14+2*M.Q_x12x22-(M.Q_x12+2*M.Q_x22)*M.Q_x12;					//dQ1_dB2
				energyDataPointer->dQdiag_dBdiag[2][i]=2*M.Q_x12x22+M.Q_x24-(2*M.Q_x12+M.Q_x22)*M.Q_x22;					//dQ2_dB1
				energyDataPointer->dQdiag_dBdiag[3][i]=M.Q_x12x22+2*M.Q_x24-(M.Q_x12+2*M.Q_x22)*M.Q_x22;					//dQ2_dB2

				calc_Q_dQdB(i);
			}
		}
	}

	energyDataPointer->energy=0;
	for(n=radius_l; n<radius_u; n++){
		for(m=theta_l; m<theta_u; m++){
			for(l=phi_l; l<phi_u; l++){
				i=Point_index(n,m,l);
				for(q_index=0; q_index<5; q_index++){
					energyDataPointer->dEnergy_dB[i+q_index*PointTotalNum]=0;
				}
			}
		}
	}
}

//通过B旋转矩阵的值和Q的特征值，计算Q的值及Q关于自变量B的偏导数的值
void calc_Q_dQdB(int i){
	double alpha=energyDataPointer->Bdiag_Angle_Elements[2*PointTotalNum+i];
	double beta=energyDataPointer->Bdiag_Angle_Elements[3*PointTotalNum+i];
	double gamma=energyDataPointer->Bdiag_Angle_Elements[4*PointTotalNum+i];
	double Q_diag[3]={energyDataPointer->Qdiag_elements[i], energyDataPointer->Qdiag_elements[PointTotalNum+i], -energyDataPointer->Qdiag_elements[i]-energyDataPointer->Qdiag_elements[PointTotalNum+i]};
	double dQdiag_dB[2][3]={{energyDataPointer->dQdiag_dBdiag[0][i], energyDataPointer->dQdiag_dBdiag[2][i], -energyDataPointer->dQdiag_dBdiag[0][i]-energyDataPointer->dQdiag_dBdiag[2][i]}, {energyDataPointer->dQdiag_dBdiag[1][i], energyDataPointer->dQdiag_dBdiag[3][i], -energyDataPointer->dQdiag_dBdiag[1][i]-energyDataPointer->dQdiag_dBdiag[3][i]}};

	double sin_alpha=sin(alpha), cos_alpha=cos(alpha);
	double sin_beta=sin(beta), cos_beta=cos(beta);
	double sin_gamma=sin(gamma), cos_gamma=cos(gamma);

	double euler[3][3]={
	{cos_alpha*cos_gamma - cos_beta*sin_alpha*sin_gamma, cos_gamma*sin_alpha + cos_alpha*cos_beta*sin_gamma, sin_beta*sin_gamma},
	{-(cos_beta*cos_gamma*sin_alpha) - cos_alpha*sin_gamma, cos_alpha*cos_beta*cos_gamma - sin_alpha*sin_gamma, cos_gamma*sin_beta},
	{sin_alpha*sin_beta, -(cos_alpha*sin_beta), cos_beta}};

	double euler_dAlpha[3][3]={
	{-(cos_gamma*sin_alpha) - cos_alpha*cos_beta*sin_gamma,cos_alpha*cos_gamma - cos_beta*sin_alpha*sin_gamma,0},
	{-(cos_alpha*cos_beta*cos_gamma) + sin_alpha*sin_gamma,-(cos_beta*cos_gamma*sin_alpha) - cos_alpha*sin_gamma,0},
	{cos_alpha*sin_beta,sin_alpha*sin_beta,0}};
	double euler_dBeta[3][3]={
	{sin_alpha*sin_beta*sin_gamma,-(cos_alpha*sin_beta*sin_gamma),cos_beta*sin_gamma},
	{cos_gamma*sin_alpha*sin_beta,-(cos_alpha*cos_gamma*sin_beta),cos_beta*cos_gamma},
	{cos_beta*sin_alpha,-(cos_alpha*cos_beta),-sin_beta}};
	double euler_dGamma[3][3]={
	{-(cos_beta*cos_gamma*sin_alpha) - cos_alpha*sin_gamma,cos_alpha*cos_beta*cos_gamma - sin_alpha*sin_gamma,cos_gamma*sin_beta},
	{-(cos_alpha*cos_gamma) + cos_beta*sin_alpha*sin_gamma,-(cos_gamma*sin_alpha) - cos_alpha*cos_beta*sin_gamma,-(sin_beta*sin_gamma)},
	{0,0,0}};

	int indexTans[2][5]={{0, 0, 0, 1, 1}, {0, 1, 2, 1, 2}};
	int k, j;
	double val, val1, val2, val3;
	for(k=0; k<5; k++){
		val1=0; val2=0; val3=0;
		for(j=0; j<3; j++){
			val=euler[indexTans[0][k]][j]*euler[indexTans[1][k]][j];
			val1+=val*Q_diag[j];
			val2+=val*dQdiag_dB[0][j];
			val3+=val*dQdiag_dB[1][j];
		}
		energyDataPointer->Q_elements[k*PointTotalNum+i]=val1;
		energyDataPointer->dQ_dBdiagAngle[i][k]=val2;
		energyDataPointer->dQ_dBdiagAngle[i][5+k]=val3;
	}

	for(k=0; k<5; k++){
		val1=0; val2=0; val3=0;
		for(j=0; j<3; j++){
			val1+=euler[indexTans[0][k]][j]*Q_diag[j]*euler_dAlpha[indexTans[1][k]][j]+euler[indexTans[1][k]][j]*Q_diag[j]*euler_dAlpha[indexTans[0][k]][j];
			val2+=euler[indexTans[0][k]][j]*Q_diag[j]*euler_dBeta[indexTans[1][k]][j]+euler[indexTans[1][k]][j]*Q_diag[j]*euler_dBeta[indexTans[0][k]][j];
			val3+=euler[indexTans[0][k]][j]*Q_diag[j]*euler_dGamma[indexTans[1][k]][j]+euler[indexTans[1][k]][j]*Q_diag[j]*euler_dGamma[indexTans[0][k]][j];
		}
		energyDataPointer->dQ_dBdiagAngle[i][2*5+k]=val1;
		energyDataPointer->dQ_dBdiagAngle[i][3*5+k]=val2;
		energyDataPointer->dQ_dBdiagAngle[i][4*5+k]=val3;
	}

}

void calc_Fbulk(int rank){
	int i;
	int n, m, l;
	int radius_l = Radius_Low(rank);
	int radius_u = Radius_Up(rank);
	int theta_l = Theta_Low(rank);
	int theta_u = Theta_Up(rank);
	int phi_l = Phi_Low(rank);
	int phi_u = Phi_Up(rank);
	double energy=0;

	//Int[Q:B - Ln(Z)] = b1*q1 + b2*q2 + (b1 + b2)*(q1 + q2) - ln(Z)
	//泛函值
	energy=0;
	for(n=radius_l; n<radius_u; n++){
		for(m=theta_l; m<theta_u; m++){
			for(l=phi_l; l<phi_u; l++){
				i=Point_index(n,m,l);
				energy += energyDataPointer->ballInteCor[i]*(energyDataPointer->Bdiag_Angle_Elements[i]*energyDataPointer->Qdiag_elements[i] + energyDataPointer->Bdiag_Angle_Elements[PointTotalNum+i]*energyDataPointer->Qdiag_elements[PointTotalNum+i] + (energyDataPointer->Bdiag_Angle_Elements[i] + energyDataPointer->Bdiag_Angle_Elements[PointTotalNum+i]) * (energyDataPointer->Qdiag_elements[i] + energyDataPointer->Qdiag_elements[PointTotalNum+i]) - log(energyDataPointer->Z_elements[i]));
			}
		}
	}
	energyDataPointer->energy += energy;
	//对自变量的导数
	for(n=radius_l; n<radius_u; n++){
		for(m=theta_l; m<theta_u; m++){
			for(l=phi_l; l<phi_u; l++){
				i=Point_index(n,m,l);
				//dEnergy_dB1
				energyDataPointer->dEnergy_dB[i] += energyDataPointer->ballInteCor[i]*( 2*energyDataPointer->Qdiag_elements[i] + energyDataPointer->Qdiag_elements[PointTotalNum+i] + (2*energyDataPointer->Bdiag_Angle_Elements[i] + energyDataPointer->Bdiag_Angle_Elements[PointTotalNum+i])*energyDataPointer->dQdiag_dBdiag[0][i] + (energyDataPointer->Bdiag_Angle_Elements[i] + 2*energyDataPointer->Bdiag_Angle_Elements[PointTotalNum+i])*energyDataPointer->dQdiag_dBdiag[2][i] - energyDataPointer->dZ_dBdiag[0][i]/energyDataPointer->Z_elements[i]);
				//dEnergy_dB2
				energyDataPointer->dEnergy_dB[PointTotalNum+i] += energyDataPointer->ballInteCor[i]*( energyDataPointer->Qdiag_elements[i] + 2*energyDataPointer->Qdiag_elements[PointTotalNum+i] + (2*energyDataPointer->Bdiag_Angle_Elements[i] + energyDataPointer->Bdiag_Angle_Elements[PointTotalNum+i])*energyDataPointer->dQdiag_dBdiag[1][i] + (energyDataPointer->Bdiag_Angle_Elements[i] + 2*energyDataPointer->Bdiag_Angle_Elements[PointTotalNum+i])*energyDataPointer->dQdiag_dBdiag[3][i] - energyDataPointer->dZ_dBdiag[1][i]/energyDataPointer->Z_elements[i]);
				//其他为0
			}
		}
	}

	//Int[|Q|^2] = q1^2 + q2^2 + (q1 + q2)^2
	//泛函值
	energy=0;
	for(n=radius_l; n<radius_u; n++){
		for(m=theta_l; m<theta_u; m++){
			for(l=phi_l; l<phi_u; l++){
				i=Point_index(n,m,l);
				energy += energyDataPointer->ballInteCor[i]*2*(energyDataPointer->Qdiag_elements[i]*energyDataPointer->Qdiag_elements[i] + energyDataPointer->Qdiag_elements[i]*energyDataPointer->Qdiag_elements[PointTotalNum+i] + energyDataPointer->Qdiag_elements[PointTotalNum+i]*energyDataPointer->Qdiag_elements[PointTotalNum+i]);
			}
		}
	}
	energyDataPointer->energy += -Alpha1*0.5 * energy;
	//对自变量的导数
	for(n=radius_l; n<radius_u; n++){
		for(m=theta_l; m<theta_u; m++){
			for(l=phi_l; l<phi_u; l++){
				i=Point_index(n,m,l);
				//dEnergy_dB1
				energyDataPointer->dEnergy_dB[i] += -Alpha1*0.5 * energyDataPointer->ballInteCor[i]*2*(energyDataPointer->Qdiag_elements[i]*(2*energyDataPointer->dQdiag_dBdiag[0][i] + energyDataPointer->dQdiag_dBdiag[2][i]) + energyDataPointer->Qdiag_elements[PointTotalNum+i]*(energyDataPointer->dQdiag_dBdiag[0][i] + 2*energyDataPointer->dQdiag_dBdiag[2][i]));
				//dEnergy_dB2
				energyDataPointer->dEnergy_dB[PointTotalNum+i] += -Alpha1*0.5 * energyDataPointer->ballInteCor[i]*2*(energyDataPointer->Qdiag_elements[i]*(2*energyDataPointer->dQdiag_dBdiag[1][i] + energyDataPointer->dQdiag_dBdiag[3][i]) + energyDataPointer->Qdiag_elements[PointTotalNum+i]*(energyDataPointer->dQdiag_dBdiag[1][i] + 2*energyDataPointer->dQdiag_dBdiag[3][i]));
				//其他为0
			}
		}
	}

}

void calc_Felas(int rank){
	int i, j, k, n, m, l, q_index;
	int indexBasis;
	int radius_l = Radius_Low(rank);
	int radius_u = Radius_Up(rank);
	int theta_l = Theta_Low(rank);
	int theta_u = Theta_Up(rank);
	int phi_l = Phi_Low(rank);
	int phi_u = Phi_Up(rank);
	double count=0, energy=0;

	for(q_index=0; q_index<5; q_index++){

		for(j=theta_l; j<theta_u; j++){
			for(k=phi_l; k<phi_u; k++){
				indexBasis=j*PointNumInPhi+k+q_index*PointTotalNum;
				for(l=radius_l; l<radius_u; l++){
					count=0;
					for(i=0; i<PointNumInRadius; i++){
						count+=energyDataPointer->gradMatrix_radius[l][i] * energyDataPointer->Q_elements[i*PointNumInTheta*PointNumInPhi+indexBasis];
					}
					energyDataPointer->dQ_drtp[0][l*PointNumInTheta*PointNumInPhi+indexBasis]=count;
				}
			}
		}

		for(i=radius_l; i<radius_u; i++){
			for(k=phi_l; k<phi_u; k++){
				indexBasis=i*PointNumInTheta*PointNumInPhi+k+q_index*PointTotalNum;
				for(l=theta_l; l<theta_u; l++){
					count=0;
					for(j=0; j<PointNumInTheta; j++){
						count+=energyDataPointer->gradMatrix_theta[l][j] * energyDataPointer->Q_elements[j*PointNumInPhi+indexBasis];
					}
					energyDataPointer->dQ_drtp[1][l*PointNumInPhi+indexBasis]=count;
				}
			}
		}

		for(i=radius_l; i<radius_u; i++){
			for(j=theta_l; j<theta_u; j++){
				indexBasis=i*PointNumInTheta*PointNumInPhi+j*PointNumInPhi+q_index*PointTotalNum;
				for(l=phi_l; l<phi_u; l++){
					count=0;
					for(k=0; k<PointNumInPhi; k++){
						count+=energyDataPointer->gradMatrix_phi[l][k] * energyDataPointer->Q_elements[k+indexBasis];
					}
					energyDataPointer->dQ_drtp[2][l+indexBasis]=count;
				}
			}
		}

	}

	double dQ_dxyz[3][5]={0};		//行dx,dy,dz, 列q1,q2,...,q5
	double x, y, z;
	double norm, norm_2, semiNorm, semiNorm_2;
	double cordTrans[3][3];

	energy=0;

	for(n=radius_l; n<radius_u; n++){
		for(m=theta_l; m<theta_u; m++){
			for(l=phi_l; l<phi_u; l++){
				i=Point_index(n,m,l);
				//将“对rtp的导数”转化为“对xyz的导数”的准备
				x=energyDataPointer->pointsInBall[0][i];
				y=energyDataPointer->pointsInBall[1][i];
				z=energyDataPointer->pointsInBall[2][i];
				semiNorm_2=x*x+y*y;
				semiNorm=sqrt(semiNorm_2);
				norm_2=semiNorm_2+z*z;
				norm=sqrt(norm_2);
				cordTrans[0][0]=x/norm; cordTrans[0][1]=x*z/semiNorm/norm_2; cordTrans[0][2]=-y/semiNorm_2;
				cordTrans[1][0]=y/norm; cordTrans[1][1]=y*z/semiNorm/norm_2; cordTrans[1][2]=x/semiNorm_2;
				cordTrans[2][0]=z/norm; cordTrans[2][1]=-semiNorm/norm_2;	 cordTrans[2][2]=0;
				//算出i点处q1-q5的对xyz的偏导数值
				for(q_index=0; q_index<5; q_index++){
					for(j=0; j<3; j++){
						dQ_dxyz[j][q_index]=cordTrans[j][0]*energyDataPointer->dQ_drtp[0][i+q_index*PointTotalNum] + cordTrans[j][1]*energyDataPointer->dQ_drtp[1][i+q_index*PointTotalNum] + cordTrans[j][2]*energyDataPointer->dQ_drtp[2][i+q_index*PointTotalNum];
					}
				}

				//计算泛函值
				count=0;
				for(q_index=0; q_index<5; q_index++){
					count+=dQ_dxyz[0][q_index]*dQ_dxyz[0][q_index]+dQ_dxyz[1][q_index]*dQ_dxyz[1][q_index]+dQ_dxyz[2][q_index]*dQ_dxyz[2][q_index];
				}
				energy += energyDataPointer->ballInteCor[i]*(count + dQ_dxyz[0][0]*dQ_dxyz[0][3]+dQ_dxyz[1][0]*dQ_dxyz[1][3]+dQ_dxyz[2][0]*dQ_dxyz[2][3]);

				//计算dEnergy_dQdrtp
				for(q_index=1; q_index<5; q_index++){		//q2, q3, q5
					if(q_index==3) continue;
					for(j=0; j<3; j++){
						energyDataPointer->dEnergy_dQdrtp[j][i+q_index*PointTotalNum]=Alpha2*energyDataPointer->ballInteCor[i]*2*(cordTrans[0][j]*dQ_dxyz[0][q_index] + cordTrans[1][j]*dQ_dxyz[1][q_index] + cordTrans[2][j]*dQ_dxyz[2][q_index]);
					}
				}
				for(j=0; j<3; j++){							//q1
					energyDataPointer->dEnergy_dQdrtp[j][i]=Alpha2*energyDataPointer->ballInteCor[i]*(cordTrans[0][j]*(2*dQ_dxyz[0][0]+dQ_dxyz[0][3]) + cordTrans[1][j]*(2*dQ_dxyz[1][0]+dQ_dxyz[1][3]) + cordTrans[2][j]*(2*dQ_dxyz[2][0]+dQ_dxyz[2][3]));
				}
				for(j=0; j<3; j++){							//q4
					energyDataPointer->dEnergy_dQdrtp[j][i+3*PointTotalNum]=Alpha2*energyDataPointer->ballInteCor[i]*(cordTrans[0][j]*(2*dQ_dxyz[0][3]+dQ_dxyz[0][0]) + cordTrans[1][j]*(2*dQ_dxyz[1][3]+dQ_dxyz[1][0]) + cordTrans[2][j]*(2*dQ_dxyz[2][3]+dQ_dxyz[2][0]));
				}

			}
		}
	}

	energyDataPointer->energy += Alpha2*energy;

	// mpi传输dEnergy_dQdrtp
	exchange_vals(rank, energyDataPointer->dEnergy_dQdrtp, energyDataPointer->dEnergy_dQdrtp, mpi_comm_space, 0);
	/*
	int n_index, m_index, l_index, block_index;
	int local_total_num = (radius_u-radius_l)*(theta_u-theta_l)*(phi_u-phi_l);
		//radius组
	double* dEnergy_dQdr_msg = (double*) malloc(5*local_total_num*sizeof(double));
	double* dEnergy_dQdr_rev = (double*) malloc(RadiusSlice*5*local_total_num*sizeof(double));

	for(n=radius_l; n<radius_u; n++){
		n_index = (n-radius_l)*(theta_u-theta_l)*(phi_u-phi_l);
		for(m=theta_l; m<theta_u; m++){
			m_index = (m-theta_l)*(phi_u-phi_l);
			for(l=phi_l; l<phi_u; l++){
				l_index = l-phi_l;
				for(q_index=0; q_index<5; q_index++){
					dEnergy_dQdr_msg[n_index+m_index+l_index+q_index*local_total_num] = energyDataPointer->dEnergy_dQdrtp[0][Point_index(n,m,l)+q_index*PointTotalNum];
				}
			}
		}
	}
	MPI_Allgather(dEnergy_dQdr_msg, 5*local_total_num, MPI_DOUBLE, dEnergy_dQdr_rev, 5*local_total_num, MPI_DOUBLE, mpi_radius_group);
	for(n=0; n<PointNumInRadius; n++){
		n_index = n%(radius_u-radius_l)*(theta_u-theta_l)*(phi_u-phi_l);
		block_index = n/(radius_u-radius_l);
		for(m=theta_l; m<theta_u; m++){
			m_index = (m-theta_l)*(phi_u-phi_l);
			for(l=phi_l; l<phi_u; l++){
				l_index = l-phi_l;
				for(q_index=0; q_index<5; q_index++){
					energyDataPointer->dEnergy_dQdrtp[0][Point_index(n,m,l)+q_index*PointTotalNum] = dEnergy_dQdr_rev[block_index*5*local_total_num+n_index+m_index+l_index+q_index*local_total_num];
				}
			}
		}
	}
	free((void *)dEnergy_dQdr_msg);
	free((void *)dEnergy_dQdr_rev);

		//theta组
	double* dEnergy_dQdt_msg = (double*) malloc(5*local_total_num*sizeof(double));
	double* dEnergy_dQdt_rev = (double*) malloc(ThetaSlice*5*local_total_num*sizeof(double));

	for(n=radius_l; n<radius_u; n++){
		n_index = (n-radius_l)*(theta_u-theta_l)*(phi_u-phi_l);
		for(m=theta_l; m<theta_u; m++){
			m_index = (m-theta_l)*(phi_u-phi_l);
			for(l=phi_l; l<phi_u; l++){
				l_index = l-phi_l;
				for(q_index=0; q_index<5; q_index++){
					dEnergy_dQdt_msg[n_index+m_index+l_index+q_index*local_total_num] = energyDataPointer->dEnergy_dQdrtp[1][Point_index(n,m,l)+q_index*PointTotalNum];
				}
			}
		}
	}
	MPI_Allgather(dEnergy_dQdt_msg, 5*local_total_num, MPI_DOUBLE, dEnergy_dQdt_rev, 5*local_total_num, MPI_DOUBLE, mpi_theta_group);
	for(m=0; m<PointNumInTheta; m++){
		m_index = m%(theta_u-theta_l)*(phi_u-phi_l);
		block_index = m/(theta_u-theta_l);
		for(n=radius_l; n<radius_u; n++){
			n_index = (n-radius_l)*(theta_u-theta_l)*(phi_u-phi_l);
			for(l=phi_l; l<phi_u; l++){
				l_index = l-phi_l;
				for(q_index=0; q_index<5; q_index++){
					energyDataPointer->dEnergy_dQdrtp[1][Point_index(n,m,l)+q_index*PointTotalNum] = dEnergy_dQdt_rev[block_index*5*local_total_num+n_index+m_index+l_index+q_index*local_total_num];
				}
			}
		}
	}
	free((void *)dEnergy_dQdt_msg);
	free((void *)dEnergy_dQdt_rev);

		//phi组
	double* dEnergy_dQdp_msg = (double*) malloc(5*local_total_num*sizeof(double));
	double* dEnergy_dQdp_rev = (double*) malloc(PhiSlice*5*local_total_num*sizeof(double));

	for(n=radius_l; n<radius_u; n++){
		n_index = (n-radius_l)*(theta_u-theta_l)*(phi_u-phi_l);
		for(m=theta_l; m<theta_u; m++){
			m_index = (m-theta_l)*(phi_u-phi_l);
			for(l=phi_l; l<phi_u; l++){
				l_index = l-phi_l;
				for(q_index=0; q_index<5; q_index++){
					dEnergy_dQdp_msg[n_index+m_index+l_index+q_index*local_total_num] = energyDataPointer->dEnergy_dQdrtp[2][Point_index(n,m,l)+q_index*PointTotalNum];
				}
			}
		}
	}
	MPI_Allgather(dEnergy_dQdp_msg, 5*local_total_num, MPI_DOUBLE, dEnergy_dQdp_rev, 5*local_total_num, MPI_DOUBLE, mpi_phi_group);
	for(l=0; l<PointNumInPhi; l++){
		l_index = l%(phi_u-phi_l);
		block_index = l/(phi_u-phi_l);
		for(n=radius_l; n<radius_u; n++){
			n_index = (n-radius_l)*(theta_u-theta_l)*(phi_u-phi_l);
			for(m=theta_l; m<theta_u; m++){
				m_index = (m-theta_l)*(phi_u-phi_l);
				for(q_index=0; q_index<5; q_index++){
					energyDataPointer->dEnergy_dQdrtp[2][Point_index(n,m,l)+q_index*PointTotalNum] = dEnergy_dQdp_rev[block_index*5*local_total_num+n_index+m_index+l_index+q_index*local_total_num];
				}
			}
		}
	}
	free((void *)dEnergy_dQdp_msg);
	free((void *)dEnergy_dQdp_rev);
	*/

	//继续计算偏导数

	for(q_index=0; q_index<5; q_index++){

		for(i=radius_l; i<radius_u; i++){
			for(j=theta_l; j<theta_u; j++){
				for(k=phi_l; k<phi_u; k++){
					count=0;
					for(l=0; l<PointNumInRadius; l++){
						count+=energyDataPointer->gradMatrix_radius[l][i] * energyDataPointer->dEnergy_dQdrtp[0][Point_index(l,j,k)+q_index*PointTotalNum];
					}
					for(l=0; l<PointNumInTheta; l++){
						count+=energyDataPointer->gradMatrix_theta[l][j] * energyDataPointer->dEnergy_dQdrtp[1][Point_index(i,l,k)+q_index*PointTotalNum];
					}
					for(l=0; l<PointNumInPhi; l++){
						count+=energyDataPointer->gradMatrix_phi[l][k] * energyDataPointer->dEnergy_dQdrtp[2][Point_index(i,j,l)+q_index*PointTotalNum];
					}
					energyDataPointer->dEnergy_dQ[Point_index(i,j,k)][q_index]=count;
				}
			}
		}

	}


	double* pointer1;
	double* pointer2;
	for(n=radius_l; n<radius_u; n++){
		for(m=theta_l; m<theta_u; m++){
			for(l=phi_l; l<phi_u; l++){
				i=Point_index(n,m,l);
				pointer1=energyDataPointer->dQ_dBdiagAngle[i];
				for(j=0; j<5; j++){
					pointer2=energyDataPointer->dEnergy_dQ[i];
					count=0;
					for(k=0; k<5; k++){
						count+=(*pointer1)*(*pointer2);
						pointer1++;
						pointer2++;
					}
					energyDataPointer->dEnergy_dB[j*PointTotalNum+i]+=count;
				}
			}
		}
	}
	pointer1=NULL;
	pointer2=NULL;

}

void calc_Fpena(int rank){

	int i, j, k, n, m, l, q_index;
	int radius_l = Radius_Low(rank);
	int radius_u = Radius_Up(rank);
	int theta_l = Theta_Low(rank);
	int theta_u = Theta_Up(rank);
	int phi_l = Phi_Low(rank);
	int phi_u = Phi_Up(rank);
	int indexBasis, index;
	double row1, row2, row3;
	double x, y, z;
	double q[5]={0};
	double q_std[5]={0};
	double p[4]={0};
	double dFpena_dq[5]={0};
	double count=0, energy=0;

	energy=0;

	for(j=theta_l; j<theta_u; j++){
		for(k=phi_l; k<phi_u; k++){
			//计算球面上某点处的函数值(q1-q5)
			for(q_index=0; q_index<5; q_index++){
				indexBasis=j*PointNumInPhi+k+q_index*PointTotalNum;
				count=0;
				for(i=0; i<PointNumInRadius; i++){
					count+=energyDataPointer->sideValVector[i] * energyDataPointer->Q_elements[i*PointNumInTheta*PointNumInPhi+indexBasis];
				}
				q[q_index]=count;
			}

			//松弛锚定边界条件
			//计算泛函值
			index=j*PointNumInPhi+k;
			x=energyDataPointer->pointsInSphere[0][index];
			y=energyDataPointer->pointsInSphere[1][index];
			z=energyDataPointer->pointsInSphere[2][index];
			p[0]=q[0]*x*y-q[1]*(x*x-1.0/3);
			p[1]=q[1]*z-q[2]*y;
			p[2]=q[3]*x*y-q[1]*(y*y-1.0/3);
			p[3]=q[1]*z-q[4]*x;

			if(radius_u==PointNumInRadius){
				count=0;
				for(i=0; i<4; i++){
					count+=p[i]*p[i];
				}
				energy += energyDataPointer->sphereInteCor[index]*count;
			}

			//计算偏导数值
			dFpena_dq[0]=Eta*energyDataPointer->sphereInteCor[index]*2*(p[0]*x*y);
			dFpena_dq[1]=Eta*energyDataPointer->sphereInteCor[index]*2*(-p[0]*(x*x-1.0/3)+p[1]*z-p[2]*(y*y-1.0/3)+p[3]*z);
			dFpena_dq[2]=Eta*energyDataPointer->sphereInteCor[index]*2*(-p[1]*y);
			dFpena_dq[3]=Eta*energyDataPointer->sphereInteCor[index]*2*(p[2]*x*y);
			dFpena_dq[4]=Eta*energyDataPointer->sphereInteCor[index]*2*(-p[3]*x);

			for(q_index=0; q_index<5; q_index++){
				for(i=radius_l; i<radius_u; i++){
					energyDataPointer->dEnergy_dQ[i*PointNumInTheta*PointNumInPhi+index][q_index]=energyDataPointer->sideValVector[i] * dFpena_dq[q_index];
				}
			}

			/*//强锚定边界条件 罚项 |Q-Q'|^2
			//计算泛函值
			index=j*PointNumInPhi+k;
			x=energyDataPointer->pointsInSphere[0][index];
			y=energyDataPointer->pointsInSphere[1][index];
			z=energyDataPointer->pointsInSphere[2][index];
			q_std[0]=S*(x*x-1.0/3);
			q_std[1]=S*x*y;
			q_std[2]=S*x*z;
			q_std[3]=S*(y*y-1.0/3);
			q_std[4]=S*y*z;

			if(radius_u==PointNumInRadius){
				count=0;
				for(i=0; i<5; i++){
					count+=(q[i]-q_std[i])*(q[i]-q_std[i]);
				}
				energy += energyDataPointer->sphereInteCor[index]*count;
			}

			//计算偏导数值
			for(i=0; i<5; i++){
				dFpena_dq[i]=Eta*energyDataPointer->sphereInteCor[index]*2*(q[i]-q_std[i]);
			}

			for(q_index=0; q_index<5; q_index++){
				for(i=radius_l; i<radius_u; i++){
					energyDataPointer->dEnergy_dQ[i*PointNumInTheta*PointNumInPhi+index][q_index]=energyDataPointer->sideValVector[i] * dFpena_dq[q_index];
				}
			}
			*/

			/*//罚项 |Q*x|^2
			//计算泛函值
			index=j*PointNumInPhi+k;
			x=energyDataPointer->pointsInSphere[0][index];
			y=energyDataPointer->pointsInSphere[1][index];
			z=energyDataPointer->pointsInSphere[2][index];
			row1=q[0]*x + q[1]*y + q[2]*z;
			row2=q[1]*x + q[3]*y + q[4]*z;
			row3=q[2]*x + q[4]*y - (q[0]+q[3])*z;

			if(radius_u==PointNumInRadius) energy += energyDataPointer->sphereInteCor[index]*(row1*row1 + row2*row2 + row3*row3);

			//计算偏导数值
			dFpena_dq[0]=Eta*energyDataPointer->sphereInteCor[index]*2*(row1*x - row3*z);
			dFpena_dq[1]=Eta*energyDataPointer->sphereInteCor[index]*2*(row1*y + row2*x);
			dFpena_dq[2]=Eta*energyDataPointer->sphereInteCor[index]*2*(row1*z + row3*x);
			dFpena_dq[3]=Eta*energyDataPointer->sphereInteCor[index]*2*(row2*y - row3*z);
			dFpena_dq[4]=Eta*energyDataPointer->sphereInteCor[index]*2*(row2*z + row3*y);

			for(q_index=0; q_index<5; q_index++){
				for(i=radius_l; i<radius_u; i++){
					energyDataPointer->dEnergy_dQ[i*PointNumInTheta*PointNumInPhi+index][q_index]=energyDataPointer->sideValVector[i] * dFpena_dq[q_index];
				}
			}
			*/
		}
	}

	if(radius_u==PointNumInRadius){
		energyDataPointer->energy += Eta*energy;
	}

	double* pointer1;
	double* pointer2;
	for(n=radius_l; n<radius_u; n++){
		for(m=theta_l; m<theta_u; m++){
			for(l=phi_l; l<phi_u; l++){
				i=Point_index(n,m,l);
				pointer1=energyDataPointer->dQ_dBdiagAngle[i];
				for(j=0; j<5; j++){
					pointer2=energyDataPointer->dEnergy_dQ[i];
					count=0;
					for(k=0; k<5; k++){
						count+=(*pointer1)*(*pointer2);
						pointer1++;
						pointer2++;
					}
					energyDataPointer->dEnergy_dB[j*PointTotalNum+i]+=count;
				}
			}
		}
	}
	pointer1=NULL;
	pointer2=NULL;

}


void setInitValue(){
	int i, j, k;

	//无序初值
	/* for(i=0; i<5*FreeVarsUnitNum; i++){
		energyDataPointer->freeVars[i]=0;
	} */

	//随机初值
	/*srand((unsigned)time(NULL));
    for(i=0;i<5*PointTotalNum;i++){
		energyDataPointer->Bdiag_Angle_Elements[i] = 0*(2.0*rand()/RAND_MAX - 1);
	}*/

	//随机初值
	/*srand((unsigned)time(NULL));
    for(i=0;i<FreeVarsTotalNum;i++){
		energyDataPointer->freeVars[i] = 2*(2.0*rand()/RAND_MAX - 1);
	}*/

	//hedgehog初值
	/*int indexBasis;
	for(i=0; i<PointNumInRadius; i++){
		energyDataPointer->freeVars[i]=-2.375*energyDataPointer->pointInRadius[i];
		for(j=0; j<PointNumInTheta; j++){
			indexBasis=i*PointNumInTheta*PointNumInPhi+j*PointNumInPhi;
			for(k=0; k<PointNumInPhi; k++){
				energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+2*PointTotalNum]=0;
				energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+3*PointTotalNum]=energyDataPointer->pointInTheta[j];
				energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+4*PointTotalNum]=PI/2-energyDataPointer->pointInPhi[k];
			}
		}
	}*/


	//hedgehog初值，ring初值
	int indexBasis;
	for(i=0; i<PointNumInRadius; i++){
		for(j=0; j<PointNumInTheta; j++){
			indexBasis=i*PointNumInTheta*PointNumInPhi+j*PointNumInPhi;
			for(k=0; k<PointNumInPhi; k++){
				energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+0*PointTotalNum]=-2.375*energyDataPointer->pointInRadius[i];
				energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+1*PointTotalNum]=-2.375*energyDataPointer->pointInRadius[i];
				energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+2*PointTotalNum]=0;
				energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+3*PointTotalNum]=energyDataPointer->pointInTheta[j];
				energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+4*PointTotalNum]=PI/2-energyDataPointer->pointInPhi[k];
			}
		}
	}

	//makeFreeVars();

	//sphereRing解
	/*int indexBasis;
	for(i=0; i<PointNumInRadius; i++){
		for(j=0; j<PointNumInTheta; j++){
			indexBasis=i*PointNumInTheta*PointNumInPhi+j*PointNumInPhi;
			for(k=0; k<PointNumInPhi; k++){
				energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+0*PointTotalNum]=-2;
				energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+1*PointTotalNum]=-2;
				energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+2*PointTotalNum]=0;
				energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+3*PointTotalNum]=0;
				energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+4*PointTotalNum]=0;
			}
		}
	}

	//makeFreeVars();*/

	//split core初值
	/*int indexBasis;
	for(i=0; i<PointNumInRadius; i++){
		for(j=0; j<PointNumInTheta; j++){
			indexBasis=i*PointNumInTheta*PointNumInPhi+j*PointNumInPhi;
			for(k=0; k<PointNumInPhi; k++){
				/*if(fabs(energyDataPointer->pointInRadius[i]*cos(energyDataPointer->pointInTheta[j]))<SplitCorePoint){
					energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+0*PointTotalNum]=-1*energyDataPointer->pointInRadius[i]-0.5;
					energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+1*PointTotalNum]=-1*energyDataPointer->pointInRadius[i]-0.5;
					energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+2*PointTotalNum]=0;
					energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+3*PointTotalNum]=PI/2;
					energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+4*PointTotalNum]=PI/2-energyDataPointer->pointInPhi[k];
				}
				else{
					energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+0*PointTotalNum]=-2*(energyDataPointer->pointInRadius[i]-SplitCorePoint);
					energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+1*PointTotalNum]=-2*(energyDataPointer->pointInRadius[i]-SplitCorePoint);
					energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+2*PointTotalNum]=0;
					if(energyDataPointer->pointInTheta[j]<PI/2){
						energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+3*PointTotalNum]=atan(energyDataPointer->pointInRadius[i]*sin(energyDataPointer->pointInTheta[j])/(energyDataPointer->pointInRadius[i]*cos(energyDataPointer->pointInTheta[j])-SplitCorePoint));
					}
					else{
						energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+3*PointTotalNum]=atan(energyDataPointer->pointInRadius[i]*sin(energyDataPointer->pointInTheta[j])/(energyDataPointer->pointInRadius[i]*cos(energyDataPointer->pointInTheta[j])+SplitCorePoint))+PI;
					}
					energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+4*PointTotalNum]=PI/2-energyDataPointer->pointInPhi[k];
				}*
				if(energyDataPointer->pointInRadius[i]<0.9){
					energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+0*PointTotalNum]=5;
					energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+1*PointTotalNum]=5;
				}
				else{
					energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+0*PointTotalNum]=-2*energyDataPointer->pointInRadius[i];
					energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+1*PointTotalNum]=-2*energyDataPointer->pointInRadius[i];
				}
				energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+2*PointTotalNum]=0;
				energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+3*PointTotalNum]=energyDataPointer->pointInTheta[j];
				energyDataPointer->Bdiag_Angle_Elements[indexBasis+k+4*PointTotalNum]=PI/2-energyDataPointer->pointInPhi[k];
			}
		}
	}

	//makeFreeVars();*/

	//读取已有结果
	/*FILE* fp;
	if((fp=fopen("result/B.txt","r"))==NULL) exit(1);
	for(i=0; i<5*PointTotalNum; i++){
		fscanf(fp, "%lf ", &energyDataPointer->Bdiag_Angle_Elements[i]);
	}
	fclose(fp);

	//makeFreeVars();*/
}
