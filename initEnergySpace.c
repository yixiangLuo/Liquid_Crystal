#include "initEnergySpace.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gaussian_nodes.h"

void energySpaceFree(energySpace* energyDataPointer){
	int i;

	if(!energyDataPointer) return;

	//释放自由变量
	//free((void *)energyDataPointer->freeVars);

	//释放自变量B
	free((void *)energyDataPointer->Bdiag_Angle_Elements);

	//释放Q的特征值
	free((void *)energyDataPointer->Qdiag_elements);

	//释放Q
	free((void *)energyDataPointer->Q_elements);

	//释放Z
	free((void *)energyDataPointer->Z_elements);

	//释放Z关于B的特征值的偏导数
	freeMatrixSpace(energyDataPointer->dZ_dBdiag, 2);

	//释放Q的特征值关于B的特征值的偏导数
	freeMatrixSpace(energyDataPointer->dQdiag_dBdiag, 4);

	//释放Q关于B的偏导数
	freeMatrixSpace(energyDataPointer->dQ_dBdiagAngle, PointTotalNum);

	//释放能量积分关于自变量B的偏导数
	free((void *)energyDataPointer->dEnergy_dB);

	//释放能量积分关于自由变量的偏导数
	//free((void *)energyDataPointer->dEnergy_freeVars);

	//释放r方向上的节点
	free((void *)energyDataPointer->pointInRadius);

	//释放theta方向上的节点
	free((void *)energyDataPointer->pointInTheta);

	//释放phi方向上的节点
	free((void *)energyDataPointer->pointInPhi);

	//释放单位球内的数值积分节点
	freeMatrixSpace(energyDataPointer->pointsInBall, 3);

	//释放单位球面上的数值积分节点
	freeMatrixSpace(energyDataPointer->pointsInSphere, 3);

	//释放单位球内的数值积分系数
	free((void *)energyDataPointer->ballInteCor);

	//释放单位球面上的数值积分系数
	free((void *)energyDataPointer->sphereInteCor);

	//多项式方法
	//释放求偏导数的转换矩阵 radius方向
	freeMatrixSpace(energyDataPointer->gradMatrix_radius, PointNumInRadius);

	//释放求偏导数的转换矩阵 theta方向
	freeMatrixSpace(energyDataPointer->gradMatrix_theta, PointNumInTheta);

	//释放求偏导数的转换矩阵 phi方向
	freeMatrixSpace(energyDataPointer->gradMatrix_phi, PointNumInPhi);

	//释放计算球面上函数值的系数向量(r方向)
	free((void *)energyDataPointer->sideValVector);

	//释放存储中间值的空间
	free((void *)energyDataPointer->dQ_drtp);
	free((void *)energyDataPointer->dEnergy_dQdrtp);
	free((void *)energyDataPointer->dEnergy_dQ);

	//释放结构体
	free((void *)energyDataPointer);

	return ;
}


energySpace* initEnergySpace(){
	double* p ;
	FILE* fp;
	int i, j, n, m, l;
	double deno, nume, temp;

	double gaussian_node[64][100] = Gaussian_Nodes;
	double gaussian_weight[64][100] = Gaussian_Weights;

	//分配指针结构体
	energySpace* energyDataPointer = (energySpace*) malloc(sizeof(energySpace));
	if(energyDataPointer==NULL){
		printf("fail to initiate energy space: menmory not enough!\n");
		exit(1);
	}

	//初始化
	//energyDataPointer->freeVars=NULL;
	energyDataPointer->Bdiag_Angle_Elements=NULL;
	energyDataPointer->Qdiag_elements=NULL;
	energyDataPointer->Q_elements=NULL;
	energyDataPointer->Z_elements=NULL;
	energyDataPointer->dZ_dBdiag=NULL;
	energyDataPointer->dQdiag_dBdiag=NULL;
	energyDataPointer->dQ_dBdiagAngle=NULL;
	energyDataPointer->dEnergy_dB=NULL;
	//energyDataPointer->dEnergy_freeVars=NULL;
	energyDataPointer->pointInRadius=NULL;
	energyDataPointer->pointInTheta=NULL;
	energyDataPointer->pointInPhi=NULL;
	energyDataPointer->pointsInBall=NULL;
	energyDataPointer->pointsInSphere=NULL;
	energyDataPointer->ballInteCor=NULL;
	energyDataPointer->sphereInteCor=NULL;
	//多项式方法
	energyDataPointer->gradMatrix_radius=NULL;
	energyDataPointer->gradMatrix_theta=NULL;
	energyDataPointer->gradMatrix_phi=NULL;
	energyDataPointer->sideValVector=NULL;
	//存储计算中间值
	energyDataPointer->dQ_drtp=NULL;
	energyDataPointer->dEnergy_dQdrtp=NULL;
	energyDataPointer->dEnergy_dQ=NULL;

	//分配自由变量的存储空间
	//energyDataPointer->freeVars = linearSpace_Energy(FreeVarsTotalNum, energyDataPointer);

	//分配自变量B的存储空间
	energyDataPointer->Bdiag_Angle_Elements = linearSpace_Energy(5*PointTotalNum, energyDataPointer);

	//分配Q的特征值的存储空间
	energyDataPointer->Qdiag_elements = linearSpace_Energy(2*PointTotalNum, energyDataPointer);

	//分配Q的存储空间
	energyDataPointer->Q_elements = linearSpace_Energy(5*PointTotalNum, energyDataPointer);

	//分配Z的存储空间
	energyDataPointer->Z_elements = linearSpace_Energy(PointTotalNum, energyDataPointer);

	//分配Z关于B的特征值的偏导数的存储空间
	energyDataPointer->dZ_dBdiag = matrixSpace_Energy(2, PointTotalNum, energyDataPointer);

	//分配Q的特征值关于B的特征值的偏导数的存储空间
	energyDataPointer->dQdiag_dBdiag = matrixSpace_Energy(4, PointTotalNum, energyDataPointer);

	//分配Q关于B的偏导数的存储空间
	energyDataPointer->dQ_dBdiagAngle = matrixSpace_Energy(PointTotalNum, 25, energyDataPointer);

	//分配能量积分关于自变量B的偏导数的存储空间
	energyDataPointer->dEnergy_dB = linearSpace_Energy(5*PointTotalNum, energyDataPointer);

	//分配能量积分关于自由变量的偏导数的存储空间
	//energyDataPointer->dEnergy_freeVars = linearSpace_Energy(FreeVarsTotalNum, energyDataPointer);

	//分配并计算r方向上的节点
	energyDataPointer->pointInRadius = linearSpace_Energy(PointNumInRadius, energyDataPointer);

		//计算
	for(i=0; i<PointNumInRadius; i++){
		energyDataPointer->pointInRadius[i] = 1.0/2*gaussian_node[Gaussian_index(PointNumInRadius)][i] + 0.5;
	}

	//分配并计算theta方向上的节点
	energyDataPointer->pointInTheta = linearSpace_Energy(PointNumInTheta, energyDataPointer);

		//计算
	for(i=0; i<PointNumInTheta; i++){
		energyDataPointer->pointInTheta[i] = PI/2.0*gaussian_node[Gaussian_index(PointNumInTheta)][i] + PI/2.0;
	}

	//分配并计算theta方向上的节点
	energyDataPointer->pointInPhi = linearSpace_Energy(PointNumInPhi, energyDataPointer);

		//计算
	for(i=0; i<PointNumInPhi; i++){
		energyDataPointer->pointInPhi[i] = PI*gaussian_node[Gaussian_index(PointNumInPhi)][i] + PI;
	}

	//分配并计算单位球中的数值积分节点
	energyDataPointer->pointsInBall = matrixSpace_Energy(3, PointTotalNum, energyDataPointer);

		//计算
	for(n=0; n<PointNumInRadius; n++){
		for(m=0; m<PointNumInTheta; m++){
			for(l=0; l<PointNumInPhi; l++){
				energyDataPointer->pointsInBall[0][Point_index(n,m,l)] = energyDataPointer->pointInRadius[n]*sin(energyDataPointer->pointInTheta[m])*cos(energyDataPointer->pointInPhi[l]);
				energyDataPointer->pointsInBall[1][Point_index(n,m,l)] = energyDataPointer->pointInRadius[n]*sin(energyDataPointer->pointInTheta[m])*sin(energyDataPointer->pointInPhi[l]);
				energyDataPointer->pointsInBall[2][Point_index(n,m,l)] = energyDataPointer->pointInRadius[n]*cos(energyDataPointer->pointInTheta[m]);
			}
		}
	}

	//分配并计算单位球面上的数值积分节点
	energyDataPointer->pointsInSphere = matrixSpace_Energy(3, PointNumInTheta*PointNumInPhi, energyDataPointer);

		//计算
	for(m=0; m<PointNumInTheta; m++){
		for(l=0; l<PointNumInPhi; l++){
			energyDataPointer->pointsInSphere[0][m*PointNumInPhi+l] = sin(energyDataPointer->pointInTheta[m])*cos(energyDataPointer->pointInPhi[l]);
			energyDataPointer->pointsInSphere[1][m*PointNumInPhi+l] = sin(energyDataPointer->pointInTheta[m])*sin(energyDataPointer->pointInPhi[l]);
			energyDataPointer->pointsInSphere[2][m*PointNumInPhi+l] = cos(energyDataPointer->pointInTheta[m]);
		}
	}

	//分配并计算单位球中的数值积分系数
	energyDataPointer->ballInteCor = linearSpace_Energy(PointTotalNum, energyDataPointer);

		//计算
	for(n=0; n<PointNumInRadius; n++){
		for(m=0; m<PointNumInTheta; m++){
			for(l=0; l<PointNumInPhi; l++){
				energyDataPointer->ballInteCor[Point_index(n,m,l)] = energyDataPointer->pointInRadius[n]*energyDataPointer->pointInRadius[n]*1.0/2*gaussian_weight[Gaussian_index(PointNumInRadius)][n] * sin(energyDataPointer->pointInTheta[m])*PI/2.0*gaussian_weight[Gaussian_index(PointNumInTheta)][m] * PI*gaussian_weight[Gaussian_index(PointNumInPhi)][l];
			}
		}
	}

	//分配并计算单位球面上的数值积分系数
	energyDataPointer->sphereInteCor = linearSpace_Energy(PointNumInTheta*PointNumInPhi, energyDataPointer);

		//计算
	for(m=0; m<PointNumInTheta; m++){
		for(l=0; l<PointNumInPhi; l++){
			energyDataPointer->sphereInteCor[m*PointNumInPhi+l] = sin(energyDataPointer->pointInTheta[m])*PI/2.0*gaussian_weight[Gaussian_index(PointNumInTheta)][m] * PI*gaussian_weight[Gaussian_index(PointNumInPhi)][l];
		}
	}

	//多项式方法
	//分配并计算求偏导数的转换矩阵 radius方向
	energyDataPointer->gradMatrix_radius = matrixSpace_Energy(PointNumInRadius, PointNumInRadius, energyDataPointer);

		//计算
	for(n=0; n<PointNumInRadius; n++){
		deno=1;
		for(i=0; i<PointNumInRadius; i++){
			if(i!=n) deno *= (energyDataPointer->pointInRadius[n]-energyDataPointer->pointInRadius[i]);
		}
		for(m=0; m<PointNumInRadius; m++){
			nume=0;
			for(l=0; l<PointNumInRadius; l++){
				if(l==n) continue;
				temp=1;
				for(i=0; i<PointNumInRadius; i++){
					if(i!=n && i!=l) temp *= (energyDataPointer->pointInRadius[m]-energyDataPointer->pointInRadius[i]);
				}
				nume+=temp;
			}
			energyDataPointer->gradMatrix_radius[m][n]=nume/deno;
		}
	}

	//分配并计算求偏导数的转换矩阵 theta方向
	energyDataPointer->gradMatrix_theta = matrixSpace_Energy(PointNumInTheta, PointNumInTheta, energyDataPointer);

		//计算
	for(n=0; n<PointNumInTheta; n++){
		deno=1;
		for(i=0; i<PointNumInTheta; i++){
			if(i!=n) deno *= (energyDataPointer->pointInTheta[n]-energyDataPointer->pointInTheta[i]);
		}
		for(m=0; m<PointNumInTheta; m++){
			nume=0;
			for(l=0; l<PointNumInTheta; l++){
				if(l==n) continue;
				temp=1;
				for(i=0; i<PointNumInTheta; i++){
					if(i!=n && i!=l) temp *= (energyDataPointer->pointInTheta[m]-energyDataPointer->pointInTheta[i]);
				}
				nume+=temp;
			}
			energyDataPointer->gradMatrix_theta[m][n]=nume/deno;
		}
	}

	//分配并计算求偏导数的转换矩阵 phi方向
	energyDataPointer->gradMatrix_phi = matrixSpace_Energy(PointNumInPhi, PointNumInPhi, energyDataPointer);

		//计算
	for(n=0; n<PointNumInPhi; n++){
		deno=1;
		for(i=0; i<PointNumInPhi; i++){
			if(i!=n) deno *= (energyDataPointer->pointInPhi[n]-energyDataPointer->pointInPhi[i]);
		}
		for(m=0; m<PointNumInPhi; m++){
			nume=0;
			for(l=0; l<PointNumInPhi; l++){
				if(l==n) continue;
				temp=1;
				for(i=0; i<PointNumInPhi; i++){
					if(i!=n && i!=l) temp *= (energyDataPointer->pointInPhi[m]-energyDataPointer->pointInPhi[i]);
				}
				nume+=temp;
			}
			energyDataPointer->gradMatrix_phi[m][n]=nume/deno;
		}
	}

	//分配并计算计算球面上函数值的系数向量(r方向)
	energyDataPointer->sideValVector = linearSpace_Energy(PointNumInRadius, energyDataPointer);

		//计算
	for(n=0; n<PointNumInRadius; n++){
		deno=1;
		nume=1;
		for(i=0; i<PointNumInRadius; i++){
			if(i==n) continue;
			deno *= (energyDataPointer->pointInRadius[n]-energyDataPointer->pointInRadius[i]);
			nume *= (1-energyDataPointer->pointInRadius[i]);
		}
		energyDataPointer->sideValVector[n]=nume/deno;
	}

	//分配存储计算中间值的空间
	energyDataPointer->dQ_drtp=matrixSpace_Energy(3, 5*PointTotalNum, energyDataPointer);
	energyDataPointer->dEnergy_dQdrtp=matrixSpace_Energy(3, 5*PointTotalNum, energyDataPointer);
	energyDataPointer->dEnergy_dQ=matrixSpace_Energy(PointTotalNum, 5, energyDataPointer);

	return energyDataPointer;
}

double* linearSpace_Energy(int size, energySpace* energyDataPointer){
	double* dataSpace = (double*) malloc(size*sizeof(double));
	if(dataSpace==NULL){
		energySpaceFree(energyDataPointer);
		printf("fail to allocate temp space: menmory not enough!\n");
		exit(1);
	}
	return dataSpace;
}

void freeMatrixSpace(double** dataSpace, int row){
	int i;
	if(dataSpace){
		for(i=0;i<row;i++){
			free((void *)dataSpace[i]);
		}
	}
	free((void *)dataSpace);
}

double** matrixSpace_Energy(int row, int col, energySpace* energyDataPointer){
	int i;
	double** dataSpace = (double **) malloc(row * sizeof(double * ));
	if(dataSpace==NULL){
		energySpaceFree(energyDataPointer);
		printf("fail to allocate temp space: menmory not enough!\n");
		exit(1);
	}
	for(i=0;i<row;i++) dataSpace[i] = NULL;
	for(i=0;i<row;i++){
		dataSpace[i] = (double*) malloc(col * sizeof(double));
		if(dataSpace[i]==NULL){
			energySpaceFree(energyDataPointer);
			printf("fail to allocate temp space: menmory not enough!\n");
			exit(1);
		}
	}
	return dataSpace;
}
