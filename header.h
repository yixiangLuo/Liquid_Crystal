#ifndef ENERGYFUNCTIONAL
#define ENERGYFUNCTIONAL

#define PI 3.14159265359

#define Alpha1 22.0	//11~-2 16~-6 22~-12
#define Alpha2 0.04
#define Eta 50.0
#define S 0.48	//强锚定边界条件的特征值参数(0,1)
#define PointNumInRadius 32
#define PointNumInTheta 32
#define PointNumInPhi 32
#define PointTotalNum (PointNumInRadius*PointNumInTheta*PointNumInPhi)
//#define FreeVarsUnitNum (PointNumInRadius*PointNumInTheta)
//#define FreeVarsTotalNum (5*FreeVarsUnitNum)

#define RadiusSlice 2
#define ThetaSlice 2
#define PhiSlice 2

#define SplitCorePoint 0.5		//0.01-0.5

#define Gaussian_index(N)  ((N) <= (64) ? (N-2) : (63))
#define Point_index(i,j,k)  ((i)*PointNumInTheta*PointNumInPhi+(j)*PointNumInPhi+(k))

#define max2(a, b)      ((a) >= (b) ? (a) : (b))
#define max3(a, b, c)   max2(max2((a), (b)), (c))


typedef struct energySpace{
	double energy;						//能量泛函的值
	//double* freeVars;					//自由变量，由此生成B，数量FreeVarsUnitNum视限制条件而定，如旋转对称/粘接球体，具体见initEnergySpace
	double* Bdiag_Angle_Elements;		//自变量，5*PointTotalNum个线性排开，前两个为B的两个特征值，后三个为将B对角化的旋转矩阵的欧拉角
	double* Qdiag_elements;				//Q的2个特征值，2*PointTotalNum个线性排开
	double* Q_elements;					//Q的5个元素，5*PointTotalNum个线性排开
	double* Z_elements;					//标准化常数Z
	double** dZ_dBdiag;					//Z对B的特征值的偏导数，二维数组
	double** dQdiag_dBdiag;				//Q的特征值对B的特征值的偏导数，4行分别为dQ1_dB1, dQ1_dB2, dQ2_dB1, dQ2_dB2, PointTotalNum列
	double** dQ_dBdiagAngle;			//Q对B的偏导数，25列，每5列是Qi对Bi(特征值or欧拉角，共5个，Qi的指标先变化)的偏导数
	double* dEnergy_dB;					//能量泛函对5个自变量的偏导数，5*PointTotalNum个线性排开
	//double* dEnergy_freeVars;			//对自由变量的偏导数
	double* pointInRadius;				//r方向上的节点
	double* pointInTheta;				//theta方向上的节点
	double* pointInPhi;					//phi方向上的节点
	double** pointsInBall;				//单位球中的数值积分节点，3行，分别表示x,y,z
	double** pointsInSphere;			//单位球面上的数值积分节点，3行，分别表示x,y,z
	double* ballInteCor;				//单位球中的数值积分系数，PointTotalNum个
	double* sphereInteCor;				//单位球面上的数值积分系数，PointNumInTheta*PointNumInPhi个
	//多项式方法
	double** gradMatrix_radius;			//在r方向由q点处值求q点处对r的偏导数值的矩阵，列为与点值相乘的系数，行为某点处的偏导数
	double** gradMatrix_theta;
	double** gradMatrix_phi;
	double* sideValVector;				//利用r方向上点值计算球面上函数值的系数向量，PointNumInRadius个
	//存储计算中间值
	double** dQ_drtp;
	double** dEnergy_dQdrtp;
	double** dEnergy_dQ;
} energySpace;

#endif
