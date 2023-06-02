/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Plate.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CPlate::CPlate()
{
	NEN_ = 4;	// Each element has 2 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 12;
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CPlate::~CPlate()
{
}

//	Read element data from stream Input
bool CPlate::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2, N3, N4;	// Left node number and right node number

	Input >> N1 >> N2 >> N3 >> N4 >> MSet;
    ElementMaterial_ = dynamic_cast<CPlateMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
    nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];

	return true;
}

//	Write element data to stream
void CPlate::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
//求某一个点的b矩阵，输入用来充当结果的B矩阵，以及该点的坐标，默认为八节点单元






void CPlate::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

//	Calculate Plate length
	double* a=new double[2];		//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	a[0]=0.5*(nodes_[1]->XYZ[0] - nodes_[0]->XYZ[0]);
    a[1] =0.5*(nodes_[2]->XYZ[1] - nodes_[1]->XYZ[1]);



//	Calculate element stiffness matrix

	CPlateMaterial* material_ = dynamic_cast<CPlateMaterial*>(ElementMaterial_);	// Pointer to material of the element
    double E = material_->E;
    double v = material_->v;
    double h = material_->h;
    

   

    int i=0,j=0,m=0,n=0;
    
    double D0=E*h*h*h/(12*(1-v*v));
    Matrix[0]=(v/2.0 + 2.0)/(a[0]*a[0]*a[0]*a[0]) - ((7.0*v)/10.0 - 7.0/10.0)/(a[0]*a[0]*a[1]*a[1]);
    Matrix[1]=4.0/(3.0*a[1]*a[1]) - (8.0*(v/2.0 - 1.0/2.0))/(15.0*a[0]*a[0]);
    Matrix[2]=(4.0*v + 11.0)/(10.0*a[0]*a[0]*a[1]);
    Matrix[3]=4.0/(3.0*a[0]*a[0]) - (8.0*(v/2.0 - 1.0/2.0))/(15.0*a[1]*a[1]);
    Matrix[4]=-v/(a[0]*a[1]);
    Matrix[5]=(v/10.0 - 1/10.0)/(a[0]*a[1]*a[1]) - (v/2.0 + 1.0)/(a[0]*a[0]*a[0]);
    Matrix[6]=(v/2.0 + 2.0)/(a[0]*a[0]*a[0]*a[0]) - ((7.0*v)/10.0 - 7.0/10.0)/(a[0]*a[0]*a[1]*a[1]);
    Matrix[7]=1.0/(a[0]*a[0]*a[0]) - (v/2.0 - 1.0/2.0)/(5.0*a[0]*a[1]*a[1]);
    Matrix[8]=-(2.0*(v - 1.0))/(5.0*a[0]*a[0]*a[1]);
    Matrix[9]=((7.0*v)/10.0 - 7.0/10.0)/(a[0]*a[0]*a[1]*a[1]) - (v/2.0 + 1.0/2.0)/(a[0]*a[0]*a[0]*a[0]);
    Matrix[10]=4.0/(3.0*a[1]*a[1]) - (8.0*(v/2.0 - 1.0/2.0))/(15.0*a[0]*a[0]);
    Matrix[11]=(4.0*v + 11.0)/(10.0*a[0]*a[0]*a[1]);
    Matrix[12]=0.0;
    Matrix[13]=(8.0*(v/2.0 - 1.0/2.0))/(15.0*a[0]*a[0]) + 2.0/(3.0*a[1]*a[1]);
    Matrix[14]=-(2.0*(v - 1.0))/(5.0*a[0]*a[0]*a[1]);
    Matrix[15]=4.0/(3.0*a[0]*a[0]) - (8.0*(v/2.0 - 1.0/2.0))/(15.0*a[1]*a[1]);
    Matrix[16]=v/(a[0]*a[1]);
    Matrix[17]=(v/2.0+ 1.0)/(a[0]*a[0]*a[0]) - (v/10.0 - 1.0/10.0)/(a[0]*a[1]*a[1]);
    Matrix[18]=(v/15.0 - 1.0/15.0)/(a[1]*a[1]) + 2.0/(3.0*a[0]*a[0]);
    Matrix[19]=0.0;
    Matrix[20]=(v/10.0 - 1.0/10.0)/(a[0]*a[1]*a[1]) - 1.0/(a[0]*a[0]*a[0]);
    Matrix[21]=(v/2.0 + 2.0)/(a[0]*a[0]*a[0]*a[0]) - ((7.0*v)/10.0 - 7.0/10.0)/(a[0]*a[0]*a[1]*a[1]);
    Matrix[22]=((a[0]*a[0] - 5.0*a[1]*a[1])*(v - 1.0))/(10.0*a[0]*a[0]*a[0]*a[1]*a[1]);
    Matrix[23]=(v - 11.0)/(10.0*a[0]*a[0]*a[1]);
    Matrix[24]=((7.0*v)/10.0 - 7.0/10.0)/(a[0]*a[0]*a[1]*a[1]) - (v/2.0 + 1.0/2.0)/(a[0]*a[0]*a[0]*a[0]);
    Matrix[25]=((v/10.0 - 1.0/10.0)*(a[0]*a[0]) + a[1]*a[1]/2.0)/(a[0]*a[0]*a[0]*a[1]*a[1]);
    Matrix[26]=-(v + 4.0)/(10.0*a[0]*a[0]*a[1]);
    Matrix[27]=(v/2.0 - 1.0)/(a[0]*a[0]*a[0]*a[0]) - ((7.0*v)/10.0 - 7.0/10.0)/(a[0]*a[0]*a[1]*a[1]);
    Matrix[28]=4.0/(3.0*a[1]*a[1]) - (8.0*(v/2.0 - 1.0/2.0))/(15.0*a[0]*a[0]);
    Matrix[29]=-(4.0*v + 11.0)/(10.0*a[0]*a[0]*a[1]);
    Matrix[30]=0.0;
    Matrix[31]=(v/15.0 - 1.0/15.0)/(a[0]*a[0]) + 2.0/(3.0*a[1]*a[1]);
    Matrix[32]=-(v - 11.0)/(10.0*a[0]*a[0]*a[1]);
    Matrix[33]=0.0;
    Matrix[34]=1.0/(3.0*a[1]*a[1]) - (v/15.0 - 1.0/15.0)/(a[0]*a[0]);
    Matrix[35]=(v + 4.0)/(10.0*a[0]*a[0]*a[1]);
    Matrix[36]=4.0/(3.0*a[0]*a[0]) - (8.0*(v/2.0 - 1.0/2.0))/(15.0*a[1]*a[1]);
    Matrix[37]=-v/(a[0]*a[1]);
    Matrix[38]=(v/2.0 + 1.0)/(a[0]*a[0]*a[0]) - (v/10.0 - 1.0/10.0)/(a[0]*a[1]*a[1]);
    Matrix[39]=(8.0*(v/2.0 - 1.0/2.0))/(15.0*a[1]*a[1]) + 2.0/(3.0*a[0]*a[0]);
    Matrix[40]=0.0;
    Matrix[41]=((a[0]*a[0] - 5.0*a[1]*a[1])*(v - 1.0))/(10.0*a[0]*a[0]*a[0]*a[1]*a[1]);
    Matrix[42]=1.0/(3.0*a[0]*a[0]) - (v/15.0 - 1.0/15.0)/(a[1]*a[1]);
    Matrix[43]=0.0;
    Matrix[44]=- 1.0/(2.0*a[0]*a[0]*a[0]) - (v/10.0 - 1.0/10.0)/(a[0]*a[1]*a[1]);
    Matrix[45]=(v/2.0 + 2.0)/(a[0]*a[0]*a[0]*a[0]) - ((7.0*v)/10.0 - 7.0/10.0)/(a[0]*a[0]*a[1]*a[1]);
    Matrix[46]=(v/10.0 - 1.0/10.0)/(a[0]*a[1]*a[1]) - 1.0/(a[0]*a[0]*a[0]);
    Matrix[47]=(2.0*(v - 1.0))/(5.0*a[0]*a[0]*a[1]);
    Matrix[48]=((7.0*v)/10.0 - 7.0/10.0)/(a[0]*a[0]*a[1]*a[1]) - (v/2.0 + 1.0/2.0)/(a[0]*a[0]*a[0]*a[0]);
    Matrix[49]=- 1.0/(2.0*a[0]*a[0]*a[0]) - (v/10.0 - 1.0/10.0)/(a[0]*a[1]*a[1]);
    Matrix[50]=-(v + 4.0)/(10.0*a[0]*a[0]*a[1]);
    Matrix[51]=(v/2.0 - 1.0)/(a[0]*a[0]*a[0]*a[0]) - ((7.0*v)/10.0 - 7.0/10.0)/(a[0]*a[0]*a[1]*a[1]);
    Matrix[52]=-((a[0]*a[0] - 5.0*a[1]*a[1])*(v - 1.0))/(10.0*a[0]*a[0]*a[0]*a[1]*a[1]);
    Matrix[53]=(v - 11.0)/(10.0*a[0]*a[0]*a[1]);
    Matrix[54]=((7.0*v)/10.0 - 7.0/10.0)/(a[0]*a[0]*a[1]*a[1]) - (v/2.0 + 1.0/2.0)/(a[0]*a[0]*a[0]*a[0]);
    Matrix[55]=4.0/(3.0*a[1]*a[1]) - (8.0*(v/2.0 - 1.0/2.0))/(15.0*a[0]*a[0]);
    Matrix[56]=-(4.0*v + 11.0)/(10.0*a[0]*a[0]*a[1]);
    Matrix[57]=0.0;
    Matrix[58]=(8.0*(v/2.0 - 1.0/2.0))/(15.0*a[0]*a[0]) + 2.0/(3.0*a[1]*a[1]);
    Matrix[59]=(2.0*(v - 1.0))/(5.0*a[0]*a[0]*a[1]);
    Matrix[60]=0.0;
    Matrix[61]=1.0/(3.0*a[1]*a[1]) - (v/15.0 - 1.0/15.0)/(a[0]*a[0]);
    Matrix[62]=(v + 4.0)/(10.0*a[0]*a[0]*a[1]);
    Matrix[63]=0.0;
    Matrix[64]=(v/15.0 - 1.0/15.0)/(a[0]*a[0]) + 2.0/(3.0*a[1]*a[1]);
    Matrix[65]=-(v - 11.0)/(10.0*a[0]*a[0]*a[1]);
    Matrix[66]=4.0/(3.0*a[0]*a[0]) - (8.0*(v/2.0 - 1.0/2.0))/(15.0*a[1]*a[1]);
    Matrix[67]=v/(a[0]*a[1]);
    Matrix[68]=(v/10.0 - 1.0/10.0)/(a[0]*a[1]*a[1]) - (v/2.0 + 1.0)/(a[0]*a[0]*a[0]);
    Matrix[69]=(v/15.0 - 1.0/15.0)/(a[1]*a[1]) + 2.0/(3.0*a[0]*a[0]);
    Matrix[70]=0.0;
    Matrix[71]=1.0/(a[0]*a[0]*a[0]) - (v/10.0 - 1.0/10.0)/(a[0]*a[1]*a[1]);
    Matrix[72]=1.0/(3.0*a[0]*a[0]) - (v/15.0 - 1.0/15.0)/(a[1]*a[1]);
    Matrix[73]=0.0;
    Matrix[74]=((v/10.0 - 1.0/10.0)*(a[0]*a[0]) + a[1]*a[1]/2.0)/(a[0]*a[0]*a[0]*a[1]*a[1]);
    Matrix[75]=(8.0*(v/2.0 - 1.0/2.0))/(15.0*a[1]*a[1]) + 2.0/(3.0*a[0]*a[0]);
    Matrix[76]=0.0;
    Matrix[77]=-((a[0]*a[0] - 5.0*a[1]*a[1])*(v - 1.0))/(10.0*a[0]*a[0]*a[0]*a[1]*a[1]);

    for (i=0; i<78; i++)
    {
        Matrix[i]=Matrix[i]*D0*h*a[0]*a[1];
    }

  








    
   
}

//	Calculate element stress 
void CPlate::ElementStress(double* stress, double* Displacement)
{
	CPlateMaterial* material_ = dynamic_cast<CPlateMaterial*>(ElementMaterial_);	// Pointer to material of the element

	double DX[3];	//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	double L2 = 0;	//	Square of Plate length (L^2)

	for (unsigned int i = 0; i < 3; i++)
	{
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
		L2 = L2 + DX[i]*DX[i];
	}

	double S[6];
	for (unsigned int i = 0; i < 3; i++)
	{
		S[i] = -DX[i] * material_->E / L2;
		S[i+3] = -S[i];
	}
	
	*stress = 0.0;
	for (unsigned int i = 0; i < 6; i++)
	{
		if (LocationMatrix_[i])
			*stress += S[i] * Displacement[LocationMatrix_[i]-1];
	}
}
