/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "LoadCaseData.h"

#include <iomanip>
#include <iostream>

using namespace std;

CLoadCaseData :: ~CLoadCaseData()
{
	delete [] node;
	delete [] dof;
	delete [] load;
}

void CLoadCaseData :: Allocate(unsigned int num)
{
	nloads = num;
	node = new unsigned int[nloads];
	dof = new unsigned int[nloads];
	load = new double[nloads];
}; 

//	Read load case data from stream Input
bool CLoadCaseData :: Read(ifstream& Input)
{
//	Load case number (LL) and number of concentrated loads in this load case(NL)
	
	unsigned int NL;

	Input >> NL;
//如果NL=0，则施加均布载荷，否则施加集中载荷，默认垂直于板方向，沿z方向为正
	if (NL == 0)
	{
		p= new double;
		Input >> *p;
		return true;
	}
	else
	{
	Allocate(NL);

	for (unsigned int i = 0; i < NL; i++)
		Input >> node[i] >> dof[i] >> load[i];

	return true;
	}
}

//	Write load case data to stream
void CLoadCaseData::Write(COutputter& output)
{
	for (unsigned int i = 0; i < nloads; i++)
		output << setw(7) << node[i] << setw(13) << dof[i]  << setw(19) << load[i] << endl;
}
