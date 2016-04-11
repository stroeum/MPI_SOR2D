/*
 *  SorSolution.h
 *  Created by Jérémy Riousset on 11/19/07.
 */

#ifndef SORSOLUTION_H
#define SORSOLUTION_H

#include "Sources.h"
#include "MPI_Input.h"

/**************************************************************************************/
enum SourceType {ChargeDistribution, PotentialDistribution};
enum PointsColor {black, red};
// Allowed type of sources to use SOR solver
/**************************************************************************************/

/**************************************************************************************/
class SorSolution
{
private:
	SourceType	type;																	// Type of source
	double		epsilon;																// Maximum tolerable error
	int			MaxStep;																// Maximum number of iterations
	CMatrix1D	a,b;																	// Laplace's Equation coefficients (cont.)
	double		c[2], d[2], e[2];														// Laplace's Equation coefficients (cont.)
	CMatrix2D	f;																		// Laplace's Equation coefficients (cont.)
	double		ErrDen;																	// Normalisation of the error
//	int			is, ie, ir;

public:	
	SorSolution(){};																	// Default constructor
	SorSolution(double eepsilon,int MMaxStep, ResGrid dd,SizeGrid NN, Charge& CC);	
	// Constructor surcharge
	~SorSolution(){};																	// Destructor
	void init(double eepsilon,int MMaxStep, ResGrid dd,SizeGrid NN, Charge& CC);
	void Update_GhostRows(CMatrix2D& local_MM, SizeGrid local_NN,PointsColor CC);		// Exchange rows
	void Solve(ResGrid dd, SizeGrid NN, CMatrix2D& pphi);
	
	CMatrix2D Gather( CMatrix2D local_MM, SizeGrid local_NN)
	{
		SizeGrid NN;
		NN.r = (local_NN.r-2)*MPI_Var::dim_sizes+2;
		NN.z = local_NN.z;
		
		CMatrix2D MM;
		if( (NN.r-2)%MPI_Var::n_processes == 0 )
		{
			int			source;
			int			dest;
			int			tag		=0;
			MPI_Status	status;
			
			MM.init(NN.r,NN.z);
			
			//COPY 1st line of local_MM to 1st line of MM on root process//
			if(MPI_Var::r_rank == MPI_Var::root)
			{
				memcpy(&MM(0,0), &local_MM(0,0), local_NN.z*sizeof(double));
			};
			
			//SEND last line of local_MM to last line of MM on last process//
			if(MPI_Var::r_rank == MPI_Var::dim_sizes-1)
			{
				dest = MPI_Var::root;
				MPI_Send(	&local_MM(local_NN.r-1,0), NN.z,	MPI_DOUBLE, dest, tag, MPI_Var::r_comm);
			}
			if(MPI_Var::r_rank == MPI_Var::root)
			{
				source = MPI_Var::dim_sizes-1;
				MPI_Recv(	&MM(NN.r-1,0),		local_NN.z,		MPI_DOUBLE, source, tag, MPI_Var::r_comm, &status);
			}
			
			//GATHER all intermediate lines//
			MPI_Gather(&local_MM(1,0), (local_NN.r-2)*local_NN.z, MPI_DOUBLE, &MM(1,0), (local_NN.r-2)*local_NN.z, MPI_DOUBLE, 0, MPI_Var::r_comm);
		}
		else
		{
			cout<<"\n***number of processes incompatible with discretization***"<<endl;
			exit(101);
		}
		return MM;
	};
};

/*
private:
	double epsilon;																		// Maximum tolerable error
	int MaxStep;																		// Maximum number of iterations
	double a,b,c,d,e,f,g;																// Laplace's Equation coefficients
	CMatrix3D h;																		// Laplace's Equation coefficients (cont.)
	SourceType type;																	// Type of source
	double ErrDen;																		// Normalisation of the error
	
public:	
	SorSolution(){};																	// Default constructor
	SorSolution(CMatrix3D& pphi, double eepsilon,int MMaxStep, ResGrid dd,SizeGrid NN, Charge& CC,		const	CMatrix3D& UUn);	
																						// Constructor surcharge
	SorSolution(CMatrix3D& pphi, double eepsilon,int MMaxStep, ResGrid dd,SizeGrid NN, Potential& PP,			CMatrix3D& UUn);	
																						// Constructor surcharge
	~SorSolution(){};																	// Destructor
	void init(CMatrix3D& pphi, double eepsilon,int MMaxStep, ResGrid dd,SizeGrid NN, Charge& CC,	const	CMatrix3D& UUn);
	void init(CMatrix3D& pphi, double eepsilon,int MMaxStep, ResGrid dd,SizeGrid NN, Potential& PP,			CMatrix3D& UUn);
	void Update_GhostPlanes(CMatrix3D& MM, SizeGrid local_N, PointsColor CC);			// Exchange rows
	void Solve(ResGrid dd, SizeGrid NN, const CMatrix3D& UUn, CMatrix3D& pphi);			// Solve solution using SOR method
*/
/**************************************************************************************/

#endif // SORSOLUTION_H
