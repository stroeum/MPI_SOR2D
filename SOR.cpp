/*
 *  SOR.cpp
 *  Created by Jeremy Riousset on 5/22/09.
 */

#include "SOR.h"


double** SOR(double **rrho, double **pphi, double *rr, double *zz, int Nr, int Nz, int argc, char** argv)
{
	/*DEFINE SIMULATION VARIABLES**************************************************/
	Var::N.r =  Nr; // WE MUST HAVE N.x and N.y even (else it unbalance the load on the processes)
	Var::N.z =  Nz;
	Var::L.init(rr[Nr-1],zz[Nz-1]);
	Var::d.init(Var::L,Var::N);
	/******************************************************************************/
	
	/*START MPI********************************************************************/
	MPI_Init(&argc, &argv);
	/******************************************************************************/
	
	//CREATE COMMUNICATORS AND TOPOLOGIES//
	CreateComm();
	CreateGridComm();
	CreateCartComm();
	
	//CREATE CHARGE DENSITY ON PROCESS ROOT//
	if (MPI_Var::world_rank == MPI_Var::root)
	{
		Var::phi.init(Var::N.r,Var::N.z);
		Var::C.rho.init(Var::N.r,Var::N.z);
		Var::C.Un.init(Var::N.r,Var::N.z);
		
		memcpy(&Var::phi.pMatrix2d[0][0],	&pphi[0][0], sizeof(double));
		memcpy(&Var::C.rho.pMatrix2d[0][0], &rrho[0][0], Nr*Nz*sizeof(double));
		for(int ii=0 ; ii<Nr-1 ; ii++) for(int kk=0 ; kk<Nz-1 ; kk++) if (Var::phi[ii][kk]!=0) Var::C.Un[ii][kk] = 1;
	}	

	//CREATE LOCAL DIMENSIONS AND GHOST VECTORS//
	InitLocalDimensions();
	CreateGhostVector();
	MPI_Var::is = (MPI_Var::r_rank != 0);
	MPI_Var::ie = (MPI_Var::r_rank == MPI_Var::n_processes-1)*(Var::local_N.r-1) +(MPI_Var::r_rank != MPI_Var::n_processes-1)*(Var::local_N.r-2);

	//SCATTER MATRICES//
	Var::local_C    = Scatter(Var::C);												// Scatter Charge
	Var::local_phi  = Scatter(Var::phi);											// Scatter Potential

	//SOLVE POISSON'S EQUATION//
	SorSolution		SOR;
	SOR.init(Var::epsilon, Var::MaxStep, Var::d, Var::local_N, Var::local_C);		
	SOR.Solve(Var::d, Var::local_N, Var::local_phi);
	
	//GATHER SCATTERED MATRICES//
	Var::phi = Gather(Var::local_phi);
	
	/*STOP MPI*********************************************************************/
	FreeComm();
	MPI_Finalize();
	/******************************************************************************/
	return Var::phi.pMatrix2d;
}