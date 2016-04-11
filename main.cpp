/*
 *  main.cpp
 *  Created by Jeremy Riousset on 11/19/07.
 */

#include "SOR.h"

int main(int argc, char** argv)
{
	//INPUT PARAMETERS//
	int Nr =  402;
	int Nz =  402;
	
	double Lr = 20;
	double Lz = 20;
	double dr = (double)(Lr/(Nr-1));
	double dz = (double)(Lz/(Nz-1));
	
	CMatrix2D phi(Nr,Nz);
	CMatrix2D rho(Nr,Nz);
	CMatrix1D r(Nr);
	CMatrix1D z(Nz);
	
	//INPUT PARAMETERS (CONT.)//
	double Zq = Lz/2;
	double Rq1= 5;
	
	double Q  = 100;
	
	double volume(0);
	for(int ii=0 ; ii<Nr-1 ; ii++)
		for(int kk=0 ; kk<Nz-1 ; kk++)
			if ( (ii*dr)*(ii*dr)+(kk*dz-Zq)*(kk*dz-Zq) <= Rq1*Rq1 )
				volume +=  (ii==1) * 2*M_PI*(dr/2)*(dr/2)*dz + (ii!=1)*2*M_PI*ii*dr*dr*dz;
	
	double rho0 = Q / volume;
	
	for(int ii=0 ; ii<Nr ; ii++) for(int kk=0 ; kk<Nz ; kk++)
		if(sqrt(pow(ii*dr,2)+pow(kk*dz-Zq,2))<=Rq1)
			rho[ii][kk] = rho0;

	for(int ii=0 ; ii<Nr ; ii++) 
		r[ii] = ii*dr;
		
	for(int kk=0 ; kk<Nz ; kk++)
		z[kk] = kk*dz;
	
	//CALCULATE SOLUTION//
	double **PHI;
	PHI = SOR(rho.pMatrix2d, phi.pMatrix2d, r.pElems, z.pElems, Nr, Nz, argc, argv);

	//STORE SOLUTION//
	if (MPI_Var::world_rank == MPI_Var::root)
	{	
		if (MPI_Var::world_rank == MPI_Var::root)
		{
			CMatrix1D phi1D(Var::N.z);
			for(int kk=0 ; kk<Var::N.z ; kk++)	phi1D[kk] = PHI[0][kk];
			phi1D.fwrite("results/phiNum.dat");
		}
	}
	return 0;
}