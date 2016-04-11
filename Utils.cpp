/*
 *  Utils.cpp
 *  Created by Jérémy Riousset on 11/19/07.
 */

#include "Utils.h"

void CreateComm(void)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &MPI_Var::world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &MPI_Var::n_processes);
};

void CreateGridComm(void)
{
	//	int processes_per_dimensions = (int)pow(MPI_Var::n_processes,1/3.);
	MPI_Var::dim_sizes		= MPI_Var::n_processes;										//processes_per_dimensions;
	MPI_Var::wrap_around	= 0;														// The domain is NOT periodic in z
	int tmp_grid_rank;
	
	MPI_Cart_create(MPI_COMM_WORLD, MPI_Var::n_dims, &MPI_Var::dim_sizes, &MPI_Var::wrap_around, MPI_Var::reorder, &MPI_Var::grid_comm);
	MPI_Comm_rank(MPI_Var::grid_comm, &tmp_grid_rank);
	MPI_Cart_coords(MPI_Var::grid_comm, tmp_grid_rank, MPI_Var::max_dims, &MPI_Var::coordinates);
	MPI_Cart_rank(MPI_Var::grid_comm, &MPI_Var::coordinates, &MPI_Var::grid_rank);
}

void CreateCartComm(void)
{
	MPI_Var::free_coords = 1;
    MPI_Cart_sub(	MPI_Var::grid_comm, &MPI_Var::free_coords, &MPI_Var::r_comm);
	MPI_Comm_rank(	MPI_Var::r_comm,	&MPI_Var::r_rank);
	MPI_Cart_shift(	MPI_Var::r_comm,	MPI_Var::direction, MPI_Var::displacement, &MPI_Var::prev_r_rank, &MPI_Var::next_r_rank);
};

void InitLocalDimensions(void)
{
	Var::local_N.r = (Var::N.r-2)/MPI_Var::dim_sizes+2;
	Var::local_N.z = Var::N.z;
}

void FreeComm(void)
{
	MPI_Comm_free( &MPI_Var::grid_comm);
	MPI_Comm_free( &MPI_Var::r_comm);
};

CMatrix2D Scatter( CMatrix2D MM)
{
	CMatrix2D local_MM;
	if( (Var::N.r-2)%MPI_Var::n_processes == 0 )
	{
		int			source;
		int			dest;
		int			tag		=0;
		MPI_Status	status;
		
		local_MM.init(Var::local_N.r,Var::local_N.z);
		
		//COPY 1st line of rho to 1st line of local_rho on root process//
		if(MPI_Var::r_rank == MPI_Var::root)
		{
			memcpy(&local_MM(0,0), &MM(0,0), Var::local_N.z*sizeof(double));
		};
		
		//SEND last line of rho to last line of local_rho on last process//
		if(MPI_Var::r_rank == MPI_Var::root)
		{
			dest = MPI_Var::dim_sizes-1;
			MPI_Send(	&MM(Var::N.r-1,0), Var::N.z,	MPI_DOUBLE, dest, tag, MPI_Var::r_comm);
		}
		if(MPI_Var::r_rank == MPI_Var::dim_sizes-1)
		{
			source = MPI_Var::root;
			MPI_Recv(	&local_MM(Var::local_N.r-1,0), Var::local_N.z,	MPI_DOUBLE, source, tag, MPI_Var::r_comm, &status);
		}
		
		//SCATTER all intermediate lines//
		MPI_Scatter(&MM(1,0), (Var::local_N.r-2)*Var::local_N.z, MPI_DOUBLE, &local_MM(1,0), (Var::local_N.r-2)*Var::local_N.z, MPI_DOUBLE, MPI_Var::root, MPI_Var::r_comm);
	}
	else
	{
		cout<<"\n***number of processes incompatible with discretization***"<<endl;
		exit(10);
	}
	return local_MM;
};


Charge Scatter(Charge CC)
{
	// WE DO NOT CARRY UNECESSARY INFORMATION ABOUT THE CHARGE TYPE TO SAVE COMMUNICATION (BCAST) TIME //
	Charge local_CC;
	local_CC.rho  = Scatter(CC.rho);
	local_CC.Un   = Scatter(CC.Un);
	return local_CC;
};

CMatrix2D Gather( CMatrix2D local_MM)
{
	CMatrix2D MM;
	if( (Var::N.r-2)%MPI_Var::n_processes == 0 )
	{
		int			source;
		int			dest;
		int			tag		=0;
		MPI_Status	status;
		
		MM.init(Var::N.r,Var::N.z);
		
		//COPY 1st line of local_MM to 1st line of MM on root process//
		if(MPI_Var::r_rank == MPI_Var::root)
			memcpy(&MM(0,0), &local_MM(0,0), Var::local_N.z*sizeof(double));
		MPI_Bcast(&MM(0,0), Var::local_N.z, MPI_DOUBLE, MPI_Var::root, MPI_Var::r_comm);  
		
		//SEND last line of local_MM to last line of MM on last process//
		if(MPI_Var::r_rank == MPI_Var::dim_sizes-1)
		{
			dest = MPI_Var::root;
			MPI_Send(	&local_MM(Var::local_N.r-1,0), Var::N.z,	MPI_DOUBLE, dest, tag, MPI_Var::r_comm);
		}
		if(MPI_Var::r_rank == MPI_Var::root)
		{
			source = MPI_Var::dim_sizes-1;
			MPI_Recv(	&MM(Var::N.r-1,0),		Var::local_N.z,		MPI_DOUBLE, source, tag, MPI_Var::r_comm, &status);
		}
		MPI_Bcast(&MM(Var::N.r-1,0), Var::local_N.z, MPI_DOUBLE, MPI_Var::root, MPI_Var::r_comm);  
		
		//GATHER all intermediate lines//
		MPI_Allgather(&local_MM(1,0), (Var::local_N.r-2)*Var::local_N.z, MPI_DOUBLE, &MM(1,0), (Var::local_N.r-2)*Var::local_N.z, MPI_DOUBLE, MPI_Var::r_comm);
	}
	else
	{
		cout<<"\n***number of processes incompatible with discretization***"<<endl;
		exit(101);
	}
	return MM;
};

void CreateGhostVector(void)
{
	// Var::local_N.z MUST BE EVEN //
	MPI_Type_vector(Var::local_N.z/2, 1, 2, MPI_DOUBLE, &MPI_Var::ghost_vector);
	MPI_Type_commit(&MPI_Var::ghost_vector);	
};