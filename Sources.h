/* Sources.h */
#ifndef SOURCES_H
#define SOURCES_H

#include <list>
#include <fstream>
#include "ResGrid.h"
#include "Matrix.h"
#include "Constants.h"
using namespace std;

/**************************************************************************************/
class Charge
{
protected:
	
public:
	/* Charge Type */
	string Type;
	/* Charge Content */
	double Q;
	/* Charge Position */
	double Zq;
	/* Charge Geometrical Parameters */
	double Rq1,Rq2;
	CMatrix2D Un;
	CMatrix2D rho;
	
	Charge();											// Default constructor
	Charge(ResGrid dd, SizeGrid NN);					// Initialize Un and rho size
	Charge(double QQ, double ZZq, double RRq1, double RRq2);
														// Initialize geometrical parameters
	bool		init(double QQ, double ZZq, double RRq1, double RRq2);
														// Initiate geometrical parameters for geometrical charges
	bool		reset(ResGrid dd, SizeGrid NN);			// Reset geometrical parameters
	bool		sphere(double QQ, double ZZq, double RR, ResGrid dd, SizeGrid NN);
														// Distribute charge assuming spherical geometry
	
	// Calculate analytical solution for a SPHERE of radius Rq1 and charge Q at Zq
	CMatrix1D	MonopoleAnalyticalSolution(	ResGrid dd, SizeGrid NN);
														// ... in free space
	CMatrix1D	DipoleAnalyticalSolution(	ResGrid dd, SizeGrid NN);
														// ... above a PEC ground plane
	CMatrix1D	MultipoleAnalyticalSolution(ResGrid dd, SizeGrid NN);
														// ... between two PEC ground planes
	CMatrix1D	getParams();							// Get Charge parameters Q, Xq,Yq,Zq, Rq1,Rq2,Rq3
	string		getType();								// Get Charge type
	void		fwrite(char * title);					// Write Charge information to a file
	bool		rotate(double a, double b, double c, double u, double v, double w, double theta, ResGrid dd, SizeGrid NN);
														// Rotate Charge distribution
		
	/*Note: All the surcharged operators below does not consider the position of the 
			centers of the charges, it is basically a shortcut to sum the total charge as
			well as the Status of the lattices. After the sum, the center is affected to
			be 0,0,0 but does not really matter since it won't be used. Our objective while
			defining those operator was to ease the use of several charges.*/
	Charge		operator+=(const Charge&);
	Charge&		operator=(const Charge&);
	Charge		operator+(const Charge&) const;
	friend ostream & operator<< (ostream &, const Charge &);
	~Charge();
};
/**************************************************************************************/

/**************************************************************************************/
class Potential
{
private:
	bool	EquiPotential;							// Type of source
	double	Vo;										// Potential value
	double	Zc;										// Center Position
	double	L,H;									// Dimensions of the electrode,
													//  * L: Length/Radius
													//  * H: Heigth
	CMatrix3D rho;									// Potential Distribution
	CMatrix3D Un;									// Status of the lattice
													// 0 is modifiable point, 1 if not.
													// Insofar as no point other than the
													// boundary points, Un is no required when
													// source is a charge distribution.
public:
	Potential();									// Defaut constructor
	~Potential();
};
/**************************************************************************************/

/**************************************************************************************/
typedef list<Charge>   ListCharge;
void write(Charge&,	char *);
void write(Potential&,	char *);
/**************************************************************************************/
#endif SOURCES_H