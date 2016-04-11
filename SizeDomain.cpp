/*
 *  SizeDomain.cpp
 *  Created by Jérémy Riousset on 10/25/07.
 */

#include "SizeDomain.h"

/**************************************************************************************/
SizeDomain::SizeDomain(){r = 1; z = 1;};
SizeDomain::SizeDomain(double rr=1,double zz=1)
{SizeDomain::init(rr,zz);};

bool SizeDomain::init(double rr, double zz)
{
	r = rr; z = zz;
	return true;
};

SizeDomain::~SizeDomain(){};
/**************************************************************************************/
