/*
 *  SOR.h
 *  Created by Jeremy Riousset on 5/22/09.
 */

#ifndef SOR_H
#define SOR_H
#include "Input.h"
#include "Utils.h"

double** SOR(double **rrho, double **pphi, double *rr, double *zz, int Nr, int Nz, int argc, char** argv);

#endif // SOR_H
