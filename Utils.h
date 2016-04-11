/*
 *  Utils.h
 *  Created by Jérémy Riousset on 11/19/07.
 */

#ifndef UTILS_H
#define UTILS_H
#include "Input.h"
#include "MPI_Input.h"

void		CreateComm(void);
void		CreateGridComm(void);
void		CreateCartComm(void);
void		InitLocalDimensions(void);

void		FreeComm(void);

CMatrix2D	Scatter(CMatrix2D	MM);
Charge		Scatter(Charge		CC);
CMatrix2D	Gather(	CMatrix2D	local_MM);
void		CreateGhostVector(void);

#endif // UTILS_H
