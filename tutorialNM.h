/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Dongsu Nam
Created          : 18-05-2021
Modified         : 29-04-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.h
----------------------------------------------------------------*/

#ifndef		_MY_NM_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NM_H

#include "myMatrix.h"

// Matrix addition
extern	Matrix	addMat(Matrix _A, Matrix _B);

// Apply back-substitution
extern	Matrix	backSub(Matrix _A, Matrix _b);

//square 구조 확인 함수
extern void chesq(int row, int col);

//vetor와 matrix가 같은 차원인지 확인하는 함수
extern void chevema(int ver, int mar);

extern void gaussElim(Matrix _A, Matrix _b, Matrix _U, Matrix _d);

extern void backsub(Matrix _U, Matrix _d, Matrix _x);

extern void makeP(Matrix P);

extern void LUdecomp(Matrix a, Matrix L, Matrix P, Matrix U);

extern Matrix solveLU(Matrix L, Matrix U, Matrix P, Matrix b);


extern Matrix transpose(Matrix a);

extern void LU_nopivot(Matrix a, Matrix L, Matrix U);

extern Matrix invinv(Matrix a);

extern double getnorm(Matrix a);

extern void QRfac(Matrix& a, Matrix& R, Matrix& Q);

extern Matrix mulmat(Matrix A, Matrix B);

extern Matrix iteration(Matrix& a, Matrix& Q, Matrix& R);

extern Matrix makeU(Matrix a);

extern Matrix minuslamda(Matrix a, Matrix lam);

extern Matrix Eigenvector(Matrix a);

extern double cond(Matrix& a);

extern Matrix linearFit(Matrix a, Matrix b);

extern double funcpre(Matrix z, double t);

extern Matrix linearInterp(Matrix a, Matrix b, Matrix x);

extern Matrix gradient(Matrix _x, Matrix _y);

extern double myFunc(const double x);

extern Matrix gradientFunc(double func(const double x), Matrix xin);

extern void gradient1D(double x[], double y[], double dydx[], int m);

#endif
#pragma once
