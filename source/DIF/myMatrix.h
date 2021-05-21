/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Dongsu Nam
Created          : 10-05-2021
Modified         : 29-04-2021
Language/ver     : C++ in MSVS2019

Description      : myMatrix.h
----------------------------------------------------------------*/

#ifndef		_MY_MATRIX_H		
#define		_MY_MATRIX_H

#include <iostream>
#include <string>
#include <fstream>

typedef struct {
	double** at;
	int rows, cols;
}Matrix;

//using namespace std;

// Create Matrix with specified size
extern	Matrix	createMat(int _rows, int _cols);

// Free a memory allocated matrix
extern	void	freeMat(Matrix _A);

// Create a matrix from a text file
extern	Matrix	txt2Mat(std::string _filePath, std::string _fileName);

//// Print matrix
extern	void	printMat(Matrix _A);

extern Matrix	arr2Mat(double* _1Darray, int _rows, int _cols);



/// It is recommended to create the following functions.

// initialization of Matrix elements
extern	void	initMat(Matrix _A, double _val);

// Create matrix of all zeros
extern	Matrix	zeros(int _rows, int _cols);

// Create matrix of all ones
//extern	Matrix	ones(int _rows, int _cols);

// Create identity 
//extern	Matrix	eye(int _rows, int _cols);

// Create Transpose matrix
//extern	Matrix	transpose(Matrix _A);

// Copy matrix
//extern	Matrix	copyMat(Matrix _A);

// Copy matrix Elements from A to B
//extern	void	copyVal(Matrix _A, Matrix _B);


#endif