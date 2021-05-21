/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Jan Park
Created          : 10-05-2021
Modified         : 10-05-2021
Language/ver     : C++ in MSVS2019

Description      : [Tutorial]Differentiation_student.cpp
-------------------------------------------------------------------------------*/

#include "myNM.h"

int main(int argc, char* argv[])
{

	// PART 1
	printf("\n**************************************************");
	printf("\n|                     PART 1.                    |");
	printf("\n**************************************************\n");

	Matrix t = txt2Mat("", "Q1_vect");
	Matrix pos = txt2Mat("", "Q1_vecx");

	Matrix vel = gradient(t, pos);
	Matrix acc = gradient(t, vel);

	printf("The velocity is \n");
	printMat(vel);

	printf("The acceleration is \n");
	printMat(acc);


	// PART 2
	printf("\n**************************************************");
	printf("\n|                     PART 2.                    |");
	printf("\n**************************************************\n");

	Matrix xin = txt2Mat("", "Q2_vecxin");
	Matrix dydx = gradientFunc(myFunc, xin);

	printf("The input is \n");
	printMat(xin);

	printf("The differentiation is \n");
	printMat(dydx);

	system("pause");
	return 0;

	printf("\n**************************************************");
	printf("\n|                     PART 3.                    |");
	printf("\n**************************************************\n");

	double T_array[] = {0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4 };
	double P_array[] = {-5.87, -4.23, -2.55, -0.89, 0.67, 2.09, 3.31, 4.31, 5.06, 5.55, 5.78, 5.77, 5.52, 5.08, 4.46, 3.72, 2.88, 2.00, 1.10, 0.23, -0.59 };
	double dydx_array[] = {0};

	int m = sizeof(T_array) / sizeof(T_array[0]);

	gradient1D(T_array, P_array, dydx_array, m);

}
