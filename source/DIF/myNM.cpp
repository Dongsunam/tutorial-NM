/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Dongsu Nam
Created          : 10-05-2021
Modified         : 29-04-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.cpp
----------------------------------------------------------------*/

#include "myNM.h"
#include "math.h"


// Matrix addition
Matrix	addMat(Matrix _A, Matrix _B)
{
	if (_A.rows != _B.rows || _A.cols != _B.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'addMat' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out = createMat(_A.rows, _B.cols);
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _B.cols; j++)
			Out.at[i][j] = _A.at[i][j] + _B.at[i][j];

	return Out;
}

// Apply back-substitution
Matrix	backSub(Matrix _A, Matrix _b)
{
	Matrix Out = createMat(_b.rows, 1);

	// error check: whether _A is a triangular matrix

	//for (int i = _A.rows; i > 0; i--) {
	//	double temp = 0;
	//	for (int j = i + 1; j <= _A.cols; j++)
	//		temp += _A.at[i - 1][j - 1] * Out.at[j - 1][0];
	//	Out.at[i - 1][0] = (_b.at[i - 1][0] - temp) / _A.at[i - 1][i - 1];
	//}

	for (int i = _A.rows - 1; i >= 0; i--) {
		double temp = 0;
		for (int j = i + 1; j < _A.cols; j++)
			temp += _A.at[i][j] * Out.at[j][0];
		Out.at[i][0] = (_b.at[i][0] - temp) / _A.at[i][i];
	}

	return Out;
}

//정사각형 행렬 구조를 체크하는 함수
void chesq(int row, int col) {
	if (row != col) {
		printf("error : not square matrix\n");
		system("pause");
	}
}

//벡터와 행렬이 같은 차원인지 확인한는 함수
void chevema(int ver, int mar) {
	if (ver != mar) {
		printf("error : the dimension of vector and matrix is not same ");
		system("pause");
	}
}

//gauss elimination 함수
void gaussElim(Matrix _A, Matrix _b, Matrix _U, Matrix _d) {
	int M = _A.rows;
	int N = _A.cols;

	for (int k = 0; k < M; k++) {
		for (int i = k + 1; i < M; i++) {
			double mi = ((_A.at[i][k]) / (_A.at[k][k]));
			_A.at[i][k] = 0;

			for (int j = k + 1; j < N; j++) {
				_A.at[i][j] = _A.at[i][j] - (mi * _A.at[k][j]);
			}
			_b.at[i][0] = _b.at[i][0] - (mi * _b.at[k][0]);
		}
	}

	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			_U.at[i][j] = _A.at[i][j];
		}
		_d.at[i] = _b.at[i];

	}
}

//back subtraction
void backsub(Matrix _U, Matrix _d, Matrix _x) {
	int M = _U.rows;
	double a;

	for (int i = M - 1; i >= 0; i--) {
		a = _d.at[i][0] / _U.at[i][i];
		_x.at[i][0] = a;
		for (int h = i; h > 0; h--) {
			_d.at[h - 1][0] -= _U.at[h - 1][i] * a;
		}
	}
}

//I 행렬 만들기
void makeP(Matrix P) {
	for (int i = 0; i < P.rows; i++) {
		P.at[i][i] = 1;
	}
}

//LU decompositon
void LUdecomp(Matrix a, Matrix L, Matrix P, Matrix U) {
	double max;
	int k;
	int _k;
	double hole;
	double holep;
	double mi = 0;


	for (int i = 0; i < a.cols; i++) {
		max = a.at[i][i];
		k = i;
		_k = i;

		//Get scale partial pivot
		for (int j = i + 1; j < a.rows; j++) {
			if (fabs(a.at[j][i]) > max) {
				max = a.at[j][i];
				_k = j;
			}
		}

		//Apply row exchange
		if (k != _k) {
			for (int i = 0; i < a.rows; i++) {
				hole = a.at[k][i];
				holep = P.at[k][i];

				a.at[k][i] = a.at[_k][i];
				a.at[_k][i] = hole;

				P.at[k][i] = P.at[_k][i];
				P.at[_k][i] = holep;
			}
		}

		//row reduction
		for (int h = i + 1; h < a.rows; h++) {
			mi = a.at[h][i] / a.at[i][i];
			L.at[h][i] = mi;
			a.at[h][i] = 0;

			for (int b = i + 1; b < a.cols; b++) {
				a.at[h][b] = a.at[h][b] - (mi * a.at[i][b]);
			}
		}
	}


	//U = a_k (when k is max)
	for (int i = 0; i < a.cols; i++) {
		for (int j = 0; j < a.rows; j++) {
			U.at[i][j] = a.at[i][j];
		}
	}

	//L = L+1
	for (int i = 0; i < a.rows; i++) {
		L.at[i][i] = 1;
	}
}

//LUx=Pb=d 풀기
Matrix solveLU(Matrix L, Matrix U, Matrix P, Matrix b) {
	Matrix out = createMat(b.rows, b.cols);
	Matrix Y = createMat(b.rows, b.cols);
	Matrix D = zeros(P.rows, b.cols);

	double sum = 0;
	double c;
	double d;

	//d = Pb
	for (int j = 0; j < b.cols; j++) {
		for (int i = 0; i < P.rows; i++) {
			for (int k = 0; k < P.cols; k++) {
				sum += P.at[i][k] * b.at[k][j];
			}
			D.at[i][j] = sum;
			sum = 0;
		}
	}

	//forward sub
	for (int i = 0; i < L.rows; i++) {
		c = D.at[i][0] / L.at[i][i];
		Y.at[i][0] = c;
		for (int j = i + 1; j < L.rows; j++) {
			D.at[j][0] -= c * L.at[j][i];
		}
	}

	//backward sub
	for (int i = L.rows - 1; i >= 0; i--) {
		d = Y.at[i][0] / U.at[i][i];
		out.at[i][0] = d;
		for (int h = i; h > 0; h--) {
			Y.at[h - 1][0] -= U.at[h - 1][i] * d;
		}
	}

	return out;
}

//역행렬 계산용 no pivoting LU decomposition
void LU_nopivot(Matrix a, Matrix L, Matrix U) {
	double mi = 0;

	for (int i = 0; i < a.rows; i++) {
		for (int j = 0; j < a.cols; j++) {
			U.at[i][j] = a.at[i][j];
		}
		L.at[i][i] = 1;
	}


	for (int k = 0; k < a.rows - 1; k++) {
		for (int i = k + 1; i < a.rows; i++) {
			mi = (U.at[i][k] / U.at[k][k]);
			L.at[i][k] = mi;
			U.at[i][k] = 0;
			for (int j = k + 1; j < a.cols; j++) {
				U.at[i][j] = U.at[i][j] - mi * U.at[k][j];
			}
		}
	}
}


//역행렬 구하기 3차시도 - 성공
Matrix invinv(Matrix a) {
	Matrix inv = zeros(a.rows, a.cols);
	Matrix U = zeros(a.rows, a.cols);
	Matrix L = zeros(a.rows, a.cols);
	Matrix P = zeros(a.rows, a.cols);

	makeP(P);
	LU_nopivot(a, L, U);

	for (int k = 0; k < a.rows; k++) {
		for (int i = 0; i < a.rows; i++) {
			double sum = 0;

			for (int j = 0; j < i; j++) {
				sum = sum + L.at[i][j] * inv.at[j][k];
			}
			inv.at[i][k] = (P.at[i][k] - sum) / L.at[i][i];
		}
	}

	//back sub
	for (int k = 0; k < a.rows; k++)
	{
		for (int i = a.rows - 1; i >= 0; i--)
		{
			double sum2 = 0;

			for (int j = a.cols - 1; j > i; j--)
			{
				sum2 = sum2 + U.at[i][j] * inv.at[j][k];
			}
			inv.at[i][k] = (inv.at[i][k] - sum2) / U.at[i][i];
		}
	}

	return inv;
}

//행렬 곱셈함수
Matrix mulmat(Matrix A, Matrix B) {
	Matrix out = zeros(A.rows, B.cols);

	double sum = 0;

	for (int i = 0; i < A.rows; i++) {
		for (int j = 0; j < B.cols; j++) {
			sum = 0;
			for (int k = 0; k < B.rows; k++) {
				sum += A.at[i][k] * B.at[k][j];
			}
			out.at[i][j] = sum;
		}
	}
	return out;
}

//행렬 transpose 함수
Matrix transpose(Matrix a) {
	Matrix out = zeros(a.cols, a.rows);

	for (int i = 0; i < a.cols; i++) {
		for (int j = 0; j < a.rows; j++) {
			out.at[i][j] = a.at[j][i];
		}
	}

	return out;
}

//norm value 구하는 함수
double getnorm(Matrix a) {
	double sum = 0;

	for (int i = 0; i < a.rows; i++) {
		sum += a.at[i][0] * a.at[i][0];
	}

	double norm = sqrt(sum);
	return norm;
}

void matcopy(Matrix& a, Matrix& b) {
	for (int i = 0; i < a.rows; i++)
	{
		for (int j = 0; j < a.cols; j++)
		{
			b.at[i][j] = a.at[i][j];
		}
	}
}

//QR factorization
void QRfac(Matrix& a, Matrix& R, Matrix& Q) {

	makeP(Q);


	Matrix I = zeros(a.rows, a.cols);
	makeP(I);

	Matrix q = zeros(Q.rows, Q.cols);
	Matrix r = zeros(a.rows, a.cols);
	makeP(q);

	//R_0 = A
	for (int i = 0; i < a.rows; i++)
	{
		for (int j = 0; j < a.cols; j++)
		{
			r.at[i][j] = a.at[i][j];
		}
	}

	Matrix C = zeros(a.rows, 1);
	Matrix e = zeros(a.rows, 1);
	Matrix v = zeros(a.rows, 1);
	Matrix H = zeros(a.rows, a.cols);


	for (int i = 0; i < a.rows - 1; i++) {

		//c에 r값 대입, e값 초기화
		for (int j = 0; j < a.cols; j++) {
			C.at[j][0] = r.at[j][i];
			e.at[j][0] = 0;
			//i가 2일때, c[1][0] = 0이므로
			if (i > j) {
				C.at[j][0] = 0;
			}
		}


		//c의 row의 첫 값이 양수면 1, 음수면 -1
		if (C.at[i][0] >= 0) {
			e.at[i][0] = 1;
		}
		else {
			e.at[i][0] = -1;
		}

		double norm = getnorm(C);

		//v 구하기
		for (int j = 0; j < a.rows; j++) {
			v.at[j][0] = C.at[j][0] + (norm * e.at[j][0]);
		}

		//v transpose, v_t * v, v * v_t 구하기
		Matrix v_T = transpose(v);
		Matrix v_Tv = mulmat(v_T, v);
		Matrix vv_T = mulmat(v, v_T);

		//이때, v_T * v는 스칼라값이므로 norm을 구해준다
		double vtv = getnorm(v_Tv);

		if (vtv == 0) {
			printf("error : v_T * v is zero\n");
			system("pause");
		}

		//H 구하기
		for (int j = 0; j < vv_T.rows; j++) {
			for (int k = 0; k < vv_T.cols; k++) {
				H.at[j][k] = I.at[j][k] - ((2 / vtv) * vv_T.at[j][k]);
			}
		}

		//Q, R 업데이트
		q = mulmat(q, H);
		r = mulmat(H, r);
	}

	for (int i = 0; i < q.rows; i++) {
		for (int j = 0; j < q.cols; j++) {
			Q.at[i][j] = q.at[i][j];
			R.at[i][j] = r.at[i][j];
		}
	}
}

Matrix iteration(Matrix& a, Matrix& Q, Matrix& R) {
	Matrix Evalue;
	for (int i = 0; i < 414; i++) {
		QRfac(a, R, Q);
		Evalue = mulmat(R, Q);
	}
	return Evalue;
}

//A_i = U 이므로 람다값을 가지는 matrix U 생성 함수
Matrix makeU(Matrix a) {
	Matrix Eigenvalue = zeros(a.rows, a.cols);
	makeP(Eigenvalue);
	for (int i = 0; i < a.rows; i++) {
		Eigenvalue.at[i][i] = a.at[i][i];
	}
	return Eigenvalue;
}

//(A-람다*I)v = 0, 이때 v를 구하기 위해 A-람다를 해준다
Matrix minuslamda(Matrix a, Matrix lam) {
	Matrix out = zeros(a.rows, a.cols);

	for (int i = 0; i < a.rows; i++) {
		for (int j = 0; j < a.cols; j++) {
			out.at[i][j] = a.at[i][j];
		}
	}

	for (int i = 0; i < a.rows; i++) {
		out.at[i][i] -= lam.at[i][i];
	}

	return out;
}

Matrix Eigenvector(Matrix a) {
	Matrix zero = zeros(a.rows, 1);
	Matrix out = zeros(a.rows, 1);
	Matrix U = zeros(a.rows, a.cols);
	Matrix d = zeros(a.rows, 1);
	int row = a.rows;
	int col = a.cols;


	gaussElim(a, zero, U, d);
	backsub(U, d, out);

	return out;
}


double cond(Matrix& a) {
	Matrix tran = transpose(a);
	Matrix A = mulmat(tran, a);


	Matrix Q = zeros(a.rows, a.cols);
	makeP(Q);

	Matrix R = zeros(a.rows, a.cols);
	for (int i = 0; i < a.rows; i++) {
		for (int j = 0; j < a.cols; j++) {
			R.at[i][j] = A.at[i][j];
		}
	}

	QRfac(A, R, Q);

	Matrix lam_A = iteration(A, Q, R);

	double max_1;
	double max_2;
	double min_1;
	double min_2;

	for (int j = 0; j < 1; j++) {
		max_1 = (sqrt(lam_A.at[j][j]));
		min_1 = (sqrt(lam_A.at[j][j]));

		for (int i = 1; i < lam_A.rows; i++) {
			max_2 = (sqrt(lam_A.at[i][i]));
			min_2 = (sqrt(lam_A.at[i][i]));

			if (max_2 > max_1) {
				max_1 = max_2;
			}

			if (min_1 > min_2) {
				min_1 = min_2;
			}
		}
	}

	double out = max_1 / min_1;
	return out;
}

//Assignment#5

//part 1
Matrix linearFit(Matrix a, Matrix b) {
	Matrix z = zeros(2, 1);
	int n = a.cols;
	double sumx = 0;
	double sumx_2 = 0;
	double sumy = 0;
	double sumxy = 0;

	for (int i = 0; i < a.cols; i++) {
		sumx = sumx + a.at[0][i];
		sumx_2 = sumx_2 + a.at[0][i] * a.at[0][i];
		sumy = sumy + b.at[0][i];
		sumxy = sumxy + a.at[0][i] * b.at[0][i];
	}

	double z_1 = (n * sumxy - sumx * sumy) / (n * sumx_2 - sumx * sumx);
	double z_2 = (sumy - z_1 * sumx) / n;

	z.at[0][0] = z_2;
	z.at[0][1] = z_1;

	return z;
}

double funcpre(Matrix z, double t) {
	double x = z.at[0][0] + z.at[0][1] * t;
	return x;
}

//part 2
Matrix linearInterp(Matrix a, Matrix b, Matrix x) {
	Matrix y = zeros(x.rows, x.cols);

	int n = x.cols;
	int m = a.cols;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m - 1; j++) {
			if (x.at[0][i] >= a.at[0][j] && x.at[0][i] < a.at[0][j + 1]) {
				y.at[0][i] = b.at[0][j] * ((x.at[0][i] - a.at[0][j + 1]) / (a.at[0][j] - a.at[0][j + 1])) + b.at[0][j + 1] * ((x.at[0][i] - a.at[0][j]) / (a.at[0][j + 1] - a.at[0][j]));
			}
		}
		//만약 t 값이 100일경우, 위 loop에 걸리지 않으므로
		if (y.at[0][i] == 0) {
			y.at[0][i] = b.at[0][m - 1] * ((x.at[0][i] - a.at[0][m]) / (a.at[0][m - 1] - a.at[0][m])) + b.at[0][m] * ((x.at[0][i] - a.at[0][m - 1]) / (a.at[0][m] - a.at[0][m - 1]));
		}
	}

	return y;
}


//#6
//three point forward
Matrix gradient(Matrix _x, Matrix _y) {
	int n = _x.rows;
	Matrix dxdy = zeros(n, 1);
	chevema(_x.rows, _y.rows);


	if (n > 2) {
		dxdy.at[0][0] = (-3 * _y.at[0][0] + 4 * _y.at[1][0] - _y.at[2][0]) / (_x.at[2][0] - _x.at[0][0]);

		for (int i = 1; i < n - 1; i++) {
			dxdy.at[i][0] = _y.at[i + 1][0] - _y.at[i - 1][0] / _x.at[i + 1][0] - _x.at[i - 1][0];
		}

		dxdy.at[n - 1][0] = (3 * _y.at[n - 1][0] - 4 * _y.at[n - 2][0] + _y.at[n - 3][0]) / (_x.at[n - 1][0] - _x.at[n - 3][0]);
	}

	else if (n==2) {
		dxdy.at[0][0] = _y.at[1][0] - _y.at[0][0] / _x.at[1][0] - _x.at[0][0];
		dxdy.at[0][1] = _y.at[1][0] - _y.at[0][0] / _x.at[1][0] - _x.at[0][0];
	}

	return dxdy;
}

// Define a function that defines the target equation.
double myFunc(const double x) {
	return  x * x * x;
}

// Return the dy/dx results for the target equation. (truncation error: O(h^2))
// Move this function to myNM.h and myNM.cpp
Matrix	gradientFunc(double func(const double x), Matrix xin) {
	int n = xin.rows;
	Matrix dxdy = zeros(n, 1);

	if (n > 2) {
		dxdy.at[0][0] = (-3 * func(xin.at[0][0]) + 4 * func(xin.at[1][0]) - func(xin.at[2][0])) / (xin.at[2][0] - xin.at[0][0]);

		for (int i = 1; i < n - 1; i++) {
			dxdy.at[i][0] = func(xin.at[i + 1][0]) - func(xin.at[i - 1][0]) / xin.at[i + 1][0] - xin.at[i - 1][0];
		}

		dxdy.at[n - 1][0] = (3 * func(xin.at[n - 1][0]) - 4 * func(xin.at[n - 2][0]) + func(xin.at[n - 3][0])) / (xin.at[n - 1][0] - xin.at[n - 3][0]);
	}

	else if (n == 2){
		dxdy.at[0][0] = func(xin.at[1][0]) - func(xin.at[0][0]) / xin.at[1][0] - xin.at[0][0];
		dxdy.at[0][1] = func(xin.at[1][0]) - func(xin.at[0][0]) / xin.at[1][0] - xin.at[0][0];
	}

	return dxdy;
}

void gradient1D(double x[], double y[], double dydx[], int m) {
	Matrix xx = arr2Mat(x, m, 1);
	Matrix yy = arr2Mat(y, m, 1);

	Matrix output = gradient(xx, yy);

	for (int i = 0; i++; i < m) {
		dydx[i] = output.at[0][i];
	}
}

