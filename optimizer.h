#pragma once
#include"fvm/eigen-3.4.0/eigen-3.4.0/Eigen/Dense"
#include<stdint.h>
#include<iostream>
#include<fstream>
#include<vector>
#include<math.h>
#include<ctime>
#include<omp.h>
#include<string>

using namespace std;
using namespace Eigen;

//Matrix<double, 2, 2> operator+(const double I, Matrix<double, 2, 2>&& A);
//Matrix<double, 2, 2> dot(Matrix<double, 2, 2>& A, Matrix<double, 2, 2>& B);

class Adam {
public://Dcenter, Ncenter, Inset, t_R, t_N;
	Adam(string zfin, string alp, string bet, string ssm, string para = "", string beta11 = "",
		string beta22 = "", string Dcenter1 = "", string Ncenter1 = "", string Inset1 = "", string t_R1 = "", string t_N1 = "") :
		Z_final(zfin), Alpha_smooth(alp), Beta_smooth(bet), S_smooth(ssm), Para(para), Beta1(beta11), Beta2(beta22), Dcenter(Dcenter1), Ncenter(Ncenter1), Inset(Inset1), t_R(t_R1), t_N(t_N1)
	{
		check();
		cout << "Thread number = " << P_ara << endl;
		cout << "beta1 = " << beta1 << endl;
		cout << "beta2 = " << beta2 << endl;
		cout << "D_center = " << D_center << endl;
		cout << "N_center = " << N_center << endl;
		cout << "inset = " << inset << endl;
		cout << "r = " << r << endl;
		cout << "n = " << n << endl;
	};
	void check();
	void error(string);
	bool ReadDate();
	void All();
	
	void partial_derivative(Matrix<double, 81, 81>&Z_S, int Z_size1, int Z_size2, double space_x, double space_y);
	vector<double> min(int D_pos_num, int N_pos_num,vector<int>(&D_pos)[2],vector<int>(&N_pos)[2],Matrix<double, 81, 81>&, Matrix<double, 81, 81>&, Matrix<double, 81, 81>&);
	template <typename T,int rows, int cols>
	double evaluation_func(Matrix<T, rows, cols>& Alpha, Matrix<T, rows, cols>& Beta, Matrix<T, rows, cols>& Ssmooth, Matrix<T, rows, cols>& z_x, Matrix<T, rows, cols>& z_y, Matrix<T, rows, cols>& z_xx,
		Matrix<T, rows, cols>& z_yy, Matrix<T, rows, cols>& z_xy, double min, double max);
	template <typename T, int rows, int cols>
	double evaluation_func0(Matrix<T, rows, cols>& Alpha, Matrix<T, rows, cols>& Beta, Matrix<T, rows, cols>& Ssmooth, Matrix<T, rows, cols>& z_x, Matrix<T, rows, cols>& z_y, Matrix<T, rows, cols>& z_xx, 
		Matrix<T, rows, cols>& z_yy, Matrix<T, rows, cols>& z_xy);
	void calculate_A(double con1, double con2, vector<int>(&D_pos)[2], vector<int>(&N_pos)[2], int min_pos1, int max_pos2);
	void calculate_A22();
	void calculate_A23();
	void calculate_A24();
	void calculate_A32();
	void calculate_A42();
	void calculate_A33(double con1,double con2, vector<int>(&D_pos)[2], vector<int>(&N_pos)[2], int min_pos1, int max_pos2);
	void calculate_A34();
	void calculate_A43();
	void calculate_A44();
	template <typename T, int rows, int cols>
	void Zeros(Matrix<T, rows, cols>& matrix);

	vector<double> Adam_solver(int D_pos_num, int N_pos_num, vector<int>(&D_pos)[2], vector<int>(&N_pos)[2], Matrix<double, 81, 81>&, double);
	void newmatrix();
	void coutresult();
	~Adam() {
		delete Z_new;
		delete H_mid; 
		delete K_mid; 
		delete D_mid; 
		delete C_mid;
		delete v;
		delete s;
		delete vt;
		delete st;
		delete gt;
	}
private:

	int Z_size1 = 80;
	int Z_size2 = 80;
	double D_center = 8;
	double N_center = 14;
	double inset = 2.5;
	double r = 3;
	double n = 1.56;
	double beta1 = 0.99155;
	double beta2 = 0.999;
	Matrix<double, 1, 81 * 81>deriva_I_Z;
	Matrix<double, 81, 81>Res;
	string Z_final, Alpha_smooth, Beta_smooth, S_smooth, Para, Beta1, Beta2, Dcenter, Ncenter, Inset, t_R, t_N;
	int P_ara = 4;
	double eta;
	Matrix<double, 81, 81>Z;
	Matrix<double, 81, 81>A;
	Matrix<double, 81, 81>B;
	Matrix<double, 81, 81>S;
	Matrix<int, 81, 81>X;
	Matrix<int, 81, 81>Y;
	Matrix<double, 81, 81>DER_INCREMENTS_Z;
	vector<double>I_u;

	Matrix<double, 81, 81>Z_x, Z_y, Z_xy, Z_xx, Z_yy;

	Matrix<double, 81, 81>*v, *s, *vt, *st, *gt;
	Matrix<double, 81, 81>*Z_new, *H_mid, *K_mid, *D_mid, *C_mid;

	void ccc(Matrix<double, 1, 81 * 81>&);
	template <typename T, int rows, int cols>
	void coutmatrix(string& out, Matrix<T, rows, cols>& C_final, Matrix<T, rows, cols>& D_final);
	int ttt = 0;
};