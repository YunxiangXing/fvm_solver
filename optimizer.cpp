#include"optimizer.h"

//Matrix<double, 2, 2> operator+(const double I, Eigen::Matrix<double, 2, 2>&& A) {
//	return I + A.array();
//}
//Matrix<double, 2, 2> dot(Matrix<double, 2, 2>& A, Matrix<double, 2, 2>& B) {
//	return A.array() * B.array();
//};

void Adam::error(string name) {
	cout << "Not found " << name << endl;
	exit(-1);
}
void Adam::check(){

	//--zfin=
	string zfin = Z_final.substr(0, 7);
	if (zfin != "--zfin=") {
		error("Z_final");
	}
	Z_final = Z_final.substr(7, Z_final.size() - 7);

	//--alp=
	string alp = Alpha_smooth.substr(0, 6);
	if (alp != "--alp=") {
		error("Alpha_smooth");
	}
	Alpha_smooth = Alpha_smooth.substr(6, Alpha_smooth.size() - 6);
	
	//--bet=
	string bet = Beta_smooth.substr(0, 6);
	if (bet != "--bet=") {
		error("Beta_smooth");
	}
	Beta_smooth = Beta_smooth.substr(6, Beta_smooth.size() - 6);

	//--ssm=
	string ssm = S_smooth.substr(0, 6);
	if (ssm != "--ssm=") {
		error("S_smooth");
	}
	S_smooth = S_smooth.substr(6, S_smooth.size() - 6);

	//--para=
	if (Para.size() != 0) {
		string para = Para.substr(0, 7);
		if (para != "--para=") {
			cout << "Thread input format: --para=num" << endl;
			exit(-1);
		}
		Para = Para.substr(7, Para.size() - 7);
		for (int i = 0; i < Para.size(); i++) {
			if (!isdigit(Para[i])) {
				cout << "Thread must be positive integer" << endl;
				exit(-1);
			}
		}
		
		P_ara = stoi(Para);
	}

	//--beta1=
	if (Beta1.size() != 0) {
		string beta11 = Beta1.substr(0, 8);
		if (beta11 != "--beta1=") {
			cout << "beta format: --beta1=(0.0, 1.0)" << endl;
			exit(-1);
		}
		Beta1 = Beta1.substr(8, Beta1.size() - 8);
		for (int i = 0; i < Beta1.size(); i++) {
			if (!isdigit(Beta1[i]) && Beta1[i] != '.') {
				cout << "Beta1 must be positive double" << endl;
				exit(-1);
			}
		}

		beta1 = stod(Beta1);
	}

	//--beta2=
	if (Beta2.size() != 0) {
		string beta22 = Beta2.substr(0, 8);
		if (beta22 != "--beta2=") {
			cout << "beta format: --beta2=(0.0, 1.0)" << endl;
			exit(-1);
		}
		Beta2 = Beta2.substr(8, Beta2.size() - 8);
		for (int i = 0; i < Beta2.size(); i++) {
			if (!isdigit(Beta2[i]) && Beta2[i] != '.') {
				cout << "Beta2 must be positive double" << endl;
				exit(-1);
			}
		}

		beta2 = stod(Beta2);
	}

	//--Dcenter=
	if (Dcenter.size() != 0) {
		string Dcenter22 = Dcenter.substr(0, 10);
		if (Dcenter22 != "--Dcenter=") {
			cout << "Dcenter format: --Dcenter=1.0" << endl;
			exit(-1);
		}
		Dcenter = Dcenter.substr(10, Dcenter.size() - 10);
		for (int i = 0; i < Dcenter.size(); i++) {
			if (!isdigit(Dcenter[i]) && Dcenter[i] != '.') {
				cout << "Dcenter must be positive double" << endl;
				exit(-1);
			}
		}

		D_center = stod(Dcenter);
	}

	//--Ncenter=
	if (Ncenter.size() != 0) {
		string Ncenter22 = Ncenter.substr(0, 10);
		if (Ncenter22 != "--Ncenter=") {
			cout << "Ncenter format: --Ncenter=1.0" << endl;
			exit(-1);
		}
		Ncenter = Ncenter.substr(10, Ncenter.size() - 10);
		for (int i = 0; i < Ncenter.size(); i++) {
			if (!isdigit(Ncenter[i]) && Ncenter[i] != '.') {
				cout << "Ncenter must be positive double" << endl;
				exit(-1);
			}
		}

		N_center = stod(Ncenter);
	}

	//--Inset=
	if (Inset.size() != 0) {
		string Inset22 = Inset.substr(0, 8);
		if (Inset22 != "--Inset=") {
			cout << "Inset format: --Inset=1.0" << endl;
			exit(-1);
		}
		Inset = Inset.substr(8, Inset.size() - 8);
		for (int i = 0; i < Inset.size(); i++) {
			if (!isdigit(Inset[i]) && Inset[i] != '.') {
				if (i == 0 && Inset[i] == '-') {
					continue;
				}
				else {
					cout << "Inset must be positive double" << endl;
					exit(-1);
				}
			}
		}

		inset = stod(Inset);
	}

	//--r=
	if (t_R.size() != 0) {
		string t_R22 = t_R.substr(0, 4);
		if (t_R22 != "--r=") {
			cout << "r format: --r=1.0" << endl;
			exit(-1);
		}
		t_R = t_R.substr(4, t_R.size() - 4);
		for (int i = 0; i < t_R.size(); i++) {
			if (!isdigit(t_R[i]) && t_R[i] != '.') {
				cout << "t_R must be positive double" << endl;
				exit(-1);
			}
		}

		r = stod(t_R);
	}

	//--n=
	if (t_N.size() != 0) {
		string t_N22 = t_N.substr(0, 4);
		if (t_N22 != "--n=") {
			cout << "n format: --n=1.56" << endl;
			exit(-1);
		}
		t_N = t_N.substr(4, t_N.size() - 4);
		for (int i = 0; i < t_N.size(); i++) {
			if (!isdigit(t_N[i]) && t_N[i] != '.') {
				cout << "t_N must be positive double" << endl;
				exit(-1);
			}
		}

		n = stod(t_N);
	}
}
bool Adam::ReadDate() {
	cout << "Start reading data" << endl;

	ifstream zfin(Z_final);
	ifstream alp(Alpha_smooth);
	ifstream bet(Beta_smooth);
	ifstream ssm(S_smooth);
	if (!zfin.is_open()) {
		error(Z_final);
	}
	if (!alp.is_open()) {
		error(Alpha_smooth);
	}
	if (!bet.is_open()) {
		error(Beta_smooth);
	}
	if (!ssm.is_open()) {
		error(S_smooth);
	}
	
	double temp = 0.0;
	for (int i = 0; i < 81; i++)
		for (int j = 0; j < 81; j++) {
			zfin >> temp;
			Z(i, j) = temp;
			alp >> temp;
			A(i, j) = temp;
			bet >> temp;
			B(i, j) = temp;
			ssm >> temp;
			S(i, j) = temp;
		}

	cout << "Data reading completed" << endl;
	return 1;
}
void Adam::partial_derivative(Matrix<double, 81, 81>& Z_S, int Z_size1, int Z_size2, double space_x, double space_y) {
	Z_size1 = 81;
	Z_size2 = 81;
	for (int i = 0; i < 81; i++) {
		for (int j = 0; j < 81; j++) {
			Z_x(i, j) = 0;
			Z_y(i, j) = 0;
			Z_xx(i, j) = 0;
			Z_yy(i, j) = 0;
			Z_xy(i, j) = 0;
		}
	}

	for (int i = 0; i < Z_size1; i++) {
		for (int j = 1; j < Z_size2 - 1; j++) {
			Z_x(i, j) = (Z_S(i, j + 1) - Z_S(i, j - 1)) / (2.0 * space_x);
			Z_xx(i, j) = (Z_S(i, j + 1) + Z_S(i, j - 1) - 2 * Z_S(i, j)) / (space_x * space_x);
		}
		Z_x(i, 0) = (Z_S(i, 1) - Z_S(i, 0)) / space_x;
		Z_xx(i, 0) = (Z_S(i, 2) + Z_S(i, 0) - 2 * Z_S(i, 1)) / (space_x * space_x);
		Z_x(i, Z_size2 - 1) = (Z_S(i, Z_size2 - 1) - Z_S(i, Z_size2 - 2)) / space_x;
		Z_xx(i, Z_size2 - 1) = (Z_S(i, Z_size2 - 1) + Z_S(i, Z_size2 - 3) - 2 * Z_S(i, Z_size2 - 2)) / (space_x * space_x);
	}

	for (int j = 0; j < Z_size2; j++) {
		for (int i = 1; i < Z_size1 - 1; i++) {
			Z_y(i, j) = (Z_S(i + 1, j) - Z_S(i - 1, j)) / (2.0 * space_y);
			Z_yy(i, j) = (Z_S(i + 1, j) + Z_S(i - 1, j) - 2 * Z_S(i, j)) / (space_y * space_y);
		}
		Z_y(0, j) = (Z_S(1, j) - Z_S(0, j)) / space_y;
		Z_yy(0, j) = (Z_S(2, j) + Z_S(0, j) - 2 * Z_S(1, j)) / (space_y * space_y);
		Z_y(Z_size1 - 1, j) = (Z_S(Z_size1 - 1, j) - Z_S(Z_size1 - 2, j)) / space_y;
		Z_yy(Z_size1 - 1, j) = (Z_S(Z_size1 - 1, j) + Z_S(Z_size1 - 3, j) - 2 * Z_S(Z_size1 - 2, j)) / (space_y * space_y);
	}

	for (int i = 1; i < Z_size1 - 1; i++) {
		for (int j = 1; j < Z_size2 - 1; j++) {
			Z_xy(i, j) = (Z_S(i + 1, j + 1) + Z_S(i - 1, j - 1) - Z_S(i + 1, j - 1) - Z_S(i - 1, j + 1)) / (4 * space_x * space_y);
		}
		Z_xy(i, 0) = (Z_S(i + 1, 2) + Z_S(i - 1, 0) - Z_S(i + 1, 0) - Z_S(i - 1, 2)) / (4 * space_x * space_y);
		Z_xy(i, Z_size2 - 1) = (Z_S(i + 1, Z_size2 - 1) + Z_S(i - 1, Z_size2 - 3) - Z_S(i + 1, Z_size2 - 3) - Z_S(i - 1, Z_size2 - 1)) / (4 * space_x * space_y);
	}

	for (int j = 1; j < Z_size2 - 1; j++) {
		Z_xy(0, j) = (Z_S(2, j + 1) + Z_S(0, j - 1) - Z_S(0, j + 1) - Z_S(2, j - 1)) / (4 * space_x * space_y);
		Z_xy(Z_size1 - 1, j) = (Z_S(Z_size1 - 1, j + 1) + Z_S(Z_size1 - 3, j - 1) - Z_S(Z_size1 - 3, j + 1) - Z_S(Z_size1 - 1, j - 1)) / (4 * space_x * space_y);
	}
	Z_xy(0, 0) = (Z_S(0, 0) + Z_S(2, 2) - Z_S(0, 2) - Z_S(2, 0)) / (4 * space_x * space_y);
	Z_xy(Z_size1 - 1, 0) = (Z_S(Z_size1 - 3, 0) + Z_S(Z_size1 - 1, 2) - Z_S(Z_size1 - 3, 2) - Z_S(Z_size1 - 1, 0)) / (4 * space_x * space_y);
	Z_xy(0, Z_size2 - 1) = (Z_S(0, Z_size2 - 3) + Z_S(2, Z_size2 - 1) - Z_S(2, Z_size2 - 3) - Z_S(0, Z_size2 - 1)) / (4 * space_x * space_y);
	Z_xy(Z_size1 - 1, Z_size2 - 1) = (Z_S(Z_size1 - 1, Z_size2 - 1) + Z_S(Z_size1 - 3, Z_size2 - 3) - Z_S(Z_size1 - 3, Z_size2 - 1) - Z_S(Z_size1 - 1, Z_size2 - 3)) / (4 * space_x * space_y);
}
vector<double> Adam::min(int D_pos_num, int N_pos_num,vector<int>(&D_pos)[2],vector<int>(&N_pos)[2],Matrix<double, 81, 81>&C_first, Matrix<double, 81, 81>&D_first, Matrix<double, 81, 81>&S_smooth) {
	int DN_pos_num = D_pos_num + N_pos_num;
	auto DN_pos_x = D_pos[0];
	auto DN_pos_y = D_pos[1];
	DN_pos_x.reserve(D_pos[0].size()+N_pos[0].size());
	DN_pos_y.reserve(D_pos[1].size() + N_pos[1].size());
	for (int i = 0; i < N_pos[0].size();i++) {
		DN_pos_x.push_back(N_pos[0][i]);
		DN_pos_y.push_back(N_pos[1][i]);
	}
	vector<double>DN_pos_C;
	double min1 = 1e9, max1 = -1e9;
	int dic_min = -1, dic_max = -1;
	for (int i = 0; i < DN_pos_num; i++) {
		//DN_pos_C.push_back(C_first(DN_pos_x[i], DN_pos_y[i]));
		if (C_first(DN_pos_x[i], DN_pos_y[i]) < min1) {
			min1 = C_first(DN_pos_x[i], DN_pos_y[i]);
			dic_min = i;
		}
		double temp = abs(D_first(DN_pos_x[i], DN_pos_y[i]) - (1 - n) * S_smooth(DN_pos_x[i], DN_pos_y[i]) * 1e3);
		if (max1 < temp) {
			max1 = temp;
			dic_max = i;
		}
	}
	vector<double>A;
	A.push_back(min1);
	A.push_back(max1);
	A.push_back(dic_min);
	A.push_back(dic_max);
	return A;
}
template <typename T,int rows,int cols>
double Adam::evaluation_func(Matrix<T, rows, cols>& Alpha, Matrix<T, rows, cols>& Beta, Matrix<T, rows, cols>& Ssmooth, Matrix<T, rows, cols>& z_x, Matrix<T, rows, cols>& z_y, Matrix<T, rows, cols>& z_xx,
	Matrix<T, rows, cols>& z_yy, Matrix<T, rows, cols>& z_xy, double min, double max) {
	double I;
	double space_x = 1.0;
	double space_y = 1.0;
	auto temp = (Alpha.array() * ((((1 + z_y.array() * z_y.array()) * z_xx.array() - 2 * z_x.array() * z_y.array() * z_xy.array() + 
		(1 + z_x.array() * z_x.array()) * z_yy.array()) / (2 * (1 + z_x.array() * z_x.array() + z_y.array() * z_y.array()).pow(1.5))).pow(2) - 
		(z_xx.array() * z_yy.array() - z_xy.array() * z_xy.array()) / (1 + z_x.array() * z_x.array() + z_y.array() * z_y.array()).pow(2)) + 
		Beta.array() * (((1 + z_y.array() * z_y.array()) * z_xx.array() - 2 * z_x.array() * z_y.array() * z_xy.array() + 
			(1 + z_x.array() * z_x.array()) * z_yy.array()) / (2 * (1 + z_x.array() * z_x.array() + z_y.array().pow(2)).pow(1.5)) - 
			Ssmooth.array()).pow(2)) * (1 + z_x.array().pow(2) + z_y.array().pow(2)).pow(0.5) * space_x * space_y;
	auto temp1 = (1.0 + z_y.array().pow(2.0)) * z_xx.array();
	auto temp2 = (temp1.array() - 2.0 * z_x.array() * z_y.array() * z_xy.array() + (1 + z_x.array().pow(2.0)) * z_yy.array());
	auto temp3 = (1.0 + z_x.array().pow(2.0) + z_y.array().pow(2.0)).array().pow(1.5);
	auto temp4 = ((temp2.array() / (2 * temp3.array())).pow(2.0) - (z_xx.array() *z_yy.array() - z_xy.array().pow(2.0)) / (1 + z_x.array().pow(2.0) + z_y.array().pow(2.0)).pow(2.0));
	auto temp5 = ((1.0 + z_y.array().pow(2)) * z_xx.array() - 2.0 * z_x.array() * z_y.array() * z_xy.array() + (1.0 + z_x.array().pow(2)) * z_yy.array());
	auto temp6 = (2.0 * (1.0 + z_x.array().pow(2.0) + z_y.array().pow(2)).array().pow(1.5));
	auto a = (Alpha.array() * temp4.array() + Beta.array() * (temp5.array() / temp6.array() - Ssmooth.array()).pow(2.0)) * sqrt(1 + z_x.array().pow(2.0) + z_y.array().pow(2)) * space_x * space_y;


	double temp7 = temp.sum();
	double temp8 = a.sum();
	double RES = temp7 + min + max;
	double RES1 = temp8 + min + max;
	return RES1;
}
template <typename T,int rows, int cols>
double Adam::evaluation_func0(Matrix<T, rows, cols>& Alpha, Matrix<T, rows, cols>& Beta, Matrix<T, rows, cols>& Ss, Matrix<T, rows, cols>& z_x, Matrix<T, rows, cols>& z_y, Matrix<T, rows, cols>& z_xx, 
	Matrix<T, rows, cols>& z_yy, Matrix<T, rows, cols>& z_xy) {
	double space_x = 1.0;
	double space_y = 1.0;
	Matrix<T, rows, cols> I_u_matrix;
	I_u_matrix = (Alpha.array()*((((1 + z_y.array().pow(2)) * z_xx.array() - 2 * z_x.array() * z_y.array() * z_xy.array() + (1 + z_x.array().pow(2)) * z_yy.array()) / (2 * (1 + z_x.array().pow(2) + z_y.array().pow(2)).pow(1.5))).pow(2) - (z_xx.array() * z_yy.array() - z_xy.array().pow(2)) /
		(1 + z_x.array().pow(2) + z_y.array().pow(2)).pow(2)) + Beta.array()*(((1 + z_y.array().pow(2)) * z_xx.array() - 2 * z_x.array() * z_y.array() * z_xy.array() + (1 + z_x.array().pow(2)) * z_yy.array()) / (2 * (1 + z_x.array().pow(2) + z_y.array().pow(2)).pow(1.5)) - Ss.array()).pow(2)) * sqrt(1 + z_x.array().pow(2) + z_y.array().pow(2)) * space_x * space_y;
	double temp = I_u_matrix.sum();
	return temp;
}
template <typename T, int rows, int cols>
void Adam::Zeros(Matrix<T, rows, cols>& matrix) {
	for (int i = 0; i < matrix.rows(); i++)
		for (int j = 0; j < matrix.cols(); j++) {
			matrix(i, j) = 0.0;
		}
}

void Adam::calculate_A22() {
	
	//z_x_add1 = Zeros(2, 2); z_y_add1 = zeros(2, 2); z_xx_add1 = zeros(2, 2); z_yy_add1 = zeros(2, 2); z_xy_add1 = zeros(2, 2);
	Matrix<double, 2, 2>z_x_add1, z_y_add1, z_xx_add1, z_yy_add1, z_xy_add1;
	Zeros(z_x_add1);
	Zeros(z_y_add1);
	Zeros(z_xx_add1);
	Zeros(z_xy_add1);
	Zeros(z_yy_add1);
	z_x_add1(0, 0) = -DER_INCREMENTS_Z(0, 0); z_x_add1(0, 1) = -DER_INCREMENTS_Z(0, 0) / 2;
	z_y_add1(0, 0) = -DER_INCREMENTS_Z(0, 0); z_y_add1(1, 0) = -DER_INCREMENTS_Z(0, 0) / 2;
	z_xx_add1(0, 0) = DER_INCREMENTS_Z(0, 0); z_xx_add1(0, 1) = DER_INCREMENTS_Z(0, 0);
	z_yy_add1(0, 0) = DER_INCREMENTS_Z(0, 0); z_yy_add1(1, 0) = DER_INCREMENTS_Z(0, 0);
	z_xy_add1(0, 0) = DER_INCREMENTS_Z(0, 0) / 4; z_xy_add1(0, 1) = DER_INCREMENTS_Z(0, 0) / 4;
	z_xy_add1(1, 0) = DER_INCREMENTS_Z(0, 0) / 4; z_xy_add1(1, 1) = DER_INCREMENTS_Z(0, 0) / 4;

	Matrix<double, 2, 2>z_x_add2, z_y_add2, z_xx_add2, z_yy_add2, z_xy_add2;
	z_x_add2 = -z_x_add1; z_y_add2 = -z_y_add1; z_xx_add2 = -z_xx_add1; z_yy_add2 = -z_yy_add1; z_xy_add2 = -z_xy_add1;
	Matrix<double, 2, 2>z_xo, z_yo, z_xxo, z_yyo, z_xyo;
	Matrix<double, 2, 2>z_x_new1, z_y_new1, z_xx_new1, z_yy_new1, z_xy_new1;
	Matrix<double, 2, 2>z_x_new2, z_y_new2, z_xx_new2, z_yy_new2, z_xy_new2;
	Matrix<double, 2, 2>Alpha_new, Beta_new, P_new;
	int ix = 0, jy = 0;
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++) {
			z_xo(i, j) = Z_x(i, j);
			z_yo(i, j) = Z_y(i, j);
			z_xxo(i, j) = Z_xx(i, j);
			z_yyo(i, j) = Z_yy(i, j);
			z_xyo(i, j) = Z_xy(i, j);
		}
	z_x_new1 = z_xo + z_x_add1; z_y_new1 = z_yo + z_y_add1; z_xx_new1 = z_xxo + z_xx_add1; z_yy_new1 = z_yyo + z_yy_add1; z_xy_new1 = z_xyo + z_xy_add1;
	z_x_new2 = z_xo + z_x_add2; z_y_new2 = z_yo + z_y_add2; z_xx_new2 = z_xxo + z_xx_add2; z_yy_new2 = z_yyo + z_yy_add2; z_xy_new2 = z_xyo + z_xy_add2;

	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++) {
			Alpha_new(i, j) = A(i, j);
			Beta_new(i, j) = B(i, j);
			P_new(i, j) = S(i, j);
		}

	double I_u_1 = evaluation_func0(Alpha_new, Beta_new, P_new, z_x_new1, z_y_new1, z_xx_new1, z_yy_new1, z_xy_new1);
	double I_u_2 = evaluation_func0(Alpha_new, Beta_new, P_new, z_x_new2, z_y_new2, z_xx_new2, z_yy_new2, z_xy_new2);
	deriva_I_Z(0) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(0, 0));
	//=============================================================================================================================//
	Zeros(z_x_add1);
	Zeros(z_y_add1);
	Zeros(z_xx_add1);
	Zeros(z_xy_add1);
	Zeros(z_yy_add1);
	z_x_add1(0, 0) = DER_INCREMENTS_Z(0, Z_size2) / 2; z_x_add1(0, 1) = DER_INCREMENTS_Z(0, Z_size2);
	z_y_add1(0, 1) = -DER_INCREMENTS_Z(0, Z_size2); z_y_add1(1, 1) = -DER_INCREMENTS_Z(0, Z_size2) / 2;
	z_xx_add1(0, 0) = DER_INCREMENTS_Z(0, Z_size2); z_xx_add1(0, 1) = DER_INCREMENTS_Z(0, Z_size2);
	z_yy_add1(0, 1) = DER_INCREMENTS_Z(0, Z_size2); z_yy_add1(1, 1) = DER_INCREMENTS_Z(0, Z_size2);
	z_xy_add1(0, 0) = -DER_INCREMENTS_Z(0, Z_size2) / 4; z_xy_add1(0, 1) = -DER_INCREMENTS_Z(0, Z_size2) / 4;
	z_xy_add1(1, 0) = -DER_INCREMENTS_Z(0, Z_size2) / 4; z_xy_add1(1, 1) = -DER_INCREMENTS_Z(0, Z_size2) / 4;
	z_x_add2 = -z_x_add1; z_y_add2 = -z_y_add1; z_xx_add2 = -z_xx_add1; z_yy_add2 = -z_yy_add1; z_xy_add2 = -z_xy_add1;

	for (int i = 0; i < 2; i++) {
		jy = 0;
		for (int j = Z_size2 - 1; j <= Z_size2; j++) {
			z_xo(i, jy) = Z_x(i, j);
			z_yo(i, jy) = Z_y(i, j);
			z_xxo(i, jy) = Z_xx(i, j);
			z_yyo(i, jy) = Z_yy(i, j);
			z_xyo(i, jy) = Z_xy(i, j);
			jy++;
		}
	}
	z_x_new1 = z_xo + z_x_add1; z_y_new1 = z_yo + z_y_add1; z_xx_new1 = z_xxo + z_xx_add1; z_yy_new1 = z_yyo + z_yy_add1; z_xy_new1 = z_xyo + z_xy_add1;
	z_x_new2 = z_xo + z_x_add2; z_y_new2 = z_yo + z_y_add2; z_xx_new2 = z_xxo + z_xx_add2; z_yy_new2 = z_yyo + z_yy_add2; z_xy_new2 = z_xyo + z_xy_add2;

	for (int i = 0; i < 2; i++) {
		jy = 0;
		for (int j = Z_size2 - 1; j <= Z_size2; j++) {
			Alpha_new(i, jy) = A(i, j);
			Beta_new(i, jy) = B(i, j);
			P_new(i, jy) = S(i, j);
			jy++;
		}
	}
	jy = 0;
	I_u_1 = evaluation_func0(Alpha_new, Beta_new, P_new, z_x_new1, z_y_new1, z_xx_new1, z_yy_new1, z_xy_new1);
	I_u_2 = evaluation_func0(Alpha_new, Beta_new, P_new, z_x_new2, z_y_new2, z_xx_new2, z_yy_new2, z_xy_new2);
	deriva_I_Z(Z_size2) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(0, Z_size2));
	//=============================================================================================================================//
	Zeros(z_x_add1);
	Zeros(z_y_add1);
	Zeros(z_xx_add1);
	Zeros(z_xy_add1);
	Zeros(z_yy_add1);
	z_x_add1(1, 0) = -DER_INCREMENTS_Z(Z_size1, 0); z_x_add1(1, 1) = -DER_INCREMENTS_Z(Z_size1, 0) / 2;
	z_y_add1(0, 0) = DER_INCREMENTS_Z(Z_size1, 0) / 2; z_y_add1(1, 0) = DER_INCREMENTS_Z(Z_size1, 0);
	z_xx_add1(1, 0) = DER_INCREMENTS_Z(Z_size1, 0); z_xx_add1(1, 1) = DER_INCREMENTS_Z(Z_size1, 0);
	z_yy_add1(0, 0) = DER_INCREMENTS_Z(Z_size1, 0); z_yy_add1(1, 0) = DER_INCREMENTS_Z(Z_size1, 0);
	z_xy_add1(0, 0) = -DER_INCREMENTS_Z(Z_size1, 0) / 4; z_xy_add1(0, 1) = -DER_INCREMENTS_Z(Z_size1, 0) / 4;
	z_xy_add1(1, 0) = -DER_INCREMENTS_Z(Z_size1, 0) / 4; z_xy_add1(1, 1) = -DER_INCREMENTS_Z(Z_size1, 0) / 4;
	z_x_add2 = -z_x_add1; z_y_add2 = -z_y_add1; z_xx_add2 = -z_xx_add1; z_yy_add2 = -z_yy_add1; z_xy_add2 = -z_xy_add1;

	for (int i = Z_size1 - 1; i <= Z_size1; i++) {
		for (int j = 0; j <= 1; j++) {
			z_xo(ix, j) = Z_x(i, j);
			z_yo(ix, j) = Z_y(i, j);
			z_xxo(ix, j) = Z_xx(i, j);
			z_yyo(ix, j) = Z_yy(i, j);
			z_xyo(ix, j) = Z_xy(i, j);
		}
		ix++;
	}
	ix = 0;
	z_x_new1 = z_xo + z_x_add1; z_y_new1 = z_yo + z_y_add1; z_xx_new1 = z_xxo + z_xx_add1; z_yy_new1 = z_yyo + z_yy_add1; z_xy_new1 = z_xyo + z_xy_add1;
	z_x_new2 = z_xo + z_x_add2; z_y_new2 = z_yo + z_y_add2; z_xx_new2 = z_xxo + z_xx_add2; z_yy_new2 = z_yyo + z_yy_add2; z_xy_new2 = z_xyo + z_xy_add2;
	for (int i = Z_size1 - 1; i <= Z_size1; i++) {
		for (int j = 0; j <= 1; j++) {
			Alpha_new(ix, j) = A(i, j);
			Beta_new(ix, j) = B(i, j);
			P_new(ix, j) = S(i, j);
		}
		ix++;
	}
	ix = 0;
	I_u_1 = evaluation_func0(Alpha_new, Beta_new, P_new, z_x_new1, z_y_new1, z_xx_new1, z_yy_new1, z_xy_new1);
	I_u_2 = evaluation_func0(Alpha_new, Beta_new, P_new, z_x_new2, z_y_new2, z_xx_new2, z_yy_new2, z_xy_new2);
	deriva_I_Z((Z_size1) * (Z_size2 + 1)) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(Z_size1, 0));
	//=============================================================================================================================//
	Zeros(z_x_add1);
	Zeros(z_y_add1);
	Zeros(z_xx_add1);
	Zeros(z_xy_add1);
	Zeros(z_yy_add1);
	z_x_add1(1, 0) = DER_INCREMENTS_Z(Z_size1, Z_size2) / 2; z_x_add1(1, 1) = DER_INCREMENTS_Z(Z_size1, Z_size2);
	z_y_add1(0, 1) = DER_INCREMENTS_Z(Z_size1, Z_size2) / 2; z_y_add1(1, 1) = DER_INCREMENTS_Z(Z_size1, Z_size2);
	z_xx_add1(1, 0) = DER_INCREMENTS_Z(Z_size1, Z_size2); z_xx_add1(1, 1) = DER_INCREMENTS_Z(Z_size1, Z_size2);
	z_yy_add1(0, 1) = DER_INCREMENTS_Z(Z_size1, Z_size2); z_yy_add1(1, 1) = DER_INCREMENTS_Z(Z_size1, Z_size2);
	z_xy_add1(0, 0) = DER_INCREMENTS_Z(Z_size1, Z_size2) / 4; z_xy_add1(0, 1) = DER_INCREMENTS_Z(Z_size1, Z_size2) / 4;
	z_xy_add1(1, 0) = DER_INCREMENTS_Z(Z_size1, Z_size2) / 4; z_xy_add1(1, 1) = DER_INCREMENTS_Z(Z_size1, Z_size2) / 4;
	z_x_add2 = -z_x_add1; z_y_add2 = -z_y_add1; z_xx_add2 = -z_xx_add1; z_yy_add2 = -z_yy_add1; z_xy_add2 = -z_xy_add1;

	for (int i = Z_size1 - 1; i <= Z_size1; i++) {
		for (int j = Z_size2 - 1; j <= Z_size2; j++) {
			z_xo(ix, jy) = Z_x(i, j);
			z_yo(ix, jy) = Z_y(i, j);
			z_xxo(ix, jy) = Z_xx(i, j);
			z_yyo(ix, jy) = Z_yy(i, j);
			z_xyo(ix, jy) = Z_xy(i, j);
			jy++;
		}
		ix++; jy = 0;
	}
	ix = 0; jy = 0;
	z_x_new1 = z_xo + z_x_add1; z_y_new1 = z_yo + z_y_add1; z_xx_new1 = z_xxo + z_xx_add1; z_yy_new1 = z_yyo + z_yy_add1; z_xy_new1 = z_xyo + z_xy_add1;
	z_x_new2 = z_xo + z_x_add2; z_y_new2 = z_yo + z_y_add2; z_xx_new2 = z_xxo + z_xx_add2; z_yy_new2 = z_yyo + z_yy_add2; z_xy_new2 = z_xyo + z_xy_add2;

	for (int i = Z_size1 - 1; i <= Z_size1; i++) {
		for (int j = Z_size2 - 1; j <= Z_size2; j++) {
			Alpha_new(ix, jy) = A(i, j);
			Beta_new(ix, jy) = B(i, j);
			P_new(ix, jy) = S(i, j);
			jy++;
		}
		ix++; jy = 0;
	}
	ix = 0; jy = 0;
	I_u_1 = evaluation_func0(Alpha_new, Beta_new, P_new, z_x_new1, z_y_new1, z_xx_new1, z_yy_new1, z_xy_new1);
	I_u_2 = evaluation_func0(Alpha_new, Beta_new, P_new, z_x_new2, z_y_new2, z_xx_new2, z_yy_new2, z_xy_new2);
	deriva_I_Z((Z_size1 + 1) * (Z_size2 + 1) - 1) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(Z_size1, Z_size2));
}
void Adam::calculate_A23() {
	//=============================================================================================================================//
	Matrix<double, 2, 3> z_x_add13, z_y_add13, z_xx_add13, z_yy_add13, z_xy_add13;
	Matrix<double, 2, 3> z_x_add23, z_y_add23, z_xx_add23, z_yy_add23, z_xy_add23;
	Matrix<double, 2, 3> z_xo3, z_yo3, z_xxo3, z_yyo3, z_xyo3;
	Matrix<double, 2, 3> z_x_new13, z_y_new13, z_xx_new13, z_yy_new13, z_xy_new13;
	Matrix<double, 2, 3> z_x_new23, z_y_new23, z_xx_new23, z_yy_new23, z_xy_new23;
	Matrix<double, 2, 3> Alpha_new3, Beta_new3, P_new3;
	int ix = 0, jy = 0;
	Zeros(z_x_add13);
	Zeros(z_y_add13);
	Zeros(z_xx_add13);
	Zeros(z_xy_add13);
	Zeros(z_yy_add13);
	z_x_add13(0, 0) = DER_INCREMENTS_Z(0, 1); z_x_add13(0, 2) = -DER_INCREMENTS_Z(0, 1) / 2;
	z_y_add13(0, 1) = -DER_INCREMENTS_Z(0, 1); z_y_add13(1, 1) = -DER_INCREMENTS_Z(0, 1) / 2;
	z_xx_add13(0, 0) = -2 * DER_INCREMENTS_Z(0, 1); z_xx_add13(0, 1) = -2 * DER_INCREMENTS_Z(0, 1); z_xx_add13(0, 2) = DER_INCREMENTS_Z(0, 1);
	z_yy_add13(0, 1) = DER_INCREMENTS_Z(0, 1); z_yy_add13(1, 1) = DER_INCREMENTS_Z(0, 1);
	z_xy_add13(0, 2) = DER_INCREMENTS_Z(0, 1) / 4; z_xy_add13(1, 2) = DER_INCREMENTS_Z(0, 1) / 4;
	z_x_add23 = -z_x_add13; z_y_add23 = -z_y_add13; z_xx_add23 = -z_xx_add13; z_yy_add23 = -z_yy_add13; z_xy_add23 = -z_xy_add13;

	for (int i = 0; i <= 1; i++)
		for (int j = 0; j <= 2; j++) {
			z_xo3(i, j) = Z_x(i, j);
			z_yo3(i, j) = Z_y(i, j);
			z_xxo3(i, j) = Z_xx(i, j);
			z_yyo3(i, j) = Z_yy(i, j);
			z_xyo3(i, j) = Z_xy(i, j);
		}
	z_x_new13 = z_xo3 + z_x_add13; z_y_new13 = z_yo3 + z_y_add13; z_xx_new13 = z_xxo3 + z_xx_add13; z_yy_new13 = z_yyo3 + z_yy_add13; z_xy_new13 = z_xyo3 + z_xy_add13;
	z_x_new23 = z_xo3 + z_x_add23; z_y_new23 = z_yo3 + z_y_add23; z_xx_new23 = z_xxo3 + z_xx_add23; z_yy_new23 = z_yyo3 + z_yy_add23; z_xy_new23 = z_xyo3 + z_xy_add23;

	for (int i = 0; i <= 1; i++)
		for (int j = 0; j <= 2; j++) {
			Alpha_new3(i, j) = A(i, j);
			Beta_new3(i, j) = B(i, j);
			P_new3(i, j) = S(i, j);
		}

	double I_u_1 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new13, z_y_new13, z_xx_new13, z_yy_new13, z_xy_new13);
	double I_u_2 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new23, z_y_new23, z_xx_new23, z_yy_new23, z_xy_new23);
	deriva_I_Z(1) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(0, 1));
	//=============================================================================================================================//
	Zeros(z_x_add13);
	Zeros(z_y_add13);
	Zeros(z_xx_add13);
	Zeros(z_xy_add13);
	Zeros(z_yy_add13);
	z_x_add13(0, 0) = DER_INCREMENTS_Z(0, Z_size2 - 1) / 2; z_x_add13(0, 2) = -DER_INCREMENTS_Z(0, Z_size2 - 1);
	z_y_add13(0, 1) = -DER_INCREMENTS_Z(0, Z_size2 - 1); z_y_add13(1, 1) = -DER_INCREMENTS_Z(0, Z_size2 - 1) / 2;
	z_xx_add13(0, 0) = DER_INCREMENTS_Z(0, Z_size2 - 1); z_xx_add13(0, 1) = -2 * DER_INCREMENTS_Z(0, Z_size2 - 1); z_xx_add13(0, 2) = -2 * DER_INCREMENTS_Z(0, Z_size2 - 1);
	z_yy_add13(0, 1) = DER_INCREMENTS_Z(0, Z_size2 - 1); z_yy_add13(1, 1) = DER_INCREMENTS_Z(0, Z_size2 - 1);
	z_xy_add13(0, 0) = -DER_INCREMENTS_Z(0, Z_size2 - 1) / 4; z_xy_add13(1, 0) = -DER_INCREMENTS_Z(0, Z_size2 - 1) / 4;
	z_x_add23 = -z_x_add13; z_y_add23 = -z_y_add13; z_xx_add23 = -z_xx_add13; z_yy_add23 = -z_yy_add13; z_xy_add23 = -z_xy_add13;
	for (int i = 0; i <= 1; i++) {
		for (int j = Z_size2 - 2; j <= Z_size2; j++) {
			z_xo3(i, jy) = Z_x(i, j);
			z_yo3(i, jy) = Z_y(i, j);
			z_xxo3(i, jy) = Z_xx(i, j);
			z_yyo3(i, jy) = Z_yy(i, j);
			z_xyo3(i, jy) = Z_xy(i, j);
			jy++;
		}
		jy = 0;
	}
	jy = 0;
	z_x_new13 = z_xo3 + z_x_add13; z_y_new13 = z_yo3 + z_y_add13; z_xx_new13 = z_xxo3 + z_xx_add13; z_yy_new13 = z_yyo3 + z_yy_add13; z_xy_new13 = z_xyo3 + z_xy_add13;
	z_x_new23 = z_xo3 + z_x_add23; z_y_new23 = z_yo3 + z_y_add23; z_xx_new23 = z_xxo3 + z_xx_add23; z_yy_new23 = z_yyo3 + z_yy_add23; z_xy_new23 = z_xyo3 + z_xy_add23;
	for (int i = 0; i <= 1; i++) {
		for (int j = Z_size2 - 2; j <= Z_size2; j++) {
			Alpha_new3(i, jy) = A(i, j);
			Beta_new3(i, jy) = B(i, j);
			P_new3(i, jy) = S(i, j);
			jy++;
		}
		jy = 0;
	}
	jy = 0;
	I_u_1 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new13, z_y_new13, z_xx_new13, z_yy_new13, z_xy_new13);
	I_u_2 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new23, z_y_new23, z_xx_new23, z_yy_new23, z_xy_new23);
	deriva_I_Z(Z_size2 - 1) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(0, Z_size2 - 1));
	//=============================================================================================================================//
	Zeros(z_x_add13);
	Zeros(z_y_add13);
	Zeros(z_xx_add13);
	Zeros(z_xy_add13);
	Zeros(z_yy_add13);
	z_x_add13(1, 0) = DER_INCREMENTS_Z(Z_size1, 1); z_x_add13(1, 2) = -DER_INCREMENTS_Z(Z_size1, 1) / 2;
	z_y_add13(0, 1) = DER_INCREMENTS_Z(Z_size1, 1) / 2; z_y_add13(1, 1) = DER_INCREMENTS_Z(Z_size1, 1);
	z_xx_add13(1, 0) = -2 * DER_INCREMENTS_Z(Z_size1, 1); z_xx_add13(1, 1) = -2 * DER_INCREMENTS_Z(Z_size1, 1); z_xx_add13(1, 2) = DER_INCREMENTS_Z(Z_size1, 1);
	z_yy_add13(0, 1) = DER_INCREMENTS_Z(Z_size1, 1); z_yy_add13(1, 1) = DER_INCREMENTS_Z(Z_size1, 1);
	z_xy_add13(0, 2) = -DER_INCREMENTS_Z(Z_size1, 1) / 4; z_xy_add13(1, 2) = -DER_INCREMENTS_Z(Z_size1, 1) / 4;
	z_x_add23 = -z_x_add13; z_y_add23 = -z_y_add13; z_xx_add23 = -z_xx_add13; z_yy_add23 = -z_yy_add13; z_xy_add23 = -z_xy_add13;
	for (int i = Z_size1 - 1; i <= Z_size1; i++) {
		for (int j = 0; j <= 2; j++) {
			z_xo3(ix, j) = Z_x(i, j);
			z_yo3(ix, j) = Z_y(i, j);
			z_xxo3(ix, j) = Z_xx(i, j);
			z_yyo3(ix, j) = Z_yy(i, j);
			z_xyo3(ix, j) = Z_xy(i, j);
		}
		ix++;
	}
	ix = 0;
	z_x_new13 = z_xo3 + z_x_add13; z_y_new13 = z_yo3 + z_y_add13; z_xx_new13 = z_xxo3 + z_xx_add13; z_yy_new13 = z_yyo3 + z_yy_add13; z_xy_new13 = z_xyo3 + z_xy_add13;
	z_x_new23 = z_xo3 + z_x_add23; z_y_new23 = z_yo3 + z_y_add23; z_xx_new23 = z_xxo3 + z_xx_add23; z_yy_new23 = z_yyo3 + z_yy_add23; z_xy_new23 = z_xyo3 + z_xy_add23;
	for (int i = Z_size1 - 1; i <= Z_size1; i++) {
		for (int j = 0; j <= 2; j++) {
			Alpha_new3(ix, j) = A(i, j);
			Beta_new3(ix, j) = B(i, j);
			P_new3(ix, j) = S(i, j);
		}
		ix++;
	}
	ix = 0;
	I_u_1 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new13, z_y_new13, z_xx_new13, z_yy_new13, z_xy_new13);
	I_u_2 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new23, z_y_new23, z_xx_new23, z_yy_new23, z_xy_new23);
	deriva_I_Z((Z_size1) * (Z_size2 + 1) + 1) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(Z_size1, 1));
	//=============================================================================================================================//
	Zeros(z_x_add13);
	Zeros(z_y_add13);
	Zeros(z_xx_add13);
	Zeros(z_xy_add13);
	Zeros(z_yy_add13);
	z_x_add13(1, 0) = DER_INCREMENTS_Z(Z_size1, Z_size2 - 1) / 2; z_x_add13(1, 2) = -DER_INCREMENTS_Z(Z_size1, Z_size2 - 1);
	z_y_add13(0, 1) = DER_INCREMENTS_Z(Z_size1, Z_size2 - 1) / 2; z_y_add13(1, 1) = DER_INCREMENTS_Z(Z_size1, Z_size2 - 1);
	z_xx_add13(1, 0) = DER_INCREMENTS_Z(Z_size1, Z_size2 - 1); z_xx_add13(1, 1) = -2 * DER_INCREMENTS_Z(Z_size1, Z_size2 - 1); z_xx_add13(1, 2) = -2 * DER_INCREMENTS_Z(Z_size1, Z_size2 - 1);
	z_yy_add13(0, 1) = DER_INCREMENTS_Z(Z_size1, Z_size2 - 1); z_yy_add13(1, 1) = DER_INCREMENTS_Z(Z_size1, Z_size2 - 1);
	z_xy_add13(0, 0) = DER_INCREMENTS_Z(Z_size1, Z_size2 - 1) / 4; z_xy_add13(1, 0) = DER_INCREMENTS_Z(Z_size1, Z_size2 - 1) / 4;
	z_x_add23 = -z_x_add13; z_y_add23 = -z_y_add13; z_xx_add23 = -z_xx_add13; z_yy_add23 = -z_yy_add13; z_xy_add23 = -z_xy_add13;
	for (int i = Z_size1 - 1; i <= Z_size1; i++) {
		for (int j = Z_size2 - 2; j <= Z_size2; j++) {
			z_xo3(ix, jy) = Z_x(i, j);
			z_yo3(ix, jy) = Z_y(i, j);
			z_xxo3(ix, jy) = Z_xx(i, j);
			z_yyo3(ix, jy) = Z_yy(i, j);
			z_xyo3(ix, jy) = Z_xy(i, j);
			jy++;
		}
		ix++; jy = 0;
	}
	ix = 0; jy = 0;
	z_x_new13 = z_xo3 + z_x_add13; z_y_new13 = z_yo3 + z_y_add13; z_xx_new13 = z_xxo3 + z_xx_add13; z_yy_new13 = z_yyo3 + z_yy_add13; z_xy_new13 = z_xyo3 + z_xy_add13;
	z_x_new23 = z_xo3 + z_x_add23; z_y_new23 = z_yo3 + z_y_add23; z_xx_new23 = z_xxo3 + z_xx_add23; z_yy_new23 = z_yyo3 + z_yy_add23; z_xy_new23 = z_xyo3 + z_xy_add23;
	for (int i = Z_size1 - 1; i <= Z_size1; i++) {
		for (int j = Z_size2 - 2; j <= Z_size2; j++) {
			Alpha_new3(ix, jy) = A(i, j);
			Beta_new3(ix, jy) = B(i, j);
			P_new3(ix, jy) = S(i, j);
			jy++;
		}
		ix++; jy = 0;
	}
	ix = 0; jy = 0;
	I_u_1 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new13, z_y_new13, z_xx_new13, z_yy_new13, z_xy_new13);
	I_u_2 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new23, z_y_new23, z_xx_new23, z_yy_new23, z_xy_new23);
	deriva_I_Z((Z_size1) * (Z_size2 + 1) + Z_size2 - 1) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(Z_size1, Z_size2 - 1));
	//=============================================================================================================================//
	for (int j = 3; j <= Z_size2 - 3; j++) {
		int i = 0;
		Zeros(z_x_add13);
		Zeros(z_y_add13);
		Zeros(z_xx_add13);
		Zeros(z_xy_add13);
		Zeros(z_yy_add13);
		z_x_add13(0, 0) = DER_INCREMENTS_Z(i, j) / 2; z_x_add13(0, 2) = -DER_INCREMENTS_Z(i, j) / 2;
		z_y_add13(0, 1) = -DER_INCREMENTS_Z(i, j); z_y_add13(1, 1) = -DER_INCREMENTS_Z(i, j) / 2;
		z_xx_add13(0, 0) = DER_INCREMENTS_Z(i, j); z_xx_add13(0, 1) = -2 * DER_INCREMENTS_Z(i, j); z_xx_add13(0, 2) = DER_INCREMENTS_Z(i, j);
		z_yy_add13(0, 1) = DER_INCREMENTS_Z(i, j); z_yy_add13(1, 1) = DER_INCREMENTS_Z(i, j);
		z_xy_add13(0, 0) = -DER_INCREMENTS_Z(i, j) / 4; z_xy_add13(0, 2) = DER_INCREMENTS_Z(i, j) / 4;
		z_xy_add13(1, 0) = -DER_INCREMENTS_Z(i, j) / 4; z_xy_add13(1, 2) = DER_INCREMENTS_Z(i, j) / 4;
		z_x_add23 = -z_x_add13; z_y_add23 = -z_y_add13; z_xx_add23 = -z_xx_add13; z_yy_add23 = -z_yy_add13; z_xy_add23 = -z_xy_add13;
		for (int ii = i; ii <= i + 1; ii++) {
			for (int jj = j - 1; jj <= j + 1; jj++) {
				z_xo3(ix, jy) = Z_x(ii, jj);
				z_yo3(ix, jy) = Z_y(ii, jj);
				z_xxo3(ix, jy) = Z_xx(ii, jj);
				z_yyo3(ix, jy) = Z_yy(ii, jj);
				z_xyo3(ix, jy) = Z_xy(ii, jj);
				jy++;
			}
			ix++; jy = 0;
		}
		ix = 0; jy = 0;
		z_x_new13 = z_xo3 + z_x_add13; z_y_new13 = z_yo3 + z_y_add13; z_xx_new13 = z_xxo3 + z_xx_add13; z_yy_new13 = z_yyo3 + z_yy_add13; z_xy_new13 = z_xyo3 + z_xy_add13;
		z_x_new23 = z_xo3 + z_x_add23; z_y_new23 = z_yo3 + z_y_add23; z_xx_new23 = z_xxo3 + z_xx_add23; z_yy_new23 = z_yyo3 + z_yy_add23; z_xy_new23 = z_xyo3 + z_xy_add23;
		for (int ii = i; ii <= i + 1; ii++) {
			for (int jj = j - 1; jj <= j + 1; jj++) {
				Alpha_new3(ix, jy) = A(ii, jj);
				Beta_new3(ix, jy) = B(ii, jj);
				P_new3(ix, jy) = S(ii, jj);
				jy++;
			}
			ix++; jy = 0;
		}
		ix = 0; jy = 0;
		I_u_1 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new13, z_y_new13, z_xx_new13, z_yy_new13, z_xy_new13);
		I_u_2 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new23, z_y_new23, z_xx_new23, z_yy_new23, z_xy_new23);
		deriva_I_Z(j) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(i, j));
	}

	for (int j = 3; j <= Z_size2 - 3; j++) {
		int i = Z_size1;
		Zeros(z_x_add13);
		Zeros(z_y_add13);
		Zeros(z_xx_add13);
		Zeros(z_xy_add13);
		Zeros(z_yy_add13);
		z_x_add13(1, 0) = DER_INCREMENTS_Z(i, j) / 2; z_x_add13(1, 2) = -DER_INCREMENTS_Z(i, j) / 2;
		z_y_add13(0, 1) = DER_INCREMENTS_Z(i, j) / 2; z_y_add13(1, 1) = DER_INCREMENTS_Z(i, j);
		z_xx_add13(1, 0) = DER_INCREMENTS_Z(i, j); z_xx_add13(1, 1) = -2 * DER_INCREMENTS_Z(i, j); z_xx_add13(1, 2) = DER_INCREMENTS_Z(i, j);
		z_yy_add13(0, 1) = DER_INCREMENTS_Z(i, j); z_yy_add13(1, 1) = DER_INCREMENTS_Z(i, j);
		z_xy_add13(0, 0) = DER_INCREMENTS_Z(i, j) / 4; z_xy_add13(0, 2) = -DER_INCREMENTS_Z(i, j) / 4;
		z_xy_add13(1, 0) = DER_INCREMENTS_Z(i, j) / 4; z_xy_add13(1, 2) = -DER_INCREMENTS_Z(i, j) / 4;
		z_x_add23 = -z_x_add13; z_y_add23 = -z_y_add13; z_xx_add23 = -z_xx_add13; z_yy_add23 = -z_yy_add13; z_xy_add23 = -z_xy_add13;
		for (int ii = i - 1; ii <= i; ii++) {
			for (int jj = j - 1; jj <= j + 1; jj++) {
				z_xo3(ix, jy) = Z_x(ii, jj);
				z_yo3(ix, jy) = Z_y(ii, jj);
				z_xxo3(ix, jy) = Z_xx(ii, jj);
				z_yyo3(ix, jy) = Z_yy(ii, jj);
				z_xyo3(ix, jy) = Z_xy(ii, jj);
				jy++;
			}
			ix++; jy = 0;
		}
		ix = 0; jy = 0;
		z_x_new13 = z_xo3 + z_x_add13; z_y_new13 = z_yo3 + z_y_add13; z_xx_new13 = z_xxo3 + z_xx_add13; z_yy_new13 = z_yyo3 + z_yy_add13; z_xy_new13 = z_xyo3 + z_xy_add13;
		z_x_new23 = z_xo3 + z_x_add23; z_y_new23 = z_yo3 + z_y_add23; z_xx_new23 = z_xxo3 + z_xx_add23; z_yy_new23 = z_yyo3 + z_yy_add23; z_xy_new23 = z_xyo3 + z_xy_add23;
		for (int ii = i - 1; ii <= i; ii++) {
			for (int jj = j - 1; jj <= j + 1; jj++) {
				Alpha_new3(ix, jy) = A(ii, jj);
				Beta_new3(ix, jy) = B(ii, jj);
				P_new3(ix, jy) = S(ii, jj);
				jy++;
			}
			ix++; jy = 0;
		}
		ix = 0; jy = 0;
		I_u_1 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new13, z_y_new13, z_xx_new13, z_yy_new13, z_xy_new13);
		I_u_2 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new23, z_y_new23, z_xx_new23, z_yy_new23, z_xy_new23);
		deriva_I_Z(i * (Z_size2+1) + j) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(i, j));
	}
}
void Adam::calculate_A24() {
	Matrix<double, 2, 4> z_x_add14, z_y_add14, z_xx_add14, z_yy_add14, z_xy_add14;
	Matrix<double, 2, 4> z_x_add24, z_y_add24, z_xx_add24, z_yy_add24, z_xy_add24;
	Matrix<double, 2, 4> z_xo4, z_yo4, z_xxo4, z_yyo4, z_xyo4;
	Matrix<double, 2, 4> z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14;
	Matrix<double, 2, 4> z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24;
	Matrix<double, 2, 4> Alpha_new4, Beta_new4, P_new4;
	int ix = 0, jy = 0;
	Zeros(z_x_add14);
	Zeros(z_y_add14);
	Zeros(z_xx_add14);
	Zeros(z_xy_add14);
	Zeros(z_yy_add14);
	z_x_add14(0, 1) = DER_INCREMENTS_Z(0, 2) / 2; z_x_add14(0, 3) = -DER_INCREMENTS_Z(0, 2) / 2;
	z_y_add14(0, 2) = -DER_INCREMENTS_Z(0, 2); z_y_add14(1, 2) = -DER_INCREMENTS_Z(0, 2) / 2;
	z_xx_add14(0, 0) = DER_INCREMENTS_Z(0, 2); z_xx_add14(0, 1) = DER_INCREMENTS_Z(0, 2); 
	z_xx_add14(0, 2) = -2 * DER_INCREMENTS_Z(0, 2); z_xx_add14(0, 3) = DER_INCREMENTS_Z(0, 2);
	z_yy_add14(0, 2) = DER_INCREMENTS_Z(0, 2); z_yy_add14(1, 2) = DER_INCREMENTS_Z(0, 2);
	z_xy_add14(0, 0) = -DER_INCREMENTS_Z(0, 2) / 4; z_xy_add14(0, 1) = -DER_INCREMENTS_Z(0, 2) / 4; 
	z_xy_add14(0, 3) = DER_INCREMENTS_Z(0, 2) / 4; z_xy_add14(1, 0) = -DER_INCREMENTS_Z(0, 2) / 4; 
	z_xy_add14(1, 1) = -DER_INCREMENTS_Z(0, 2) / 4; z_xy_add14(1, 3) = DER_INCREMENTS_Z(0, 2) / 4;
	z_x_add24 = -z_x_add14; z_y_add24 = -z_y_add14; z_xx_add24 = -z_xx_add14; z_yy_add24 = -z_yy_add14; z_xy_add24 = -z_xy_add14;

	for (int i = 0; i <= 1; i++)
		for (int j = 0; j <= 3; j++) {
			z_xo4(i, j) = Z_x(i, j);
			z_yo4(i, j) = Z_y(i, j);
			z_xxo4(i, j) = Z_xx(i, j);
			z_yyo4(i, j) = Z_yy(i, j);
			z_xyo4(i, j) = Z_xy(i, j);
		}
	z_x_new14 = z_xo4 + z_x_add14; z_y_new14 = z_yo4 + z_y_add14; z_xx_new14 = z_xxo4 + z_xx_add14; z_yy_new14 = z_yyo4 + z_yy_add14; z_xy_new14 = z_xyo4 + z_xy_add14;
	z_x_new24 = z_xo4 + z_x_add24; z_y_new24 = z_yo4 + z_y_add24; z_xx_new24 = z_xxo4 + z_xx_add24; z_yy_new24 = z_yyo4 + z_yy_add24; z_xy_new24 = z_xyo4 + z_xy_add24;

	for (int i = 0; i <= 1; i++)
		for (int j = 0; j <= 3; j++) {
			Alpha_new4(i, j) = A(i, j);
			Beta_new4(i, j) = B(i, j);
			P_new4(i, j) = S(i, j);
		}

	double I_u_1 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14);
	double I_u_2 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24);
	deriva_I_Z(2) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(0, 2));
	//=============================================================================================================================//
	Zeros(z_x_add14);
	Zeros(z_y_add14);
	Zeros(z_xx_add14);
	Zeros(z_xy_add14);
	Zeros(z_yy_add14);
	z_x_add14(0, 0) = DER_INCREMENTS_Z(0, Z_size2 - 2) / 2; z_x_add14(0, 2) = -DER_INCREMENTS_Z(0, Z_size2 - 2) / 2;
	z_y_add14(0, 1) = -DER_INCREMENTS_Z(0, Z_size2 - 2); z_y_add14(1, 1) = -DER_INCREMENTS_Z(0, Z_size2 - 2) / 2;
	z_xx_add14(0, 0) = DER_INCREMENTS_Z(0, Z_size2 - 2); z_xx_add14(0, 1) = -2 * DER_INCREMENTS_Z(0, Z_size2 - 2); 
	z_xx_add14(0, 2) = DER_INCREMENTS_Z(0, Z_size2 - 2); z_xx_add14(0, 3) = DER_INCREMENTS_Z(0, Z_size2 - 2);
	z_yy_add14(0, 1) = DER_INCREMENTS_Z(0, Z_size2 - 2); z_yy_add14(1, 1) = DER_INCREMENTS_Z(0, Z_size2 - 2);
	z_xy_add14(0, 0) = -DER_INCREMENTS_Z(0, Z_size2 - 2) / 4; z_xy_add14(0, 2) = DER_INCREMENTS_Z(0, Z_size2 - 2) / 4; z_xy_add14(0, 3) = DER_INCREMENTS_Z(0, Z_size2 - 2) / 4;
	z_xy_add14(1, 0) = -DER_INCREMENTS_Z(0, Z_size2 - 2) / 4; z_xy_add14(1, 2) = DER_INCREMENTS_Z(0, Z_size2 - 2) / 4; z_xy_add14(1, 3) = DER_INCREMENTS_Z(0, Z_size2 - 2) / 4;
	z_x_add24 = -z_x_add14; z_y_add24 = -z_y_add14; z_xx_add24 = -z_xx_add14; z_yy_add24 = -z_yy_add14; z_xy_add24 = -z_xy_add14;
	
	for (int i = 0; i <= 1; i++) {
		for (int j = Z_size2 - 3; j <= Z_size2; j++) {
			z_xo4(i, jy) = Z_x(i, j);
			z_yo4(i, jy) = Z_y(i, j);
			z_xxo4(i, jy) = Z_xx(i, j);
			z_yyo4(i, jy) = Z_yy(i, j);
			z_xyo4(i, jy) = Z_xy(i, j);
			jy++;
		}
		jy = 0;
	}
	jy = 0;
	z_x_new14 = z_xo4 + z_x_add14; z_y_new14 = z_yo4 + z_y_add14; z_xx_new14 = z_xxo4 + z_xx_add14; z_yy_new14 = z_yyo4 + z_yy_add14; z_xy_new14 = z_xyo4 + z_xy_add14;
	z_x_new24 = z_xo4 + z_x_add24; z_y_new24 = z_yo4 + z_y_add24; z_xx_new24 = z_xxo4 + z_xx_add24; z_yy_new24 = z_yyo4 + z_yy_add24; z_xy_new24 = z_xyo4 + z_xy_add24;
	for (int i = 0; i <= 1; i++) {
		for (int j = Z_size2 - 3; j <= Z_size2; j++) {
			Alpha_new4(i, jy) = A(i, j);
			Beta_new4(i, jy) = B(i, j);
			P_new4(i, jy) = S(i, j);
			jy++;
		}
		jy = 0;
	}
	jy = 0;
	I_u_1 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14);
	I_u_2 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24);
	deriva_I_Z(Z_size2 - 2) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(0, Z_size2 - 2));
	//=============================================================================================================================//
	Zeros(z_x_add14);
	Zeros(z_y_add14);
	Zeros(z_xx_add14);
	Zeros(z_xy_add14);
	Zeros(z_yy_add14);
	z_x_add14(1, 1) = DER_INCREMENTS_Z(Z_size1, 2) / 2; z_x_add14(1, 3) = -DER_INCREMENTS_Z(Z_size1, 2) / 2;
	z_y_add14(0, 2) = DER_INCREMENTS_Z(Z_size1, 2) / 2; z_y_add14(1, 2) = DER_INCREMENTS_Z(Z_size1, 2);
	z_xx_add14(1, 0) = DER_INCREMENTS_Z(Z_size1, 2); z_xx_add14(1, 1) = DER_INCREMENTS_Z(Z_size1, 2); z_xx_add14(1, 2) = -2 * DER_INCREMENTS_Z(Z_size1, 2); z_xx_add14(1, 3) = DER_INCREMENTS_Z(Z_size1, 2);
	z_yy_add14(0, 2) = DER_INCREMENTS_Z(Z_size1, 2); z_yy_add14(1, 2) = DER_INCREMENTS_Z(Z_size1, 2);
	z_xy_add14(0, 0) = DER_INCREMENTS_Z(Z_size1, 2) / 4; z_xy_add14(0, 1) = DER_INCREMENTS_Z(Z_size1, 2) / 4; z_xy_add14(0, 3) = -DER_INCREMENTS_Z(Z_size1, 2) / 4;
	z_xy_add14(1, 0) = DER_INCREMENTS_Z(Z_size1, 2) / 4; z_xy_add14(1, 1) = DER_INCREMENTS_Z(Z_size1, 2) / 4; z_xy_add14(1, 3) = -DER_INCREMENTS_Z(Z_size1, 2) / 4;
	z_x_add24 = -z_x_add14; z_y_add24 = -z_y_add14; z_xx_add24 = -z_xx_add14; z_yy_add24 = -z_yy_add14; z_xy_add24 = -z_xy_add14;
	for (int i = Z_size1 - 1; i <= Z_size1; i++) {
		for (int j = 0; j <= 3; j++) {
			z_xo4(ix, j) = Z_x(i, j);
			z_yo4(ix, j) = Z_y(i, j);
			z_xxo4(ix, j) = Z_xx(i, j);
			z_yyo4(ix, j) = Z_yy(i, j);
			z_xyo4(ix, j) = Z_xy(i, j);
		}
		ix++;
	}
	ix = 0;
	z_x_new14 = z_xo4 + z_x_add14; z_y_new14 = z_yo4 + z_y_add14; z_xx_new14 = z_xxo4 + z_xx_add14; z_yy_new14 = z_yyo4 + z_yy_add14; z_xy_new14 = z_xyo4 + z_xy_add14;
	z_x_new24 = z_xo4 + z_x_add24; z_y_new24 = z_yo4 + z_y_add24; z_xx_new24 = z_xxo4 + z_xx_add24; z_yy_new24 = z_yyo4 + z_yy_add24; z_xy_new24 = z_xyo4 + z_xy_add24;
	for (int i = Z_size1 - 1; i <= Z_size1; i++) {
		for (int j = 0; j <= 3; j++) {
			Alpha_new4(ix, j) = A(i, j);
			Beta_new4(ix, j) = B(i, j);
			P_new4(ix, j) = S(i, j);
		}
		ix++;
	}
	ix = 0;
	I_u_1 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14);
	I_u_2 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24);
	deriva_I_Z((Z_size1) * (Z_size2 + 1) + 2) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(Z_size1, 2));
	//=============================================================================================================================//
	Zeros(z_x_add14);
	Zeros(z_y_add14);
	Zeros(z_xx_add14);
	Zeros(z_xy_add14);
	Zeros(z_yy_add14);
	z_x_add14(1, 0) = DER_INCREMENTS_Z(Z_size1, Z_size2 - 2) / 2; z_x_add14(1, 2) = -DER_INCREMENTS_Z(Z_size1, Z_size2 - 2) / 2;
	z_y_add14(0, 1) = DER_INCREMENTS_Z(Z_size1, Z_size2 - 2) / 2; z_y_add14(1, 1) = DER_INCREMENTS_Z(Z_size1, Z_size2 - 2);
	z_xx_add14(1, 0) = DER_INCREMENTS_Z(Z_size1, Z_size2 - 2); z_xx_add14(1, 1) = -2 * DER_INCREMENTS_Z(Z_size1, Z_size2 - 2); 
	z_xx_add14(1, 2) = DER_INCREMENTS_Z(Z_size1, Z_size2 - 2); z_xx_add14(1, 3) = DER_INCREMENTS_Z(Z_size1, Z_size2 - 2);
	z_yy_add14(0, 1) = DER_INCREMENTS_Z(Z_size1, Z_size2 - 2); z_yy_add14(1, 1) = DER_INCREMENTS_Z(Z_size1, Z_size2 - 2);
	z_xy_add14(0, 0) = DER_INCREMENTS_Z(Z_size1, Z_size2 - 2) / 4; z_xy_add14(0, 2) = -DER_INCREMENTS_Z(Z_size1, Z_size2 - 2) / 4; z_xy_add14(0, 3) = -DER_INCREMENTS_Z(Z_size1, Z_size2 - 2) / 4;
	z_xy_add14(1, 0) = DER_INCREMENTS_Z(Z_size1, Z_size2 - 2) / 4; z_xy_add14(1, 2) = -DER_INCREMENTS_Z(Z_size1, Z_size2 - 2) / 4; z_xy_add14(1, 3) = -DER_INCREMENTS_Z(Z_size1, Z_size2 - 2) / 4;
	z_x_add24 = -z_x_add14; z_y_add24 = -z_y_add14; z_xx_add24 = -z_xx_add14; z_yy_add24 = -z_yy_add14; z_xy_add24 = -z_xy_add14;
	//z_xo = z_x0(Z_size1 - 1:Z_size1, Z_size2 - 3 : Z_size2); z_yo = z_y0(Z_size1 - 1:Z_size1, Z_size2 - 3 : Z_size2); z_xxo = z_xx0(Z_size1 - 1:Z_size1, Z_size2 - 3 : Z_size2); z_yyo = z_yy0(Z_size1 - 1:Z_size1, Z_size2 - 3 : Z_size2); z_xyo = z_xy0(Z_size1 - 1:Z_size1, Z_size2 - 3 : Z_size2);
	for (int i = Z_size1 - 1; i <= Z_size1; i++) {
		for (int j = Z_size2 - 3; j <= Z_size2; j++) {
			z_xo4(ix, jy) = Z_x(i, j);
			z_yo4(ix, jy) = Z_y(i, j);
			z_xxo4(ix, jy) = Z_xx(i, j);
			z_yyo4(ix, jy) = Z_yy(i, j);
			z_xyo4(ix, jy) = Z_xy(i, j);
			jy++;
		}
		ix++; jy = 0;
	}
	ix = 0; jy = 0;
	z_x_new14 = z_xo4 + z_x_add14; z_y_new14 = z_yo4 + z_y_add14; z_xx_new14 = z_xxo4 + z_xx_add14; z_yy_new14 = z_yyo4 + z_yy_add14; z_xy_new14 = z_xyo4 + z_xy_add14;
	z_x_new24 = z_xo4 + z_x_add24; z_y_new24 = z_yo4 + z_y_add24; z_xx_new24 = z_xxo4 + z_xx_add24; z_yy_new24 = z_yyo4 + z_yy_add24; z_xy_new24 = z_xyo4 + z_xy_add24;
	for (int i = Z_size1 - 1; i <= Z_size1; i++) {
		for (int j = Z_size2 - 3; j <= Z_size2; j++) {
			Alpha_new4(ix, jy) = A(i, j);
			Beta_new4(ix, jy) = B(i, j);
			P_new4(ix, jy) = S(i, j);
			jy++;
		}
		ix++; jy = 0;
	}
	ix = 0; jy = 0;
	I_u_1 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14);
	I_u_2 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24);
	deriva_I_Z((Z_size1) * (Z_size2+1) + Z_size2 - 2) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(Z_size1, Z_size2 - 2));
}
void Adam::calculate_A32() {
	Matrix<double, 3, 2> z_x_add13, z_y_add13, z_xx_add13, z_yy_add13, z_xy_add13;
	Matrix<double, 3, 2> z_x_add23, z_y_add23, z_xx_add23, z_yy_add23, z_xy_add23;
	Matrix<double, 3, 2> z_xo3, z_yo3, z_xxo3, z_yyo3, z_xyo3;
	Matrix<double, 3, 2> z_x_new13, z_y_new13, z_xx_new13, z_yy_new13, z_xy_new13;
	Matrix<double, 3, 2> z_x_new23, z_y_new23, z_xx_new23, z_yy_new23, z_xy_new23;
	Matrix<double, 3, 2> Alpha_new3, Beta_new3, P_new3;
	int ix = 0, jy = 0;
	Zeros(z_x_add13);
	Zeros(z_y_add13);
	Zeros(z_xx_add13);
	Zeros(z_xy_add13);
	Zeros(z_yy_add13);
	z_x_add13(1, 0) = -DER_INCREMENTS_Z(1, 0); z_x_add13(1, 1) = -DER_INCREMENTS_Z(1, 0) / 2;
	z_y_add13(0, 0) = DER_INCREMENTS_Z(1, 0); z_y_add13(2, 0) = -DER_INCREMENTS_Z(1, 0) / 2;
	z_xx_add13(1, 0) = DER_INCREMENTS_Z(1, 0); z_xx_add13(1, 1) = DER_INCREMENTS_Z(1, 0);
	z_yy_add13(0, 0) = -2 * DER_INCREMENTS_Z(1, 0); z_yy_add13(1, 0) = -2 * DER_INCREMENTS_Z(1, 0); z_yy_add13(2, 0) = DER_INCREMENTS_Z(1, 0);
	z_xy_add13(2, 0) = DER_INCREMENTS_Z(1, 0) / 4; z_xy_add13(2, 1) = DER_INCREMENTS_Z(1, 0) / 4;
	z_x_add23 = -z_x_add13; z_y_add23 = -z_y_add13; z_xx_add23 = -z_xx_add13; z_yy_add23 = -z_yy_add13; z_xy_add23 = -z_xy_add13;
	for (int i = 0; i <= 2; i++)
		for (int j = 0; j <= 1; j++) {
			z_xo3(i, j) = Z_x(i, j);
			z_yo3(i, j) = Z_y(i, j);
			z_xxo3(i, j) = Z_xx(i, j);
			z_yyo3(i, j) = Z_yy(i, j);
			z_xyo3(i, j) = Z_xy(i, j);
		}
	z_x_new13 = z_xo3 + z_x_add13; z_y_new13 = z_yo3 + z_y_add13; z_xx_new13 = z_xxo3 + z_xx_add13; z_yy_new13 = z_yyo3 + z_yy_add13; z_xy_new13 = z_xyo3 + z_xy_add13;
	z_x_new23 = z_xo3 + z_x_add23; z_y_new23 = z_yo3 + z_y_add23; z_xx_new23 = z_xxo3 + z_xx_add23; z_yy_new23 = z_yyo3 + z_yy_add23; z_xy_new23 = z_xyo3 + z_xy_add23;
	for (int i = 0; i <= 2; i++)
		for (int j = 0; j <= 1; j++) {
			Alpha_new3(i, j) = A(i, j);
			Beta_new3(i, j) = B(i, j);
			P_new3(i, j) = S(i, j);
		}
	double I_u_1 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new13, z_y_new13, z_xx_new13, z_yy_new13, z_xy_new13);
	double I_u_2 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new23, z_y_new23, z_xx_new23, z_yy_new23, z_xy_new23);
	deriva_I_Z(Z_size2 + 1) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(1, 0));
	//=============================================================================================================================//
	Zeros(z_x_add13);
	Zeros(z_y_add13);
	Zeros(z_xx_add13);
	Zeros(z_xy_add13);
	Zeros(z_yy_add13);
	z_x_add13(1, 0) = -DER_INCREMENTS_Z(Z_size1 - 1, 0); z_x_add13(1, 1) = -DER_INCREMENTS_Z(Z_size1 - 1, 0) / 2;
	z_y_add13(0, 0) = DER_INCREMENTS_Z(Z_size1 - 1, 0) / 2; z_y_add13(2, 0) = -DER_INCREMENTS_Z(Z_size1 - 1, 0);
	z_xx_add13(1, 0) = DER_INCREMENTS_Z(Z_size1 - 1, 0); z_xx_add13(1, 1) = DER_INCREMENTS_Z(Z_size1 - 1, 0);
	z_yy_add13(0, 0) = DER_INCREMENTS_Z(Z_size1 - 1, 0); z_yy_add13(1, 0) = -2 * DER_INCREMENTS_Z(Z_size1 - 1, 0); z_yy_add13(2, 0) = -2 * DER_INCREMENTS_Z(Z_size1 - 1, 0);
	z_xy_add13(0, 0) = -DER_INCREMENTS_Z(Z_size1 - 1, 0) / 4; z_xy_add13(0, 1) = -DER_INCREMENTS_Z(Z_size1 - 1, 0) / 4;
	z_x_add23 = -z_x_add13; z_y_add23 = -z_y_add13; z_xx_add23 = -z_xx_add13; z_yy_add23 = -z_yy_add13; z_xy_add23 = -z_xy_add13;
	//z_xo = z_x0(Z_size1 - 2:Z_size1, 1 : 2); z_yo = z_y0(Z_size1 - 2:Z_size1, 1 : 2); z_xxo = z_xx0(Z_size1 - 2:Z_size1, 1 : 2); z_yyo = z_yy0(Z_size1 - 2:Z_size1, 1 : 2); z_xyo = z_xy0(Z_size1 - 2:Z_size1, 1 : 2);
	for (int i = Z_size1 - 2; i <= Z_size1; i++) {
		for (int j = 0; j <= 1; j++) {
			z_xo3(ix, j) = Z_x(i, j);
			z_yo3(ix, j) = Z_y(i, j);
			z_xxo3(ix, j) = Z_xx(i, j);
			z_yyo3(ix, j) = Z_yy(i, j);
			z_xyo3(ix, j) = Z_xy(i, j);
		}
		ix++;
	}
	ix = 0;
	z_x_new13 = z_xo3 + z_x_add13; z_y_new13 = z_yo3 + z_y_add13; z_xx_new13 = z_xxo3 + z_xx_add13; z_yy_new13 = z_yyo3 + z_yy_add13; z_xy_new13 = z_xyo3 + z_xy_add13;
	z_x_new23 = z_xo3 + z_x_add23; z_y_new23 = z_yo3 + z_y_add23; z_xx_new23 = z_xxo3 + z_xx_add23; z_yy_new23 = z_yyo3 + z_yy_add23; z_xy_new23 = z_xyo3 + z_xy_add23;
	for (int i = Z_size1 - 2; i <= Z_size1; i++) {
		for (int j = 0; j <= 1; j++) {
			Alpha_new3(ix, j) = A(i, j);
			Beta_new3(ix, j) = B(i, j);
			P_new3(ix, j) = S(i, j);
		}
		ix++;
	}
	ix = 0;
	I_u_1 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new13, z_y_new13, z_xx_new13, z_yy_new13, z_xy_new13);
	I_u_2 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new23, z_y_new23, z_xx_new23, z_yy_new23, z_xy_new23);
	deriva_I_Z((Z_size1 - 1) * (Z_size2+1)) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(Z_size1 - 1, 0));
	//=============================================================================================================================//
	Zeros(z_x_add13);
	Zeros(z_y_add13);
	Zeros(z_xx_add13);
	Zeros(z_xy_add13);
	Zeros(z_yy_add13);
	z_x_add13(1, 0) = DER_INCREMENTS_Z(1, Z_size2) / 2; z_x_add13(1, 1) = DER_INCREMENTS_Z(1, Z_size2);
	z_y_add13(0, 1) = DER_INCREMENTS_Z(1, Z_size2); z_y_add13(2, 1) = -DER_INCREMENTS_Z(1, Z_size2) / 2;
	z_xx_add13(1, 0) = DER_INCREMENTS_Z(1, Z_size2); z_xx_add13(1, 1) = DER_INCREMENTS_Z(1, Z_size2);
	z_yy_add13(0, 1) = -2 * DER_INCREMENTS_Z(1, Z_size2); z_yy_add13(1, 1) = -2 * DER_INCREMENTS_Z(1, Z_size2); z_yy_add13(2, 1) = DER_INCREMENTS_Z(1, Z_size2);
	z_xy_add13(2, 0) = -DER_INCREMENTS_Z(1, Z_size2) / 4; z_xy_add13(2, 1) = -DER_INCREMENTS_Z(1, Z_size2) / 4;
	z_x_add23 = -z_x_add13; z_y_add23 = -z_y_add13; z_xx_add23 = -z_xx_add13; z_yy_add23 = -z_yy_add13; z_xy_add23 = -z_xy_add13;
	for (int i = 0; i <= 2; i++) {
		for (int j = Z_size2 - 1; j <= Z_size2; j++) {
			z_xo3(i, jy) = Z_x(i, j);
			z_yo3(i, jy) = Z_y(i, j);
			z_xxo3(i, jy) = Z_xx(i, j);
			z_yyo3(i, jy) = Z_yy(i, j);
			z_xyo3(i, jy) = Z_xy(i, j);
			jy++;
		}
		jy = 0;
	}
	jy = 0;
	z_x_new13 = z_xo3 + z_x_add13; z_y_new13 = z_yo3 + z_y_add13; z_xx_new13 = z_xxo3 + z_xx_add13; z_yy_new13 = z_yyo3 + z_yy_add13; z_xy_new13 = z_xyo3 + z_xy_add13;
	z_x_new23 = z_xo3 + z_x_add23; z_y_new23 = z_yo3 + z_y_add23; z_xx_new23 = z_xxo3 + z_xx_add23; z_yy_new23 = z_yyo3 + z_yy_add23; z_xy_new23 = z_xyo3 + z_xy_add23;
	for (int i = 0; i <= 2; i++) {
		for (int j = Z_size2 - 1; j <= Z_size2; j++) {
			Alpha_new3(i, jy) = A(i, j);
			Beta_new3(i, jy) = B(i, j);
			P_new3(i, jy) = S(i, j);
			jy++;
		}
		jy = 0;
	}
	jy = 0;
	I_u_1 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new13, z_y_new13, z_xx_new13, z_yy_new13, z_xy_new13);
	I_u_2 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new23, z_y_new23, z_xx_new23, z_yy_new23, z_xy_new23);
	deriva_I_Z(2 * (Z_size2 + 1) - 1) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(2, Z_size2));
	//=============================================================================================================================//
	Zeros(z_x_add13);
	Zeros(z_y_add13);
	Zeros(z_xx_add13);
	Zeros(z_xy_add13);
	Zeros(z_yy_add13);
	z_x_add13(1, 0) = DER_INCREMENTS_Z(Z_size1 - 1, Z_size2) / 2; z_x_add13(1, 1) = DER_INCREMENTS_Z(Z_size1 - 1, Z_size2);
	z_y_add13(0, 1) = DER_INCREMENTS_Z(Z_size1 - 1, Z_size2) / 2; z_y_add13(2, 1) = -DER_INCREMENTS_Z(Z_size1 - 1, Z_size2);
	z_xx_add13(1, 0) = DER_INCREMENTS_Z(Z_size1 - 1, Z_size2); z_xx_add13(1, 1) = DER_INCREMENTS_Z(Z_size1 - 1, Z_size2);
	z_yy_add13(0, 1) = DER_INCREMENTS_Z(Z_size1 - 1, Z_size2); z_yy_add13(1, 1) = -2 * DER_INCREMENTS_Z(Z_size1 - 1, Z_size2); z_yy_add13(2, 1) = -2 * DER_INCREMENTS_Z(Z_size1 - 1, Z_size2);
	z_xy_add13(0, 0) = DER_INCREMENTS_Z(Z_size1 - 1, Z_size2) / 4; z_xy_add13(0, 1) = DER_INCREMENTS_Z(Z_size1 - 1, Z_size2) / 4;
	z_x_add23 = -z_x_add13; z_y_add23 = -z_y_add13; z_xx_add23 = -z_xx_add13; z_yy_add23 = -z_yy_add13; z_xy_add23 = -z_xy_add13;
	//z_xo = z_x0(Z_size1 - 2:Z_size1, Z_size2 - 1 : Z_size2); z_yo = z_y0(Z_size1 - 2:Z_size1, Z_size2 - 1 : Z_size2); z_xxo = z_xx0(Z_size1 - 2:Z_size1, Z_size2 - 1 : Z_size2); z_yyo = z_yy0(Z_size1 - 2:Z_size1, Z_size2 - 1 : Z_size2); z_xyo = z_xy0(Z_size1 - 2:Z_size1, Z_size2 - 1 : Z_size2);
	for (int i = Z_size1 - 2; i <= Z_size1; i++) {
		for (int j = Z_size2 - 1; j <= Z_size2; j++) {
			z_xo3(ix, jy) = Z_x(i, j);
			z_yo3(ix, jy) = Z_y(i, j);
			z_xxo3(ix, jy) = Z_xx(i, j);
			z_yyo3(ix, jy) = Z_yy(i, j);
			z_xyo3(ix, jy) = Z_xy(i, j);
			jy++;
		}
		ix++; jy = 0;
	}
	ix = 0; jy = 0;
	z_x_new13 = z_xo3 + z_x_add13; z_y_new13 = z_yo3 + z_y_add13; z_xx_new13 = z_xxo3 + z_xx_add13; z_yy_new13 = z_yyo3 + z_yy_add13; z_xy_new13 = z_xyo3 + z_xy_add13;
	z_x_new23 = z_xo3 + z_x_add23; z_y_new23 = z_yo3 + z_y_add23; z_xx_new23 = z_xxo3 + z_xx_add23; z_yy_new23 = z_yyo3 + z_yy_add23; z_xy_new23 = z_xyo3 + z_xy_add23;
	for (int i = Z_size1 - 2; i <= Z_size1; i++) {
		for (int j = Z_size2 - 1; j <= Z_size2; j++) {
			Alpha_new3(ix, jy) = A(i, j);
			Beta_new3(ix, jy) = B(i, j);
			P_new3(ix, jy) = S(i, j);
			jy++;
		}
		ix++; jy = 0;
	}
	ix = 0; jy = 0;
	I_u_1 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new13, z_y_new13, z_xx_new13, z_yy_new13, z_xy_new13);
	I_u_2 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new23, z_y_new23, z_xx_new23, z_yy_new23, z_xy_new23);
	deriva_I_Z((Z_size1) * (Z_size2+1)-1) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(Z_size1 - 1, Z_size2));
	//=============================================================================================================================//
	for (int i = 3; i <= Z_size1 - 3; i++) {
		int j = 0;
		Zeros(z_x_add13);
		Zeros(z_y_add13);
		Zeros(z_xx_add13);
		Zeros(z_xy_add13);
		Zeros(z_yy_add13);
		z_x_add13(1, 0) = -DER_INCREMENTS_Z(i, j); z_x_add13(1, 1) = -DER_INCREMENTS_Z(i, j) / 2;
		z_y_add13(0, 0) = DER_INCREMENTS_Z(i, j) / 2; z_y_add13(2, 0) = -DER_INCREMENTS_Z(i, j) / 2;
		z_xx_add13(1, 0) = DER_INCREMENTS_Z(i, j); z_xx_add13(1, 1) = DER_INCREMENTS_Z(i, j);
		z_yy_add13(0, 0) = DER_INCREMENTS_Z(i, j); z_yy_add13(1, 0) = -2 * DER_INCREMENTS_Z(i, j); z_yy_add13(2, 0) = DER_INCREMENTS_Z(i, j);
		z_xy_add13(0, 0) = -DER_INCREMENTS_Z(i, j) / 4; z_xy_add13(0, 1) = -DER_INCREMENTS_Z(i, j) / 4;
		z_xy_add13(2, 0) = DER_INCREMENTS_Z(i, j) / 4; z_xy_add13(2, 1) = DER_INCREMENTS_Z(i, j) / 4;
		z_x_add23 = -z_x_add13; z_y_add23 = -z_y_add13; z_xx_add23 = -z_xx_add13; z_yy_add23 = -z_yy_add13; z_xy_add23 = -z_xy_add13;
		for (int ii = i - 1; ii <= i + 1; ii++) {
			for (int jj = j; jj <= j + 1; jj++) {
				z_xo3(ix, jy) = Z_x(ii, jj);
				z_yo3(ix, jy) = Z_y(ii, jj);
				z_xxo3(ix, jy) = Z_xx(ii, jj);
				z_yyo3(ix, jy) = Z_yy(ii, jj);
				z_xyo3(ix, jy) = Z_xy(ii, jj);
				jy++;
			}
			ix++; jy = 0;
		}
		ix = 0; jy = 0;
		z_x_new13 = z_xo3 + z_x_add13; z_y_new13 = z_yo3 + z_y_add13; z_xx_new13 = z_xxo3 + z_xx_add13; z_yy_new13 = z_yyo3 + z_yy_add13; z_xy_new13 = z_xyo3 + z_xy_add13;
		z_x_new23 = z_xo3 + z_x_add23; z_y_new23 = z_yo3 + z_y_add23; z_xx_new23 = z_xxo3 + z_xx_add23; z_yy_new23 = z_yyo3 + z_yy_add23; z_xy_new23 = z_xyo3 + z_xy_add23;
		for (int ii = i - 1; ii <= i + 1; ii++) {
			for (int jj = j; jj <= j + 1; jj++) {
				Alpha_new3(ix, jy) = A(ii, jj);
				Beta_new3(ix, jy) = B(ii, jj);
				P_new3(ix, jy) = S(ii, jj);
				jy++;
			}
			ix++; jy = 0;
		}
		ix = 0; jy = 0;
		I_u_1 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new13, z_y_new13, z_xx_new13, z_yy_new13, z_xy_new13);
		I_u_2 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new23, z_y_new23, z_xx_new23, z_yy_new23, z_xy_new23);
		deriva_I_Z((i) * (Z_size2+1) + j) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(i, j));
	}
	//=============================================================================================================================//
	for (int i = 3; i <= Z_size1 - 3; i++) {
		int j = Z_size2;
		Zeros(z_x_add13);
		Zeros(z_y_add13);
		Zeros(z_xx_add13);
		Zeros(z_xy_add13);
		Zeros(z_yy_add13);
		z_x_add13(1, 0) = DER_INCREMENTS_Z(i, j) / 2; z_x_add13(1, 1) = DER_INCREMENTS_Z(i, j);
		z_y_add13(0, 1) = DER_INCREMENTS_Z(i, j) / 2; z_y_add13(2, 1) = -DER_INCREMENTS_Z(i, j) / 2;
		z_xx_add13(1, 0) = DER_INCREMENTS_Z(i, j); z_xx_add13(1, 1) = DER_INCREMENTS_Z(i, j);
		z_yy_add13(0, 1) = DER_INCREMENTS_Z(i, j); z_yy_add13(1, 1) = -2 * DER_INCREMENTS_Z(i, j); z_yy_add13(2, 1) = DER_INCREMENTS_Z(i, j);
		z_xy_add13(0, 0) = DER_INCREMENTS_Z(i, j) / 4; z_xy_add13(0, 1) = DER_INCREMENTS_Z(i, j) / 4;
		z_xy_add13(2, 0) = -DER_INCREMENTS_Z(i, j) / 4; z_xy_add13(2, 1) = -DER_INCREMENTS_Z(i, j) / 4;
		z_x_add23 = -z_x_add13; z_y_add23 = -z_y_add13; z_xx_add23 = -z_xx_add13; z_yy_add23 = -z_yy_add13; z_xy_add23 = -z_xy_add13;
		for (int ii = i - 1; ii <= i + 1; ii++) {
			for (int jj = j - 1; jj <= j; jj++) {
				z_xo3(ix, jy) = Z_x(ii, jj);
				z_yo3(ix, jy) = Z_y(ii, jj);
				z_xxo3(ix, jy) = Z_xx(ii, jj);
				z_yyo3(ix, jy) = Z_yy(ii, jj);
				z_xyo3(ix, jy) = Z_xy(ii, jj);
				jy++;
			}
			ix++; jy = 0;
		}
		ix = 0; jy = 0;
		z_x_new13 = z_xo3 + z_x_add13; z_y_new13 = z_yo3 + z_y_add13; z_xx_new13 = z_xxo3 + z_xx_add13; z_yy_new13 = z_yyo3 + z_yy_add13; z_xy_new13 = z_xyo3 + z_xy_add13;
		z_x_new23 = z_xo3 + z_x_add23; z_y_new23 = z_yo3 + z_y_add23; z_xx_new23 = z_xxo3 + z_xx_add23; z_yy_new23 = z_yyo3 + z_yy_add23; z_xy_new23 = z_xyo3 + z_xy_add23;
		for (int ii = i - 1; ii <= i + 1; ii++) {
			for (int jj = j - 1; jj <= j; jj++) {
				Alpha_new3(ix, jy) = A(ii, jj);
				Beta_new3(ix, jy) = B(ii, jj);
				P_new3(ix, jy) = S(ii, jj);
				jy++;
			}
			ix++; jy = 0;
		}
		ix = 0; jy = 0;
		I_u_1 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new13, z_y_new13, z_xx_new13, z_yy_new13, z_xy_new13);
		I_u_2 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new23, z_y_new23, z_xx_new23, z_yy_new23, z_xy_new23);
		deriva_I_Z((i) * (Z_size2+1) + j) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(i, j));
	}
}
void Adam::calculate_A42() {
	Matrix<double, 4, 2> z_x_add14, z_y_add14, z_xx_add14, z_yy_add14, z_xy_add14;
	Matrix<double, 4, 2> z_x_add24, z_y_add24, z_xx_add24, z_yy_add24, z_xy_add24;
	Matrix<double, 4, 2> z_xo4, z_yo4, z_xxo4, z_yyo4, z_xyo4;
	Matrix<double, 4, 2> z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14;
	Matrix<double, 4, 2> z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24;
	Matrix<double, 4, 2> Alpha_new4, Beta_new4, P_new4;
	int ix = 0, jy = 0;
	Zeros(z_x_add14);
	Zeros(z_y_add14);
	Zeros(z_xx_add14);
	Zeros(z_xy_add14);
	Zeros(z_yy_add14);
	z_x_add14(2, 0) = -DER_INCREMENTS_Z(2, 0); z_x_add14(2, 1) = -DER_INCREMENTS_Z(2, 0) / 2;
	z_y_add14(1, 0) = DER_INCREMENTS_Z(2, 0) / 2; z_y_add14(3, 0) = -DER_INCREMENTS_Z(2, 0) / 2;
	z_xx_add14(2, 0) = DER_INCREMENTS_Z(2, 0); z_xx_add14(2, 1) = DER_INCREMENTS_Z(2, 0);
	z_yy_add14(0, 0) = DER_INCREMENTS_Z(2, 0); z_yy_add14(1, 0) = DER_INCREMENTS_Z(2, 0); 
	z_yy_add14(2, 0) = -2 * DER_INCREMENTS_Z(2, 0); z_yy_add14(3, 0) = DER_INCREMENTS_Z(2, 0);
	z_xy_add14(0, 0) = -DER_INCREMENTS_Z(2, 0) / 4; z_xy_add14(1, 0) = -DER_INCREMENTS_Z(2, 0) / 4; z_xy_add14(3, 0) = DER_INCREMENTS_Z(2, 0) / 4;
	z_xy_add14(0, 1) = -DER_INCREMENTS_Z(2, 0) / 4; z_xy_add14(1, 1) = -DER_INCREMENTS_Z(2, 0) / 4; z_xy_add14(3, 1) = DER_INCREMENTS_Z(2, 0) / 4;
	z_x_add24 = -z_x_add14; z_y_add24 = -z_y_add14; z_xx_add24 = -z_xx_add14; z_yy_add24 = -z_yy_add14; z_xy_add24 = -z_xy_add14;
	for (int i = 0; i <= 3; i++)
		for (int j = 0; j <= 1; j++) {
			z_xo4(i, j) = Z_x(i, j);
			z_yo4(i, j) = Z_y(i, j);
			z_xxo4(i, j) = Z_xx(i, j);
			z_yyo4(i, j) = Z_yy(i, j);
			z_xyo4(i, j) = Z_xy(i, j);
		}
	z_x_new14 = z_xo4 + z_x_add14; z_y_new14 = z_yo4 + z_y_add14; z_xx_new14 = z_xxo4 + z_xx_add14; z_yy_new14 = z_yyo4 + z_yy_add14; z_xy_new14 = z_xyo4 + z_xy_add14;
	z_x_new24 = z_xo4 + z_x_add24; z_y_new24 = z_yo4 + z_y_add24; z_xx_new24 = z_xxo4 + z_xx_add24; z_yy_new24 = z_yyo4 + z_yy_add24; z_xy_new24 = z_xyo4 + z_xy_add24;
	//Alpha_new = Alpha(1:4, 1 : 2); Beta_new = Beta(1:4, 1 : 2); P_new = P(1:4, 1 : 2);
	for (int i = 0; i <= 3; i++)
		for (int j = 0; j <= 1; j++) {
			Alpha_new4(i, j) = A(i, j);
			Beta_new4(i, j) = B(i, j);
			P_new4(i, j) = S(i, j);
		}
	double I_u_1 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14);
	double I_u_2 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24);
	deriva_I_Z(2 * (Z_size2+1)) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(2, 0));
	//=============================================================================================================================//
	Zeros(z_x_add14);
	Zeros(z_y_add14);
	Zeros(z_xx_add14);
	Zeros(z_xy_add14);
	Zeros(z_yy_add14);
	z_x_add14(1, 0) = -DER_INCREMENTS_Z(Z_size1 - 2, 0); z_x_add14(1, 1) = -DER_INCREMENTS_Z(Z_size1 - 2, 0) / 2;
	z_y_add14(0, 0) = DER_INCREMENTS_Z(Z_size1 - 2, 0) / 2; z_y_add14(2, 0) = -DER_INCREMENTS_Z(Z_size1 - 2, 0) / 2;
	z_xx_add14(1, 0) = DER_INCREMENTS_Z(Z_size1 - 2, 0); z_xx_add14(1, 1) = DER_INCREMENTS_Z(Z_size1 - 2, 0);
	z_yy_add14(0, 0) = DER_INCREMENTS_Z(Z_size1 - 2, 0); z_yy_add14(1, 0) = -2 * DER_INCREMENTS_Z(Z_size1 - 2, 0); 
	z_yy_add14(2, 0) = DER_INCREMENTS_Z(Z_size1 - 2, 0); z_yy_add14(3, 0) = DER_INCREMENTS_Z(Z_size1 - 2, 0);
	z_xy_add14(0, 0) = -DER_INCREMENTS_Z(Z_size1 - 2, 0) / 4; z_xy_add14(2, 0) = DER_INCREMENTS_Z(Z_size1 - 2, 0) / 4; z_xy_add14(3, 0) = DER_INCREMENTS_Z(Z_size1 - 2, 0) / 4;
	z_xy_add14(0, 1) = -DER_INCREMENTS_Z(Z_size1 - 2, 0) / 4; z_xy_add14(2, 1) = DER_INCREMENTS_Z(Z_size1 - 2, 0) / 4; z_xy_add14(3, 1) = DER_INCREMENTS_Z(Z_size1 - 2, 0) / 4;
	z_x_add24 = -z_x_add14; z_y_add24 = -z_y_add14; z_xx_add24 = -z_xx_add14; z_yy_add24 = -z_yy_add14; z_xy_add24 = -z_xy_add14;

	for (int i = Z_size1 - 3; i <= Z_size1; i++) {
		for (int j = 0; j <= 1; j++) {
			z_xo4(ix, j) = Z_x(i, j);
			z_yo4(ix, j) = Z_y(i, j);
			z_xxo4(ix, j) = Z_xx(i, j);
			z_yyo4(ix, j) = Z_yy(i, j);
			z_xyo4(ix, j) = Z_xy(i, j);
		}
		ix++;
	}
	ix = 0;
	z_x_new14 = z_xo4 + z_x_add14; z_y_new14 = z_yo4 + z_y_add14; z_xx_new14 = z_xxo4 + z_xx_add14; z_yy_new14 = z_yyo4 + z_yy_add14; z_xy_new14 = z_xyo4 + z_xy_add14;
	z_x_new24 = z_xo4 + z_x_add24; z_y_new24 = z_yo4 + z_y_add24; z_xx_new24 = z_xxo4 + z_xx_add24; z_yy_new24 = z_yyo4 + z_yy_add24; z_xy_new24 = z_xyo4 + z_xy_add24;
	for (int i = Z_size1 - 3; i <= Z_size1; i++) {
		for (int j = 0; j <= 1; j++) {
			Alpha_new4(ix, j) = A(i, j);
			Beta_new4(ix, j) = B(i, j);
			P_new4(ix, j) = S(i, j);
		}
		ix++;
	}
	ix = 0;
	I_u_1 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14);
	I_u_2 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24);
	deriva_I_Z((Z_size1 - 2) * (Z_size2+1)) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(Z_size1 - 2, 0));
	//=============================================================================================================================//
	Zeros(z_x_add14);
	Zeros(z_y_add14);
	Zeros(z_xx_add14);
	Zeros(z_xy_add14);
	Zeros(z_yy_add14);
	z_x_add14(2, 0) = DER_INCREMENTS_Z(2, Z_size2) / 2; z_x_add14(2, 1) = DER_INCREMENTS_Z(2, Z_size2);
	z_y_add14(1, 1) = DER_INCREMENTS_Z(2, Z_size2) / 2; z_y_add14(3, 1) = -DER_INCREMENTS_Z(2, Z_size2) / 2;
	z_xx_add14(2, 0) = DER_INCREMENTS_Z(2, Z_size2); z_xx_add14(2, 1) = DER_INCREMENTS_Z(2, Z_size2);
	z_yy_add14(0, 1) = DER_INCREMENTS_Z(2, Z_size2); z_yy_add14(1, 1) = DER_INCREMENTS_Z(2, Z_size2); 
	z_yy_add14(2, 1) = -2 * DER_INCREMENTS_Z(2, Z_size2); z_yy_add14(3, 1) = DER_INCREMENTS_Z(2, Z_size2);
	z_xy_add14(0, 0) = DER_INCREMENTS_Z(2, Z_size2) / 4; z_xy_add14(1, 0) = DER_INCREMENTS_Z(2, Z_size2) / 4; z_xy_add14(3, 0) = -DER_INCREMENTS_Z(2, Z_size2) / 4;
	z_xy_add14(0, 1) = DER_INCREMENTS_Z(2, Z_size2) / 4; z_xy_add14(1, 1) = DER_INCREMENTS_Z(2, Z_size2) / 4; z_xy_add14(3, 1) = -DER_INCREMENTS_Z(2, Z_size2) / 4;
	z_x_add24 = -z_x_add14; z_y_add24 = -z_y_add14; z_xx_add24 = -z_xx_add14; z_yy_add24 = -z_yy_add14; z_xy_add24 = -z_xy_add14;
	for (int i = 0; i <= 3; i++) {
		for (int j = Z_size2 - 1; j <= Z_size2; j++) {
			z_xo4(i, jy) = Z_x(i, j);
			z_yo4(i, jy) = Z_y(i, j);
			z_xxo4(i, jy) = Z_xx(i, j);
			z_yyo4(i, jy) = Z_yy(i, j);
			z_xyo4(i, jy) = Z_xy(i, j);
			jy++;
		}
		jy = 0;
	}
	jy = 0;
	z_x_new14 = z_xo4 + z_x_add14; z_y_new14 = z_yo4 + z_y_add14; z_xx_new14 = z_xxo4 + z_xx_add14; z_yy_new14 = z_yyo4 + z_yy_add14; z_xy_new14 = z_xyo4 + z_xy_add14;
	z_x_new24 = z_xo4 + z_x_add24; z_y_new24 = z_yo4 + z_y_add24; z_xx_new24 = z_xxo4 + z_xx_add24; z_yy_new24 = z_yyo4 + z_yy_add24; z_xy_new24 = z_xyo4 + z_xy_add24;
	for (int i = 0; i <= 3; i++) {
		for (int j = Z_size2 - 1; j <= Z_size2; j++) {
			Alpha_new4(i, jy) = A(i, j);
			Beta_new4(i, jy) = B(i, j);
			P_new4(i, jy) = S(i, j);
			jy++;
		}
		jy = 0;
	}
	jy = 0;
	I_u_1 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14);
	I_u_2 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24);
	deriva_I_Z(3 * (Z_size2 + 1) - 1) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(2, Z_size2));
	//=============================================================================================================================//
	Zeros(z_x_add14);
	Zeros(z_y_add14);
	Zeros(z_xx_add14);
	Zeros(z_xy_add14);
	Zeros(z_yy_add14);
	z_x_add14(1, 0) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2) / 2; z_x_add14(1, 1) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2);
	z_y_add14(0, 1) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2) / 2; z_y_add14(2, 1) = -DER_INCREMENTS_Z(Z_size1 - 2, Z_size2) / 2;
	z_xx_add14(1, 0) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2); z_xx_add14(1, 1) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2);
	z_yy_add14(0, 1) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2); z_yy_add14(1, 1) = -2 * DER_INCREMENTS_Z(Z_size1 - 2, Z_size2); 
	z_yy_add14(2, 1) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2); z_yy_add14(3, 1) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2);
	z_xy_add14(0, 0) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2) / 4; z_xy_add14(2, 0) = -DER_INCREMENTS_Z(Z_size1 - 2, Z_size2) / 4; z_xy_add14(3, 0) = -DER_INCREMENTS_Z(Z_size1 - 2, Z_size2) / 4;
	z_xy_add14(0, 1) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2) / 4; z_xy_add14(2, 1) = -DER_INCREMENTS_Z(Z_size1 - 2, Z_size2) / 4; z_xy_add14(3, 1) = -DER_INCREMENTS_Z(Z_size1 - 2, Z_size2) / 4;
	z_x_add24 = -z_x_add14; z_y_add24 = -z_y_add14; z_xx_add24 = -z_xx_add14; z_yy_add24 = -z_yy_add14; z_xy_add24 = -z_xy_add14;
	//z_xo = z_x0(Z_size1 - 3:Z_size1, Z_size2 - 1 : Z_size2); z_yo = z_y0(Z_size1 - 3:Z_size1, Z_size2 - 1 : Z_size2); z_xxo = z_xx0(Z_size1 - 3:Z_size1, Z_size2 - 1 : Z_size2); z_yyo = z_yy0(Z_size1 - 3:Z_size1, Z_size2 - 1 : Z_size2); z_xyo = z_xy0(Z_size1 - 3:Z_size1, Z_size2 - 1 : Z_size2);
	for (int i = Z_size1 - 3; i <= Z_size1; i++) {
		for (int j = Z_size2 - 1; j <= Z_size2; j++) {
			z_xo4(ix, jy) = Z_x(i, j);
			z_yo4(ix, jy) = Z_y(i, j);
			z_xxo4(ix, jy) = Z_xx(i, j);
			z_yyo4(ix, jy) = Z_yy(i, j);
			z_xyo4(ix, jy) = Z_xy(i, j);
			jy++;
		}
		ix++; jy = 0;
	}
	ix = 0; jy = 0;
	z_x_new14 = z_xo4 + z_x_add14; z_y_new14 = z_yo4 + z_y_add14; z_xx_new14 = z_xxo4 + z_xx_add14; z_yy_new14 = z_yyo4 + z_yy_add14; z_xy_new14 = z_xyo4 + z_xy_add14;
	z_x_new24 = z_xo4 + z_x_add24; z_y_new24 = z_yo4 + z_y_add24; z_xx_new24 = z_xxo4 + z_xx_add24; z_yy_new24 = z_yyo4 + z_yy_add24; z_xy_new24 = z_xyo4 + z_xy_add24;
	for (int i = Z_size1 - 3; i <= Z_size1; i++) {
		for (int j = Z_size2 - 1; j <= Z_size2; j++) {
			Alpha_new4(ix, jy) = A(i, j);
			Beta_new4(ix, jy) = B(i, j);
			P_new4(ix, jy) = S(i, j);
			jy++;
		}
		ix++; jy = 0;
	}
	ix = 0; jy = 0;
	I_u_1 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14);
	I_u_2 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24);
	deriva_I_Z((Z_size1 - 1)* (Z_size2 + 1) - 1) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(Z_size1 - 2, Z_size2));
}
void Adam::calculate_A33(double con1, double con2, vector<int>(&D_pos)[2], vector<int>(&N_pos)[2], int min_pos1, int max_pos2) {
	Matrix<double, 3, 3> z_x_add13, z_y_add13, z_xx_add13, z_yy_add13, z_xy_add13;
	Matrix<double, 3, 3> z_x_add23, z_y_add23, z_xx_add23, z_yy_add23, z_xy_add23;
	Matrix<double, 3, 3> z_xo3, z_yo3, z_xxo3, z_yyo3, z_xyo3;
	Matrix<double, 3, 3> z_x_new13, z_y_new13, z_xx_new13, z_yy_new13, z_xy_new13;
	Matrix<double, 3, 3> z_x_new23, z_y_new23, z_xx_new23, z_yy_new23, z_xy_new23;
	Matrix<double, 3, 3> Alpha_new3, Beta_new3, P_new3;
	int ix = 0, jy = 0;
	Zeros(z_x_add13);
	Zeros(z_y_add13);
	Zeros(z_xx_add13);
	Zeros(z_xy_add13);
	Zeros(z_yy_add13);
	z_x_add13(1, 0) = DER_INCREMENTS_Z(1, 1); z_x_add13(1, 2) = -DER_INCREMENTS_Z(1, 1) / 2;
	z_y_add13(0, 1) = DER_INCREMENTS_Z(1, 1); z_y_add13(2, 1) = -DER_INCREMENTS_Z(1, 1) / 2;
	z_xx_add13(1, 0) = -2 * DER_INCREMENTS_Z(1, 1); z_xx_add13(1, 1) = -2 * DER_INCREMENTS_Z(1, 1); z_xx_add13(1, 2) = DER_INCREMENTS_Z(1, 1);
	z_yy_add13(0, 1) = -2 * DER_INCREMENTS_Z(1, 1); z_yy_add13(1, 1) = -2 * DER_INCREMENTS_Z(1, 1); z_yy_add13(2, 1) = DER_INCREMENTS_Z(1, 1);
	z_xy_add13(2, 2) = DER_INCREMENTS_Z(1, 1) / 4;
	z_x_add23 = -z_x_add13; z_y_add23 = -z_y_add13; z_xx_add23 = -z_xx_add13; z_yy_add23 = -z_yy_add13; z_xy_add23 = -z_xy_add13;
	for (int i = 0; i <= 2; i++)
		for (int j = 0; j <= 2; j++) {
			z_xo3(i, j) = Z_x(i, j);
			z_yo3(i, j) = Z_y(i, j);
			z_xxo3(i, j) = Z_xx(i, j);
			z_yyo3(i, j) = Z_yy(i, j);
			z_xyo3(i, j) = Z_xy(i, j);
		}
	z_x_new13 = z_xo3 + z_x_add13; z_y_new13 = z_yo3 + z_y_add13; z_xx_new13 = z_xxo3 + z_xx_add13; z_yy_new13 = z_yyo3 + z_yy_add13; z_xy_new13 = z_xyo3 + z_xy_add13;
	z_x_new23 = z_xo3 + z_x_add23; z_y_new23 = z_yo3 + z_y_add23; z_xx_new23 = z_xxo3 + z_xx_add23; z_yy_new23 = z_yyo3 + z_yy_add23; z_xy_new23 = z_xyo3 + z_xy_add23;
	for (int i = 0; i <= 2; i++)
		for (int j = 0; j <= 2; j++) {
			Alpha_new3(i, j) = A(i, j);
			Beta_new3(i, j) = B(i, j);
			P_new3(i, j) = S(i, j);
		}
	double I_u_1 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new13, z_y_new13, z_xx_new13, z_yy_new13, z_xy_new13);
	double I_u_2 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new23, z_y_new23, z_xx_new23, z_yy_new23, z_xy_new23);
	deriva_I_Z(Z_size2 + 2) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(1, 1));
	//=============================================================================================================================//
	Zeros(z_x_add13);
	Zeros(z_y_add13);
	Zeros(z_xx_add13);
	Zeros(z_xy_add13);
	Zeros(z_yy_add13);
	z_x_add13(1, 0) = DER_INCREMENTS_Z(1, Z_size2 - 1) / 2; z_x_add13(1, 2) = -DER_INCREMENTS_Z(1, Z_size2 - 1);
	z_y_add13(0, 1) = DER_INCREMENTS_Z(1, Z_size2 - 1); z_y_add13(2, 1) = -DER_INCREMENTS_Z(1, Z_size2 - 1) / 2;
	z_xx_add13(1, 0) = DER_INCREMENTS_Z(1, Z_size2 - 1); z_xx_add13(1, 1) = -2 * DER_INCREMENTS_Z(1, Z_size2 - 1); z_xx_add13(1, 2) = -2 * DER_INCREMENTS_Z(1, Z_size2 - 1);
	z_yy_add13(0, 1) = -2 * DER_INCREMENTS_Z(1, Z_size2 - 1); z_yy_add13(1, 1) = -2 * DER_INCREMENTS_Z(1, Z_size2 - 1); z_yy_add13(2, 1) = DER_INCREMENTS_Z(1, Z_size2 - 1);
	z_xy_add13(2, 0) = -DER_INCREMENTS_Z(1, Z_size2 - 1) / 4;
	z_x_add23 = -z_x_add13; z_y_add23 = -z_y_add13; z_xx_add23 = -z_xx_add13; z_yy_add23 = -z_yy_add13; z_xy_add23 = -z_xy_add13;
	for (int i = 0; i <= 2; i++) {
		for (int j = Z_size2 - 2; j <= Z_size2; j++) {
			z_xo3(i, jy) = Z_x(i, j);
			z_yo3(i, jy) = Z_y(i, j);
			z_xxo3(i, jy) = Z_xx(i, j);
			z_yyo3(i, jy) = Z_yy(i, j);
			z_xyo3(i, jy) = Z_xy(i, j);
			jy++;
		}
		jy = 0;
	}
	jy = 0;
	z_x_new13 = z_xo3 + z_x_add13; z_y_new13 = z_yo3 + z_y_add13; z_xx_new13 = z_xxo3 + z_xx_add13; z_yy_new13 = z_yyo3 + z_yy_add13; z_xy_new13 = z_xyo3 + z_xy_add13;
	z_x_new23 = z_xo3 + z_x_add23; z_y_new23 = z_yo3 + z_y_add23; z_xx_new23 = z_xxo3 + z_xx_add23; z_yy_new23 = z_yyo3 + z_yy_add23; z_xy_new23 = z_xyo3 + z_xy_add23;
	for (int i = 0; i <= 2; i++) {
		for (int j = Z_size2 - 2; j <= Z_size2; j++) {
			Alpha_new3(i, jy) = A(i, j);
			Beta_new3(i, jy) = B(i, j);
			P_new3(i, jy) = S(i, j);
			jy++;
		}
		jy = 0;
	}
	jy = 0;
	I_u_1 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new13, z_y_new13, z_xx_new13, z_yy_new13, z_xy_new13);
	I_u_2 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new23, z_y_new23, z_xx_new23, z_yy_new23, z_xy_new23);
	deriva_I_Z(Z_size2 + Z_size2) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(1, Z_size2 - 1));
	//=============================================================================================================================//
	Zeros(z_x_add13);
	Zeros(z_y_add13);
	Zeros(z_xx_add13);
	Zeros(z_xy_add13);
	Zeros(z_yy_add13);
	z_x_add13(1, 0) = DER_INCREMENTS_Z(Z_size1 - 1, 1); z_x_add13(1, 2) = -DER_INCREMENTS_Z(Z_size1 - 1, 1) / 2;
	z_y_add13(0, 1) = DER_INCREMENTS_Z(Z_size1 - 1, 1) / 2; z_y_add13(2, 1) = -DER_INCREMENTS_Z(Z_size1 - 1, 1);
	z_xx_add13(1, 0) = -2 * DER_INCREMENTS_Z(Z_size1 - 1, 1); z_xx_add13(1, 1) = -2 * DER_INCREMENTS_Z(Z_size1 - 1, 1); z_xx_add13(1, 2) = DER_INCREMENTS_Z(Z_size1 - 1, 1);
	z_yy_add13(0, 1) = DER_INCREMENTS_Z(Z_size1 - 1, 1); z_yy_add13(1, 1) = -2 * DER_INCREMENTS_Z(Z_size1 - 1, 1); z_yy_add13(2, 1) = -2 * DER_INCREMENTS_Z(Z_size1 - 1, 1);
	z_xy_add13(0, 2) = -DER_INCREMENTS_Z(Z_size1 - 1, 1) / 4;
	z_x_add23 = -z_x_add13; z_y_add23 = -z_y_add13; z_xx_add23 = -z_xx_add13; z_yy_add23 = -z_yy_add13; z_xy_add23 = -z_xy_add13;
	for (int i = Z_size1 - 2; i <= Z_size1; i++) {
		for (int j = 0; j <= 2; j++) {
			z_xo3(ix, j) = Z_x(i, j);
			z_yo3(ix, j) = Z_y(i, j);
			z_xxo3(ix, j) = Z_xx(i, j);
			z_yyo3(ix, j) = Z_yy(i, j);
			z_xyo3(ix, j) = Z_xy(i, j);
		}
		ix++;
	}
	ix = 0;
	z_x_new13 = z_xo3 + z_x_add13; z_y_new13 = z_yo3 + z_y_add13; z_xx_new13 = z_xxo3 + z_xx_add13; z_yy_new13 = z_yyo3 + z_yy_add13; z_xy_new13 = z_xyo3 + z_xy_add13;
	z_x_new23 = z_xo3 + z_x_add23; z_y_new23 = z_yo3 + z_y_add23; z_xx_new23 = z_xxo3 + z_xx_add23; z_yy_new23 = z_yyo3 + z_yy_add23; z_xy_new23 = z_xyo3 + z_xy_add23;
	for (int i = Z_size1 - 2; i <= Z_size1; i++) {
		for (int j = 0; j <= 2; j++) {
			Alpha_new3(ix, j) = A(i, j);
			Beta_new3(ix, j) = B(i, j);
			P_new3(ix, j) = S(i, j);
		}
		ix++;
	}
	ix = 0;
	I_u_1 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new13, z_y_new13, z_xx_new13, z_yy_new13, z_xy_new13);
	I_u_2 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new23, z_y_new23, z_xx_new23, z_yy_new23, z_xy_new23);
	deriva_I_Z((Z_size1 - 1) * (Z_size2 + 1) + 1) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(Z_size1 - 1, 1));
	//=============================================================================================================================//
	Zeros(z_x_add13);
	Zeros(z_y_add13);
	Zeros(z_xx_add13);
	Zeros(z_xy_add13);
	Zeros(z_yy_add13);
	z_x_add13(1, 0) = DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 1) / 2; z_x_add13(1, 2) = -DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 1);
	z_y_add13(0, 1) = DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 1) / 2; z_y_add13(2, 1) = -DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 1);
	z_xx_add13(1, 0) = DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 1); z_xx_add13(1, 1) = -2 * DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 1); z_xx_add13(1, 2) = -2 * DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 1);
	z_yy_add13(0, 1) = DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 1); z_yy_add13(1, 1) = -2 * DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 1); z_yy_add13(2, 1) = -2 * DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 1);
	z_xy_add13(0, 0) = DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 1) / 4;
	z_x_add23 = -z_x_add13; z_y_add23 = -z_y_add13; z_xx_add23 = -z_xx_add13; z_yy_add23 = -z_yy_add13; z_xy_add23 = -z_xy_add13;
	for (int i = Z_size1 - 2; i <= Z_size1; i++) {
		for (int j = Z_size2 - 2; j <= Z_size2; j++) {
			z_xo3(ix, jy) = Z_x(i, j);
			z_yo3(ix, jy) = Z_y(i, j);
			z_xxo3(ix, jy) = Z_xx(i, j);
			z_yyo3(ix, jy) = Z_yy(i, j);
			z_xyo3(ix, jy) = Z_xy(i, j);
			jy++;
		}
		ix++; jy = 0;
	}
	ix = 0; jy = 0;
	z_x_new13 = z_xo3 + z_x_add13; z_y_new13 = z_yo3 + z_y_add13; z_xx_new13 = z_xxo3 + z_xx_add13; z_yy_new13 = z_yyo3 + z_yy_add13; z_xy_new13 = z_xyo3 + z_xy_add13;
	z_x_new23 = z_xo3 + z_x_add23; z_y_new23 = z_yo3 + z_y_add23; z_xx_new23 = z_xxo3 + z_xx_add23; z_yy_new23 = z_yyo3 + z_yy_add23; z_xy_new23 = z_xyo3 + z_xy_add23;
	for (int i = Z_size1 - 2; i <= Z_size1; i++) {
		for (int j = Z_size2 - 2; j <= Z_size2; j++) {
			Alpha_new3(ix, jy) = A(i, j);
			Beta_new3(ix, jy) = B(i, j);
			P_new3(ix, jy) = S(i, j);
			jy++;
		}
		ix++; jy = 0;
	}
	ix = 0; jy = 0;
	I_u_1 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new13, z_y_new13, z_xx_new13, z_yy_new13, z_xy_new13);
	I_u_2 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new23, z_y_new23, z_xx_new23, z_yy_new23, z_xy_new23);
	deriva_I_Z((Z_size1 - 1) * (Z_size2 + 1) + Z_size2 - 1) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 1));
	//=============================================================================================================================//
	for (int j = 3; j <= Z_size2 - 3; j++) {
		int i = 1;
		Zeros(z_x_add13);
		Zeros(z_y_add13);
		Zeros(z_xx_add13);
		Zeros(z_xy_add13);
		Zeros(z_yy_add13);
		z_x_add13(1, 0) = DER_INCREMENTS_Z(i, j) / 2; z_x_add13(1, 2) = -DER_INCREMENTS_Z(i, j) / 2;
		z_y_add13(0, 1) = DER_INCREMENTS_Z(i, j); z_y_add13(2, 1) = -DER_INCREMENTS_Z(i, j) / 2;
		z_xx_add13(1, 0) = DER_INCREMENTS_Z(i, j); z_xx_add13(1, 1) = -2 * DER_INCREMENTS_Z(i, j); z_xx_add13(1, 2) = DER_INCREMENTS_Z(i, j);
		z_yy_add13(0, 1) = -2 * DER_INCREMENTS_Z(i, j); z_yy_add13(1, 1) = -2 * DER_INCREMENTS_Z(i, j); z_yy_add13(2, 1) = DER_INCREMENTS_Z(i, j);
		z_xy_add13(2, 0) = -DER_INCREMENTS_Z(i, j) / 4; z_xy_add13(2, 2) = DER_INCREMENTS_Z(i, j) / 4;
		z_x_add23 = -z_x_add13; z_y_add23 = -z_y_add13; z_xx_add23 = -z_xx_add13; z_yy_add23 = -z_yy_add13; z_xy_add23 = -z_xy_add13;
		for (int ii = i - 1; ii <= i + 1; ii++) {
			for (int jj = j - 1; jj <= j + 1; jj++) {
				z_xo3(ix, jy) = Z_x(ii, jj);
				z_yo3(ix, jy) = Z_y(ii, jj);
				z_xxo3(ix, jy) = Z_xx(ii, jj);
				z_yyo3(ix, jy) = Z_yy(ii, jj);
				z_xyo3(ix, jy) = Z_xy(ii, jj);
				jy++;
			}
			ix++; jy = 0;
		}
		ix = 0; jy = 0;
		z_x_new13 = z_xo3 + z_x_add13; z_y_new13 = z_yo3 + z_y_add13; z_xx_new13 = z_xxo3 + z_xx_add13; z_yy_new13 = z_yyo3 + z_yy_add13; z_xy_new13 = z_xyo3 + z_xy_add13;
		z_x_new23 = z_xo3 + z_x_add23; z_y_new23 = z_yo3 + z_y_add23; z_xx_new23 = z_xxo3 + z_xx_add23; z_yy_new23 = z_yyo3 + z_yy_add23; z_xy_new23 = z_xyo3 + z_xy_add23;
		for (int ii = i - 1; ii <= i + 1; ii++) {
			for (int jj = j - 1; jj <= j + 1; jj++) {
				Alpha_new3(ix, jy) = A(ii, jj);
				Beta_new3(ix, jy) = B(ii, jj);
				P_new3(ix, jy) = S(ii, jj);
				jy++;
			}
			ix++; jy = 0;
		}
		ix = 0; jy = 0;
		I_u_1 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new13, z_y_new13, z_xx_new13, z_yy_new13, z_xy_new13);
		I_u_2 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new23, z_y_new23, z_xx_new23, z_yy_new23, z_xy_new23);
		deriva_I_Z((i) * (Z_size2 + 1) + j) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(i, j));
	}
	//=============================================================================================================================//
	for (int j = 3; j <= Z_size2 - 3; j++) {
		int i = Z_size1 - 1;
		Zeros(z_x_add13);
		Zeros(z_y_add13);
		Zeros(z_xx_add13);
		Zeros(z_xy_add13);
		Zeros(z_yy_add13);
		z_x_add13(1, 0) = DER_INCREMENTS_Z(i, j) / 2; z_x_add13(1, 2) = -DER_INCREMENTS_Z(i, j) / 2;
		z_y_add13(0, 1) = DER_INCREMENTS_Z(i, j) / 2; z_y_add13(2, 1) = -DER_INCREMENTS_Z(i, j);
		z_xx_add13(1, 0) = DER_INCREMENTS_Z(i, j); z_xx_add13(1, 1) = -2 * DER_INCREMENTS_Z(i, j); z_xx_add13(1, 2) = DER_INCREMENTS_Z(i, j);
		z_yy_add13(0, 1) = DER_INCREMENTS_Z(i, j); z_yy_add13(1, 1) = -2 * DER_INCREMENTS_Z(i, j); z_yy_add13(2, 1) = -2 * DER_INCREMENTS_Z(i, j);
		z_xy_add13(0, 0) = DER_INCREMENTS_Z(i, j) / 4; z_xy_add13(0, 2) = -DER_INCREMENTS_Z(i, j) / 4;
		z_x_add23 = -z_x_add13; z_y_add23 = -z_y_add13; z_xx_add23 = -z_xx_add13; z_yy_add23 = -z_yy_add13; z_xy_add23 = -z_xy_add13;
		for (int ii = i - 1; ii <= i + 1; ii++) {
			for (int jj = j - 1; jj <= j + 1; jj++) {
				z_xo3(ix, jy) = Z_x(ii, jj);
				z_yo3(ix, jy) = Z_y(ii, jj);
				z_xxo3(ix, jy) = Z_xx(ii, jj);
				z_yyo3(ix, jy) = Z_yy(ii, jj);
				z_xyo3(ix, jy) = Z_xy(ii, jj);
				jy++;
			}
			ix++; jy = 0;
		}
		ix = 0; jy = 0;
		z_x_new13 = z_xo3 + z_x_add13; z_y_new13 = z_yo3 + z_y_add13; z_xx_new13 = z_xxo3 + z_xx_add13; z_yy_new13 = z_yyo3 + z_yy_add13; z_xy_new13 = z_xyo3 + z_xy_add13;
		z_x_new23 = z_xo3 + z_x_add23; z_y_new23 = z_yo3 + z_y_add23; z_xx_new23 = z_xxo3 + z_xx_add23; z_yy_new23 = z_yyo3 + z_yy_add23; z_xy_new23 = z_xyo3 + z_xy_add23;
		for (int ii = i - 1; ii <= i + 1; ii++) {
			for (int jj = j - 1; jj <= j + 1; jj++) {
				Alpha_new3(ix, jy) = A(ii, jj);
				Beta_new3(ix, jy) = B(ii, jj);
				P_new3(ix, jy) = S(ii, jj);
				jy++;
			}
			ix++; jy = 0;
		}
		ix = 0; jy = 0;
		I_u_1 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new13, z_y_new13, z_xx_new13, z_yy_new13, z_xy_new13);
		I_u_2 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new23, z_y_new23, z_xx_new23, z_yy_new23, z_xy_new23);
		deriva_I_Z((i) * (Z_size2 + 1) + j) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(i, j));
	}
	//=============================================================================================================================//
	for (int i = 3; i <= Z_size1 - 3; i++) {
		int j = 1;
		Zeros(z_x_add13);
		Zeros(z_y_add13);
		Zeros(z_xx_add13);
		Zeros(z_xy_add13);
		Zeros(z_yy_add13);
		z_x_add13(1, 0) = DER_INCREMENTS_Z(i, j); z_x_add13(1, 2) = -DER_INCREMENTS_Z(i, j) / 2;
		z_y_add13(0, 1) = DER_INCREMENTS_Z(i, j) / 2; z_y_add13(2, 1) = -DER_INCREMENTS_Z(i, j) / 2;
		z_xx_add13(1, 0) = -2 * DER_INCREMENTS_Z(i, j); z_xx_add13(1, 1) = -2 * DER_INCREMENTS_Z(i, j); z_xx_add13(1, 2) = DER_INCREMENTS_Z(i, j);
		z_yy_add13(0, 1) = DER_INCREMENTS_Z(i, j); z_yy_add13(1, 1) = -2 * DER_INCREMENTS_Z(i, j); z_yy_add13(2, 1) = DER_INCREMENTS_Z(i, j);
		z_xy_add13(0, 2) = -DER_INCREMENTS_Z(i, j) / 4; z_xy_add13(2, 2) = DER_INCREMENTS_Z(i, j) / 4;
		z_x_add23 = -z_x_add13; z_y_add23 = -z_y_add13; z_xx_add23 = -z_xx_add13; z_yy_add23 = -z_yy_add13; z_xy_add23 = -z_xy_add13;
		for (int ii = i - 1; ii <= i + 1; ii++) {
			for (int jj = j - 1; jj <= j + 1; jj++) {
				z_xo3(ix, jy) = Z_x(ii, jj);
				z_yo3(ix, jy) = Z_y(ii, jj);
				z_xxo3(ix, jy) = Z_xx(ii, jj);
				z_yyo3(ix, jy) = Z_yy(ii, jj);
				z_xyo3(ix, jy) = Z_xy(ii, jj);
				jy++;
			}
			ix++; jy = 0;
		}
		ix = 0; jy = 0;
		z_x_new13 = z_xo3 + z_x_add13; z_y_new13 = z_yo3 + z_y_add13; z_xx_new13 = z_xxo3 + z_xx_add13; z_yy_new13 = z_yyo3 + z_yy_add13; z_xy_new13 = z_xyo3 + z_xy_add13;
		z_x_new23 = z_xo3 + z_x_add23; z_y_new23 = z_yo3 + z_y_add23; z_xx_new23 = z_xxo3 + z_xx_add23; z_yy_new23 = z_yyo3 + z_yy_add23; z_xy_new23 = z_xyo3 + z_xy_add23;
		for (int ii = i - 1; ii <= i + 1; ii++) {
			for (int jj = j - 1; jj <= j + 1; jj++) {
				Alpha_new3(ix, jy) = A(ii, jj);
				Beta_new3(ix, jy) = B(ii, jj);
				P_new3(ix, jy) = S(ii, jj);
				jy++;
			}
			ix++; jy = 0;
		}
		ix = 0; jy = 0;
		I_u_1 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new13, z_y_new13, z_xx_new13, z_yy_new13, z_xy_new13);
		I_u_2 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new23, z_y_new23, z_xx_new23, z_yy_new23, z_xy_new23);
		deriva_I_Z((i) * (Z_size2 + 1) + j) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(i, j));
	}
	//=============================================================================================================================//
	for (int i = 3; i <= Z_size1 - 3; i++) {
		int j = Z_size2 - 1;
		Zeros(z_x_add13);
		Zeros(z_y_add13);
		Zeros(z_xx_add13);
		Zeros(z_xy_add13);
		Zeros(z_yy_add13);
		z_x_add13(1, 0) = DER_INCREMENTS_Z(i, j) / 2; z_x_add13(1, 2) = -DER_INCREMENTS_Z(i, j);
		z_y_add13(0, 1) = DER_INCREMENTS_Z(i, j) / 2; z_y_add13(2, 1) = -DER_INCREMENTS_Z(i, j) / 2;
		z_xx_add13(1, 0) = DER_INCREMENTS_Z(i, j); z_xx_add13(1, 1) = -2 * DER_INCREMENTS_Z(i, j); z_xx_add13(1, 2) = -2 * DER_INCREMENTS_Z(i, j);
		z_yy_add13(0, 1) = DER_INCREMENTS_Z(i, j); z_yy_add13(1, 1) = -2 * DER_INCREMENTS_Z(i, j); z_yy_add13(2, 1) = DER_INCREMENTS_Z(i, j);
		z_xy_add13(0, 0) = DER_INCREMENTS_Z(i, j) / 4; z_xy_add13(2, 0) = -DER_INCREMENTS_Z(i, j) / 4;
		z_x_add23 = -z_x_add13; z_y_add23 = -z_y_add13; z_xx_add23 = -z_xx_add13; z_yy_add23 = -z_yy_add13; z_xy_add23 = -z_xy_add13;
		for (int ii = i - 1; ii <= i + 1; ii++) {
			for (int jj = j - 1; jj <= j + 1; jj++) {
				z_xo3(ix, jy) = Z_x(ii, jj);
				z_yo3(ix, jy) = Z_y(ii, jj);
				z_xxo3(ix, jy) = Z_xx(ii, jj);
				z_yyo3(ix, jy) = Z_yy(ii, jj);
				z_xyo3(ix, jy) = Z_xy(ii, jj);
				jy++;
			}
			ix++; jy = 0;
		}
		ix = 0; jy = 0;
		z_x_new13 = z_xo3 + z_x_add13; z_y_new13 = z_yo3 + z_y_add13; z_xx_new13 = z_xxo3 + z_xx_add13; z_yy_new13 = z_yyo3 + z_yy_add13; z_xy_new13 = z_xyo3 + z_xy_add13;
		z_x_new23 = z_xo3 + z_x_add23; z_y_new23 = z_yo3 + z_y_add23; z_xx_new23 = z_xxo3 + z_xx_add23; z_yy_new23 = z_yyo3 + z_yy_add23; z_xy_new23 = z_xyo3 + z_xy_add23;
		for (int ii = i - 1; ii <= i + 1; ii++) {
			for (int jj = j - 1; jj <= j + 1; jj++) {
				Alpha_new3(ix, jy) = A(ii, jj);
				Beta_new3(ix, jy) = B(ii, jj);
				P_new3(ix, jy) = S(ii, jj);
				jy++;
			}
			ix++; jy = 0;
		}
		ix = 0; jy = 0;
		I_u_1 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new13, z_y_new13, z_xx_new13, z_yy_new13, z_xy_new13);
		I_u_2 = evaluation_func0(Alpha_new3, Beta_new3, P_new3, z_x_new23, z_y_new23, z_xx_new23, z_yy_new23, z_xy_new23);
		deriva_I_Z((i) * (Z_size2 + 1) + j) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(i, j));
	}
	//=============================================================================================================================//
	double con_temp1 = 0.0;
	double con_temp2 = 0.0;
	double con_temp3 = 0.0;
	double con_temp4 = 0.0;
	for (int i = 3; i <= Z_size1 - 3; i++) {
		for (int j = 3; j <= Z_size2 - 3; j++) {
			Zeros(z_x_add13);
			Zeros(z_y_add13);
			Zeros(z_xx_add13);
			Zeros(z_xy_add13);
			Zeros(z_yy_add13);
			z_x_add13(1, 0) = DER_INCREMENTS_Z(i, j) / 2; z_x_add13(1, 2) = -DER_INCREMENTS_Z(i, j) / 2;
			z_y_add13(0, 1) = DER_INCREMENTS_Z(i, j) / 2; z_y_add13(2, 1) = -DER_INCREMENTS_Z(i, j) / 2;
			z_xx_add13(1, 0) = DER_INCREMENTS_Z(i, j); z_xx_add13(1, 1) = -2 * DER_INCREMENTS_Z(i, j); z_xx_add13(1, 2) = DER_INCREMENTS_Z(i, j);
			z_yy_add13(0, 1) = DER_INCREMENTS_Z(i, j); z_yy_add13(1, 1) = -2 * DER_INCREMENTS_Z(i, j); z_yy_add13(2, 1) = DER_INCREMENTS_Z(i, j);
			z_xy_add13(0, 0) = DER_INCREMENTS_Z(i, j) / 4; z_xy_add13(0, 2) = -DER_INCREMENTS_Z(i, j) / 4;
			z_xy_add13(2, 0) = -DER_INCREMENTS_Z(i, j) / 4; z_xy_add13(2, 2) = DER_INCREMENTS_Z(i, j) / 4;
			z_x_add23 = -z_x_add13; z_y_add23 = -z_y_add13; z_xx_add23 = -z_xx_add13; z_yy_add23 = -z_yy_add13; z_xy_add23 = -z_xy_add13;
			for (int ii = i - 1; ii <= i + 1; ii++) {
				for (int jj = j - 1; jj <= j + 1; jj++) {
					z_xo3(ix, jy) = Z_x(ii, jj);
					z_yo3(ix, jy) = Z_y(ii, jj);
					z_xxo3(ix, jy) = Z_xx(ii, jj);
					z_yyo3(ix, jy) = Z_yy(ii, jj);
					z_xyo3(ix, jy) = Z_xy(ii, jj);
					jy++;
				}
				ix++; jy = 0;
			}
			ix = 0; jy = 0;
			z_x_new13 = z_xo3 + z_x_add13; z_y_new13 = z_yo3 + z_y_add13; z_xx_new13 = z_xxo3 + z_xx_add13; z_yy_new13 = z_yyo3 + z_yy_add13; z_xy_new13 = z_xyo3 + z_xy_add13;
			z_x_new23 = z_xo3 + z_x_add23; z_y_new23 = z_yo3 + z_y_add23; z_xx_new23 = z_xxo3 + z_xx_add23; z_yy_new23 = z_yyo3 + z_yy_add23; z_xy_new23 = z_xyo3 + z_xy_add23;
			for (int ii = i - 1; ii <= i + 1; ii++) {
				for (int jj = j - 1; jj <= j + 1; jj++) {
					Alpha_new3(ix, jy) = A(ii, jj);
					Beta_new3(ix, jy) = B(ii, jj);
					P_new3(ix, jy) = S(ii, jj);
					jy++;
				}
				ix++; jy = 0;
			}
			ix = 0; jy = 0;
			if (con1 == 0) {
				con_temp1 = 0.0;
				con_temp2 = 0.0;
			}
			else {
				int temp1, temp2;
				if (min_pos1 >= D_pos[0].size()) {
					temp1 = N_pos[0][min_pos1 - D_pos[0].size()];
					temp2 = N_pos[1][min_pos1 - D_pos[0].size()];
				}
				else {
					temp1 = D_pos[0][min_pos1];
					temp2 = D_pos[1][min_pos1];
				}
				if (temp1 >= i - 1 && temp1 <= i + 1 && temp2 >= j - 1 && temp2 <= j + 1) {
					temp1 = temp1 - i + 1;
					temp2 = temp2 - j + 1;
					double temp3 = z_y_new13(temp1, temp2) * z_y_new13(temp1, temp2);
					double temp33 = z_x_new13(temp1, temp2) * z_x_new13(temp1, temp2);
					double temp4 = pow(1 + temp33 + temp3, 1.5);
					double temp5 = pow((((1 + temp3) * z_xx_new13(temp1, temp2) - 2 * z_x_new13(temp1, temp2) * z_y_new13(temp1, temp2) * z_xy_new13(temp1, temp2)
						+ (1 + z_x_new13(temp1, temp2) * z_x_new13(temp1, temp2)) * z_yy_new13(temp1, temp2)) / (2 * temp4)) * 1e3, 2.0);
					double temp6 = z_xx_new13(temp1, temp2) * z_yy_new13(temp1, temp2) - z_xy_new13(temp1, temp2) * z_xy_new13(temp1, temp2);
					con_temp1 = abs((1 - n) * 2 * sqrt(temp5 - (temp6 / ((1 + temp33 + temp3) * (1 + temp33 + temp3))) * 1e6)) - 0.12;
					double tempp3 = z_y_new23(temp1, temp2) * z_y_new23(temp1, temp2);
					double tempp4 = z_x_new23(temp1, temp2) * z_x_new23(temp1, temp2);
					double tempp5 = pow(1 + tempp4 + tempp3, 1.5);
					double tempp6 = pow((((1 + tempp3) * z_xx_new23(temp1, temp2) - 2 * z_x_new23(temp1, temp2) * z_y_new23(temp1, temp2) * z_xy_new23(temp1, temp2) +
						(1 + tempp4) * z_yy_new23(temp1, temp2)) / (2 * tempp5)) * 1e3, 2.0);
					con_temp2 = abs((1 - n) * 2 * sqrt(tempp6 - ((z_xx_new23(temp1, temp2) * z_yy_new23(temp1, temp2) - z_xy_new23(temp1, temp2) * z_xy_new23(temp1, temp2)) /
						((1 + tempp4 + tempp3) * (1 + tempp4 + tempp3))) * 1e6)) - 0.12;
				}
				else {
					con_temp1 = con1;
					con_temp2 = con1;
				}
			}
			if (con2 == 0) {
				con_temp3 = 0.0;
				con_temp4 = 0.0;
			}
			else {
				int temp1, temp2;
				if (max_pos2 >= D_pos[0].size()) {
					temp1 = N_pos[0][max_pos2 - D_pos[0].size()];
					temp2 = N_pos[1][max_pos2 - D_pos[0].size()];
				}
				else {
					temp1 = D_pos[0][max_pos2];
					temp2 = D_pos[1][max_pos2];
				}
				if (temp1 >= i - 1 && temp1 <= i + 1 && temp2 >= j - 1 && temp2 <= j + 1) {
					temp1 = temp1 - i + 1;
					temp2 = temp2 - j + 1;
					double temp3 = z_y_new13(temp1, temp2) * z_y_new13(temp1, temp2);
					double temp4 = z_x_new13(temp1, temp2) * z_x_new13(temp1, temp2);
					con_temp3 = abs((1 - n) * (((1 + temp3) * z_xx_new13(temp1, temp2) - 2 * z_x_new13(temp1, temp2) * z_y_new13(temp1, temp2) * z_xy_new13(temp1, temp2) + (1 + temp4) * z_yy_new13(temp1, temp2)) / (2 * pow(1 + temp4 + temp3, 1.5))) * 1e3 - (1 - n) * S(temp1 + i - 1, temp2 + j - 1) * 1e3) - 0.06;
					con_temp4 = abs((1 - n) * (((1 + z_y_new23(temp1, temp2) * z_y_new23(temp1, temp2)) * z_xx_new23(temp1, temp2) - 2 * z_x_new23(temp1, temp2) * z_y_new23(temp1, temp2) * z_xy_new23(temp1, temp2) +
						(1 + z_x_new23(temp1, temp2) * z_x_new23(temp1, temp2)) * z_yy_new23(temp1, temp2)) / (2 * pow(1 + z_x_new23(temp1, temp2) * z_x_new23(temp1, temp2) + z_y_new23(temp1, temp2) * z_y_new23(temp1, temp2), 1.5))) * 1e3 - (1 - n) * S(temp1 + i - 1, temp2 + j - 1) * 1e3) - 0.06;
				}
				else {
					con_temp3 = con2;
					con_temp4 = con2;
				}
			}
			I_u_1 = evaluation_func(Alpha_new3, Beta_new3, P_new3, z_x_new13, z_y_new13, z_xx_new13, z_yy_new13, z_xy_new13, 0.0, 0.0);
			I_u_2 = evaluation_func(Alpha_new3, Beta_new3, P_new3, z_x_new23, z_y_new23, z_xx_new23, z_yy_new23, z_xy_new23, 0.0, 0.0);
			double I_u_3 = con_temp2 - con_temp1 + con_temp4 - con_temp3;
			deriva_I_Z((i) * (Z_size2 + 1) + j) = (I_u_2 - I_u_1 + I_u_3) / (-2 * DER_INCREMENTS_Z(i, j));
		}
	}
}
void Adam::calculate_A34() {
	Matrix<double, 3, 4> z_x_add14, z_y_add14, z_xx_add14, z_yy_add14, z_xy_add14;
	Matrix<double, 3, 4> z_x_add24, z_y_add24, z_xx_add24, z_yy_add24, z_xy_add24;
	Matrix<double, 3, 4> z_xo4, z_yo4, z_xxo4, z_yyo4, z_xyo4;
	Matrix<double, 3, 4> z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14;
	Matrix<double, 3, 4> z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24;
	Matrix<double, 3, 4> Alpha_new4, Beta_new4, P_new4;
	int ix = 0, jy = 0;
	Zeros(z_x_add14);
	Zeros(z_y_add14);
	Zeros(z_xx_add14);
	Zeros(z_xy_add14);
	Zeros(z_yy_add14);
	z_x_add14(1, 1) = DER_INCREMENTS_Z(1, 2) / 2; z_x_add14(1, 3) = -DER_INCREMENTS_Z(1, 2) / 2;
	z_y_add14(0, 2) = DER_INCREMENTS_Z(1, 2); z_y_add14(2, 2) = -DER_INCREMENTS_Z(1, 2) / 2;
	z_xx_add14(1, 0) = DER_INCREMENTS_Z(1, 2); z_xx_add14(1, 1) = DER_INCREMENTS_Z(1, 2); z_xx_add14(1, 2) = -2 * DER_INCREMENTS_Z(1, 2); z_xx_add14(1, 3) = DER_INCREMENTS_Z(1, 2);
	z_yy_add14(0, 2) = -2 * DER_INCREMENTS_Z(1, 2); z_yy_add14(1, 2) = -2 * DER_INCREMENTS_Z(1, 2); z_yy_add14(2, 2) = DER_INCREMENTS_Z(1, 2);
	z_xy_add14(2, 0) = -DER_INCREMENTS_Z(1, 2) / 4; z_xy_add14(2, 1) = -DER_INCREMENTS_Z(1, 2) / 4; z_xy_add14(2, 3) = DER_INCREMENTS_Z(1, 2) / 4;
	z_x_add24 = -z_x_add14; z_y_add24 = -z_y_add14; z_xx_add24 = -z_xx_add14; z_yy_add24 = -z_yy_add14; z_xy_add24 = -z_xy_add14;
	for (int i = 0; i <= 2; i++)
		for (int j = 0; j <= 3; j++) {
			z_xo4(i, j) = Z_x(i, j);
			z_yo4(i, j) = Z_y(i, j);
			z_xxo4(i, j) = Z_xx(i, j);
			z_yyo4(i, j) = Z_yy(i, j);
			z_xyo4(i, j) = Z_xy(i, j);
		}
	z_x_new14 = z_xo4 + z_x_add14; z_y_new14 = z_yo4 + z_y_add14; z_xx_new14 = z_xxo4 + z_xx_add14; z_yy_new14 = z_yyo4 + z_yy_add14; z_xy_new14 = z_xyo4 + z_xy_add14;
	z_x_new24 = z_xo4 + z_x_add24; z_y_new24 = z_yo4 + z_y_add24; z_xx_new24 = z_xxo4 + z_xx_add24; z_yy_new24 = z_yyo4 + z_yy_add24; z_xy_new24 = z_xyo4 + z_xy_add24;
	//Alpha_new = Alpha(1:4, 1 : 2); Beta_new = Beta(1:4, 1 : 2); P_new = P(1:4, 1 : 2);
	for (int i = 0; i <= 2; i++)
		for (int j = 0; j <= 3; j++) {
			Alpha_new4(i, j) = A(i, j);
			Beta_new4(i, j) = B(i, j);
			P_new4(i, j) = S(i, j);
		}
	double I_u_1 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14);
	double I_u_2 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24);
	deriva_I_Z(Z_size2 + 3) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(1, 2));
	//=============================================================================================================================//
	Zeros(z_x_add14);
	Zeros(z_y_add14);
	Zeros(z_xx_add14);
	Zeros(z_xy_add14);
	Zeros(z_yy_add14);
	z_x_add14(1, 0) = DER_INCREMENTS_Z(1, Z_size2 - 2) / 2; z_x_add14(1, 2) = -DER_INCREMENTS_Z(1, Z_size2 - 2) / 2;
	z_y_add14(0, 1) = DER_INCREMENTS_Z(1, Z_size2 - 2); z_y_add14(2, 1) = -DER_INCREMENTS_Z(1, Z_size2 - 2) / 2;
	z_xx_add14(1, 0) = DER_INCREMENTS_Z(1, Z_size2 - 2); z_xx_add14(1, 1) = -2 * DER_INCREMENTS_Z(1, Z_size2 - 2); 
	z_xx_add14(1, 2) = DER_INCREMENTS_Z(1, Z_size2 - 2); z_xx_add14(1, 3) = DER_INCREMENTS_Z(1, Z_size2 - 2);
	z_yy_add14(0, 1) = -2 * DER_INCREMENTS_Z(1, Z_size2 - 2); z_yy_add14(1, 1) = -2 * DER_INCREMENTS_Z(1, Z_size2 - 2); z_yy_add14(2, 1) = DER_INCREMENTS_Z(1, Z_size2 - 2);
	z_xy_add14(2, 0) = -DER_INCREMENTS_Z(1, Z_size2 - 2) / 4; z_xy_add14(2, 2) = DER_INCREMENTS_Z(1, Z_size2 - 2) / 4; z_xy_add14(2, 3) = DER_INCREMENTS_Z(1, Z_size2 - 2) / 4;
	z_x_add24 = -z_x_add14; z_y_add24 = -z_y_add14; z_xx_add24 = -z_xx_add14; z_yy_add24 = -z_yy_add14; z_xy_add24 = -z_xy_add14;
	for (int i = 0; i <= 2; i++) {
		for (int j = Z_size2 - 3; j <= Z_size2; j++) {
			z_xo4(i, jy) = Z_x(i, j);
			z_yo4(i, jy) = Z_y(i, j);
			z_xxo4(i, jy) = Z_xx(i, j);
			z_yyo4(i, jy) = Z_yy(i, j);
			z_xyo4(i, jy) = Z_xy(i, j);
			jy++;
		}
		jy = 0;
	}
	jy = 0;
	z_x_new14 = z_xo4 + z_x_add14; z_y_new14 = z_yo4 + z_y_add14; z_xx_new14 = z_xxo4 + z_xx_add14; z_yy_new14 = z_yyo4 + z_yy_add14; z_xy_new14 = z_xyo4 + z_xy_add14;
	z_x_new24 = z_xo4 + z_x_add24; z_y_new24 = z_yo4 + z_y_add24; z_xx_new24 = z_xxo4 + z_xx_add24; z_yy_new24 = z_yyo4 + z_yy_add24; z_xy_new24 = z_xyo4 + z_xy_add24;
	//Alpha_new = Alpha(1:4, 1 : 2); Beta_new = Beta(1:4, 1 : 2); P_new = P(1:4, 1 : 2);
	for (int i = 0; i <= 2; i++) {
		for (int j = Z_size2 - 3; j <= Z_size2; j++) {
			Alpha_new4(i, jy) = A(i, j);
			Beta_new4(i, jy) = B(i, j);
			P_new4(i, jy) = S(i, j);
			jy++;
		}
		jy = 0;
	}
	jy = 0;
	I_u_1 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14);
	I_u_2 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24);
	deriva_I_Z(Z_size2 + Z_size2 - 1) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(1, Z_size2 - 2));
	//=============================================================================================================================//
	Zeros(z_x_add14);
	Zeros(z_y_add14);
	Zeros(z_xx_add14);
	Zeros(z_xy_add14);
	Zeros(z_yy_add14);
	z_x_add14(1, 1) = DER_INCREMENTS_Z(Z_size1 - 1, 2) / 2; z_x_add14(1, 3) = -DER_INCREMENTS_Z(Z_size1 - 1, 2) / 2;
	z_y_add14(0, 2) = DER_INCREMENTS_Z(Z_size1 - 1, 2) / 2; z_y_add14(2, 2) = -DER_INCREMENTS_Z(Z_size1 - 1, 2);
	z_xx_add14(1, 0) = DER_INCREMENTS_Z(Z_size1 - 1, 2); z_xx_add14(1, 1) = DER_INCREMENTS_Z(Z_size1 - 1, 2); z_xx_add14(1, 2) = -2 * DER_INCREMENTS_Z(Z_size1 - 1, 2); z_xx_add14(1, 3) = DER_INCREMENTS_Z(Z_size1 - 1, 2);
	z_yy_add14(0, 2) = DER_INCREMENTS_Z(Z_size1 - 1, 2); z_yy_add14(1, 2) = -2 * DER_INCREMENTS_Z(Z_size1 - 1, 2); z_yy_add14(2, 2) = -2 * DER_INCREMENTS_Z(Z_size1 - 1, 2);
	z_xy_add14(0, 0) = DER_INCREMENTS_Z(Z_size1 - 1, 2) / 4; z_xy_add14(0, 1) = DER_INCREMENTS_Z(Z_size1 - 1, 2) / 4; z_xy_add14(0, 3) = -DER_INCREMENTS_Z(Z_size1 - 1, 2) / 4;
	z_x_add24 = -z_x_add14; z_y_add24 = -z_y_add14; z_xx_add24 = -z_xx_add14; z_yy_add24 = -z_yy_add14; z_xy_add24 = -z_xy_add14;
	for (int i = Z_size1 - 2; i <= Z_size1; i++) {
		for (int j = 0; j <= 3; j++) {
			z_xo4(ix, j) = Z_x(i, j);
			z_yo4(ix, j) = Z_y(i, j);
			z_xxo4(ix, j) = Z_xx(i, j);
			z_yyo4(ix, j) = Z_yy(i, j);
			z_xyo4(ix, j) = Z_xy(i, j);
		}
		ix++;
	}
	ix = 0;
	z_x_new14 = z_xo4 + z_x_add14; z_y_new14 = z_yo4 + z_y_add14; z_xx_new14 = z_xxo4 + z_xx_add14; z_yy_new14 = z_yyo4 + z_yy_add14; z_xy_new14 = z_xyo4 + z_xy_add14;
	z_x_new24 = z_xo4 + z_x_add24; z_y_new24 = z_yo4 + z_y_add24; z_xx_new24 = z_xxo4 + z_xx_add24; z_yy_new24 = z_yyo4 + z_yy_add24; z_xy_new24 = z_xyo4 + z_xy_add24;
	//Alpha_new = Alpha(1:4, 1 : 2); Beta_new = Beta(1:4, 1 : 2); P_new = P(1:4, 1 : 2);
	for (int i = Z_size1 - 2; i <= Z_size1; i++) {
		for (int j = 0; j <= 3; j++) {
			Alpha_new4(ix, j) = A(i, j);
			Beta_new4(ix, j) = B(i, j);
			P_new4(ix, j) = S(i, j);
		}
		ix++;
	}
	ix = 0;
	I_u_1 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14);
	I_u_2 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24);
	deriva_I_Z((Z_size1 - 1) * (Z_size2+1) + 2) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(Z_size1 - 1, 2));
	//=============================================================================================================================//
	Zeros(z_x_add14);
	Zeros(z_y_add14);
	Zeros(z_xx_add14);
	Zeros(z_xy_add14);
	Zeros(z_yy_add14);
	z_x_add14(1, 0) = DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 2) / 2; z_x_add14(1, 2) = -DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 2) / 2;
	z_y_add14(0, 1) = DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 2) / 2; z_y_add14(2, 1) = -DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 2);
	z_xx_add14(1, 0) = DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 2); z_xx_add14(1, 1) = -2 * DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 2); z_xx_add14(1, 2) = DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 2); 
	z_xx_add14(1, 3) = DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 2);
	z_yy_add14(0, 1) = DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 2); z_yy_add14(1, 1) = -2 * DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 2); z_yy_add14(2, 1) = -2 * DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 2);
	z_xy_add14(0, 0) = DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 2) / 4; z_xy_add14(0, 2) = -DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 2) / 4; z_xy_add14(0, 3) = -DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 2) / 4;
	z_x_add24 = -z_x_add14; z_y_add24 = -z_y_add14; z_xx_add24 = -z_xx_add14; z_yy_add24 = -z_yy_add14; z_xy_add24 = -z_xy_add14;
	for (int i = Z_size1 - 2; i <= Z_size1; i++) {
		for (int j = Z_size2 - 3; j <= Z_size2; j++) {
			z_xo4(ix, jy) = Z_x(i, j);
			z_yo4(ix, jy) = Z_y(i, j);
			z_xxo4(ix, jy) = Z_xx(i, j);
			z_yyo4(ix, jy) = Z_yy(i, j);
			z_xyo4(ix, jy) = Z_xy(i, j);
			jy++;
		}
		ix++; jy = 0;
	}
	ix = 0; jy = 0;
	z_x_new14 = z_xo4 + z_x_add14; z_y_new14 = z_yo4 + z_y_add14; z_xx_new14 = z_xxo4 + z_xx_add14; z_yy_new14 = z_yyo4 + z_yy_add14; z_xy_new14 = z_xyo4 + z_xy_add14;
	z_x_new24 = z_xo4 + z_x_add24; z_y_new24 = z_yo4 + z_y_add24; z_xx_new24 = z_xxo4 + z_xx_add24; z_yy_new24 = z_yyo4 + z_yy_add24; z_xy_new24 = z_xyo4 + z_xy_add24;
	//Alpha_new = Alpha(1:4, 1 : 2); Beta_new = Beta(1:4, 1 : 2); P_new = P(1:4, 1 : 2);
	for (int i = Z_size1 - 2; i <= Z_size1; i++) {
		for (int j = Z_size2 - 3; j <= Z_size2; j++) {
			Alpha_new4(ix, jy) = A(i, j);
			Beta_new4(ix, jy) = B(i, j);
			P_new4(ix, jy) = S(i, j);
			jy++;
		}
		ix++; jy = 0;
	}
	ix = 0; jy = 0;
	I_u_1 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14);
	I_u_2 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24);
	deriva_I_Z((Z_size1 - 1) * (Z_size2+1) + Z_size2 - 2) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(Z_size1 - 1, Z_size2 - 2));
	//=============================================================================================================================//
	for (int i = 3; i <= Z_size1 - 3; i++) {
		int j = 2;
		Zeros(z_x_add14);
		Zeros(z_y_add14);
		Zeros(z_xx_add14);
		Zeros(z_xy_add14);
		Zeros(z_yy_add14);
		z_x_add14(1, 1) = DER_INCREMENTS_Z(i, j) / 2; z_x_add14(1, 3) = -DER_INCREMENTS_Z(i, j) / 2;
		z_y_add14(0, 2) = DER_INCREMENTS_Z(i, j) / 2; z_y_add14(2, 2) = -DER_INCREMENTS_Z(i, j) / 2;
		z_xx_add14(1, 0) = DER_INCREMENTS_Z(i, j); z_xx_add14(1, 1) = DER_INCREMENTS_Z(i, j); z_xx_add14(1, 2) = -2 * DER_INCREMENTS_Z(i, j); z_xx_add14(1, 3) = DER_INCREMENTS_Z(i, j);
		z_yy_add14(0, 2) = DER_INCREMENTS_Z(i, j); z_yy_add14(1, 2) = -2 * DER_INCREMENTS_Z(i, j); z_yy_add14(2, 2) = DER_INCREMENTS_Z(i, j);
		z_xy_add14(0, 0) = DER_INCREMENTS_Z(i, j) / 4; z_xy_add14(0, 1) = DER_INCREMENTS_Z(i, j) / 4; z_xy_add14(0, 3) = -DER_INCREMENTS_Z(i, j) / 4;
		z_xy_add14(2, 0) = -DER_INCREMENTS_Z(i, j) / 4; z_xy_add14(2, 1) = -DER_INCREMENTS_Z(i, j) / 4; z_xy_add14(2, 3) = DER_INCREMENTS_Z(i, j) / 4;
		z_x_add24 = -z_x_add14; z_y_add24 = -z_y_add14; z_xx_add24 = -z_xx_add14; z_yy_add24 = -z_yy_add14; z_xy_add24 = -z_xy_add14;
		for (int ii = i - 1; ii <= i + 1; ii++) {
			for (int jj = j - 2; jj <= j + 1; jj++) {
				z_xo4(ix, jy) = Z_x(ii, jj);
				z_yo4(ix, jy) = Z_y(ii, jj);
				z_xxo4(ix, jy) = Z_xx(ii, jj);
				z_yyo4(ix, jy) = Z_yy(ii, jj);
				z_xyo4(ix, jy) = Z_xy(ii, jj);
				jy++;
			}
			ix++; jy = 0;
		}
		ix = 0; jy = 0;
		z_x_new14 = z_xo4 + z_x_add14; z_y_new14 = z_yo4 + z_y_add14; z_xx_new14 = z_xxo4 + z_xx_add14; z_yy_new14 = z_yyo4 + z_yy_add14; z_xy_new14 = z_xyo4 + z_xy_add14;
		z_x_new24 = z_xo4 + z_x_add24; z_y_new24 = z_yo4 + z_y_add24; z_xx_new24 = z_xxo4 + z_xx_add24; z_yy_new24 = z_yyo4 + z_yy_add24; z_xy_new24 = z_xyo4 + z_xy_add24;
		for (int ii = i - 1; ii <= i + 1; ii++) {
			for (int jj = j - 2; jj <= j + 1; jj++) {
				Alpha_new4(ix, jy) = A(ii, jj);
				Beta_new4(ix, jy) = B(ii, jj);
				P_new4(ix, jy) = S(ii, jj);
				jy++;
			}
			ix++; jy = 0;
		}
		ix = 0; jy = 0;
		I_u_1 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14);
		I_u_2 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24);
		deriva_I_Z((i) * (Z_size2+1) + j) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(i, j));
	}
	//=============================================================================================================================//
	for (int i = 3; i <= Z_size1 - 3; i++) {
		int j = Z_size2 - 2;
		Zeros(z_x_add14);
		Zeros(z_y_add14);
		Zeros(z_xx_add14);
		Zeros(z_xy_add14);
		Zeros(z_yy_add14);
		z_x_add14(1, 0) = DER_INCREMENTS_Z(i, j) / 2; z_x_add14(1, 2) = -DER_INCREMENTS_Z(i, j) / 2;
		z_y_add14(0, 1) = DER_INCREMENTS_Z(i, j) / 2; z_y_add14(2, 1) = -DER_INCREMENTS_Z(i, j) / 2;
		z_xx_add14(1, 0) = DER_INCREMENTS_Z(i, j); z_xx_add14(1, 1) = -2 * DER_INCREMENTS_Z(i, j); 
		z_xx_add14(1, 2) = DER_INCREMENTS_Z(i, j); z_xx_add14(1, 3) = DER_INCREMENTS_Z(i, j);
		z_yy_add14(0, 1) = DER_INCREMENTS_Z(i, j); z_yy_add14(1, 1) = -2 * DER_INCREMENTS_Z(i, j); z_yy_add14(2, 1) = DER_INCREMENTS_Z(i, j);
		z_xy_add14(0, 0) = DER_INCREMENTS_Z(i, j) / 4; z_xy_add14(0, 2) = -DER_INCREMENTS_Z(i, j) / 4; z_xy_add14(0, 3) = -DER_INCREMENTS_Z(i, j) / 4;
		z_xy_add14(2, 0) = -DER_INCREMENTS_Z(i, j) / 4; z_xy_add14(2, 2) = DER_INCREMENTS_Z(i, j) / 4; z_xy_add14(2, 3) = DER_INCREMENTS_Z(i, j) / 4;
		z_x_add24 = -z_x_add14; z_y_add24 = -z_y_add14; z_xx_add24 = -z_xx_add14; z_yy_add24 = -z_yy_add14; z_xy_add24 = -z_xy_add14;
		for (int ii = i - 1; ii <= i + 1; ii++) {
			for (int jj = j - 1; jj <= j + 2; jj++) {
				z_xo4(ix, jy) = Z_x(ii, jj);
				z_yo4(ix, jy) = Z_y(ii, jj);
				z_xxo4(ix, jy) = Z_xx(ii, jj);
				z_yyo4(ix, jy) = Z_yy(ii, jj);
				z_xyo4(ix, jy) = Z_xy(ii, jj);
				jy++;
			}
			ix++; jy = 0;
		}
		ix = 0; jy = 0;
		z_x_new14 = z_xo4 + z_x_add14; z_y_new14 = z_yo4 + z_y_add14; z_xx_new14 = z_xxo4 + z_xx_add14; z_yy_new14 = z_yyo4 + z_yy_add14; z_xy_new14 = z_xyo4 + z_xy_add14;
		z_x_new24 = z_xo4 + z_x_add24; z_y_new24 = z_yo4 + z_y_add24; z_xx_new24 = z_xxo4 + z_xx_add24; z_yy_new24 = z_yyo4 + z_yy_add24; z_xy_new24 = z_xyo4 + z_xy_add24;
		for (int ii = i - 1; ii <= i + 1; ii++) {
			for (int jj = j - 1; jj <= j + 2; jj++) {
				Alpha_new4(ix, jy) = A(ii, jj);
				Beta_new4(ix, jy) = B(ii, jj);
				P_new4(ix, jy) = S(ii, jj);
				jy++;
			}
			ix++; jy = 0;
		}
		ix = 0; jy = 0;
		I_u_1 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14);
		I_u_2 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24);
		deriva_I_Z((i) * (Z_size2+1) + j) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(i, j));
	}
}
void Adam::calculate_A43() {
	Matrix<double, 4, 3> z_x_add14, z_y_add14, z_xx_add14, z_yy_add14, z_xy_add14;
	Matrix<double, 4, 3> z_x_add24, z_y_add24, z_xx_add24, z_yy_add24, z_xy_add24;
	Matrix<double, 4, 3> z_xo4, z_yo4, z_xxo4, z_yyo4, z_xyo4;
	Matrix<double, 4, 3> z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14;
	Matrix<double, 4, 3> z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24;
	Matrix<double, 4, 3> Alpha_new4, Beta_new4, P_new4;
	int ix = 0, jy = 0;
	Zeros(z_x_add14);
	Zeros(z_y_add14);
	Zeros(z_xx_add14);
	Zeros(z_xy_add14);
	Zeros(z_yy_add14);
	z_x_add14(2, 0)= DER_INCREMENTS_Z(2, 1); z_x_add14(2, 2) = -DER_INCREMENTS_Z(2, 1) / 2;
	z_y_add14(1, 1) = DER_INCREMENTS_Z(2, 1) / 2; z_y_add14(3, 1) = -DER_INCREMENTS_Z(2, 1) / 2;
	z_xx_add14(2, 0) = -2 * DER_INCREMENTS_Z(2, 1); z_xx_add14(2, 1) = -2 * DER_INCREMENTS_Z(2, 1); z_xx_add14(2, 2) = DER_INCREMENTS_Z(2, 1);
	z_yy_add14(0, 1) = DER_INCREMENTS_Z(2, 1); z_yy_add14(1, 1) = DER_INCREMENTS_Z(2, 1); z_yy_add14(2, 1) = -2 * DER_INCREMENTS_Z(2, 1); z_yy_add14(3, 1) = DER_INCREMENTS_Z(2, 1);
	z_xy_add14(0, 2) = -DER_INCREMENTS_Z(2, 1) / 4; z_xy_add14(1, 2) = -DER_INCREMENTS_Z(2, 1) / 4; z_xy_add14(3, 2) = DER_INCREMENTS_Z(2, 1) / 4;
	z_x_add24 = -z_x_add14; z_y_add24 = -z_y_add14; z_xx_add24 = -z_xx_add14; z_yy_add24 = -z_yy_add14; z_xy_add24 = -z_xy_add14;
	for (int i = 0; i <= 3; i++) {
		for (int j = 0; j <= 2; j++) {
			z_xo4(i, j) = Z_x(i, j);
			z_yo4(i, j) = Z_y(i, j);
			z_xxo4(i, j) = Z_xx(i, j);
			z_yyo4(i, j) = Z_yy(i, j);
			z_xyo4(i, j) = Z_xy(i, j);
		}
	}
	z_x_new14 = z_xo4 + z_x_add14; z_y_new14 = z_yo4 + z_y_add14; z_xx_new14 = z_xxo4 + z_xx_add14; z_yy_new14 = z_yyo4 + z_yy_add14; z_xy_new14 = z_xyo4 + z_xy_add14;
	z_x_new24 = z_xo4 + z_x_add24; z_y_new24 = z_yo4 + z_y_add24; z_xx_new24 = z_xxo4 + z_xx_add24; z_yy_new24 = z_yyo4 + z_yy_add24; z_xy_new24 = z_xyo4 + z_xy_add24;
	//Alpha_new = Alpha(1:4, 1 : 2); Beta_new = Beta(1:4, 1 : 2); P_new = P(1:4, 1 : 2);
	for (int i = 0; i <= 3; i++) {
		for (int j = 0; j <= 2; j++) {
			Alpha_new4(i, j) = A(i, j);
			Beta_new4(i, j) = B(i, j);
			P_new4(i, j) = S(i, j);
		}
	}
	double I_u_1 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14);
	double I_u_2 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24);
	deriva_I_Z(2 * (Z_size2+1) + 1) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(2, 1));
	//=============================================================================================================================//
	Zeros(z_x_add14);
	Zeros(z_y_add14);
	Zeros(z_xx_add14);
	Zeros(z_xy_add14);
	Zeros(z_yy_add14);
	z_x_add14(1, 0) = DER_INCREMENTS_Z(Z_size1 - 2, 1); z_x_add14(1, 2) = -DER_INCREMENTS_Z(Z_size1 - 2, 1) / 2;
	z_y_add14(0, 1) = DER_INCREMENTS_Z(Z_size1 - 2, 1) / 2; z_y_add14(2, 1) = -DER_INCREMENTS_Z(Z_size1 - 2, 1) / 2;
	z_xx_add14(1, 0) = -2 * DER_INCREMENTS_Z(Z_size1 - 2, 1); z_xx_add14(1, 1) = -2 * DER_INCREMENTS_Z(Z_size1 - 2, 1); z_xx_add14(1, 2) = DER_INCREMENTS_Z(Z_size1 - 2, 1);
	z_yy_add14(0, 1) = DER_INCREMENTS_Z(Z_size1 - 2, 1); z_yy_add14(1, 1) = -2 * DER_INCREMENTS_Z(Z_size1 - 2, 1); z_yy_add14(2, 1) = DER_INCREMENTS_Z(Z_size1 - 2, 1); 
	z_yy_add14(3, 1) = DER_INCREMENTS_Z(Z_size1 - 2, 1);
	z_xy_add14(0, 2) = -DER_INCREMENTS_Z(Z_size1 - 2, 1) / 4; z_xy_add14(2, 2) = DER_INCREMENTS_Z(Z_size1 - 2, 1) / 4; z_xy_add14(3, 2) = DER_INCREMENTS_Z(Z_size1 - 2, 1) / 4;
	z_x_add24 = -z_x_add14; z_y_add24 = -z_y_add14; z_xx_add24 = -z_xx_add14; z_yy_add24 = -z_yy_add14; z_xy_add24 = -z_xy_add14;
	for (int i = Z_size1 - 3; i <= Z_size1; i++) {
		for (int j = 0; j <= 2; j++) {
			z_xo4(ix, j) = Z_x(i, j);
			z_yo4(ix, j) = Z_y(i, j);
			z_xxo4(ix, j) = Z_xx(i, j);
			z_yyo4(ix, j) = Z_yy(i, j);
			z_xyo4(ix, j) = Z_xy(i, j);
		}
		ix++;
	}
	ix = 0;
	z_x_new14 = z_xo4 + z_x_add14; z_y_new14 = z_yo4 + z_y_add14; z_xx_new14 = z_xxo4 + z_xx_add14; z_yy_new14 = z_yyo4 + z_yy_add14; z_xy_new14 = z_xyo4 + z_xy_add14;
	z_x_new24 = z_xo4 + z_x_add24; z_y_new24 = z_yo4 + z_y_add24; z_xx_new24 = z_xxo4 + z_xx_add24; z_yy_new24 = z_yyo4 + z_yy_add24; z_xy_new24 = z_xyo4 + z_xy_add24;
	//Alpha_new = Alpha(1:4, 1 : 2); Beta_new = Beta(1:4, 1 : 2); P_new = P(1:4, 1 : 2);
	for (int i = Z_size1 - 3; i <= Z_size1; i++) {
		for (int j = 0; j <= 2; j++) {
			Alpha_new4(ix, j) = A(i, j);
			Beta_new4(ix, j) = B(i, j);
			P_new4(ix, j) = S(i, j);
		}
		ix++;
	}
	ix = 0;
	I_u_1 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14);
	I_u_2 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24);
	deriva_I_Z((Z_size1 - 2) * (Z_size2+1) + 1) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(Z_size1 - 2, 1));
	//=============================================================================================================================//
	Zeros(z_x_add14);
	Zeros(z_y_add14);
	Zeros(z_xx_add14);
	Zeros(z_xy_add14);
	Zeros(z_yy_add14);
	z_x_add14(2, 0) = DER_INCREMENTS_Z(2, Z_size2 - 1) / 2; z_x_add14(2, 2) = -DER_INCREMENTS_Z(2, Z_size2 - 1);
	z_y_add14(1, 1) = DER_INCREMENTS_Z(2, Z_size2 - 1) / 2; z_y_add14(3, 1) = -DER_INCREMENTS_Z(2, Z_size2 - 1) / 2;
	z_xx_add14(2, 0) = DER_INCREMENTS_Z(2, Z_size2 - 1); z_xx_add14(2, 1) = -2 * DER_INCREMENTS_Z(2, Z_size2 - 1); z_xx_add14(2, 2) = -2 * DER_INCREMENTS_Z(2, Z_size2 - 1);
	z_yy_add14(0, 1) = DER_INCREMENTS_Z(2, Z_size2 - 1); z_yy_add14(1, 1) = DER_INCREMENTS_Z(2, Z_size2 - 1); z_yy_add14(2, 1) = -2 * DER_INCREMENTS_Z(2, Z_size2 - 1);
	z_yy_add14(3, 1) = DER_INCREMENTS_Z(2, Z_size2 - 1);
	z_xy_add14(0, 0) = DER_INCREMENTS_Z(2, Z_size2 - 1) / 4; z_xy_add14(1, 0) = DER_INCREMENTS_Z(2, Z_size2 - 1) / 4; z_xy_add14(3, 0) = -DER_INCREMENTS_Z(2, Z_size2 - 1) / 4;
	z_x_add24 = -z_x_add14; z_y_add24 = -z_y_add14; z_xx_add24 = -z_xx_add14; z_yy_add24 = -z_yy_add14; z_xy_add24 = -z_xy_add14;
	for (int i = 0; i <= 3; i++) {
		for (int j = Z_size2 - 2; j <= Z_size2; j++) {
			z_xo4(i, jy) = Z_x(i, j);
			z_yo4(i, jy) = Z_y(i, j);
			z_xxo4(i, jy) = Z_xx(i, j);
			z_yyo4(i, jy) = Z_yy(i, j);
			z_xyo4(i, jy) = Z_xy(i, j);
			jy++;
		}
		jy = 0;
	}
	jy = 0;
	z_x_new14 = z_xo4 + z_x_add14; z_y_new14 = z_yo4 + z_y_add14; z_xx_new14 = z_xxo4 + z_xx_add14; z_yy_new14 = z_yyo4 + z_yy_add14; z_xy_new14 = z_xyo4 + z_xy_add14;
	z_x_new24 = z_xo4 + z_x_add24; z_y_new24 = z_yo4 + z_y_add24; z_xx_new24 = z_xxo4 + z_xx_add24; z_yy_new24 = z_yyo4 + z_yy_add24; z_xy_new24 = z_xyo4 + z_xy_add24;
	for (int i = 0; i <= 3; i++) {
		for (int j = Z_size2 - 2; j <= Z_size2; j++) {
			Alpha_new4(i, jy) = A(i, j);
			Beta_new4(i, jy) = B(i, j);
			P_new4(i, jy) = S(i, j);
			jy++;
		}
		jy = 0;
	}
	jy = 0;
	I_u_1 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14);
	I_u_2 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24);
	deriva_I_Z(2 * (Z_size2+1) + Z_size2 - 1) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(2, Z_size2 - 1));
	//=============================================================================================================================//
	Zeros(z_x_add14);
	Zeros(z_y_add14);
	Zeros(z_xx_add14);
	Zeros(z_xy_add14);
	Zeros(z_yy_add14);
	z_x_add14(1, 0) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 1) / 2; z_x_add14(1, 2) = -DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 1);
	z_y_add14(0, 1) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 1) / 2; z_y_add14(2, 1) = -DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 1) / 2;
	z_xx_add14(1, 0) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 1); z_xx_add14(1, 1) = -2 * DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 1); z_xx_add14(1, 2) = -2 * DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 1);
	z_yy_add14(0, 1) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 1); z_yy_add14(1, 1) = -2 * DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 1); 
	z_yy_add14(2, 1) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 1); z_yy_add14(3, 1) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 1);
	z_xy_add14(0, 0) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 1) / 4; z_xy_add14(2, 0) = -DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 1) / 4; z_xy_add14(3, 0) = -DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 1) / 4;
	z_x_add24 = -z_x_add14; z_y_add24 = -z_y_add14; z_xx_add24 = -z_xx_add14; z_yy_add24 = -z_yy_add14; z_xy_add24 = -z_xy_add14;
	for (int i = Z_size1 - 3; i <= Z_size1; i++) {
		for (int j = Z_size2 - 2; j <= Z_size2; j++) {
			z_xo4(ix, jy) = Z_x(i, j);
			z_yo4(ix, jy) = Z_y(i, j);
			z_xxo4(ix, jy) = Z_xx(i, j);
			z_yyo4(ix, jy) = Z_yy(i, j);
			z_xyo4(ix, jy) = Z_xy(i, j);
			jy++;
		}
		ix++; jy = 0;
	}
	ix = 0; jy = 0;
	z_x_new14 = z_xo4 + z_x_add14; z_y_new14 = z_yo4 + z_y_add14; z_xx_new14 = z_xxo4 + z_xx_add14; z_yy_new14 = z_yyo4 + z_yy_add14; z_xy_new14 = z_xyo4 + z_xy_add14;
	z_x_new24 = z_xo4 + z_x_add24; z_y_new24 = z_yo4 + z_y_add24; z_xx_new24 = z_xxo4 + z_xx_add24; z_yy_new24 = z_yyo4 + z_yy_add24; z_xy_new24 = z_xyo4 + z_xy_add24;
	for (int i = Z_size1 - 3; i <= Z_size1; i++) {
		for (int j = Z_size2 - 2; j <= Z_size2; j++) {
			Alpha_new4(ix, jy) = A(i, j);
			Beta_new4(ix, jy) = B(i, j);
			P_new4(ix, jy) = S(i, j);
			jy++;
		}
		ix++; jy = 0;
	}
	ix = 0; jy = 0;
	I_u_1 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14);
	I_u_2 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24);
	deriva_I_Z((Z_size1 - 2) * (Z_size2+1) + Z_size2 - 1) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 1));
	//=============================================================================================================================//
	for (int j = 3; j <= Z_size2 - 3; j++) {
		int i = 2;
		Zeros(z_x_add14);
		Zeros(z_y_add14);
		Zeros(z_xx_add14);
		Zeros(z_xy_add14);
		Zeros(z_yy_add14);
		z_x_add14(2, 0) = DER_INCREMENTS_Z(i, j) / 2; z_x_add14(2, 2) = -DER_INCREMENTS_Z(i, j) / 2;
		z_y_add14(1, 1) = DER_INCREMENTS_Z(i, j) / 2; z_y_add14(3, 1) = -DER_INCREMENTS_Z(i, j) / 2;
		z_xx_add14(2, 0) = DER_INCREMENTS_Z(i, j); z_xx_add14(2, 1) = -2 * DER_INCREMENTS_Z(i, j); z_xx_add14(2, 2) = DER_INCREMENTS_Z(i, j);
		z_yy_add14(0, 1) = DER_INCREMENTS_Z(i, j); z_yy_add14(1, 1) = DER_INCREMENTS_Z(i, j); z_yy_add14(2, 1) = -2 * DER_INCREMENTS_Z(i, j); z_yy_add14(3, 1) = DER_INCREMENTS_Z(i, j);
		z_xy_add14(0, 0) = DER_INCREMENTS_Z(i, j) / 4; z_xy_add14(1, 0) = DER_INCREMENTS_Z(i, j) / 4; z_xy_add14(3, 0) = -DER_INCREMENTS_Z(i, j) / 4;
		z_xy_add14(0, 2) = -DER_INCREMENTS_Z(i, j) / 4; z_xy_add14(1, 2) = -DER_INCREMENTS_Z(i, j) / 4; z_xy_add14(3, 2) = DER_INCREMENTS_Z(i, j) / 4;
		z_x_add24 = -z_x_add14; z_y_add24 = -z_y_add14; z_xx_add24 = -z_xx_add14; z_yy_add24 = -z_yy_add14; z_xy_add24 = -z_xy_add14;
		for (int ii = i - 2; ii <= i + 1; ii++) {
			for (int jj = j - 1; jj <= j + 1; jj++) {
				z_xo4(ix, jy) = Z_x(ii, jj);
				z_yo4(ix, jy) = Z_y(ii, jj);
				z_xxo4(ix, jy) = Z_xx(ii, jj);
				z_yyo4(ix, jy) = Z_yy(ii, jj);
				z_xyo4(ix, jy) = Z_xy(ii, jj);
				jy++;
			}
			ix++; jy = 0;
		}
		ix = 0; jy = 0;
		z_x_new14 = z_xo4 + z_x_add14; z_y_new14 = z_yo4 + z_y_add14; z_xx_new14 = z_xxo4 + z_xx_add14; z_yy_new14 = z_yyo4 + z_yy_add14; z_xy_new14 = z_xyo4 + z_xy_add14;
		z_x_new24 = z_xo4 + z_x_add24; z_y_new24 = z_yo4 + z_y_add24; z_xx_new24 = z_xxo4 + z_xx_add24; z_yy_new24 = z_yyo4 + z_yy_add24; z_xy_new24 = z_xyo4 + z_xy_add24;
		for (int ii = i - 2; ii <= i + 1; ii++) {
			for (int jj = j - 1; jj <= j + 1; jj++) {
				Alpha_new4(ix, jy) = A(ii, jj);
				Beta_new4(ix, jy) = B(ii, jj);
				P_new4(ix, jy) = S(ii, jj);
				jy++;
			}
			ix++; jy = 0;
		}
		ix = 0; jy = 0;
		I_u_1 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14);
		I_u_2 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24);
		deriva_I_Z((i) * (Z_size2+1) + j) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(i, j));
	}
	//=============================================================================================================================//
	for (int j = 3; j <= Z_size2 - 3; j++) {
		int i = Z_size1 - 2;
		Zeros(z_x_add14);
		Zeros(z_y_add14);
		Zeros(z_xx_add14);
		Zeros(z_xy_add14);
		Zeros(z_yy_add14);
		z_x_add14(1, 0) = DER_INCREMENTS_Z(i, j) / 2; z_x_add14(1, 2) = -DER_INCREMENTS_Z(i, j) / 2;
		z_y_add14(0, 1) = DER_INCREMENTS_Z(i, j) / 2; z_y_add14(2, 1) = -DER_INCREMENTS_Z(i, j) / 2;
		z_xx_add14(1, 0) = DER_INCREMENTS_Z(i, j); z_xx_add14(1, 1) = -2 * DER_INCREMENTS_Z(i, j); z_xx_add14(1, 2) = DER_INCREMENTS_Z(i, j);
		z_yy_add14(0, 1) = DER_INCREMENTS_Z(i, j); z_yy_add14(1, 1) = -2 * DER_INCREMENTS_Z(i, j); z_yy_add14(2, 1) = DER_INCREMENTS_Z(i, j); z_yy_add14(3, 1) = DER_INCREMENTS_Z(i, j);
		z_xy_add14(0, 0) = DER_INCREMENTS_Z(i, j) / 4; z_xy_add14(2, 0) = -DER_INCREMENTS_Z(i, j) / 4; z_xy_add14(3, 0) = -DER_INCREMENTS_Z(i, j) / 4;
		z_xy_add14(0, 2) = -DER_INCREMENTS_Z(i, j) / 4; z_xy_add14(2, 2) = DER_INCREMENTS_Z(i, j) / 4; z_xy_add14(3, 2) = DER_INCREMENTS_Z(i, j) / 4;
		z_x_add24 = -z_x_add14; z_y_add24 = -z_y_add14; z_xx_add24 = -z_xx_add14; z_yy_add24 = -z_yy_add14; z_xy_add24 = -z_xy_add14;
		for (int ii = i - 1; ii <= i + 2; ii++) {
			for (int jj = j - 1; jj <= j + 1; jj++) {
				z_xo4(ix, jy) = Z_x(ii, jj);
				z_yo4(ix, jy) = Z_y(ii, jj);
				z_xxo4(ix, jy) = Z_xx(ii, jj);
				z_yyo4(ix, jy) = Z_yy(ii, jj);
				z_xyo4(ix, jy) = Z_xy(ii, jj);
				jy++;
			}
			ix++; jy = 0;
		}
		ix = 0; jy = 0;
		z_x_new14 = z_xo4 + z_x_add14; z_y_new14 = z_yo4 + z_y_add14; z_xx_new14 = z_xxo4 + z_xx_add14; z_yy_new14 = z_yyo4 + z_yy_add14; z_xy_new14 = z_xyo4 + z_xy_add14;
		z_x_new24 = z_xo4 + z_x_add24; z_y_new24 = z_yo4 + z_y_add24; z_xx_new24 = z_xxo4 + z_xx_add24; z_yy_new24 = z_yyo4 + z_yy_add24; z_xy_new24 = z_xyo4 + z_xy_add24;
		for (int ii = i - 1; ii <= i + 2; ii++) {
			for (int jj = j - 1; jj <= j + 1; jj++) {
				Alpha_new4(ix, jy) = A(ii, jj);
				Beta_new4(ix, jy) = B(ii, jj);
				P_new4(ix, jy) = S(ii, jj);
				jy++;
			}
			ix++; jy = 0;
		}
		ix = 0; jy = 0;
		I_u_1 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14);
		I_u_2 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24);
		deriva_I_Z((i) * (Z_size2+1) + j) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(i, j));
	}
}
void Adam::calculate_A44() {
	Matrix<double, 4, 4> z_x_add14, z_y_add14, z_xx_add14, z_yy_add14, z_xy_add14;
	Matrix<double, 4, 4> z_x_add24, z_y_add24, z_xx_add24, z_yy_add24, z_xy_add24;
	Matrix<double, 4, 4> z_xo4, z_yo4, z_xxo4, z_yyo4, z_xyo4;
	Matrix<double, 4, 4> z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14;
	Matrix<double, 4, 4> z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24;
	Matrix<double, 4, 4> Alpha_new4, Beta_new4, P_new4;
	int ix = 0, jy = 0;
	Zeros(z_x_add14);
	Zeros(z_y_add14);
	Zeros(z_xx_add14);
	Zeros(z_xy_add14);
	Zeros(z_yy_add14);
	z_x_add14(2, 1) = DER_INCREMENTS_Z(2, 2) / 2; z_x_add14(2, 3) = -DER_INCREMENTS_Z(2, 2) / 2;
	z_y_add14(1, 2) = DER_INCREMENTS_Z(2, 2) / 2; z_y_add14(3, 2) = -DER_INCREMENTS_Z(2, 2) / 2;
	z_xx_add14(2, 0) = DER_INCREMENTS_Z(2, 2); z_xx_add14(2, 1) = DER_INCREMENTS_Z(2, 2); z_xx_add14(2, 2) = -2 * DER_INCREMENTS_Z(2, 2); z_xx_add14(2, 3) = DER_INCREMENTS_Z(2, 2);
	z_yy_add14(0, 2) = DER_INCREMENTS_Z(2, 2); z_yy_add14(1, 2) = DER_INCREMENTS_Z(2, 2); z_yy_add14(2, 2) = -2 * DER_INCREMENTS_Z(2, 2); z_yy_add14(3, 2) = DER_INCREMENTS_Z(2, 2);
	z_xy_add14(0, 0) = DER_INCREMENTS_Z(2, 2) / 4; z_xy_add14(0, 1) = DER_INCREMENTS_Z(2, 2) / 4; z_xy_add14(0, 3) = -DER_INCREMENTS_Z(2, 2) / 4;
	z_xy_add14(1, 0) = DER_INCREMENTS_Z(2, 2) / 4; z_xy_add14(1, 1) = DER_INCREMENTS_Z(2, 2) / 4; z_xy_add14(1, 3) = -DER_INCREMENTS_Z(2, 2) / 4;
	z_xy_add14(3, 0) = -DER_INCREMENTS_Z(2, 2) / 4; z_xy_add14(3, 1) = -DER_INCREMENTS_Z(2, 2) / 4; z_xy_add14(3, 3) = DER_INCREMENTS_Z(2, 2) / 4;
	z_x_add24 = -z_x_add14; z_y_add24 = -z_y_add14; z_xx_add24 = -z_xx_add14; z_yy_add24 = -z_yy_add14; z_xy_add24 = -z_xy_add14;
	for (int i = 0; i <= 3; i++) {
		for (int j = 0; j <= 3; j++) {
			z_xo4(i, j) = Z_x(i, j);
			z_yo4(i, j) = Z_y(i, j);
			z_xxo4(i, j) = Z_xx(i, j);
			z_yyo4(i, j) = Z_yy(i, j);
			z_xyo4(i, j) = Z_xy(i, j);
		}
	}
	z_x_new14 = z_xo4 + z_x_add14; z_y_new14 = z_yo4 + z_y_add14; z_xx_new14 = z_xxo4 + z_xx_add14; z_yy_new14 = z_yyo4 + z_yy_add14; z_xy_new14 = z_xyo4 + z_xy_add14;
	z_x_new24 = z_xo4 + z_x_add24; z_y_new24 = z_yo4 + z_y_add24; z_xx_new24 = z_xxo4 + z_xx_add24; z_yy_new24 = z_yyo4 + z_yy_add24; z_xy_new24 = z_xyo4 + z_xy_add24;
	for (int i = 0; i <= 3; i++) {
		for (int j = 0; j <= 3; j++) {
			Alpha_new4(i, j) = A(i, j);
			Beta_new4(i, j) = B(i, j);
			P_new4(i, j) = S(i, j);
		}
	}
	double I_u_1 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14);
	double I_u_2 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24);
	deriva_I_Z(2 * (Z_size2+1) + 2) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(2, 2));
	//=============================================================================================================================//
	Zeros(z_x_add14);
	Zeros(z_y_add14);
	Zeros(z_xx_add14);
	Zeros(z_xy_add14);
	Zeros(z_yy_add14);
	z_x_add14(2, 0) = DER_INCREMENTS_Z(2, Z_size2 - 2) / 2; z_x_add14(2, 2) = -DER_INCREMENTS_Z(2, Z_size2 - 2) / 2;
	z_y_add14(1, 1) = DER_INCREMENTS_Z(2, Z_size2 - 2) / 2; z_y_add14(3, 1) = -DER_INCREMENTS_Z(2, Z_size2 - 2) / 2;
	z_xx_add14(2, 0) = DER_INCREMENTS_Z(2, Z_size2 - 2); z_xx_add14(2, 1) = -2 * DER_INCREMENTS_Z(2, Z_size2 - 2); 
	z_xx_add14(2, 2) = DER_INCREMENTS_Z(2, Z_size2 - 2); z_xx_add14(2, 3) = DER_INCREMENTS_Z(2, Z_size2 - 2);
	z_yy_add14(0, 1) = DER_INCREMENTS_Z(2, Z_size2 - 2); z_yy_add14(1, 1) = DER_INCREMENTS_Z(2, Z_size2 - 2); 
	z_yy_add14(2, 1) = -2 * DER_INCREMENTS_Z(2, Z_size2 - 2); z_yy_add14(3, 1) = DER_INCREMENTS_Z(2, Z_size2 - 2);
	z_xy_add14(0, 0) = DER_INCREMENTS_Z(2, Z_size2 - 2) / 4; z_xy_add14(0, 2) = -DER_INCREMENTS_Z(2, Z_size2 - 2) / 4; z_xy_add14(0, 3) = -DER_INCREMENTS_Z(2, Z_size2 - 2) / 4;
	z_xy_add14(1, 0) = DER_INCREMENTS_Z(2, Z_size2 - 2) / 4; z_xy_add14(1, 2) = -DER_INCREMENTS_Z(2, Z_size2 - 2) / 4; z_xy_add14(1, 3) = -DER_INCREMENTS_Z(2, Z_size2 - 2) / 4;
	z_xy_add14(3, 0) = -DER_INCREMENTS_Z(2, Z_size2 - 2) / 4; z_xy_add14(3, 2) = DER_INCREMENTS_Z(2, Z_size2 - 2) / 4; z_xy_add14(3, 3) = DER_INCREMENTS_Z(2, Z_size2 - 2) / 4;
	z_x_add24 = -z_x_add14; z_y_add24 = -z_y_add14; z_xx_add24 = -z_xx_add14; z_yy_add24 = -z_yy_add14; z_xy_add24 = -z_xy_add14;
	for (int i = 0; i <= 3; i++) {
		for (int j = Z_size2 - 3; j <= Z_size2; j++) {
			z_xo4(i, jy) = Z_x(i, j);
			z_yo4(i, jy) = Z_y(i, j);
			z_xxo4(i, jy) = Z_xx(i, j);
			z_yyo4(i, jy) = Z_yy(i, j);
			z_xyo4(i, jy) = Z_xy(i, j);
			jy++;
		}
		jy = 0;
	}
	jy = 0;
	z_x_new14 = z_xo4 + z_x_add14; z_y_new14 = z_yo4 + z_y_add14; z_xx_new14 = z_xxo4 + z_xx_add14; z_yy_new14 = z_yyo4 + z_yy_add14; z_xy_new14 = z_xyo4 + z_xy_add14;
	z_x_new24 = z_xo4 + z_x_add24; z_y_new24 = z_yo4 + z_y_add24; z_xx_new24 = z_xxo4 + z_xx_add24; z_yy_new24 = z_yyo4 + z_yy_add24; z_xy_new24 = z_xyo4 + z_xy_add24;
	for (int i = 0; i <= 3; i++) {
		for (int j = Z_size2 - 3; j <= Z_size2; j++) {
			Alpha_new4(i, jy) = A(i, j);
			Beta_new4(i, jy) = B(i, j);
			P_new4(i, jy) = S(i, j);
			jy++;
		}
		jy = 0;
	}
	jy = 0;
	I_u_1 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14);
	I_u_2 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24);
	deriva_I_Z(2 * (Z_size2+1) + Z_size2 - 2) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(2, Z_size2 - 2));
	//=============================================================================================================================//
	Zeros(z_x_add14);
	Zeros(z_y_add14);
	Zeros(z_xx_add14);
	Zeros(z_xy_add14);
	Zeros(z_yy_add14);
	z_x_add14(1, 1) = DER_INCREMENTS_Z(Z_size1 - 2, 2) / 2; z_x_add14(1, 3) = -DER_INCREMENTS_Z(Z_size1 - 2, 2) / 2;
	z_y_add14(0, 2) = DER_INCREMENTS_Z(Z_size1 - 2, 2) / 2; z_y_add14(2, 2) = -DER_INCREMENTS_Z(Z_size1 - 2, 2) / 2;
	z_xx_add14(1, 0) = DER_INCREMENTS_Z(Z_size1 - 2, 2); z_xx_add14(1, 1) = DER_INCREMENTS_Z(Z_size1 - 2, 2);
	z_xx_add14(1, 2) = -2 * DER_INCREMENTS_Z(Z_size1 - 2, 2); z_xx_add14(1, 3) = DER_INCREMENTS_Z(Z_size1 - 2, 2);
	z_yy_add14(0, 2) = DER_INCREMENTS_Z(Z_size1 - 2, 2); z_yy_add14(1, 2) = -2 * DER_INCREMENTS_Z(Z_size1 - 2, 2);
	z_yy_add14(2, 2) = DER_INCREMENTS_Z(Z_size1 - 2, 2); z_yy_add14(3, 2) = DER_INCREMENTS_Z(Z_size1 - 2, 2);
	z_xy_add14(0, 0) = DER_INCREMENTS_Z(Z_size1 - 2, 2) / 4; z_xy_add14(0, 1) = DER_INCREMENTS_Z(Z_size1 - 2, 2) / 4; z_xy_add14(0, 3) = -DER_INCREMENTS_Z(Z_size1 - 2, 2) / 4;
	z_xy_add14(2, 0) = -DER_INCREMENTS_Z(Z_size1 - 2, 2) / 4; z_xy_add14(2, 1) = -DER_INCREMENTS_Z(Z_size1 - 2, 2) / 4; z_xy_add14(2, 3) = DER_INCREMENTS_Z(Z_size1 - 2, 2) / 4;
	z_xy_add14(3, 0) = -DER_INCREMENTS_Z(Z_size1 - 2, 2) / 4; z_xy_add14(3, 1) = -DER_INCREMENTS_Z(Z_size1 - 2, 2) / 4; z_xy_add14(3, 3) = DER_INCREMENTS_Z(Z_size1 - 2, 2) / 4;
	z_x_add24 = -z_x_add14; z_y_add24 = -z_y_add14; z_xx_add24 = -z_xx_add14; z_yy_add24 = -z_yy_add14; z_xy_add24 = -z_xy_add14;
	for (int i = Z_size1 - 3; i <= Z_size1; i++) {
		for (int j = 0; j <= 3; j++) {
			z_xo4(ix, j) = Z_x(i, j);
			z_yo4(ix, j) = Z_y(i, j);
			z_xxo4(ix, j) = Z_xx(i, j);
			z_yyo4(ix, j) = Z_yy(i, j);
			z_xyo4(ix, j) = Z_xy(i, j);
		}
		ix++;
	}
	ix = 0;
	z_x_new14 = z_xo4 + z_x_add14; z_y_new14 = z_yo4 + z_y_add14; z_xx_new14 = z_xxo4 + z_xx_add14; z_yy_new14 = z_yyo4 + z_yy_add14; z_xy_new14 = z_xyo4 + z_xy_add14;
	z_x_new24 = z_xo4 + z_x_add24; z_y_new24 = z_yo4 + z_y_add24; z_xx_new24 = z_xxo4 + z_xx_add24; z_yy_new24 = z_yyo4 + z_yy_add24; z_xy_new24 = z_xyo4 + z_xy_add24;
	for (int i = Z_size1 - 3; i <= Z_size1; i++) {
		for (int j = 0; j <= 3; j++) {
			Alpha_new4(ix, j) = A(i, j);
			Beta_new4(ix, j) = B(i, j);
			P_new4(ix, j) = S(i, j);
		}
		ix++;
	}
	ix = 0;
	I_u_1 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14);
	I_u_2 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24);
	deriva_I_Z((Z_size1 - 2) * (Z_size2+1) + 2) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(Z_size1 - 2, 2));
	//=============================================================================================================================//
	Zeros(z_x_add14);
	Zeros(z_y_add14);
	Zeros(z_xx_add14);
	Zeros(z_xy_add14);
	Zeros(z_yy_add14);
	z_x_add14(1, 0) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 2) / 2; z_x_add14(1, 2) = -DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 2) / 2;
	z_y_add14(0, 1) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 2) / 2; z_y_add14(2, 1) = -DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 2) / 2;
	z_xx_add14(1, 0) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 2); z_xx_add14(1, 1) = -2 * DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 2); 
	z_xx_add14(1, 2) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 2); z_xx_add14(1, 3) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 2);
	z_yy_add14(0, 1) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 2); z_yy_add14(1, 1) = -2 * DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 2);
	z_yy_add14(2, 1) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 2); z_yy_add14(3, 1) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 2);
	z_xy_add14(0, 0) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 2) / 4; z_xy_add14(0, 2) = -DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 2) / 4; z_xy_add14(0, 3) = -DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 2) / 4;
	z_xy_add14(2, 0) = -DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 2) / 4; z_xy_add14(2, 2) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 2) / 4; z_xy_add14(2, 3) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 2) / 4;
	z_xy_add14(3, 0) = -DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 2) / 4; z_xy_add14(3, 2) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 2) / 4; z_xy_add14(3, 3) = DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 2) / 4;
	z_x_add24 = -z_x_add14; z_y_add24 = -z_y_add14; z_xx_add24 = -z_xx_add14; z_yy_add24 = -z_yy_add14; z_xy_add24 = -z_xy_add14;
	for (int i = Z_size1 - 3; i <= Z_size1; i++) {
		for (int j = Z_size2 - 3; j <= Z_size2; j++) {
			z_xo4(ix, jy) = Z_x(i, j);
			z_yo4(ix, jy) = Z_y(i, j);
			z_xxo4(ix, jy) = Z_xx(i, j);
			z_yyo4(ix, jy) = Z_yy(i, j);
			z_xyo4(ix, jy) = Z_xy(i, j);
			jy++;
		}
		ix++; jy = 0;
	}
	ix = 0; jy = 0;
	z_x_new14 = z_xo4 + z_x_add14; z_y_new14 = z_yo4 + z_y_add14; z_xx_new14 = z_xxo4 + z_xx_add14; z_yy_new14 = z_yyo4 + z_yy_add14; z_xy_new14 = z_xyo4 + z_xy_add14;
	z_x_new24 = z_xo4 + z_x_add24; z_y_new24 = z_yo4 + z_y_add24; z_xx_new24 = z_xxo4 + z_xx_add24; z_yy_new24 = z_yyo4 + z_yy_add24; z_xy_new24 = z_xyo4 + z_xy_add24;
	for (int i = Z_size1 - 3; i <= Z_size1; i++) {
		for (int j = Z_size2 - 3; j <= Z_size2; j++) {
			Alpha_new4(ix, jy) = A(i, j);
			Beta_new4(ix, jy) = B(i, j);
			P_new4(ix, jy) = S(i, j);
			jy++;
		}
		ix++; jy = 0;
	}
	ix = 0; jy = 0;
	I_u_1 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new14, z_y_new14, z_xx_new14, z_yy_new14, z_xy_new14);
	I_u_2 = evaluation_func0(Alpha_new4, Beta_new4, P_new4, z_x_new24, z_y_new24, z_xx_new24, z_yy_new24, z_xy_new24);
	deriva_I_Z((Z_size1 - 2) * (Z_size2+1) + Z_size2 - 2) = (I_u_2 - I_u_1) / (-2 * DER_INCREMENTS_Z(Z_size1 - 2, Z_size2 - 2));
}
void Adam::calculate_A(double con1, double con2, vector<int>(&D_pos)[2], vector<int>(&N_pos)[2], int min_pos1, int max_pos2) {
	
	omp_set_num_threads(P_ara);
	#pragma omp parallel sections
	{
		#pragma omp section
		{
			calculate_A22();
		}

		#pragma omp section
		{
			calculate_A23();
		}

		#pragma omp section
		{
			calculate_A24();
		}

		#pragma omp section
		{
			calculate_A32();
		}

		#pragma omp section
		{
			calculate_A42();
		}

		#pragma omp section
		{
			calculate_A33(con1, con2, D_pos, N_pos, min_pos1, max_pos2);
		}

		#pragma omp section
		{
			calculate_A34();
		}

		#pragma omp section
		{
			calculate_A43();
		}

		#pragma omp section
		{
			calculate_A44();
		}
	}
	/*calculate_A22();
	calculate_A23();
	calculate_A24();
	calculate_A32();
	calculate_A42();
	calculate_A33(con1, con2, D_pos, N_pos, min_pos1, max_pos2);
	calculate_A34();
	calculate_A43();
	calculate_A44();*/
	//if(ttt<=2) ccc(deriva_I_Z);
	ttt++;
	double temp_sum = 0.0;
	for (int i = 0; i < 81; i++) {
		for (int j = 0; j < 81; j++) {
			Res(i, j) = deriva_I_Z(i * 81 + j);
			temp_sum += Res(i, j);
		}
	}

	cout << "*********************************************" << endl;
	cout << "差分矩阵之和：" << temp_sum << endl;
}
void Adam::ccc(Matrix<double,1,81*81>& de) {
	ifstream file("C:\\Users\\yunxiang.xing\\Desktop\\LBM\\fvm_solver\\test\\cc" + to_string(ttt) + ".txt");
	double temp; 
	for (int i = 0; i < 6561; i++) {
		file >> temp;
		if (abs(temp - de(i)) > 1e-7) {
			cout << " ";
		}
	}

}
vector<double> Adam::Adam_solver(int D_pos_num, int N_pos_num, vector<int>(&D_pos)[2], vector<int>(&N_pos)[2], Matrix<double, 81, 81>&Ssmooth, double I_u_old) {
	int maxt = 1;
	double I_u_new = 0.0;
	eta = 0.1;
	vector<double>min_max;
	while (1) {
		*Z_new = Z.array() - eta * (*gt).array();
		partial_derivative(*Z_new, Z_size1, Z_size2, 1.0, 1.0);

		*H_mid = (((1 + Z_y.array().pow(2)) * Z_xx.array() - 2 * Z_x.array() * Z_y.array() * Z_xy.array() + (1 + Z_x.array().pow(2)) * Z_yy.array()) / (2 * (1 + Z_x.array().pow(2) + Z_y.array().pow(2)).pow(1.5))) * 1e3;
		*K_mid = ((Z_xx.array() * Z_yy.array() - Z_xy.array().pow(2)) / (1 + Z_x.array().pow(2) + Z_y.array().pow(2)).pow(2)) * 1e6;
		*D_mid = (1 - n) * (*H_mid).array();
		*C_mid = (1 - n) * 2 * ((*H_mid).array().pow(2) - (*K_mid).array()).pow(0.5);

		min_max = Adam::min(D_pos_num, N_pos_num, D_pos, N_pos, *C_mid, *D_mid, Ssmooth);
		if (abs(min_max[0]) - 0.12 < 0)
			min_max[0] = 0;
		else
			min_max[0] = abs(min_max[0]) - 0.12;
		if (min_max[1] - 0.06 < 0)
			min_max[1] = 0;
		else
			min_max[1] = min_max[1] - 0.06;

		I_u_new = evaluation_func(A, B, S, Z_x, Z_y, Z_xx, Z_yy, Z_xy, min_max[0], min_max[1]);

		if (I_u_old > I_u_new || maxt > 20) {
			break;
		}
		eta = eta * 0.5;
		maxt = maxt + 1;
	}
	Z = *Z_new;
	
	cout << "约束1：" << min_max[0] << endl;
	cout << "约束2：" << min_max[1] << endl;
	min_max.push_back(I_u_new);
	return min_max;
}
void Adam::All() {
	ReadDate();
	//x,y
	for (int i = 0; i < 81; i++)
		for (int j = 0; j < 81; j++) {
			X(i, j) = -40 + j;
			Y(i, j) = -40 + i;
			DER_INCREMENTS_Z(i, j) = 0;
		}


	

	vector<int>D_pos[2];
	vector<int>N_pos[2];

	clock_t start = clock();

	for (int i = D_center - r + 41; i <= D_center + r + 41; i++)
		for (int j = -r + 41; j <= r + 41; j++) {
			if (X(i-1, j-1) * X(i-1, j-1) + pow(Y(i-1, j-1) - D_center, 2.0) <= 1.0*r * r) {
				D_pos[0].push_back(i-1);
				D_pos[1].push_back(j-1);
			}
		}
	for(int i = -N_center - r + 41; i <= -N_center + r + 41; i++)
		for (int j = floor(inset) - r + 41; j <= ceil(inset) + r + 41; j++) {
			if (pow(X(i-1, j-1) - inset, 2) + pow(Y(i-1, j-1) + N_center, 2.0) <= 1.0*r * r) {
				N_pos[0].push_back(i-1);
				N_pos[1].push_back(j-1);
			}
		}
	S = S * 1e-3;

	int opttimes = 1;
	int space_x = 1;
	int space_y = 1;

	partial_derivative(Z, Z_size1, Z_size2, space_x, space_y);

	Matrix<double, 81, 81>H_first, K_first, D_first, C_first;
	H_first = (((1.0 + Z_y.array() * Z_y.array()) * Z_xx.array() - 2.0 * Z_x.array() * Z_y.array() * Z_xy.array() + (1 + Z_x.array() * Z_x.array()) * Z_yy.array()) / (2.0 * (1.0 + Z_x.array() * Z_x.array() + Z_y.array() * Z_y.array()).array().pow(1.5))) * 1e3;
	K_first = ((Z_xx.array() * Z_yy.array() - Z_xy.array() * Z_xy.array()) / (1 + Z_x.array() * Z_x.array() + Z_y.array() * Z_y.array()).array().pow(2)) * 1e6;
	D_first = (1 - n) * H_first.array();
	C_first = 2 * (1 - n) * (H_first.array() * H_first.array() - K_first.array()).array().pow(0.5);

	vector<double>Adam_m1 = Adam::min(D_pos[0].size(), N_pos[0].size(), D_pos, N_pos, C_first, D_first, S);
	if (abs(Adam_m1[0]) - 0.12 < 0)
		Adam_m1[0] = 0;
	else
		Adam_m1[0] = abs(Adam_m1[0]) - 0.12;
	if (Adam_m1[1] - 0.06 < 0)
		Adam_m1[1] = 0;
	else
		Adam_m1[1] = Adam_m1[1] - 0.06;

	I_u.push_back(evaluation_func(A, B, S, Z_x, Z_y, Z_xx, Z_yy, Z_xy, Adam_m1[0], Adam_m1[1]));

	cout << "迭代次数：" << "0" << endl;
	cout << "损失函数：" << I_u[0] * 1000 << endl;

	for (int i = 0; i < 81; i++) {
		for (int j = 0; j < 81; j++) {
			DER_INCREMENTS_Z(i, j) = 1e-6;
		}
	}
	
	
	newmatrix();
	
	int num = 1;
	while (Adam_m1[0] > 0 || Adam_m1[1] > 0) {

		calculate_A(Adam_m1[0], Adam_m1[1], D_pos, N_pos, Adam_m1[2], Adam_m1[3]);
		num++;

		*v = beta1 * (*v).array() + (1 - beta1) * Res.array();
		*s = beta2 * (*s).array() + (1 - beta2) * Res.array() * Res.array();
		*vt = (*v).array() / (1 - pow(beta1, num - 1));
		*st = (*s).array() / (1 - pow(beta2, num - 1));
		*gt = (*vt).array() / ((*st).array().pow(0.5) + 1e-6);

		Adam_m1 = Adam_solver(D_pos[0].size(), N_pos[0].size(), D_pos, N_pos, S, I_u[I_u.size() - 1]);
		I_u.push_back(Adam_m1[Adam_m1.size() - 1]);

		cout << "迭代次数：" << num - 1 << endl;
		
		cout << "损失函数：" << I_u[I_u.size() - 1] * 1e3 << endl;
		cout << "*********************************************" << endl;
	}

	clock_t end = clock();
	double du = double(end - start) / CLOCKS_PER_SEC;
	cout << "Time consumption: " << du << "s" << endl;

	coutresult();
	
}
void Adam::coutresult() {
	partial_derivative(Z, Z_size1, Z_size2, 1.0, 1.0);
	*H_mid = (((1.0 + Z_y.array() * Z_y.array()) * Z_xx.array() - 2.0 * Z_x.array() * Z_y.array() * Z_xy.array() + (1 + Z_x.array() * Z_x.array()) * Z_yy.array()) / (2.0 * (1.0 + Z_x.array() * Z_x.array() + Z_y.array() * Z_y.array()).array().pow(1.5))) * 1e3;
	*K_mid = ((Z_xx.array() * Z_yy.array() - Z_xy.array() * Z_xy.array()) / (1 + Z_x.array() * Z_x.array() + Z_y.array() * Z_y.array()).array().pow(2)) * 1e6;
	*D_mid = (1 - n) * (*H_mid).array();
	*C_mid = 2 * (1 - n) * ((*H_mid).array() * (*H_mid).array() - (*K_mid).array()).array().pow(0.5);

	coutmatrix(Z_final, *C_mid, *D_mid);
}
template <typename T, int rows, int cols>
void Adam::coutmatrix(string& out, Matrix<T, rows, cols>& C_final, Matrix<T, rows, cols>& D_final) {
	std::size_t pos = out.find_last_of("/\\");
	string directory;
	// 取路径的部分（不包含文件名）
	if (pos != std::string::npos) {
		directory = out.substr(0, pos);
	}
	else {
		directory = "";
	}

	string temp = directory + "Z_final.csv";
	string temp1 = directory + "C_final.csv";
	string temp2 = directory + "D_final.csv";
	string temp3 = directory + "I_u_all.csv";
	ofstream file(temp);
	ofstream file1(temp1);
	ofstream file2(temp2);
	ofstream file3(temp3);

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			file << Z(i, j) << ",";
		}
		file << endl;
	}

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			file1 << C_final(i, j) << ",";
		}
		file1 << endl;
	}

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			file2 << D_final(i, j) << ",";
		}
		file2 << endl;
	}

	for (int i = 0; i < I_u.size(); i++) {
		file3 << i << "," << I_u[i] << endl;
	}

	//to_sdf

}
void Adam::newmatrix() {
	v = new Matrix<double, 81, 81>();
	s = new Matrix<double, 81, 81>();
	vt = new Matrix<double, 81, 81>();
	st = new Matrix<double, 81, 81>();
	gt = new Matrix<double, 81, 81>();
	Z_new = new Matrix<double, 81, 81>();
	H_mid = new Matrix<double, 81, 81>();
	C_mid = new Matrix<double, 81, 81>();
	D_mid = new Matrix<double, 81, 81>();
	K_mid = new Matrix<double, 81, 81>();
	Zeros(*v);
	Zeros(*s);
	Zeros(*gt);
}