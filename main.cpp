//#include "cvfem/Mesh_cvfem.h"
//#include "fvm/FVM.h"
#include"optimizer.h"

int main(int argc, char* argv[]) {

	//Mesh_cvfem test;
	/*vector<double>b;
	b.resize(3);
	b[0] = 20;
	b[1] = 33;
	b[2] = 36;
	vector<vector<double>>A;
	A.resize(3);
	for (int i = 0; i < 3; i++) {
		A[i].resize(3);
	}
	A[0][0] = 8;
	A[0][1] = -3;
	A[0][2] = 2;
	A[1][0] = 4;
	A[1][1] = 11;
	A[1][2] = -1;
	A[2][0] = 6;
	A[2][1] = 3;
	A[2][2] = 12;
	auto x = test.solver_equtionGaussSeidel(A, b);
	auto y = test.solver_equtionJacobi(A, b);
	auto z = test.solver_equtionSOR(A, b, 0.7);*/


	//vector<string>file1;
	//file1.resize(3);
	//file1[0] = "C:\\Users\\yunxiang.xing\\Desktop\\test\\dingzi_model1.msh";
	//file1[1] = "C:\\Users\\yunxiang.xing\\Desktop\\test\\zhuanzi_model1.msh";
	//file1[2] = "C:\\Users\\yunxiang.xing\\Desktop\\test\\zhijia_model.msh";
	//FVM::Fvm test;
	//test.marge_msh(file1);
	//test.cal_Diff("C:\\Users\\freedom\\Desktop\\×Ô±àÐ´Çó½âÆ÷\\fvm_solver\\reference\\test.rmsh");
	//test.cal("C:\\Users\\yunxiang.xing\\Desktop\\LBM\\fvm_solver\\reference\\test.rmsh");

	if (argc < 5 || argc > 13) {
		cout << "--zfin=path to Z_final matrix" << endl;
		cout << "--alp=path to alpha_smooth matrix" << endl;
		cout << "--bet=pathh to beta_smooth matrix" << endl;
		cout << "--ssm=path to S_smooth matrix" << endl;
		cout << "--para=thread number" << endl;
		cout << "--beta1=number(0.0, 1.0)" << endl;
		cout << "--beta2=number(0.0,1.0)" << endl;
		cout << "--Dcenter=number" << endl;
		cout << "--Ncenter=number" << endl;
		cout << "--Inset=number" << endl;
		cout << "--r=number" << endl;
		cout << "--n=number" << endl;
		return 0;
	}
	
	string temp5 = "";
	string temp6 = "";
	string temp7 = "";
	string temp8 = "";
	string temp9 = "";
	string temp10 = "";
	string temp11 = "";
	string temp12 = "";

	if (argc >= 6) {
		temp5 = argv[5];
	}
	if (argc >= 7) {
		temp6 = argv[6];
	}
	if (argc >= 8) {
		temp7 = argv[7];
	}
	if (argc >= 9) {
		temp8 = argv[8];
	}
	if (argc >= 10) {
		temp9 = argv[9];
	}
	if (argc >= 11) {
		temp10 = argv[10];
	}
	if (argc >= 12) {
		temp11 = argv[11];
	}
	if (argc >= 13) {
		temp12 = argv[12];
	}
	Adam test(argv[1], argv[2], argv[3], argv[4], temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12);
	test.All();
	
	return 0;
}