#include "cvfem/Mesh_cvfem.h"
#include "fvm/FVM.h"
//int main(int argc, char *argv[])
//{
//    QApplication a(argc, argv);
//    QtWidgetsApplication1 w;
//    w.show();
//    return a.exec();
//}

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
	FVM::Fvm test;
	//test.marge_msh(file1);
	test.cal_Diff("C:\\Users\\freedom\\Desktop\\�Ա�д�����\\fvm_solver\\reference\\test.rmsh");
	return 0;
}