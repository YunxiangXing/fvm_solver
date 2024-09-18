#ifndef CVFEM_H
#define CVFEM_H
#include"../cvfem/Mesh_cvfem.h"
#include<map>
#include <algorithm>
//1.over-relaxed: Ef = Sf * Sf / (ecf * Sf) * ecf
//
namespace FVM {
	enum bir {
		inlet = 1,
		outlet = 3
	};
	struct gradvalue {
		Point C, f1, f2, f3, f4;
		int idf1 = -1, idf2 = -1, idf3 = -1, idf4 = -1;
	};
	struct gradcon {
		Point data;
		int id;
	};
	struct Msh_Physical_Group {
		vector<int>eletype;
		vector<int>phyid;
		vector<string>phyname;
		vector<int>newphyid;
		int num_node = 0;
		int num_ele = 0;
	};
	class Fvm :public Mesh_cvfem {
	public:
		void findsharedface();
		void init(string cwd);
		void cal(string cwd);
		void cal_Diff(string cwd);
		void cal_Convertion(string cwd);
		void cal_gradient1(vector<double>&, vector<double>&);
		void cal_gradient(vector<double>&, vector<double>&);
		void insert_fai(double, uint64_t);
		//ѹ���ٶ�����㷨
		void simple();
		void read(string cwd);
		void marge_msh(vector<string>filepath);
		double error_Fvm(vector<double>&, vector<double>&);
		void Write(vector<double>& x, string cwd);
		inline uint64_t faceid(uint64_t, uint64_t, uint64_t);
		//virtual void write() override;
		//
	private:
		vector<double> solver_equtionGaussSeidel(vector<vector<double>>& A, vector<double>& b) override;//��˹���¶�����
		void pgrad(Point&, int, int, Point&, Point&, double);
		int birin = 1;
		int birout = 2;
		//vector<Element_cvfem>mesh_eles;//��Ԫ
		//vector<Point>mesh_nodes;//�ڵ�
		//vector<vector<int>>mesh_nei;//�ھ�
		//vector<vector<int>>mesh_neihaxi;//���ټ����ھӽڵ��ϣ��
		//vector<vector<int>>mesh_bj;//�߽�����������
		//vector<Point>mesh_gradient;//��Ԫ���ݶ�
		map<uint64_t, vector<int>>sharedface;
		map<uint64_t, double>sharedfai;
		vector<Point>Grad;
		vector<double>diff_k;

		vector<vector<gradcon>>Gradelements;
		vector<gradvalue>grad;

		map<uint64_t, int>fvm_bj;
		vector<double>res;
	};
}
#endif