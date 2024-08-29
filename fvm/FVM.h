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
	class Fvm :public Mesh_cvfem {
	public:
		void findsharedface();
		void init(string cwd);
		void cal_Diff(string cwd);
		void cal_gradient1(vector<double>&);
		void cal_gradient(vector<double>&);
		//ѹ���ٶ�����㷨
		void simple();
		void read(string cwd);
		double error_Fvm(vector<double>&, vector<double>&);
		void Write(vector<double>& x, string cwd);
		inline uint64_t faceid(uint64_t, uint64_t, uint64_t);
		//virtual void write() override;
		//
	private:
		void pgrad(Point&, int, int, Point&, Point&, double);
		//vector<Element_cvfem>mesh_eles;//��Ԫ
		//vector<Point>mesh_nodes;//�ڵ�
		//vector<vector<int>>mesh_nei;//�ھ�
		//vector<vector<int>>mesh_neihaxi;//���ټ����ھӽڵ��ϣ��
		//vector<vector<int>>mesh_bj;//�߽�����������
		//vector<Point>mesh_gradient;//��Ԫ���ݶ�
		map<uint64_t, vector<int>>sharedface;
		vector<Point>Grad;
		vector<double>diff_k;

		vector<vector<gradcon>>Gradelements;
		vector<gradvalue>grad;

		map<uint64_t, int>fvm_bj;
	};
}
#endif