#ifndef MESH_CVFEM_H
#define MESH_CVFEM_H
#include"../base/base.h"

using namespace std;

class Point {
public:
	Point(double x = 0.0, double y = 0.0, double z = 0.0) :x(x), y(y), z(z) {};
	const double& getx();
	const double& gety();
	const double& getz();
	void setx(double&);
	void sety(double&);
	void setz(double&);
	void setpoint(Point*);
	Point cross(Point&);
private:
	double x, y, z;
};

Point cross(Point a, Point b);
double norm(Point a);
Point operator+(Point a, Point b);
Point operator-(Point a, Point b);
double operator*(Point a, Point b);
Point operator*(double a, Point b);

//四面体单元
class  Element_cvfem {
public:
	//Element_cvfem() {};
	Point& getA();
	Point& getB();
	Point& getC();
	Point& getD();
	void setA(Point*);
	void setB(Point*);
	void setC(Point*);
	void setD(Point*);
	void setelement(Element_cvfem*);
	int& getdict();
	void setdict(int&);
	double& getvol();
	double& getS_abc();
	double& getS_abd();
	double& getS_acd();
	double& getS_bcd();
	Point getface_ABC();
	Point getface_ABD();
	Point getface_ACD();
	Point getface_BCD();
	Point getcent();
	vector<Point> getinterfacenormal();
	vector<Point> getelementsurfacenormal();
	int numA, numB, numC, numD;
private:
	void calvol();
	void calS(int);
	Point calsurface_normal(Point& a, Point& b, Point& c);
	bool cal_normal = 0;
	Point a, b, c, d;
	Point abc_normal, abd_normal, acd_normal, bcd_normal;
	int dict = 0;
	double vol = 0, S_abc = 0, S_abd = 0, S_acd = 0, S_bcd = 0;
};
class Mesh_cvfem {
public:
	bool readRmsh(string);
	const double max(const double&, const double&);
	void findneighbor();
	void init();
	void cal(string cwd);
	double error(vector<vector<double>>&, vector<double>&, vector<double>&, vector<double>&);
	vector<double> solver_equtionJacobi(vector<vector<double>>& A, vector<double>& b);//雅克比迭代
	vector<double> solver_equtionGaussSeidel(vector<vector<double>>& A, vector<double>& b);//高斯赛德尔迭代
	vector<double> solver_equtionSOR(vector<vector<double>>& A, vector<double>& b, double w);//SOR迭代
	void setdiffk(int,double[]);
	void yanzheng(vector<vector<double>>&, vector<double>&);
	//void WritingMeshRmsh();
	//void WritingDat();
	void WritingRes(vector<double>&);
	//virtual void write() = 0;
//private:
	vector<Element_cvfem>mesh_eles;//记录单元矩阵
	vector<Point>mesh_nodes;//记录节点矩阵
	vector<vector<int>>mesh_nei;//记录邻居节点矩阵
	vector<vector<int>>mesh_neihaxi;//快速检索邻居节点哈希表
	vector<vector<int>>mesh_bj;//边界三角形网格
	vector<double>vol;//记录控制体积
	vector<Point>vel;//记录速度
	vector<double>ap;//记录方程组apx-anbpx=bp
	vector<vector<double>>anb;
	vector<double>bp;
	vector<vector<double>>A;//系数矩阵
	double NA, NB;
	int maxnnb = 0;
	//int Physical;
	//vector<string>PhysicalNames;
};

#endif