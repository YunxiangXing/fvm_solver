#include"Mesh_cvfem.h"

Point cross(Point a, Point b) {
	double x, y, z;
	x = a.gety() * b.getz() - a.getz() * b.gety();
	y = a.getz() * b.getx() - a.getx() * b.getz();
	z = a.getx() * b.gety() - a.gety() * b.getx();
	Point c(x, y, z);
	return c;
}
double norm(Point a) {
	double b;
	b = sqrt(a.getx() * a.getx() + a.gety() * a.gety() + a.getz() * a.getz());
	return b;
}
Point operator+(Point a, Point b) {
	Point c;
	double x, y, z;
	x = a.getx() + b.getx();
	y = a.gety() + b.gety();
	z = a.getz() + b.getz();
	c.setx(x);
	c.sety(y);
	c.setz(z);
	return c;
}
Point operator-(Point a, Point b) {
	Point c;
	double x, y, z;
	x = a.getx() - b.getx();
	y = a.gety() - b.gety();
	z = a.getz() - b.getz();
	c.setx(x);
	c.sety(y);
	c.setz(z);
	return c;
}
double operator*(Point a, Point b) {
	double c;
	c = a.getx() * b.getx() + a.gety() * b.gety() + a.getz() * b.getz();
	return c;
}
Point operator*(double a, Point b) {
	Point c;
	double x, y, z;
	x = a * b.getx();
	y = a * b.gety();
	z = a * b.getz();
	c.setx(x);
	c.sety(y);
	c.setz(z);
	return c;
}
const double& Point::getx() {
	return this->x;
}
const double& Point::gety() {
	return this->y;
}
const double& Point::getz() {
	return this->z;
}
void Point::setx(double& x) {
	this->x = x;
}
void Point::sety(double& y) {
	this->y = y;
}
void Point::setz(double& z) {
	this->z = z;
}
void Point::setpoint(Point* poi) {
	this->x = poi->getx();
	this->y = poi->gety();
	this->z = poi->getz();
}
Point Point::cross(Point& b) {
	double x, y, z;
	x = this->x * b.getz() - this->z * b.gety();
	y = this->z * b.getx() - this->x * b.getz();
	z = this->x * b.gety() - this->y * b.getx();
	Point c(x, y, z);
	return c;
}
Point& Element_cvfem::getA() {
	return this->a;
}
Point& Element_cvfem::getB() {
	return this->b;
}
Point& Element_cvfem::getC() {
	return this->c;
}
Point& Element_cvfem::getD() {
	return this->d;
}
void Element_cvfem::setA(Point* A) {
	this->a.setpoint(A);
}
void Element_cvfem::setB(Point* B) {
	this->b.setpoint(B);
}
void Element_cvfem::setC(Point* C) {
	this->c.setpoint(C);
}
void Element_cvfem::setD(Point* D) {
	this->d.setpoint(D);
}
int& Element_cvfem::getdict() {
	return this->dict;
}
void Element_cvfem::setdict(int& dict) {
	this->dict = dict;
}
void Element_cvfem::setelement(Element_cvfem* ele) {
	this->a = ele->a;
	this->b = ele->b;
	this->c = ele->c;
	this->d = ele->d;
}
void Element_cvfem::calvol() {
	this->vol = abs(1.0 / 6.0 * (this->d - this->a) * cross(this->b - this->a, this->c - this->a));
}
double& Element_cvfem::getvol() {
	if (this->vol < Min) {
		calvol();
	}
	return this->vol;
}
void Element_cvfem::calS(int t) {
	if (t == 1) S_abc = 0.5 * norm(cross(b - a, c - a));
	else if (t == 2) S_abd = 0.5 * norm(cross(b - a, d - a));
	else if (t == 3) S_acd = 0.5 * norm(cross(c - a, d - a));
	else if (t == 4) S_bcd = 0.5 * norm(cross(c - b, d - b));
}
double& Element_cvfem::getS_abc() {
	if (this->S_abc < Min) {
		calS(1);
	}
	return this->S_abc;
}
double& Element_cvfem::getS_abd() {
	if (this->S_abd < Min) {
		calS(2);
	}
	return this->S_abd;
}
double& Element_cvfem::getS_acd() {
	if (this->S_acd < Min) {
		calS(3);
	}
	return this->S_acd;
}
double& Element_cvfem::getS_bcd() {
	if (this->S_bcd < Min) {
		calS(4);
	}
	return this->S_bcd;
}
Point Element_cvfem::getface_ABC() {
	return 1.0 / 3.0 * (getA() + getB() + getC());
}
Point Element_cvfem::getface_ABD() {
	return 1.0 / 3.0 * (getA() + getB() + getD());
}
Point Element_cvfem::getface_ACD() {
	return 1.0 / 3.0 * (getA() + getC() + getD());
}
Point Element_cvfem::getface_BCD() {
	return 1.0 / 3.0 * (getB() + getC() + getD());
}
Point Element_cvfem::getcent() {
	Point temp;
	temp = 0.25 * (a + b + c + d);
	return temp;
}
vector<Point> Element_cvfem::getinterfacenormal() {
	vector<Point>temp;
	Point ac = cross(a, c);
	Point ab = cross(a, b);
	Point ad = cross(a, d);
	Point bc = cross(b, c);
	Point bd = cross(b, d);
	Point cd = cross(c, d);
	temp.resize(6);
	temp[0] = 1.0 / 24.0 * ac - 1.0 / 24.0 * ad + 1.0 / 24.0 * bc - 1.0 / 24.0 * bd + 1.0 / 12.0 * cd;
	temp[1] = -1.0 / 24.0 * ab + 1.0 / 24.0 * ad + 1.0 / 24.0 * bc - 1.0 / 12.0 * bd + 1.0 / 24.0 * cd;
	temp[2] = 1.0 / 24.0 * ab - 1.0 / 24.0 * ac + 1.0 / 12.0 * bc - 1.0 / 24.0 * bd + 1.0 / 24.0 * cd;
	temp[3] = -1.0 / 24.0 * ab - 1.0 / 24.0 * ac + 1.0 / 12.0 * ad - 1.0 / 24.0 * bd - 1.0 / 24.0 * cd;
	temp[4] = 1.0 / 24.0 * ab - 1.0 / 12.0 * ac + 1.0 / 24.0 * ad + 1.0 / 24.0 * bc - 1.0 / 24.0 * cd;
	temp[5] = 1.0 / 12.0 * ab - 1.0 / 24.0 * ac - 1.0 / 24.0 * ad + 1.0 / 24.0 * bc + 1.0 / 24.0 * bd;
	return temp;
}
vector<Point> Element_cvfem::getelementsurfacenormal() {
	if (cal_normal == false) {
		this->abc_normal = calsurface_normal(a, b, c);
		this->abd_normal = calsurface_normal(a, b, d);
		this->acd_normal = calsurface_normal(a, c, d);
		this->bcd_normal = calsurface_normal(b, c, d);
		cal_normal = true;
	}
	vector<Point>temp;
	temp.push_back(this->abc_normal);
	temp.push_back(this->abd_normal);
	temp.push_back(this->acd_normal);
	temp.push_back(this->bcd_normal);
	return temp;
}
Point Element_cvfem::calsurface_normal(Point& a, Point& b, Point& c) {
	Point cen = getcent();
	Point surf_abc = 1.0 / 3.0 * (a + b + c);
	Point diction = surf_abc - cen;
	Point ab = b - a;
	Point ac = c - a;
	Point normal = 0.5 * cross(ab, ac);
	if (diction * normal < 0) normal = -1 * normal;
	return normal;
}
bool Mesh_cvfem::readRmsh(string cwd) {
	ifstream file1;
	file1.open(cwd);
	if (!file1.is_open()) {
		return 0;
	}
	string line;
	vector<int> Phy_name;
	while (!file1.eof()) {
		getline(file1, line);
		if (line == "$MeshFormat") {
			getline(file1, line);
			if (line != "4 0 8") {
				return 0;
			}
		}
		if (line == "$Entities") {
			int geo[4];//point,line,face,comp
			file1 >> geo[0] >> geo[1] >> geo[2] >> geo[3];
			double coor[6];
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < geo[i]; j++) {
					int num;
					file1 >> num;//编号
					for (int k = 0; k < 6; k++) {
						file1 >> coor[k];//坐标
					}
					if (i >= 2) {
						int Phy_num;//物理标签个数
						int sign;
						file1 >> Phy_num;
						if (Phy_num != 0) {
							int Phy_name1;
							for (int l = 0; l < Phy_num; l++) {
								file1 >> sign;
								Phy_name.push_back(sign);
							}
						}
						else {
							Phy_name.push_back(-1);
						}
					}
					getline(file1, line);
				}
			}
		}
		if (line == "$Nodes") {
			int num_nodes, num_geo;
			file1 >> num_geo >> num_nodes;
			mesh_nodes.resize(num_nodes);
			for (int i = 0; i < num_geo; i++) {
				int num1, num2, num3, num4;
				file1 >> num1 >> num2 >> num3 >> num4;
				for (int j = 0; j < num4; j++) {
					int num;
					double coorx, coory, coorz;
					file1 >> num >> coorx >> coory >> coorz;
					mesh_nodes[num - 1].setx(coorx);
					mesh_nodes[num - 1].sety(coory);
					mesh_nodes[num - 1].setz(coorz);
				}
			}
		}
		if (line == "$Elements") {
			int num_eles, num_geo, temp = 0;
			file1 >> num_geo >> num_eles;
			for (int i = 0; i < num_geo; i++) {
				int num1, num2, num3, num4;//编号，维度，网格类型，元素个数
				file1 >> num1 >> num2 >> num3 >> num4;
				if (num2 == 3) mesh_eles.resize(num4);
				getline(file1, line);
				for (int j = 0; j < num4; j++) {
					if (num2 == 2 && Phy_name[num1] != -1) {
						int num, dicta, dictb, dictc;
						file1 >> num >> dicta >> dictb >> dictc;
						vector<int>tempv;
						tempv.push_back(dicta - 1);
						tempv.push_back(dictb - 1);
						tempv.push_back(dictc - 1);
						tempv.push_back(Phy_name[num1 - 1]);
						mesh_bj.push_back(tempv);
					}
					if (num2 == 3) {
						int num, dicta, dictb, dictc, dictd;
						file1 >> num >> dicta >> dictb >> dictc >> dictd;
						mesh_eles[temp].setA(&mesh_nodes[dicta - 1]);
						mesh_eles[temp].setB(&mesh_nodes[dictb - 1]);
						mesh_eles[temp].setC(&mesh_nodes[dictc - 1]);
						mesh_eles[temp].setD(&mesh_nodes[dictd - 1]);
						mesh_eles[temp].numA = dicta - 1;
						mesh_eles[temp].numB = dictb - 1;
						mesh_eles[temp].numC = dictc - 1;
						mesh_eles[temp].numD = dictd - 1;
						temp++;
					}
					else {
						getline(file1, line);
					}
				}
			}
		}
	}
	return 1;
}
const double Mesh_cvfem::max(const double& qf, const double& temp) {
	if (qf > temp) return qf;
	else return temp;
}
void Mesh_cvfem::findneighbor() {
	//寻找邻居点
	//遍历单元：
	int num_nodes = mesh_nodes.size();
	int num_ele = mesh_eles.size();
	mesh_neihaxi.resize(num_nodes);
	mesh_nei.resize(num_nodes);
	for (int i = 0; i < num_nodes; i++) {
		mesh_neihaxi[i].resize(num_nodes);
		for (int j = 0; j < num_nodes; j++)
			mesh_neihaxi[i][j] = -1;
	}
	for (int i = 0; i < num_ele; i++) {
		int temp[4] = { mesh_eles[i].numA,mesh_eles[i].numB,mesh_eles[i].numC,mesh_eles[i].numD };
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				if (j != k) {
					mesh_neihaxi[temp[j]][temp[k]] = 1;
				}
			}
		}
	}
	maxnnb = 0;
	int tempnnb = 0;
	for (int i = 0; i < num_nodes; i++) {
		for (int j = 0; j < num_nodes; j++) {
			if (mesh_neihaxi[i][j] == 1) {
				tempnnb++;
			}
		}
		if (tempnnb > maxnnb) maxnnb = tempnnb;
		tempnnb = 0;
	}
	tempnnb = 0;
	for (int i = 0; i < num_nodes; i++) {
		mesh_nei[i].resize(maxnnb);
		for (int j = 0; j < num_nodes; j++) {
			if (mesh_neihaxi[i][j] == 1) {
				mesh_nei[i][tempnnb] = j;
				tempnnb++;
			}
		}
		tempnnb = 0;
	}
	//对哈希表重新排序为联系anb的矩阵
	tempnnb = -1;
	for (int i = 0; i < num_nodes; i++) {
		for (int j = 0; j < num_nodes; j++) {
			if (mesh_neihaxi[i][j] == 1) {
				mesh_neihaxi[i][j] = mesh_neihaxi[i][j] + tempnnb;
				tempnnb++;
			}
		}
		tempnnb = -1;
	}
}
void Mesh_cvfem::init() {
	//初始化速度
	NA = 17.0 / 48.0;
	NB = 7.0 / 48.0;
	vel.resize(mesh_nodes.size());
	ap.resize(mesh_nodes.size());
	anb.resize(mesh_nodes.size());
	vol.resize(mesh_nodes.size());
	bp.resize(mesh_nodes.size());
	for (int i = 0; i < mesh_nodes.size(); i++) {
		anb[i].resize(maxnnb);
		ap[i] = 0.0;
		bp[i] = 0.0;
		//v=1/r
		double r = mesh_nodes[i].getx() * mesh_nodes[i].getx() + mesh_nodes[i].getz() * mesh_nodes[i].getz();
		double vx = mesh_nodes[i].getx() / r;
		double vz = mesh_nodes[i].getz() / r;
		vel[i].setx(vx);
		vel[i].setz(vz);
		for (int j = 0; j < maxnnb; j++) {
			anb[i][j] = 0.0;
		}
	}
	for (int i = 0; i < mesh_nodes.size(); i++) {
		vol[i] = 0.0;
	}
}
//获得系数矩阵
void Mesh_cvfem::setdiffk(int i,double v[]) {
	double temp_v[4];
	temp_v[0] = 1.0 / sqrt(mesh_eles[i].getA().getx() * mesh_eles[i].getA().getx() + mesh_eles[i].getA().getz() * mesh_eles[i].getA().getz());
	temp_v[1] = 1.0 / sqrt(mesh_eles[i].getB().getx() * mesh_eles[i].getB().getx() + mesh_eles[i].getB().getz() * mesh_eles[i].getB().getz());
	temp_v[2] = 1.0 / sqrt(mesh_eles[i].getC().getx() * mesh_eles[i].getC().getx() + mesh_eles[i].getC().getz() * mesh_eles[i].getC().getz());
	temp_v[3] = 1.0 / sqrt(mesh_eles[i].getD().getx() * mesh_eles[i].getD().getx() + mesh_eles[i].getD().getz() * mesh_eles[i].getD().getz());
	v[0] = NA * temp_v[0] + NA * temp_v[1] + NB * temp_v[2] + NB * temp_v[3];
	v[1] = NA * temp_v[0] + NB * temp_v[1] + NA * temp_v[2] + NB * temp_v[3];
	v[2] = NA * temp_v[0] + NB * temp_v[1] + NB * temp_v[2] + NA * temp_v[3];
	v[3] = NB * temp_v[0] + NA * temp_v[1] + NA * temp_v[2] + NB * temp_v[3];
	v[4] = NB * temp_v[0] + NA * temp_v[1] + NB * temp_v[2] + NA * temp_v[3];
	v[5] = NB * temp_v[0] + NB * temp_v[1] + NA * temp_v[2] + NA * temp_v[3];
}
void Mesh_cvfem::cal(string cwd) {
	//diffusive flux/advective flux
	readRmsh(cwd);
	findneighbor();
	init();
	int num_ele = mesh_eles.size();
	int num_nodes = mesh_nodes.size();
	double df[6][4];
	double v[6];
	vector<Point>f;
	for (int i = 0; i < num_ele; i++) {
		f = mesh_eles[i].getinterfacenormal();
		Point gradN[4];
		gradN[0] = 1.0 / 6.0 * (1.0 / mesh_eles[i].getvol()) * cross(mesh_eles[i].getD() - mesh_eles[i].getB(), mesh_eles[i].getC() - mesh_eles[i].getB());
		gradN[1] = 1.0 / 6.0 * (1.0 / mesh_eles[i].getvol()) * cross(mesh_eles[i].getC() - mesh_eles[i].getA(), mesh_eles[i].getD() - mesh_eles[i].getA());
		gradN[2] = 1.0 / 6.0 * (1.0 / mesh_eles[i].getvol()) * cross(mesh_eles[i].getD() - mesh_eles[i].getA(), mesh_eles[i].getB() - mesh_eles[i].getA());
		gradN[3] = 1.0 / 6.0 * (1.0 / mesh_eles[i].getvol()) * cross(mesh_eles[i].getB() - mesh_eles[i].getA(), mesh_eles[i].getC() - mesh_eles[i].getA());
		//v=1/r
		for (int j = 0; j < 6; j++)
			v[j] = 1.0;
		setdiffk(i, v);
		for (int j = 0; j < 6; j++) {
			for (int k = 0; k < 4; k++) {
				df[j][k] = v[j] * gradN[k] * f[j];
			}
		}
		//qfi=vfi*nfi
		Point vf[6];
		double qf[6];
		vf[0] = NA * vel[mesh_eles[i].numA] + NA * vel[mesh_eles[i].numB] + NB * vel[mesh_eles[i].numC] + NB * vel[mesh_eles[i].numD];
		vf[1] = NA * vel[mesh_eles[i].numA] + NB * vel[mesh_eles[i].numB] + NA * vel[mesh_eles[i].numC] + NB * vel[mesh_eles[i].numD];
		vf[2] = NA * vel[mesh_eles[i].numA] + NB * vel[mesh_eles[i].numB] + NB * vel[mesh_eles[i].numC] + NA * vel[mesh_eles[i].numD];
		vf[3] = NB * vel[mesh_eles[i].numA] + NA * vel[mesh_eles[i].numB] + NA * vel[mesh_eles[i].numC] + NB * vel[mesh_eles[i].numD];
		vf[4] = NB * vel[mesh_eles[i].numA] + NA * vel[mesh_eles[i].numB] + NB * vel[mesh_eles[i].numC] + NA * vel[mesh_eles[i].numD];
		vf[5] = NB * vel[mesh_eles[i].numA] + NB * vel[mesh_eles[i].numB] + NA * vel[mesh_eles[i].numC] + NA * vel[mesh_eles[i].numD];
		for (int j = 0; j < 6; j++)
			qf[j] = vf[j] * f[j];
		//Node A
		ap[mesh_eles[i].numA] = ap[mesh_eles[i].numA] - df[0][0] - df[1][0] - df[2][0] + max(qf[0], 0.0) + max(qf[1], 0.0) + max(qf[2], 0.0);
		anb[mesh_eles[i].numA][mesh_neihaxi[mesh_eles[i].numA][mesh_eles[i].numB]] = anb[mesh_eles[i].numA][mesh_neihaxi[mesh_eles[i].numA][mesh_eles[i].numB]] + df[0][1] + df[1][1] +
			df[2][1] + max(-qf[0], 0.0);
		anb[mesh_eles[i].numA][mesh_neihaxi[mesh_eles[i].numA][mesh_eles[i].numC]] = anb[mesh_eles[i].numA][mesh_neihaxi[mesh_eles[i].numA][mesh_eles[i].numC]] + df[0][2] + df[1][2] +
			df[2][2] + max(-qf[1], 0.0);
		anb[mesh_eles[i].numA][mesh_neihaxi[mesh_eles[i].numA][mesh_eles[i].numD]] = anb[mesh_eles[i].numA][mesh_neihaxi[mesh_eles[i].numA][mesh_eles[i].numD]] + df[0][3] + df[1][3] +
			df[2][3] + max(-qf[2], 0.0);
		vol[mesh_eles[i].numA] += 0.25 * mesh_eles[i].getvol();
		//Node B
		ap[mesh_eles[i].numB] = ap[mesh_eles[i].numB] + df[0][1] - df[3][1] - df[4][1] + max(-qf[0], 0.0) + max(qf[3], 0.0) + max(qf[4], 0.0);
		anb[mesh_eles[i].numB][mesh_neihaxi[mesh_eles[i].numB][mesh_eles[i].numA]] = anb[mesh_eles[i].numB][mesh_neihaxi[mesh_eles[i].numB][mesh_eles[i].numA]] - df[0][0] + df[3][0] +
			df[4][0] + max(qf[0], 0.0);
		anb[mesh_eles[i].numB][mesh_neihaxi[mesh_eles[i].numB][mesh_eles[i].numC]] = anb[mesh_eles[i].numB][mesh_neihaxi[mesh_eles[i].numB][mesh_eles[i].numC]] - df[0][2] + df[3][2] +
			df[4][2] + max(-qf[3], 0.0);
		anb[mesh_eles[i].numB][mesh_neihaxi[mesh_eles[i].numB][mesh_eles[i].numD]] = anb[mesh_eles[i].numB][mesh_neihaxi[mesh_eles[i].numB][mesh_eles[i].numD]] - df[0][3] + df[3][3] +
			df[4][3] + max(-qf[4], 0.0);
		vol[mesh_eles[i].numB] += 0.25 * mesh_eles[i].getvol();
		//Node C
		ap[mesh_eles[i].numC] = ap[mesh_eles[i].numC] + df[1][2] + df[3][2] - df[5][2] + max(-qf[1], 0.0) + max(-qf[3], 0.0) + max(qf[5], 0.0);
		anb[mesh_eles[i].numC][mesh_neihaxi[mesh_eles[i].numC][mesh_eles[i].numA]] = anb[mesh_eles[i].numC][mesh_neihaxi[mesh_eles[i].numC][mesh_eles[i].numA]] - df[1][0] - df[3][0] +
			df[5][0] + max(qf[1], 0.0);
		anb[mesh_eles[i].numC][mesh_neihaxi[mesh_eles[i].numC][mesh_eles[i].numB]] = anb[mesh_eles[i].numC][mesh_neihaxi[mesh_eles[i].numC][mesh_eles[i].numB]] - df[1][1] - df[3][1] +
			df[5][1] + max(qf[3], 0.0);
		anb[mesh_eles[i].numC][mesh_neihaxi[mesh_eles[i].numC][mesh_eles[i].numD]] = anb[mesh_eles[i].numC][mesh_neihaxi[mesh_eles[i].numC][mesh_eles[i].numD]] - df[1][3] - df[3][3] +
			df[5][3] + max(-qf[5], 0.0);
		vol[mesh_eles[i].numC] += 0.25 * mesh_eles[i].getvol();
		//Node D
		ap[mesh_eles[i].numD] = ap[mesh_eles[i].numD] + df[2][3] + df[4][3] + df[5][3] + max(-qf[2], 0.0) + max(-qf[4], 0.0) + max(-qf[5], 0.0);
		anb[mesh_eles[i].numD][mesh_neihaxi[mesh_eles[i].numD][mesh_eles[i].numA]] = anb[mesh_eles[i].numD][mesh_neihaxi[mesh_eles[i].numD][mesh_eles[i].numA]] - df[2][0] - df[4][0] -
			df[5][0] + max(qf[2], 0.0);
		anb[mesh_eles[i].numD][mesh_neihaxi[mesh_eles[i].numD][mesh_eles[i].numB]] = anb[mesh_eles[i].numD][mesh_neihaxi[mesh_eles[i].numD][mesh_eles[i].numB]] - df[2][1] - df[4][1] -
			df[5][1] + max(qf[4], 0.0);
		anb[mesh_eles[i].numD][mesh_neihaxi[mesh_eles[i].numD][mesh_eles[i].numC]] = anb[mesh_eles[i].numD][mesh_neihaxi[mesh_eles[i].numD][mesh_eles[i].numC]] - df[2][2] - df[4][2] -
			df[5][2] + max(qf[5], 0.0);
		vol[mesh_eles[i].numD] += 0.25 * mesh_eles[i].getvol();
	}
	//加边界条件
	double h = 1e16;
	for (int i = 0; i < mesh_bj.size(); i++) {
		//inlet
		double Sabc = 0.5 * norm(cross(mesh_nodes[mesh_bj[i][1]] - mesh_nodes[mesh_bj[i][0]], mesh_nodes[mesh_bj[i][2]] - mesh_nodes[mesh_bj[i][0]]));
		if (mesh_bj[i][3] == 1) {
			for (int j = 0; j < 3; j++) {
				ap[mesh_bj[i][j]] += 1.0 / 3.0 * Sabc * h;
				bp[mesh_bj[i][j]] += 1.0 / 3.0 * Sabc * h * 1.0;
			}
		}
		//outlet
		else if (mesh_bj[i][3] == 2) {
			for (int j = 0; j < 3; j++) {
				ap[mesh_bj[i][j]] += 1.0 / 3.0 * Sabc * h;
				bp[mesh_bj[i][j]] += 1.0 / 3.0 * Sabc * h * 0.0;
			}
		}
		//df/dy=0,h=0,f=finite value
	}
	//加源场:Q = a/r
	/*double Q;
	for (int i = 0; i < mesh_nodes.size(); i++) {
		double r = sqrt(mesh_nodes[i].getx() * mesh_nodes[i].getx() + mesh_nodes[i].getz() * mesh_nodes[i].getz());
		bp[i] += 2.0 / r * vol[i];
	}*/
	//拼装系数矩阵
	A.resize(num_nodes);
	for (int i = 0; i < num_nodes; i++) {
		A[i].resize(num_nodes);
		for (int j = 0; j < num_nodes; j++) {
			A[i][j] = 0.0;
			if (i == j) {
				A[i][j] += ap[i];
			}
		}
	}
	for (int i = 0; i < num_nodes; i++) {
		for (int j = 0; j < maxnnb; j++) {
			if (mesh_nei[i][j] != 0) {
				if (abs(A[i][mesh_nei[i][j]]) > Min)
					return;
				A[i][mesh_nei[i][j]] -= anb[i][j];
				if (mesh_neihaxi[i][mesh_nei[i][j]] != j) {
					return;
				}
			}
		}
	}
	//yanzheng(A, bp);
	//solvering equation
	auto x = solver_equtionGaussSeidel(A, bp);

	WritingRes(x);
}
//迭代求解方程组
double Mesh_cvfem::error(vector<vector<double>>& A, vector<double>& b, vector<double>& x, vector<double>& x_old) {
	double error1 = 0.0;
	double max_err = 0.0;
	for (int i = 0; i < x.size(); i++) {
		error1 = abs(x[i] - x_old[i]);
		if (max_err < error1) max_err = error1;
	}
	/*double error2 = 0.0, max_err1 = 0.0;
	for (int i = 0; i < x.size(); i++) {
		for (int j = 0; j < x.size(); j++) {
			error2 += A[i][j] * x[j];
		}
		error2 = abs(error2 - b[i]);
		if (max_err1 < error2) max_err1 = error2;
	}*/
	cout << "解变量最大误差：" << max_err << endl;
	//cout << "实际最大误差：" << max_err1 << endl;
	return max_err;
}
vector<double> Mesh_cvfem::solver_equtionGaussSeidel(vector<vector<double>>& A, vector<double>& b) {
	//1.调整矩阵形式x1=...,,,x2=....
	//2.初始值0
	vector<double>x;
	vector<double>x_old;
	if (A.size() != A[0].size()) {
		cout << "必须方阵" << endl;
		return x;
	}
#if(1)
	for (int i = 0; i < A.size(); i++) {
		if (A[i][i] < Min) {
			cout << "第" << i << "行对角元素小于0" << endl;
			//return x;
		}
	}
#endif
	x.resize(b.size());
	x_old.resize(b.size());
	for (int i = 0; i < x.size(); i++) {
		x[i] = 0.0;
		x_old[i] = 100.0;
	}
	int num = 0;
	double err = 100.0;
	while (num < 10000 && err > 1e-4) {
		x_old = x;
		for (int i = 0; i < A.size(); i++) {
			x[i] = 0.0;
			for (int j = 0; j < A[i].size(); j++) {
				if (j != i) {
					x[i] = x[i] - A[i][j] * x[j];
				}
			}
			x[i] += b[i];
			x[i] /= A[i][i];
		}
		num++;
		if (num % 10 == 0) {
			cout << "GaussSeidel迭代次数：" << num << endl;
			err = error(A, b, x, x_old);
		}
	}
	for (int i = 0; i < x.size(); i++) {
		if (x[i] < Min) x[i] = 0.0;
	}
	return x;
}
vector<double> Mesh_cvfem::solver_equtionJacobi(vector<vector<double>>& A, vector<double>& b) {
	vector<double>x;
	vector<double>x_old;
	if (A.size() != A[0].size()) {
		cout << "必须方阵" << endl;
		return x;
	}
#if(1)
	for (int i = 0; i < A.size(); i++) {
		if (A[i][i] < Min) {
			cout << "第" << i << "行对角元素小于0" << endl;
			//return x;
		}
	}
#endif
	x.resize(b.size());
	x_old.resize(b.size());
	for (int i = 0; i < x.size(); i++) {
		x[i] = 0.0;
		x_old[i] = 100.0;
	}
	int num = 0;
	while (num < 10000 && error(A, b, x, x_old)>1e-5) {
		x_old = x;
		for (int i = 0; i < A.size(); i++) {
			x[i] = 0.0;
			for (int j = 0; j < A[i].size(); j++) {
				if (j != i) {
					x[i] = x[i] - A[i][j] * x_old[j];
				}
			}
			x[i] += b[i];
			x[i] /= A[i][i];
		}
		num++;
		cout << "Jacobi迭代次数：" << num << endl;
	}
	return x;
}
vector<double> Mesh_cvfem::solver_equtionSOR(vector<vector<double>>& A, vector<double>& b, double w) {
	vector<double>x;
	vector<double>x_old;
	if (A.size() != A[0].size()) {
		cout << "必须方阵" << endl;
		return x;
	}
	for (int i = 0; i < A.size(); i++) {
		if (A[i][i] < Min) {
			cout << "对角必须大于0" << endl;
			return x;
		}
	}
	x.resize(b.size());
	x_old.resize(b.size());
	for (int i = 0; i < x.size(); i++) {
		x[i] = 0.0;
		x_old[i] = 100.0;
	}
	int num = 0;
	while (num < 10000 && error(A, b, x, x_old)>1e-5) {
		x_old = x;
		for (int i = 0; i < A.size(); i++) {
			x[i] = 0.0;
			for (int j = 0; j < A[i].size(); j++) {
				if (j != i) {
					x[i] = x[i] - A[i][j] * x_old[j] * w;
				}
			}
			x[i] += b[i] * w;
			x[i] /= A[i][i];
			x[i] += x_old[i] * (1 - w);
		}
		num++;
		cout << "SOR迭代次数：" << num << endl;
	}
	return x;
}
void Mesh_cvfem::WritingRes(vector<double>& x) {
	string cwd = "C:\\Users\\freedom\\Desktop\\FVM\\test\\test1.txt";
	ofstream file(cwd);
	if (!file.is_open()) return;
	for (int i = 0; i < x.size(); i++) {
		double dis = sqrt(mesh_nodes[i].getx() * mesh_nodes[i].getx() + mesh_nodes[i].getz() * mesh_nodes[i].getz());
		file << dis << " " << x[i] << endl;
	}
}
void Mesh_cvfem::yanzheng(vector<vector<double>>& A, vector<double>& b) {
	for (int i = 0; i < mesh_eles.size(); i++) {
		Point ad = mesh_eles[i].getD() - mesh_eles[i].getA();
		Point ab = mesh_eles[i].getB() - mesh_eles[i].getA();
		Point bc = mesh_eles[i].getC() - mesh_eles[i].getB();
		double t = cross(ab, bc) * ad;
		if (t < 0)
			cout << "m";
	}
	for (int i = 0; i < mesh_bj.size(); i++) {
		double temp;
		if (mesh_bj[i][3] == 1) {
			for (int j = 0; j < 3; j++) {
				temp = sqrt(mesh_nodes[mesh_bj[i][j]].getx() * mesh_nodes[mesh_bj[i][j]].getx() + mesh_nodes[mesh_bj[i][j]].getz() * mesh_nodes[mesh_bj[i][j]].getz());
				if (abs(temp - 1.0) > 1e-5) {
					cout << "m";
				}
			}
		}
		else if (mesh_bj[i][3] == 2) {
			for (int j = 0; j < 3; j++) {
				temp = sqrt(mesh_nodes[mesh_bj[i][j]].getx() * mesh_nodes[mesh_bj[i][j]].getx() + mesh_nodes[mesh_bj[i][j]].getz() * mesh_nodes[mesh_bj[i][j]].getz());
				if (abs(temp - 2.0) > 1e-5) {
					cout << "m";
				}
			}
		}
		else if (mesh_bj[i][3] == 3) {
			for (int j = 0; j < 3; j++)
				if (abs(mesh_nodes[mesh_bj[i][j]].getz()) > 1e-5) {
					cout << "m";
				}
		}
		else if (mesh_bj[i][3] == 4) {
			for (int j = 0; j < 3; j++)
				if (abs(mesh_nodes[mesh_bj[i][j]].getx()) > 1e-5) {
					cout << "m";
				}
		}
		else if (mesh_bj[i][3] == 5) {
			for (int j = 0; j < 3; j++)
				if (abs(mesh_nodes[mesh_bj[i][j]].gety()) - 1.0 > 1e-5) {
					cout << "m";
				}
		}
		else if (mesh_bj[i][3] == 6) {
			for (int j = 0; j < 3; j++)
				if (abs(mesh_nodes[mesh_bj[i][j]].gety()) > 1e-5) {
					cout << "m";
				}
		}
	}
	vector<double>x;
	double error1 = 0.0;
	x.resize(mesh_nodes.size());
	for (int i = 0; i < x.size(); i++) {
		double r = sqrt(mesh_nodes[i].getx() * mesh_nodes[i].getx() + mesh_nodes[i].getz() * mesh_nodes[i].getz());
		x[i] = 1.0 - log(r) / log(2.0);
	}
	for (int i = 0; i < x.size(); i++) {
		for (int j = 0; j < x.size(); j++) {
			error1 += A[i][j] * x[j];
		}
		error1 -= b[i];
		cout << i << "节点误差为：" << error1 << endl;
	}
	auto temp_x = solver_equtionGaussSeidel(A, bp);
	for (int i = 0; i < mesh_nodes.size(); i++) {
		double sum = ap[i] * temp_x[i];
		for (int j = 0; j < maxnnb; j++) {
			sum += anb[i][j] * temp_x[mesh_nei[i][j]];
		}
		cout << sum - bp[i] << endl;
	}

}
//瞬态计算：多一个时间步的迭代
//1. 完全显示，时间间隔必须小于 vol_i * f_inew = vol_i * f_i + dt * ( sum(ap_ij * f_g_ij) - ap_i * f_i );
//2. 完全隐式，求解线性方程组即可