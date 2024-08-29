#include"FVM.h"
using namespace FVM;
void Fvm::pgrad(Point& g, int id, int i, Point& df, Point& Sf,double y) {
	if (id != -1) {
		gradcon temp;
		temp.data = 1.0 / mesh_eles[i].getvol() * (0.5 * df * g + y) * Sf;
		temp.id = id;
		Gradelements[i].push_back(temp);
	}
}
void Fvm::cal_gradient1(vector<double>& fai) {
	//step1:计算全局梯度使用现有的fai
	Grad.resize(mesh_eles.size());
	for (int i = 0; i < mesh_eles.size(); i++) {
		vector<Point>Sf = mesh_eles[i].getelementsurfacenormal();
		uint64_t id[4];
		id[0] = faceid(mesh_eles[i].numA, mesh_eles[i].numB, mesh_eles[i].numC);
		id[1] = faceid(mesh_eles[i].numA, mesh_eles[i].numB, mesh_eles[i].numD);
		id[2] = faceid(mesh_eles[i].numA, mesh_eles[i].numC, mesh_eles[i].numD);
		id[3] = faceid(mesh_eles[i].numB, mesh_eles[i].numC, mesh_eles[i].numD);
		Point f[4];
		f[0] = mesh_eles[i].getface_ABC();
		f[1] = mesh_eles[i].getface_ABD();
		f[2] = mesh_eles[i].getface_ACD();
		f[3] = mesh_eles[i].getface_BCD();
		//step1: 判断是内部面还是外部面
		double gDiff[4] = { 0.0,0.0,0.0,0.0 };
		Point grad_fai_c;
		for (int j = 0; j < 4; j++) {
			if (sharedface[id[j]][1] != -1) {
				int share;
				if (sharedface[id[j]][0] != i) share = sharedface[id[j]][0];
				else share = sharedface[id[j]][1];
				//gradfai_f * Tf
				//gradfai_f = g * fai_c + (1 - g) * fai_F
				Point dCfi, dfiF;
				dCfi = f[j] - mesh_eles[i].getcent();
				dfiF = mesh_eles[share].getcent() - f[j];
				float gfi = norm(dfiF) / (norm(dCfi) + norm(dfiF));
				//fai_fj
				double fai_fj = gfi * fai[i] + (1.0 - gfi) * fai[share];
				grad_fai_c = grad_fai_c + fai_fj * Sf[j];
			}
			else {
				//边界
				if (fvm_bj[id[j]] == 1) grad_fai_c = grad_fai_c + 1.0 * Sf[j];
				else if (fvm_bj[id[j]] == 2) grad_fai_c = grad_fai_c + 0.0 * Sf[j];
				else grad_fai_c = grad_fai_c + fai[i] * Sf[j];
			}
		}
		//grad_faiC
		grad_fai_c = 1.0 / mesh_eles[i].getvol() * grad_fai_c;
		Grad[i] = grad_fai_c;
	}

	//step2:update
	vector<Point> Grad_new;
	Grad_new.resize(mesh_eles.size());
	for (int i = 0; i < mesh_eles.size(); i++) {
		vector<Point>Sf = mesh_eles[i].getelementsurfacenormal();
		uint64_t id[4];
		id[0] = faceid(mesh_eles[i].numA, mesh_eles[i].numB, mesh_eles[i].numC);
		id[1] = faceid(mesh_eles[i].numA, mesh_eles[i].numB, mesh_eles[i].numD);
		id[2] = faceid(mesh_eles[i].numA, mesh_eles[i].numC, mesh_eles[i].numD);
		id[3] = faceid(mesh_eles[i].numB, mesh_eles[i].numC, mesh_eles[i].numD);
		Point f[4];
		f[0] = mesh_eles[i].getface_ABC();
		f[1] = mesh_eles[i].getface_ABD();
		f[2] = mesh_eles[i].getface_ACD();
		f[3] = mesh_eles[i].getface_BCD();
		for (int j = 0; j < 4; j++) {
			if (sharedface[id[j]][1] != -1) {
				int share;
				if (sharedface[id[j]][0] != i) share = sharedface[id[j]][0];
				else share = sharedface[id[j]][1];
				Point dCfi, dfiF;
				dCfi = f[j] - mesh_eles[i].getcent();
				dfiF = mesh_eles[share].getcent() - f[j];
				float gfi = norm(dfiF) / (norm(dCfi) + norm(dfiF));
				double fai_f = gfi * fai[i] + (1.0 - gfi) * fai[share];
				double fai_fnew = fai_f + gfi * Grad[i] * (f[j] - mesh_eles[i].getcent()) + (1.0 - gfi) * Grad[share] * (f[j] - mesh_eles[share].getcent());
				Grad_new[i] = Grad_new[i] + fai_fnew * Sf[j];
			}
			else {
				//边界
				if (fvm_bj[id[j]] == 1) Grad_new[i] = Grad_new[i] + 1.0 * Sf[j];
				else if (fvm_bj[id[j]] == 2) Grad_new[i] = Grad_new[i] + 0.0 * Sf[j];
				else Grad_new[i] = Grad_new[i] + fai[i] * Sf[j];
			}
		
		}
		Grad_new[i] = 1.0 / mesh_eles[i].getvol() * Grad_new[i];
	}

	Grad = move(Grad_new);


	for (int i = 0; i < mesh_eles.size(); i++) {
		vector<Point>Sf = mesh_eles[i].getelementsurfacenormal();
		uint64_t id[4];
		id[0] = faceid(mesh_eles[i].numA, mesh_eles[i].numB, mesh_eles[i].numC);
		id[1] = faceid(mesh_eles[i].numA, mesh_eles[i].numB, mesh_eles[i].numD);
		id[2] = faceid(mesh_eles[i].numA, mesh_eles[i].numC, mesh_eles[i].numD);
		id[3] = faceid(mesh_eles[i].numB, mesh_eles[i].numC, mesh_eles[i].numD);
		Point f[4];
		f[0] = mesh_eles[i].getface_ABC();
		f[1] = mesh_eles[i].getface_ABD();
		f[2] = mesh_eles[i].getface_ACD();
		f[3] = mesh_eles[i].getface_BCD();
		//step1: 判断是内部面还是外部面
		double gDiff[4] = { 0.0,0.0,0.0,0.0 };
		for (int j = 0; j < 4; j++) {
			if (sharedface[id[j]][1] != -1) {
				int share;
				if (sharedface[id[j]][0] != i) share = sharedface[id[j]][0];
				else share = sharedface[id[j]][1];
				Point CF = mesh_eles[share].getcent() - mesh_eles[i].getcent();
				double Dcf = norm(CF);
				Point ecf = 1.0 / Dcf * CF;
				Point Ef = (Sf[j] * Sf[j]) / (ecf * Sf[j]) * ecf;
				Point Tf;
				Point dCfi, dfiF;
				dCfi = f[j] - mesh_eles[i].getcent();
				dfiF = mesh_eles[share].getcent() - f[j];
				float gfi = norm(dfiF) / (norm(dCfi) + norm(dfiF));
				//fai_fj
				Point grad_fai_fj = gfi * Grad[i] + (1.0 - gfi) * Grad[share];
				Tf = Sf[j] - Ef;
				bp[i] += grad_fai_fj * Tf * (gfi * diff_k[i] + (1.0 - gfi) * diff_k[share]);
			}
			else {
				//边界条件中的交叉扩散项   grad_fai_b * kiff * Tb
				/*Point ecb = f[j] - mesh_eles[i].getcent();
				double dcb = norm(ecb);
				ecb = 1.0 / dcb * ecb;
				Point Eb = (Sf[j] * Sf[j]) / (ecb * Sf[j]) * ecb;
				gDiff[j] = norm(Eb) / dcb;
				Point Tb = Sf[j] - Eb;
				Point ef = 1.0 / norm(Sf[j]) * Sf[j];*/
				//指定边界通过向内差分计算梯度 grad_fai_b
				/*if (fvm_bj[id[j]] == 1) {
					Point d = f[j] - mesh_eles[i].getcent();
					bp[i] += (1.0 - fai[i]) / norm(d) * ef * Tb * diff_k[i];
				}
				else if (fvm_bj[id[j]] == 2) {
					Point d = f[j] - mesh_eles[i].getcent();
					bp[i] += (0.0 - fai[i]) / norm(d) * ef * Tb * diff_k[i];
				}*/
				//bp[i] += Grad[i] * Tb * diff_k[i];
			}
		}
	}
}
void Fvm::cal_gradient(vector<double>& fai) {
	//step1:计算全局梯度使用现有的fai
	Grad.resize(mesh_eles.size());
	for (int i = 0; i < mesh_eles.size(); i++) {
		vector<Point>Sf = mesh_eles[i].getelementsurfacenormal();
		uint64_t id[4];
		id[0] = faceid(mesh_eles[i].numA, mesh_eles[i].numB, mesh_eles[i].numC);
		id[1] = faceid(mesh_eles[i].numA, mesh_eles[i].numB, mesh_eles[i].numD);
		id[2] = faceid(mesh_eles[i].numA, mesh_eles[i].numC, mesh_eles[i].numD);
		id[3] = faceid(mesh_eles[i].numB, mesh_eles[i].numC, mesh_eles[i].numD);
		Point f[4];
		f[0] = mesh_eles[i].getface_ABC();
		f[1] = mesh_eles[i].getface_ABD();
		f[2] = mesh_eles[i].getface_ACD();
		f[3] = mesh_eles[i].getface_BCD();
		//step1: 判断是内部面还是外部面
		double gDiff[4] = { 0.0,0.0,0.0,0.0 };
		Point grad_fai_c;
		for (int j = 0; j < 4; j++) {
			if (sharedface[id[j]][1] != -1) {
				int share;
				if (sharedface[id[j]][0] != i) share = sharedface[id[j]][0];
				else share = sharedface[id[j]][1];
				//gradfai_f * Tf
				//gradfai_f = g * fai_c + (1 - g) * fai_F
				Point dCfi, dfiF;
				dCfi = f[j] - mesh_eles[i].getcent();
				dfiF = mesh_eles[share].getcent() - f[j];
				float gfi = norm(dfiF) / (norm(dCfi) + norm(dfiF));
				//fai_fj
				double fai_fj = gfi * fai[i] + (1.0 - gfi) * fai[share];
				grad_fai_c = grad_fai_c + fai_fj * Sf[j];
			}
			else {
				//边界
				if(fvm_bj[id[j]] == 1) grad_fai_c = grad_fai_c + 1.0 * Sf[j];
				else if(fvm_bj[id[j]] == 2) grad_fai_c = grad_fai_c + 0.0 * Sf[j];
				else grad_fai_c = grad_fai_c + fai[i] * Sf[j];
			}
		}
		//grad_faiC
		grad_fai_c = 1.0 / mesh_eles[i].getvol() * grad_fai_c;
		Grad[i] = grad_fai_c;
	}


	for (int i = 0; i < mesh_eles.size(); i++) {
		vector<Point>Sf = mesh_eles[i].getelementsurfacenormal();
		uint64_t id[4];
		id[0] = faceid(mesh_eles[i].numA, mesh_eles[i].numB, mesh_eles[i].numC);
		id[1] = faceid(mesh_eles[i].numA, mesh_eles[i].numB, mesh_eles[i].numD);
		id[2] = faceid(mesh_eles[i].numA, mesh_eles[i].numC, mesh_eles[i].numD);
		id[3] = faceid(mesh_eles[i].numB, mesh_eles[i].numC, mesh_eles[i].numD);
		Point f[4];
		f[0] = mesh_eles[i].getface_ABC();
		f[1] = mesh_eles[i].getface_ABD();
		f[2] = mesh_eles[i].getface_ACD();
		f[3] = mesh_eles[i].getface_BCD();
		//step1: 判断是内部面还是外部面
		double gDiff[4] = { 0.0,0.0,0.0,0.0 };
		for (int j = 0; j < 4; j++) {
			if (sharedface[id[j]][1] != -1) {
				int share;
				if (sharedface[id[j]][0] != i) share = sharedface[id[j]][0];
				else share = sharedface[id[j]][1];
				Point CF = mesh_eles[share].getcent() - mesh_eles[i].getcent();
				double Dcf = norm(CF);
				Point ecf = 1.0 / Dcf * CF;
				Point Ef = (Sf[j] * Sf[j]) / (ecf * Sf[j]) * ecf;
				Point Tf;
				Point dCfi, dfiF;
				dCfi = f[j] - mesh_eles[i].getcent();
				dfiF = mesh_eles[share].getcent() - f[j];
				float gfi = norm(dfiF) / (norm(dCfi) + norm(dfiF));
				//fai_fj
				Point grad_fai_fj = gfi * Grad[i] + (1.0 - gfi) * Grad[share];
				Tf = Sf[j] - Ef;
				bp[i] += grad_fai_fj * Tf * (gfi * diff_k[i] + (1.0 - gfi) * diff_k[share]);
			}
			else {
				//边界条件中的交叉扩散项   grad_fai_b * kiff * Tb
				Point ecb = f[j] - mesh_eles[i].getcent();
				double dcb = norm(ecb);
				ecb = 1.0 / dcb * ecb;
				Point Eb = (Sf[j] * Sf[j]) / (ecb * Sf[j]) * ecb;
				gDiff[j] = norm(Eb) / dcb;
				Point Tb = Sf[j] - Eb;
				Point ef = 1.0 / norm(Sf[j]) * Sf[j];
				//指定边界通过向内差分计算梯度 grad_fai_b
				/*if (fvm_bj[id[j]] == 1) {
					Point d = f[j] - mesh_eles[i].getcent();
					bp[i] += (1.0 - fai[i]) / norm(d) * ef * Tb * diff_k[i];
				}
				else if (fvm_bj[id[j]] == 2) {
					Point d = f[j] - mesh_eles[i].getcent();
					bp[i] += (0.0 - fai[i]) / norm(d) * ef * Tb * diff_k[i];
				}*/
				//bp[i] += Grad[i] * Tb * diff_k[i];
			}
		}
	}
}
void Fvm::findsharedface() {
	uint64_t num = mesh_nodes.size();
	for (int i = 0; i < mesh_eles.size(); i++) {
		vector<uint64_t>temp;
		vector<int>data;
		data.resize(2);
		data[1] = -1;
		temp.push_back(mesh_eles[i].numA);
		temp.push_back(mesh_eles[i].numB);
		temp.push_back(mesh_eles[i].numC);
		temp.push_back(mesh_eles[i].numD);
		sort(temp.begin(), temp.end());
		uint64_t id[4];
		uint64_t m = temp[3] * num * num;
		id[0] = temp[0] + temp[1] * num + temp[2] * num * num;
		id[1] = temp[0] + temp[1] * num + temp[3] * num * num;
		id[2] = temp[0] + temp[2] * num + temp[3] * num * num;
		id[3] = temp[1] + temp[2] * num + temp[3] * num * num;
		data[0] = i;
		for (int j = 0; j < 4; j++) {
			auto result = sharedface.insert({ id[j],data });
			if (!result.second) {
				//已经存放过
				auto temp = sharedface[id[j]];
				temp[1] = data[0];
				sharedface[id[j]] = temp;
			}
		}
	}
	for (int i = 0; i < mesh_bj.size(); i++) {
		uint64_t id = faceid(mesh_bj[i][0], mesh_bj[i][1], mesh_bj[i][2]);
		auto r = sharedface.find(id);
		if (r == sharedface.end()) {
			cout << "Error in " << __LINE__ << endl;
		}
		if (sharedface[id][1] != -1) {
			cout << "Error in " << __LINE__ << endl;
		}
	}
}
inline uint64_t Fvm::faceid(uint64_t a, uint64_t b, uint64_t c) {
	vector<uint64_t>temp;
	temp.push_back(a);
	temp.push_back(b);
	temp.push_back(c);
	sort(temp.begin(), temp.end());
	uint64_t num = mesh_nodes.size();
	uint64_t id = temp[0] + temp[1] * num + temp[2] * num * num;
	return id;
}
void Fvm::init(string cwd) {

	readRmsh(cwd);
	//findneighbor();
	findsharedface();
	//cal_gradient();

	ap.clear();
	anb.clear();
	bp.clear();
	ap.resize(mesh_eles.size());
	anb.resize(ap.size());
	A.clear();
	A.resize(ap.size());
	bp.resize(ap.size());
	for (int i = 0; i < ap.size(); i++) {
		ap[i] = 0.0;
		bp[i] = 0.0;
		anb[i].resize(ap.size());
		A[i].resize(ap.size());
		for (int j = 0; j < ap.size(); j++) {
			anb[i][j] = 0.0;
			A[i][j] = 0.0;
		}
	}

	for (int i = 0; i < mesh_bj.size(); i++) {
		uint64_t id = faceid(mesh_bj[i][0], mesh_bj[i][1], mesh_bj[i][2]);
		fvm_bj.insert({ id, mesh_bj[i][3] });
	}

	diff_k.resize(mesh_eles.size());
	for (int i = 0; i < mesh_eles.size(); i++) {
		double dis = sqrt(mesh_eles[i].getcent().getx() * mesh_eles[i].getcent().getx() +
			mesh_eles[i].getcent().getz() * mesh_eles[i].getcent().getz());
		diff_k[i] = 1.0 / dis;

		bp[i] += mesh_eles[i].getvol() * 1.0 / dis;
	}
}
void Fvm::cal_Diff(string cwd) {
	//计算单元面与哪个单元连接：
	//1. 建立单元面哈希表

	init(cwd);
	for (int i = 0; i < mesh_eles.size(); i++) {
		vector<Point>Sf = mesh_eles[i].getelementsurfacenormal();
		uint64_t id[4];
		id[0] = faceid(mesh_eles[i].numA, mesh_eles[i].numB, mesh_eles[i].numC);
		id[1] = faceid(mesh_eles[i].numA, mesh_eles[i].numB, mesh_eles[i].numD);
		id[2] = faceid(mesh_eles[i].numA, mesh_eles[i].numC, mesh_eles[i].numD);
		id[3] = faceid(mesh_eles[i].numB, mesh_eles[i].numC, mesh_eles[i].numD);
		Point f[4];
		f[0] = mesh_eles[i].getface_ABC();
		f[1] = mesh_eles[i].getface_ABD();
		f[2] = mesh_eles[i].getface_ACD();
		f[3] = mesh_eles[i].getface_BCD();
		//step1: 判断是内部面还是外部面
		double gDiff[4] = { 0.0,0.0,0.0,0.0 };
		float gfi = 0.0;
		for (int j = 0; j < 4; j++) {
			if (sharedface[id[j]][1] == -1) {
				//边界条件 id[j]
				if (fvm_bj[id[j]] == 1 || fvm_bj[id[j]] == 2) {
					Point ecb = f[j] - mesh_eles[i].getcent();
					double dcb = norm(ecb);
					ecb = 1.0 / dcb * ecb;
					Point Eb = (Sf[j] * Sf[j]) / (ecb * Sf[j]) * ecb;
					gDiff[j] = norm(Eb) / dcb;
					Point Tb = Sf[j] - Eb;
					if (fvm_bj[id[j]] == 1) {
						bp[i] += gDiff[j] * 1.0 * diff_k[i];

						/*for (int k = 0; k < Gradelements[i].size(); k++) {
							anb[i][Gradelements[i][k].id] -= Gradelements[i][k].data * Tb;
						}*/

					}
					else { 
						bp[i] += gDiff[j] * 0.0 * diff_k[i]; 

						/*for (int k = 0; k < Gradelements[i].size(); k++) {
							anb[i][Gradelements[i][k].id] -= Gradelements[i][k].data * Tb;
						}*/
					}
				}
			}
			else {
				int share;
				if (sharedface[id[j]][0] != i) share = sharedface[id[j]][0];
				else share = sharedface[id[j]][1];
				Point CF = mesh_eles[share].getcent() - mesh_eles[i].getcent();
				double Dcf = norm(CF);
				Point ecf = 1.0 / Dcf * CF;
				Point Ef = (Sf[j] * Sf[j]) / (ecf * Sf[j]) * ecf;
				gDiff[j] = norm(Ef) / Dcf;
				//gradfai_f * Tf
				//gradfai_f = g * fai_c + (1 - g) * fai_F
				Point Tf;
				Point dCfi, dfiF;
				Tf = Sf[j] - Ef;
				dCfi = f[j] - mesh_eles[i].getcent();
				dfiF = mesh_eles[share].getcent() - f[j];
				gfi = norm(dCfi) / (norm(dCfi) + norm(dfiF));
				gDiff[j] = gDiff[j] * (gfi * diff_k[share] + (1.0 - gfi) * diff_k[i]);
				anb[i][share] += -gDiff[j];
				//grad_ffi = gfi * grad_fF + (1 - gfi) * grad_fC
				//插值因子也需要修正
				//Point grad_ffi;
				//for (int k = 0; k < Gradelements[share].size(); k++) {
				//	anb[i][Gradelements[share][k].id] -= (gfi * Gradelements[share][k].data * Tf);// -gfi * Gradelements[share][k].data * ecf * ecf * Tf;
				//}
				//for (int k = 0; k < Gradelements[i].size(); k++) {
				//	anb[i][Gradelements[i][k].id] -= ((1.0 - gfi) * Gradelements[i][k].data * Tf);// -(1.0 - gfi) * Gradelements[i][k].data * ecf * ecf * Tf;
				//}
				//anb[i][share] -= ecf * Tf / Dcf;
				//anb[i][i] -= -1.0 * ecf * Tf / Dcf;
			}
		}
		ap[i] += gDiff[0] + gDiff[1] + gDiff[2] + gDiff[3];
	}

	//AX = b
	for (int i = 0; i < ap.size(); i++) {
		for (int j = 0; j < ap.size(); j++) {
			if (i == j) {
				A[i][j] = anb[i][j] + ap[i];
			}
			else {
				A[i][j] = anb[i][j];
			}
		}
	}

	//Write(bp);

	auto res = solver_equtionGaussSeidel(A, bp);
	
	string cwd1 = "C:\\Users\\yunxiang.xing\\Desktop\\LBM\\fvm_solver\\test_init";
	string cwd11 = ".txt";
	Write(res, cwd1 + cwd11);

	cout << endl;
	cout << "首次迭代已结束" << endl;
	cout << endl;

	//迭代求解
	//step1:插值计算全局梯度
	//step2:插值计算全局面梯度
	//step3:更新bp,重新求解
	vector<double> res_old = res;
	int num = 0;
	double RES = 100.0;

	double la = 1.0;
	while (RES > 1e-3) {
		num++;
		cout << "第" << num << "次迭代开始" << endl << endl;
		cal_gradient1(res);
		for (int i = 0; i < mesh_eles.size(); i++) {
			ap[i] = ap[i] / la;
			bp[i] = bp[i] + (1 - la) * ap[i] / la * res[i];
			A[i][i] = anb[i][i] + ap[i];
		}
		res_old = res;
		res = solver_equtionGaussSeidel(A, bp);
		RES = error_Fvm(res, res_old);
		cout << "误差为" << RES << endl;
		Write(res, cwd1 + '_' + to_string(num) + cwd11);
	}

	string cwd2 = "C:\\Users\\yunxiang.xing\\Desktop\\LBM\\fvm_solver\\test_res.txt";
	Write(res, cwd2);
}
double Fvm::error_Fvm(vector<double>& res, vector<double>& res_old) {
	if (res.size() != res_old.size()) {
		cout << "Error in" << __LINE__ << endl;
		exit(1);
	}
	double RES = 0.0;
	for (int i = 0; i < res.size(); i++) {
		double temp = abs(res[i] - res_old[i]);
		RES += temp;
	}
	return RES;
}
void Fvm::Write(vector<double>& x, string cwd) {
	ofstream file(cwd);
	if (!file.is_open()) return;
	for (int i = 0; i < x.size(); i++) {
		double dis = sqrt(mesh_eles[i].getcent().getx() * mesh_eles[i].getcent().getx() + 
			 mesh_eles[i].getcent().getz() * mesh_eles[i].getcent().getz());
		file << dis << " " << x[i] << endl;
	}
}