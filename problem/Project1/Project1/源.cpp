#include"标头.h"
#include<memory>
#include<algorithm>
using namespace std;
extern int nn;
constexpr int max_size = 100;//编译阶段已知
void ssum(const int& a) {

}
class A{
public:
	void sum() {};
	A() {};
	int a;
	void getb() {};
protected:
	int b;
};

class B :public A {
public:
	int geta() { return a; };
	int getb() { return b; };
protected:
	int getaa() { return a; };
};


int find(vector<int> n, int num) {
	int left = 0;
	int right = n.size() - 1;
	while (left <= right) {
		int temp = (left + right) / 2;
		if (n[temp] < num) {
			left = temp+1;
		}
		else if (n[temp] > num) {
			right = temp-1;
		}
		else {
			return temp;
		}
	}
	return -1;
}

//void sort(vector<int>& nums) {
//
//	for (int j = 1; j < nums.size(); j++) {
//		for (int i = 0; i < nums.size() - j; i++)
//		{
//			if (nums[i] > nums[i + 1]) {
//				int temp = nums[i];
//				nums[i] = nums[i + 1];
//				nums[i + 1] = temp;
//			}
//		}
//	}
//}

void sum(int* a) {

}
void sum(const int* a) {

}

class tree{
public:
	tree(double Val) :left(nullptr), right(nullptr), val(Val) {};
	unique_ptr<tree>left;
	unique_ptr<tree>right;
	double val;

	
};
unique_ptr<tree> seet(unique_ptr<tree> temp, double v) {
	if (temp == nullptr) {
		return make_unique<tree>(v);
	}

	if (temp->val<v) {
		temp->right = seet(move(temp->right), v);
	}
	else {
		temp->left = seet(move(temp->left), v);
	}

	return temp;
};

void print(unique_ptr<tree>temp, int id) {
	if (temp == nullptr) {
		return;
	}
	if(id == 0)
		cout << temp->val << "left";
	else
		cout << temp->val << "right";
	print(move(temp->left), 0);
	print(move(temp->right), 1);
}

int fei(int k) {
	vector<int>dp;
	dp.resize(k+1);
	dp[0] = 1;
	dp[1] = 1;
	for (int i = 2; i <= k; i++) {
		dp[i] = dp[i - 1] + dp[i - 2];
	}

	return dp[k];
}


int main() {
	unordered_map<int, int>p;
	vector<int>nums;
	nums.push_back(1);
	p[1] = 0;
	p[3] = 2;
	cout<<p.begin()->first;
	cout << p.begin()->second;
	cout << (++p.begin())->first;
	for (auto it = p.begin(); it != p.end(); it++) {
		it->first;
	}
	sort(nums.begin(),nums.end());

	auto pp2 = p.begin();
	auto pp = p.find(1);
	auto pp1 = p.find(2);
	if (pp1 == p.end()) {
		cout << "not found";
	}
	p.erase(1);
	return 0;
}