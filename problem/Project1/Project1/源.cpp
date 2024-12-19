#include"标头.h"
#include<memory>
using namespace std;
extern int nn;
constexpr int max_size = 100;//编译阶段已知
void ssum(const int& a) {

}
class A{
	void sum() {};
	A() {};
	int a;
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

void sort(vector<int>& nums) {

	for (int j = 1; j < nums.size(); j++) {
		for (int i = 0; i < nums.size() - j; i++)
		{
			if (nums[i] > nums[i + 1]) {
				int temp = nums[i];
				nums[i] = nums[i + 1];
				nums[i + 1] = temp;
			}
		}
	}
}

void sum(int* a) {

}
void sum(const int* a) {

}

int main() {
	//find
	int target = 9;
	sum(&target);
	vector<int>a = {7,6,5,4,32,1,-1,99,111};
	sort(a);
	for (auto& aa : a) {
		cout << aa << " ";
	}
}