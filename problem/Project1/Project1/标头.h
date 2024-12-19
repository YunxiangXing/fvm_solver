#pragma once
#include<iostream>
#include<map>
#include<unordered_map>
static int k = 0;
static void sum() {
	std::cout << "sum" << std::endl;
}
class test {
public:
	test() { std::cout << "num_s: " << num << std::endl; };
	~test() { std::cout << "num_x: " << num << std::endl; };
	static int const A = 10;
	static const int num;
	const int B = 5;
	const int getB() { return B; };
	bool operator<(const test& a) const {
		if (this->numm < a.numm) {
			return true;
		}
		return false;
	}
private:
	int numm = 0;
};

const int test::num = 0;//可以共同使用

extern int nn;
void P() { nn++; };