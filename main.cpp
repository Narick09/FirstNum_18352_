#include<iostream>
#include"RNK.h"


int main() {
	setlocale(LC_ALL, "");
	RNK rna1;
	size_t size = 25;
	RNK rna2(size, C);
	RNK rna3(rna2);

	std::cout << rna2;

	rna3[56] = G;
	std::cout << rna3;

	system("pause");
	return 0;
}
