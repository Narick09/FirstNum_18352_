#include<iostream>
#include"RNK.h"


int main() {
	setlocale(LC_ALL, "");
	RNK rna1;
	size_t size = 25;
	RNK rna2(size, C);
	RNK rna3(rna2);
	std::cout << rna2;

//	rna1[1000001] = T;
//	rna1.trim(10);
	int s = 0;
	rna1[2] = A;
	rna1[3] = C;
	rna1[4] = G;
	rna1[5] = T;//
	rna1[6] = C;//not correct
	rna1[7] = G;//
	rna1[8] = T;//
	rna1[9] = C;//
	rna1[10] = T;//
	rna1[11] = G;//
	rna1[12] = C;//
	rna1[13] = T;//
	rna1[14] = G;//
	rna1[15] = C;//
	rna1[16] = A;//
	rna1[17] = G;
	rna1[18] = T;
	rna1[19] = C;
	rna1[20] = T;
	rna1[21] = T;
	rna1[22] = T;
	rna1[23] = C;
	rna1[24] = G;
	rna1[25] = C;//
	rna1[26] = T;//not correct
	rna1[27] = T;//
	rna1[28] = G;



//	std::cin >> s;
//	!rna3;
	std::cout << rna1;

//	rna1 = rna2 + rna3;
//	std::cout << "rna1: \n" << rna3 + rna2;

	system("pause");
	return 0;
}
