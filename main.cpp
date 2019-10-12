#include<iostream>
#include"RNK.h"


int main() {
	setlocale(LC_ALL, "");
	RNK rna1;
	size_t size = 25;
	RNK rna2(size, C);
	RNK rna3(rna2);
	std::cout << rna2 <<  rna3;

	int s = 0;
	rna1[0] = T;
	rna1[1] = C;
	rna1[2] = G;
//	std::cout << rna1.getPSize() << "\n";
	rna1[3] = T;
	rna1[4] = C;
	rna1[5] = G;
//	std::cout << rna1.getPSize() << "\n";
	rna1[6] = C;
	rna1[7] = G;
	rna1[8] = T;
	rna1[9] = C;
	rna1[10] = C;
	rna1[11] = G;
	rna1[12] = T;
	rna1[13] = T;
//	std::cout << "-------------------------------------------------\n";
	rna1[14] = G;//1
	rna1[15] = C;//2
//	std::cout << "------------------------------------------------\n";
	rna1[16] = T;
//	std::cout << rna1.getPSize() << "\n";
	rna1[17] = G;
//	std::cout << rna1.getPSize() << "\n";
	rna1[18] = C;
	rna1[19] = C;
	rna1[20] = T;
	rna1[21] = T;
	rna1[22] = T;
	rna1[23] = C;
	rna1[24] = G;
	rna1[25] = C;
	rna1[26] = T;
	rna1[27] = T;
	rna1[28] = G;
	rna1[29] = C;
//	std::cout << "-------------------------------------------------\n";
	rna1[30] = T;
	rna1[31] = G;
//	std::cout << "-------------------------------------------------\n";
	rna1[32] = G;
/*
	std::cout << "\"op[]()\":\n";
	for (size_t i = 0; i < rna1.length(); i++) {
		std::cout << (Nucleotide)rna1[i] << " ";
	}
	std::cout << "\n";
	std::cout << "\"tmpGetNuc\":\n";
	for (size_t i = 0; i < rna1.length(); i++) {
		std::cout << rna1.TmpGeterNuc(i)<<" ";
	}
	std::cout << "\n";
	*/
//	std::cin >> s;
//	!rna3;
//	std::cout << "\"cout<<\":\n";
	std::cout << rna1;

	std::cout << "rna2 + rna1: \n" << rna2 + rna1;
	rna2 = rna1.split(5);
//	rna1.trim(5);
	std::cout << "rna1 trimmed: \n" << rna1 << "rna2 - splited rna1: \n" << rna2;

	

	system("pause");
	return 0;
}
