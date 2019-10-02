#include<iostream>
#include"RNK.h"

int main() {
	setlocale(LC_ALL, "");
	RNK rna1;
	size_t size = 25;
	RNK rna2(size, C);
	RNK rna3(rna2);
	for (size_t i = 0; i < size; i++)
		rna3[i];
	std::cout << "\n";

	!rna3;
	for (size_t i = 0; i < size; i++)
		rna3[i];
	std::cout << "\n";

	std::cout << rna2.isComplementary(rna3) << "\n";
	
	rna3 = rna2;
//	rna3[0] = T;
//	rna3[5] = G;
	rna3[1] = A;
	rna3[250] = T;
	for (size_t i = 0; i < rna3.getSize() + 1; i++)
		std::cout << rna3[i] << " ";
	std::cout << "\n";
	
	Nucleotide n = rna3[249];			//alloc mem
	std::cout << rna3[250] << "\n";


	rna3.trim(48);
	rna1 = rna3;
	!rna3;
	for (size_t i = 0; i < rna3.getSize() + 1; i++)
		std::cout << rna3[i] << " ";
	std::cout << "\n";
	for (size_t i = 0; i < rna1.getSize() + 1; i++)
		std::cout << rna1[i] << " ";
	std::cout << "\n" << rna1.isComplementary(rna3);

//	std::cout << "\n" << rna3[2];

//	rna1 = rna2; // rna1.operator=(rna2);
//	Nucleotide n = rna3.TmpGeterNuc(0);//(rna3.operator[](ind)).operator Nucleotide(&n);

//	reference r = rna3.operator[](0);//именно для доп класса реализован =.
//	r.operator=(T);

//	rna3[1] = A; //TA
//	rna3[10] = A; //ERROR

//	rna3[0] = rna2[10];
//	int a = 9;
//	int & i = a, &j = a;
//	i = j;
//	std::cout << "Hello, World!" << std::endl;
	
	return 0;
}
