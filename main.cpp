#include<iostream>
#include"RNK.h"

int main() {
	setlocale(LC_ALL, "");
	RNK rna1;
	size_t size = 25;
	RNK rna2(size, C);
	RNK rna3(rna2);
	for (size_t i = 0; i < rna3.length(); i++)
		std::cout << rna3[i] << " ";
	std::cout << "\n";

//	rna3 = rna2;
	rna3[250] = T;//alloc mem
	for (size_t i = 0; i < rna3.length(); i++)
		std::cout << rna3[i] << " ";
	std::cout << "\n";

	rna1[4] = T;//alloc mem
	std::cout << rna1.length() << "\n";
	std::cout << "array:\n";
	for (size_t i = 0; i < rna1.length(); i++)
		std::cout << (Nucleotide)rna1[i] << " ";
	std::cout << "end:\n";
//	Nucleotide n = rna3[249];			//correct
//	std::cout << rna3[250] << "\n";

	rna3.trim(45);
	
	!rna3;
	for (size_t i = 0; i < rna3.length(); i++)
		std::cout << rna3[i] << " ";
	std::cout << "\n";

	std::cout << "\n" << rna1.isComplementary(rna3);

//	std::cout << "\n" << rna3[2];


	system("pause");
	return 0;
}
