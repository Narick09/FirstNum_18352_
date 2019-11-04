#include"pch.h"
#include<iostream>
#include"RNK.h"
#include<cassert>
#include<ctime>

int main(int argc, char *argv[]) {
	setlocale(LC_ALL, "");
	/*RNK rnk(1000, A);
	size_t allocLength = rnk.capacity();
	assert(allocLength >= 1000 * 2 / 8 / sizeof(size_t));
	size_t startSetClass = clock();
	RNK set;

	size_t startSetBig = clock();
	set[1000000] = A;
	size_t endSetBig = clock();

	assert(allocLength < set.length());

	std::cout << "runtime copy-constructor: " << startSetBig - startSetClass;
	std::cout << "\nruntime memory allocate: " << endSetBig - startSetBig << "\n";

	RNK tmp;
	size_t startSplit = clock();
	tmp = set.split(100000);
	size_t endSplit = clock();
	std::cout << "\nruntime split: " << ((endSplit - startSplit)) << "ms = " << (double)((endSplit - startSplit)) / 60000 << "minutes\n";
	//	std::cout << set;
		//need check split 
//	system("pause");
	*/
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}