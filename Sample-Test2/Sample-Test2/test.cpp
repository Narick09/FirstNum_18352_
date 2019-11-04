#include "pch.h"
#include"RNK.h"

//test constructors
TEST(TestInit, TestRnkConstructor) {
	Nucleotide n = C;
	const RNK rn(1000, n);
	RNK rn2;
	bool tmp = true;
	for (size_t i = 0; i < rn.length(); i++) {
		if (rn[i] != n) {
			tmp = false;
			break;
		}
	}

	EXPECT_EQ(1000, rn.length());
	EXPECT_EQ(0, rn2.length());
	EXPECT_TRUE(tmp);
}
//copy and "=="
TEST(TestCopy, TestCopyRNK) {
	Nucleotide n = C;
	const RNK rn1(1000, n);					//chenged to const operator[] - need to check
	const RNK rn2(rn1);
	EXPECT_EQ(rn1.length(), rn2.length());
	EXPECT_EQ(rn1.capacity(), rn2.capacity());
	bool tmp = true;
	bool IsEqual = (rn1 == rn2);
	for (size_t i = 0; i < rn1.length(); i++) {
		if (rn1[i] != rn2[i]) {
			tmp = false;
			break;
		}
	}
	EXPECT_TRUE(tmp);
	EXPECT_TRUE(IsEqual);
}
//operator=
TEST(TestAssigning, TestOperatorAssigningRNK) {
	Nucleotide n = C;
	RNK rn1(1000, n);
	RNK rn2;
	rn2 = rn1;
	EXPECT_EQ(rn1.length(), rn2.length());
	EXPECT_EQ(rn1.capacity(), rn2.capacity());
	bool tmp = true;
	bool IsEqual = (rn1 == rn2);
	for (size_t i = 0; i < rn1.length(); i++) {
		if (rn1[i] != rn2[i]) {
			tmp = false;
			break;
		}
	}
	EXPECT_TRUE(tmp);
	EXPECT_TRUE(IsEqual);
}
TEST(TestAssigning, TestAssigningThisToThis) {
	Nucleotide n = C;
	RNK rn1(1000, n);
	RNK rn2(rn1);
	rn1 = rn1;
	bool tmp = true;
	for (size_t i = 0; i < rn1.length(); i++) {
		if (rn1[i] != rn2[i]) {
			tmp = false;
			break;
		}
	}
	EXPECT_TRUE(rn1.length() == rn2.length());
	EXPECT_TRUE(tmp);
}

//operator!
TEST(TestNot, TestNotRnk) {
	Nucleotide n = C;
	RNK rn1(1000, n);
	RNK rn2;
	rn2 = rn1;
	EXPECT_EQ(rn1.length(), rn2.length());
	EXPECT_EQ(rn1.capacity(), rn2.capacity());
	bool tmp = true;
	bool IsEqual = (rn1 == rn2);
	for (size_t i = 0; i < rn1.length(); i++) {
		if (rn1[i] != rn2[i]) {
			tmp = false;
			break;
		}
	}
	EXPECT_TRUE(tmp);
	EXPECT_TRUE(IsEqual);
}

//void trim();
TEST(TrimTest, TestTrimRNK) {
	RNK rn2(1000, G);
	rn2.trim(0);
	EXPECT_EQ(rn2.length(), 0);
	EXPECT_EQ(rn2.capacity(), 0);
}

//run time
TEST(TestTime, RunTimeSmall) {
	RNK rn;
	int size = 10000;
	//	size_t start = clock();
	for (int i = 0; i < size; i++) {
		rn[i] = G;//(Nucleotide)((i*i) % 4);
	}
//	size = 1000000000;
	//	size_t end = clock();
	//	std::cout << "runTime: " << end - start << "\n";
	EXPECT_TRUE(true);
}
//run time with 100000000
TEST(TestTime, RunTimeBig) {
	RNK rn;
	int size = 100000000;
	//	size_t start = clock();
	for (int i = 0; i < size; i++) {
//		rn[i] = G;//(Nucleotide)((i*i) % 4);
	}
	EXPECT_TRUE(true);
}

TEST(TestOperatorSum, TestCorrect) {
	size_t size1 = 1000;
	size_t size2 = 1000;
	Nucleotide n1 = C;
	Nucleotide n2 = G;

	RNK rn1(size1, n1);
	RNK rn2(size2, n2);
	const RNK rn3 = rn1 + rn2;
	bool isEqual = true;
	for (size_t i = 0; i < size1; i++) {
		if (rn3[i] != n1) {
			isEqual = false;
			break;
		}
	}
	if (isEqual == true) {
		for (size_t i = size1; i < rn3.length(); i++) {
			if (rn3[i] != n2) {
				isEqual = false;
				break;
			}
		}
	}

	EXPECT_TRUE(isEqual);
}
TEST(TestOperatorSum, TestBoundaryCauseSizeOf) {
	size_t size1 = sizeof(size_t);
	size_t size2 = 3 * sizeof(size_t);
	Nucleotide n1 = C;
	Nucleotide n2 = G;

	RNK rn1(size1, n1);
	RNK rn2(size2, n2);
	const RNK rn3 = rn1 + rn2;
	bool isEqual = true;
	for (size_t i = 0; i < size1; i++) {
		if (rn3[i] != n1) {
			isEqual = false;
			break;
		}
	}
	if (isEqual == true) {
		for (size_t i = size1; i < rn3.length(); i++) {
			if (rn3[i] != n2) {
				isEqual = false;
				break;
			}
		}
	}
	EXPECT_TRUE(isEqual);
}
TEST(TestOperatorSum, TestBoundaryCauseSizeOfMinusOne) {
	size_t size1 = (2 * sizeof(size_t)) - 1;
	size_t size2 = 3 * sizeof(size_t);
	Nucleotide n1 = C;
	Nucleotide n2 = G;

	RNK rn1(size1, n1);
	RNK rn2(size2, n2);
	const RNK rn3 = rn1 + rn2;
	bool isEqual = true;
	for (size_t i = 0; i < size1; i++) {
		if (rn3[i] != n1) {
			isEqual = false;
			break;
		}
	}
	if (isEqual == true) {
		for (size_t i = size1; i < rn3.length(); i++) {
			if (rn3[i] != n2) {
				isEqual = false;
				break;
			}
		}
	}
	EXPECT_TRUE(isEqual);
}
//operator!=
TEST(TestNotEqual, TestNotEqualRnk) {
	Nucleotide n = C;
	Nucleotide N = G;
	RNK rn1(1000, n);
	RNK rn2(1000, N);
	RNK rn3;
//	rn2 = rn1;
//	!rn1;
//	EXPECT_FALSE(rn1.length() == rn2.length());
//	EXPECT_EQ(rn1.capacity(), rn2.capacity());
	bool IsEqual = (rn1 == rn2);
	EXPECT_FALSE(IsEqual);
	IsEqual = (rn1 == rn3);
	EXPECT_FALSE(IsEqual);
	IsEqual = (rn2 == rn3);
	EXPECT_FALSE(IsEqual);
}

TEST(TestIsCompl, TestIsComplRnk) {
	size_t size = 10000;
	Nucleotide n = C;
	RNK rn1(size, n);
	RNK rn2 = rn1;
	!rn1;
	EXPECT_TRUE(rn1.isComplementary(rn2));
}

TEST(TestSplit, TestSplitRunTime) {
	size_t size = 10000;// 000;
	size_t splitNum = 100;//00;
	Nucleotide n = C;
	RNK rn1(size, n);
	const RNK rn2 = rn1.split(splitNum);
	!rn1;
	size_t tmpSize = rn1.length();
	bool tmp = true;
	for (size_t i = 0; i < tmpSize; i++) {
		if (((Nucleotide)rn1[i]) != G) {
			tmp = false;
			break;
		}
		if (rn2[i + tmpSize] != C) {
			tmp = false;
			break;
		}
	}
	EXPECT_TRUE(tmp);
}