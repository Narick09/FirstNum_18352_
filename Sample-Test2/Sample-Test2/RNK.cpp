#include"pch.h"
#include"RNK.h"

RNK::reference::reference(RNK * const ptr, size_t ind) :rna(ptr), ind(ind) {}

//not correct
RNK::reference & RNK::reference::operator=(Nucleotide nucl) {
	if (ind + 1 <= rna->size_) {
		//error 
		std::cout << "ERROR. RNA don't can be changed in a middle\n";
	}
	else
	{
		size_t newPhysSize = (ind + 1) / (4 * sizeof(size_t));

		if ((ind + 1) % (4 * sizeof(size_t)) != 0)
			newPhysSize++;
		//проверить, нужно ли выделять память

		if (newPhysSize <= rna->physSize) {											//changed
			//don't need to alloc
			//but we need to change old mem, in with index > rna->size, to 0
			// 1) << 2) >>

			size_t tmpShift = (sizeof(size_t) * 8 - ((rna->size_ * 2) % (sizeof(size_t) * 8)));//in bits
			if (rna->physSize != 0) {
				rna->Nuc_[rna->physSize - 1] = (rna->Nuc_[rna->physSize - 1] << tmpShift) >> tmpShift;
			}//чтобы в случае 0й размерности не вылетал/ тут мб даже и не проверять
			size_t tmpMusk = rna->getMusk(ind, nucl);
			rna->Nuc_[newPhysSize - 1] = rna->Nuc_[newPhysSize - 1] | tmpMusk;
			rna->size_ = ind + 1;
		}
		else
		{
			//need to alloc memory																//correct vrode
			//but we need to change old mem, in with index > rna->size, to 0 else
			newPhysSize = 2 * rna->physSize;											//changed
			size_t *tmpArray = new size_t[newPhysSize];
			for (size_t i = 0; i < rna->physSize; i++) {
				tmpArray[i] = rna->Nuc_[i];
			}
			size_t tmpShift = (sizeof(size_t) * 8 - ((rna->size_ * 2) % (sizeof(size_t) * 8)));//in bits
			if (rna->physSize != 0) {
				tmpArray[rna->physSize - 1] = (tmpArray[rna->physSize - 1] << tmpShift) >> tmpShift;
			}//чтобы в случае 0й размерности не вылетал

			for (size_t i = rna->physSize; i < newPhysSize; i++)
				tmpArray[i] = 0;
			//самому последнему нуклеотиду присвоить принимаемое значение - бубен с маской
			tmpArray[newPhysSize - 1] = tmpArray[newPhysSize - 1] | rna->getMusk(ind, nucl);//просто инд, а не инд+1, пушо в маске 0я позиция -
																							//это 1я физическая позиция
			delete[] rna->Nuc_;
			rna->Nuc_ = tmpArray;
			rna->physSize = newPhysSize;
			rna->size_ = ind + 1;
		}
	}
	return *this;
}

RNK::reference::operator Nucleotide() const {
	if (ind > rna->size_)
		return A;
	//	size_t Musk = rna->getMusk(i, T);
	size_t countNucInType = (4 * sizeof(size_t));
	size_t d = ind % countNucInType;
	size_t c = ind / countNucInType;

	size_t tmp = rna->Nuc_[c];
	tmp = tmp << 2 * (countNucInType - 1 - d);
	tmp = tmp >> 2 * (countNucInType - 1);
	//	std::cout << "tmp " << tmp << "\n";

	Nucleotide n = A;
	switch (tmp)
	{
	case A:
		n = A;
		//		std::cout << "A ";
		break;
	case C:
		n = C;
		//		std::cout << "C ";
		break;
	case G:
		n = G;
		//		std::cout << "G ";
		break;
	case T:
		n = T;
		//		std::cout << "T ";
		break;
	default:
		break;
	}
	return n;
}
//------------------------------------------------------------------------------------

size_t RNK::getMusk(size_t i, size_t Nuc) const {
	if (Nuc == A) return 0;
	size_t d = i % (4 * sizeof(size_t));
	size_t musk = Nuc;
	musk = musk << (d * 2);
	return musk;
}

//operator[]
RNK::reference RNK::operator[](size_t ind) {
	return reference(this, ind);
}
Nucleotide RNK::operator[](size_t ind) const {
	RNK tmp(*this);
	const reference r(&tmp, ind);
	return r;
}
//constructor
RNK::RNK(size_t s, Nucleotide n) :size_(s)
{
	size_t countNucInType = (4 * sizeof(size_t));
	if (size_ % countNucInType != 0)
		physSize = ((size_ / countNucInType) + 1);
	else
		physSize = size_ / countNucInType;
	Nuc_ = new size_t[physSize];
	size_t musk = n;
	for (size_t i = 0; i < countNucInType - 1; i++) {
		musk = musk << 2;
		musk = musk | n;
	}
	//	std::cout << "musk in constructor: " << musk << std::endl;
	for (size_t i = 0; i < physSize; i++) {
		Nuc_[i] = 0;
		Nuc_[i] = musk;
	}
}

//destructor
RNK::~RNK() {
	delete[] Nuc_;
	Nuc_ = nullptr;
}

size_t RNK::length() const {
	return this->size_;
}

//copy constructor 
RNK::RNK(const RNK& other) {
	this->size_ = other.size_;
	this->physSize = other.physSize;
	Nuc_ = new size_t[physSize];
	for (size_t i = 0; i < physSize; i++) {
		Nuc_[i] = other.Nuc_[i];
	}
}

RNK RNK::operator !() {
	for (size_t i = 0; i < physSize; i++)
		Nuc_[i] = ~Nuc_[i];
	return *this;
}

bool RNK::isComplementary(const RNK& rnk) const {
	if (size_ != rnk.size_) return false;
	bool tmp = true;
	size_t countNucInType = 4 * sizeof(size_t);
	size_t tmpSize = size_ / countNucInType;
	size_t surplus = size_ % countNucInType;
	size_t i = 0;
	while (tmp && (i < tmpSize)) {
		size_t tmp1 = ~Nuc_[i];
		if (tmp1 != rnk.Nuc_[i]) {
			tmp = false;
			break;
		}
		i++;
	}
	if (i == tmpSize) {
		for (size_t j = 0; j < surplus; j++) {
			size_t tmp1 = ~Nuc_[i];
			size_t tmp2 = rnk.Nuc_[i];
		}
	}
	return tmp;
}

RNK & RNK::operator=(const RNK & rna) {
	if (&rna == this)
		return *this;
	if (this->physSize == rna.physSize) {
		for (size_t i = 0; i < physSize; i++)
			this->Nuc_[i] = rna.Nuc_[i];
	}
	else
	{
		delete this->Nuc_;
		this->Nuc_ = new size_t[rna.physSize];
		for (size_t i = 0; i < rna.physSize; i++)
			this->Nuc_[i] = rna.Nuc_[i];
	}
	this->physSize = rna.physSize;
	this->size_ = rna.size_;

	return *this;
}

//it was so hard
//hard indexes again
RNK RNK::operator+(const RNK& rnk) const {
	if (this->size_ == 0) {
		RNK returnRnk(rnk);
		return returnRnk;
	}
	if (rnk.size_ == 0) {
		RNK returnRnk(*this);
		return returnRnk;
	}

	RNK returnRnk(this->size_ + rnk.size_, A);
	returnRnk.Nuc_ = new size_t[returnRnk.physSize];
	//from 1st RNK
	for (size_t i = 0; i < this->physSize; i++) {
		returnRnk.Nuc_[i] = this->Nuc_[i];
	}
	//from 2nd RNK
	size_t countNucInType = 4 * sizeof(size_t);
	size_t shiftLeft = 2 * ((this->size_) % countNucInType);//in bits
	size_t shiftRight = 2 * countNucInType - shiftLeft;//in bits
	returnRnk.Nuc_[this->physSize - 1] = returnRnk.Nuc_[this->physSize - 1] << shiftRight;//занулить ненужные нуклеотиды
	returnRnk.Nuc_[this->physSize - 1] = returnRnk.Nuc_[this->physSize - 1] >> shiftRight;//left-right
	//это нужно, чтобы записать со сдвигом:
	for (size_t i = this->physSize; i < returnRnk.physSize; i++) {
		returnRnk.Nuc_[i] = rnk.Nuc_[i - this->physSize] >> shiftRight;
		size_t tmp = rnk.Nuc_[i - this->physSize] << shiftLeft;
		returnRnk.Nuc_[i - 1] = returnRnk.Nuc_[i - 1] | tmp;
	}
	//----------------------------------------------------------------------------------------------------there are trables may be.(or not)
	size_t tmp = rnk.Nuc_[rnk.physSize - 1] << shiftLeft;
	returnRnk.Nuc_[returnRnk.physSize - 1] = returnRnk.Nuc_[returnRnk.physSize - 1] | tmp;
	//----------------------------------------------------------------------------------------------------
	return returnRnk;
}

void RNK::trim(size_t lastIndex) {
	if (lastIndex < this->size_) {
		this->size_ = lastIndex;
		this->physSize = (lastIndex) / (4 * sizeof(size_t));
		if ((lastIndex) % (4 * sizeof(size_t)) != 0)
			this->physSize++;
		size_t * NucNew = new size_t[this->physSize];
		for (size_t i = 0; i < this->physSize; i++) {
			NucNew[i] = Nuc_[i];
		}
		delete[] Nuc_;
		Nuc_ = NucNew;
	}
}

//need to check index
bool RNK::operator==(const RNK& rna) const {
	//	RNK tmpRnk(rna);
	//	bool tmp = true;
	if (this->size_ != rna.size_) return false;
	for (size_t i = 0; i < rna.physSize - 1; i++) {
		if (this->Nuc_[i] != rna.Nuc_[i]) {
			return false;
		}
	}
	for (size_t i = rna.size_ - 1 - ((rna.physSize - 1) * sizeof(size_t)); i < rna.size_; i++) {
		if (this->operator[](i) != rna[i]) {
			return false;
		}
	}
	return true;
}

bool RNK::operator!=(const RNK& rnk) const {
	if (*this == rnk)
		return false;
	else
		return true;
}

//need to test
RNK RNK::split(size_t ind) {
	RNK rnkTmp(this->size_ - ind, A);
	//	RNK rnkTmpThis(*this);
	size_t countNucInType = 4 * sizeof(size_t);
	size_t shiftLeft = 2 * ((this->size_) % countNucInType);//in bits
	size_t shiftRight = 2 * countNucInType - shiftLeft;//in bits
	size_t indInPhys = ind / countNucInType;
	for (size_t i = 0; i < this->physSize - indInPhys; i++) {
		rnkTmp.Nuc_[i] = this->Nuc_[i + indInPhys] >> shiftRight;
		size_t tmpShift = this->Nuc_[i + indInPhys + 1] << shiftLeft;
		rnkTmp.Nuc_[i] = rnkTmp.Nuc_[i] | tmpShift;
		//		rnkTmp[i] = (Nucleotide)rnkTmpThis[i + ind];
	}
	if (ind != 0)
		this->trim(ind);
	else
	{
		delete[]this->Nuc_;
		this->size_ = 0;
		this->physSize = 0;
	}
	return rnkTmp;
}
//without checking:
size_t RNK::cardinality(Nucleotide value) const {
	size_t returnNum = 0;
	for (size_t i = 0; i < this->size_; i++) {
		if (this->operator[](i) == value) {
			returnNum++;
		}
	}
	return returnNum;
}

std::unordered_map < Nucleotide, int, std::hash<int> > RNK::cardinality() const {
	std::unordered_map < Nucleotide, int, std::hash<int> > returnNums;
	for (size_t i = 0; i < this->size_; i++) {
		switch (this->operator[](i))
		{
		case A:
			returnNums[A]++;
			break;
		case C:
			returnNums[C]++;
			break;
		case G:
			returnNums[G]++;
			break;
		case T:
			returnNums[T]++;
			break;
		default:
			break;
		}
	}
	return returnNums;
}

std::ostream& operator << (std::ostream & out, const RNK& other) {
	for (size_t i = 0; i < other.length(); i++) {
		switch (other[i])
		{
		case A:
			out << "A ";
			break;
		case C:
			out << "C ";
			break;
		case G:
			out << "G ";
			break;
		case T:
			out << "T ";
			break;
		default:
			break;
		}
		//		out << other[i] << " ";
	}
	out << "\n";
	return out;
}
//for debug:
/*----------------------------------------------------------------------------------------------for debug
Nucleotide RNK::TmpGeterNuc(size_t i) const {
	if (i > this->size_)
		return A;
	size_t Musk = this->getMusk(i, T);
	size_t tmp = 0;
	size_t countNucInType = (4 * sizeof(size_t));
	size_t d = i % countNucInType;
	size_t c = i / countNucInType;
	tmp = Nuc_[c];
	tmp = tmp << 2 * (countNucInType - 1 - d);
	tmp = tmp >> 2 * (countNucInType - 1);
	Nucleotide n = A;
	switch (tmp)
	{
	case A:
		n = A;
		std::cout << "A ";
		break;
	case C:
		n = C;
		std::cout << "C ";
		break;
	case G:
		n = G;
		std::cout << "G ";
		break;
	case T:
		n = T;
		std::cout << "T ";
		break;
	default:
		break;
	}
	return n;
}
*/