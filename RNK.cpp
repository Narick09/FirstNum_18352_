#include"RNK.h"

RNK::reference::reference(RNK * const ptr, size_t ind):rna(ptr), ind(ind) {}

//correct vrode. It was veryvery hard чтобы не запутаться в индексах. Ух. Нужны еще тесты...
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
		if (newPhysSize == rna->physSize) {
			//don't need to alloc
			//but we need to change old mem, in with index > rna->size, to 0
			// 1) << 2) >>

			size_t tmpShift = 8 * (sizeof(size_t) - (rna->physSize % sizeof(size_t))) - 2;
			if (rna->physSize != 0) {
				rna->Nuc_[rna->physSize - 1] = (rna->Nuc_[rna->physSize - 1] << tmpShift) >> tmpShift;
			}//чтобы в случае 0й размерности не вылетал/ тут мб даже и не проверять

			rna->Nuc_[newPhysSize - 1] = rna->Nuc_[newPhysSize - 1] | rna->getMusk(ind, nucl);
			rna->size_= ind + 1;
		}
		else {
			//need to alloc memory																//correct vrode
			//but we need to change old mem, in with index > rna->size, to 0 else
			size_t *tmpArray = new size_t[newPhysSize];
			for (size_t i = 0; i < rna->physSize; i++) {
				tmpArray[i] = rna->Nuc_[i];
			}
//		std::cout << "smt";
			size_t tmpShift = 8 * (sizeof(size_t) - (rna->physSize % sizeof(size_t))) - 2;
			if (rna->physSize != 0) {//чтобы в случае 0й размерности рнк не вылетал
				tmpArray[rna->physSize - 1] = (tmpArray[rna->physSize - 1] << tmpShift) >> tmpShift;
			}
			
			for (size_t i = rna->physSize; i < newPhysSize; i++)
				tmpArray[i] = 0;
			//самому последнему нуклеотиду присвоить принимаемое значение - бубен с маской
			tmpArray[newPhysSize - 1] = tmpArray[newPhysSize - 1] | rna->getMusk(ind, nucl);
			//инд, а не инд+1 пушо в маске счет от нуля. но фактически 1й нуклеотид на 1м месте																				//это 1я физическая позиция
			delete[] rna->Nuc_;
			rna->Nuc_ = tmpArray;
			rna->physSize = newPhysSize;
			rna->size_ = ind + 1;
		}

	}
	return *this;
}

RNK::reference::operator Nucleotide() {
	if (ind > rna->size_)
		return A;
	size_t countNucInType = (4 * sizeof(size_t));
	size_t d = ind % countNucInType;
	size_t c = ind / countNucInType;

	size_t tmp = rna->Nuc_[c]; 
	tmp = tmp << 2 * (countNucInType - 1 - d);
	tmp = tmp >> 2 * (countNucInType - 1);

	Nucleotide n = A;
	switch (tmp)
	{
	case A:
		n = A;
		break;
	case C:
		n = C;
		break;
	case G:
		n = G;
		break;
	case T:
		n = T;
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
RNK::reference RNK::operator[](size_t ind){
	return reference(this, ind);
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
	Nuc_ = new unsigned int[physSize];
	for (size_t i = 0; i < physSize; i++) {
		Nuc_[i] = other.Nuc_[i];
	}
}

RNK RNK::operator !(){
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

//not correct in ind==0
void RNK::trim(size_t lastIndex) {
	if (lastIndex < this->size_) {		//< / <= ?
		this->size_ = lastIndex;
		this->physSize = lastIndex / (4 * sizeof(size_t));
		if (lastIndex % (4 * sizeof(size_t)) != 0)
			this->physSize++;
		unsigned int * NucNew = new unsigned int[this->physSize];
		for (size_t i = 0; i < this->physSize; i++) {
			NucNew[i] = Nuc_[i];
		}
		delete[] Nuc_;
		Nuc_ = NucNew;
	}
}


//may be faled in efficiency factor(KPD)
bool RNK::operator==(const RNK& rna) const {
	RNK tmpRnk(rna);
	bool tmp = this->isComplementary(!tmpRnk);
	return tmp;
}


bool RNK::operator!=(const RNK& rnk) const {
	if (*this == rnk)
		return false;
	else
		return true;
}

//need to check all things in this func
RNK RNK::split(size_t ind) {
	RNK rnkTmp;
	RNK rnkTmpThis(*this);
	for (size_t i = 0; i < this->size_ - ind; i++) {
		rnkTmp[i] = (Nucleotide)rnkTmpThis[i + ind];
	}
	this->trim(ind);
	return rnkTmp;
}
