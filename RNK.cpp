#include"RNK.h"

RNK::reference::reference(RNK * const ptr, size_t ind):rna(ptr), ind(ind) {}

//need to check may be:
RNK::reference & RNK::reference::operator=(Nucleotide nucl) {
	if (ind < rna->size_) {
		//error 
		std::cout << "ERROR. RNA don't can be changed in a middle\n";
	}
	else {						//we can alloc memory if we need. - correct vrode
		//realloc
		//put in data
		//add data
		size_t newPhysSize = (ind + 1) / (4 * sizeof(size_t));
		if ((ind + 1) % (4 * sizeof(size_t)) != 0)
			newPhysSize++;

		if (newPhysSize != rna->physSize) {					//we need to alloc mem
			size_t * tmpArr = new size_t[newPhysSize];
			for (size_t i = 0; i < rna->physSize; i++)
				tmpArr[i] = rna->Nuc_[i];

			size_t tmpShift = 4 * (sizeof(size_t) - (rna->physSize % sizeof(size_t)));
			size_t tmp = tmpArr[rna->physSize - 1] << (tmpShift * 2 - 2);
			tmpArr[rna->physSize - 1] = tmp >> (tmpShift * 2 - 2);


			for (size_t i = rna->physSize; i < newPhysSize; i++)
				tmpArr[i] = 0;
			
			size_t Musk = rna->getMusk(ind, nucl);
			std::cout << Musk << " - Musk\n";
			tmpArr[newPhysSize - 1] = tmpArr[newPhysSize - 1] | Musk;
			/*add data in last size_t*/
			//movie with musk
			std::cout << "Musc from alloced mem:" << Musk <<"\n";

			rna->physSize = newPhysSize;
			rna->size_ = ind;
			delete[] rna->Nuc_;
			rna->Nuc_ = tmpArr;
		}
		else {									//we don't need to alloc memory - correct vrode.
			/*add data in last size_t*/
			//movies with musks
			size_t surplus = ind % (4 * sizeof(size_t));
			for (size_t i = rna->size_; i < ind; i++) {
				size_t clearMusk = rna->getMusk(i, T);
				clearMusk = ~clearMusk;
				std::cout << "clearMusk from \"dont need musk\" " << clearMusk;
				std::cout << "\n";
				rna->Nuc_[newPhysSize - 1] = rna->Nuc_[newPhysSize - 1] & clearMusk;
			}
			size_t Musk = rna->getMusk(ind, nucl);
			rna->Nuc_[newPhysSize - 1] = rna->Nuc_[newPhysSize - 1] | Musk;
			rna->size_ = ind;
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

size_t RNK::getSize() const {
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

//correct vrode
void RNK::trim(size_t lastIndex) {
	if (lastIndex < this->size_) {		//< / <= ?
//		size_t size = 0;
		this->size_ = lastIndex - 1;
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
