#include<iostream>
#include<unordered_map>
enum Nucleotide { A, G, C, T };


class RNK {
private:
	friend class reference;
	class reference
	{
		RNK *rna;
		size_t ind;
	public:
		reference(RNK * const ptr, size_t ind);
		reference & operator=(Nucleotide nucleotide);
		operator Nucleotide() const;//just "toType" - operator overload
	};
	size_t size_;
	size_t *Nuc_;
	size_t physSize; // actual physical Size 
public:
	RNK() : Nuc_(nullptr), size_(0), physSize(0) {}
	RNK(size_t s, Nucleotide n);

	RNK(const RNK& other); //copy-constr: memory allocation

	size_t getMusk(size_t i, size_t Nuc) const;
	size_t length() const;

	RNK operator !();
	bool operator==(const RNK& rnk) const;
	bool operator!=(const RNK& rnk) const;
	RNK operator+(const RNK& rnk) const;

	bool isComplementary(const RNK& rnk) const;

	void trim(size_t ind);
	RNK split(size_t ibdex);

	//для нуклеотида  - число значений
	size_t cardinality(Nucleotide value) const;
	//аналогично но сразу для всех типов тритов
	std::unordered_map < Nucleotide, int, std::hash<int> > cardinality() const;

	RNK & operator=(const RNK & rna);

	reference operator[](size_t ind);
	Nucleotide operator[](size_t ind) const;
	friend std::ostream& operator << (std::ostream & out, const RNK& other);
//	Nucleotide TmpGeterNuc(size_t i) const;						//for debag							//for debag
//	size_t getPSize() const{ return this->physSize;}
	~RNK();
};
