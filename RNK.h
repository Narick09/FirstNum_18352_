#include<iostream>
#include<unordered_map>
enum Nucleotide { A, G, C, T };


class RNK {
private:
	//   sizeof(size_t)
	friend class reference;
	class reference
	{
		RNK *rna;				//why not "&"?
		size_t ind;
	public:
		reference(RNK * const ptr, size_t ind);
		reference & operator=(Nucleotide nucleotide);
		operator Nucleotide();//just "toType" - operator overload
		//Can we write operator overload in "RNK", but not in "ref."? What will changed after that?
		//It's will be problem with rna[i] = Nucl;
	};
	size_t size_;
	size_t *Nuc_;
	size_t physSize; // actual physSize 
public:
	RNK() : Nuc_(nullptr), size_(0), physSize(0){}
	RNK(size_t s, Nucleotide n);

	RNK(const RNK& other); //copy-constr: memory allocation

	size_t getMusk(size_t i, size_t Nuc) const;
	size_t getSize() const;

	RNK operator !();
	bool operator==(const RNK& rna) const;
	bool operator!=(const RNK& rna) const;
	RNK operator+(const RNK& rna);


	bool isComplementary(const RNK& rnk) const;

	void trim(size_t ind);
//------------------------------------------------------------------------
//	Nucleotide TmpGeterNuc(size_t i) const;
	RNK & operator=(const RNK & rna);

	reference operator[](size_t ind);

	~RNK();



};
