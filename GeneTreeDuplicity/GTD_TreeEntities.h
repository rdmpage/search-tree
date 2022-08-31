#ifndef __GTD_TREE_ENTITIES__
#define __GTD_TREE_ENTITIES__

#include"GTD_common.h"

class MulTreeNode{
public:
	const static int INFO_EQUAL = 0;
	const static int INFO_GREATER = 1;
	const static int INFO_LESSER = -1;
	const static int INFO_NA = INT_MAX;//comparison not applicable

	const static int DO_COLLAPSE = 1;
	MulTreeNode()
	{
		reset();
	}
	void reset()
	{
		firstChInd = UINT_MAX;
		Sdn = Sup = 0;
		tmpVal =0;
	}

	unsigned short numChild;// will also act as taxaId if firstChInd=UINT_MAX
	unsigned short Sdn;
	unsigned short Sup;
	unsigned short tmpVal;//to store some temporary flag attributes
	unsigned int firstChInd;
	unsigned int lastChInd;// will also act as gi_id if firstChInd=UINT_MAX, i.e., a leaf node
	unsigned int nextSibInd;
	unsigned int prevSibInd;
	//unsigned int parentInd;
	
	
	void addChild(unsigned int parentI, unsigned int childInd);
	void collapseChild(unsigned int chInd);
	void deleteChild(unsigned int chInd);
	void markForCollapse();
	void chkN_Collapse(unsigned int chInd);
	bool isInformative();
	int compInfWidChld(unsigned int ndI);
	int compInfWidSib(unsigned int ndI);
};

class MulTree{
public:
	unsigned int rootInd;
	unsigned short numNodes;	
};

//wrapper for taxa for convenient insertion/deletion in the set
class TaxaLabel{
public:
	unsigned short ind;//points to its position in taxaBuf
};

bool  operator <(const TaxaLabel & lhs, const TaxaLabel & rhs);
bool  operator >(const TaxaLabel & lhs, const TaxaLabel & rhs);
bool  operator ==(const TaxaLabel & lhs, const TaxaLabel & rhs);



#endif
