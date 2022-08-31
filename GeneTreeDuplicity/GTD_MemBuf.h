#ifndef __GTD_MEMBUF__
#define __GTD_MEMBUF__

#include"GTD_TreeEntities.h"
#include"GTD_util.h"

class MemBuf
{
public:
	MemBuf();
	~MemBuf();
	static const unsigned int TAXA_MAX = USHRT_MAX/2;
	static const unsigned int TAXA_LEN_MAX = 100;

	void reset();
	void init();

	CustArrContainer<MulTreeNode> nodeBuf;
	unsigned int nodeBufInd;//index to the next free slot in nodeBuf	

	CustArrContainer<MulTree> treeBuf;
	unsigned int treeBufInd;//index to the next free slot in treeBuf

	char (*taxaBuf)[TAXA_LEN_MAX];
	unsigned short taxaBufInd;//index to the next free slot in taxaBuf, also gives num of uniq taxa
	size_t *taxaLenBuf;//to store string len of each taxa so that at printing, time could be saved

	set<TaxaLabel> taxaIdSet;
	bool *isTaxaValidArr;
	

	//tmp bufs
	unsigned int *uintStack;
	unsigned int *uintStack2;
	unsigned int *uintStack3;
	unsigned int *printStack;
	unsigned short *ushrtStack;	
	unsigned short *ushrtStack2;
	bool *boolStack;
	bool *boolStack2;
	char *chBuf;
};


#endif //__GTD_MEMBUF__
