#include"GTD_MemBuf.h"

MemBuf::MemBuf()
{
	init();
}

void MemBuf::init()
{
	reset();

	taxaBuf = new char[TAXA_MAX][TAXA_LEN_MAX];
	FIXED_MEM_SZ+= sizeof(char)*TAXA_MAX*TAXA_LEN_MAX;

	taxaLenBuf = new size_t[TAXA_MAX];
	FIXED_MEM_SZ+= sizeof(size_t)*TAXA_MAX;

	isTaxaValidArr = new bool[TAXA_MAX];
	FIXED_MEM_SZ+= sizeof(bool)*TAXA_MAX;

	//tmp bufs
	uintStack = new unsigned int[2*TAXA_MAX];
	FIXED_MEM_SZ+= sizeof(unsigned int)*2*TAXA_MAX;

	uintStack2 = new unsigned int[2*TAXA_MAX];
	FIXED_MEM_SZ+= sizeof(unsigned int)*2*TAXA_MAX;

	uintStack3 = new unsigned int[2*TAXA_MAX];
	FIXED_MEM_SZ+= sizeof(unsigned int)*2*TAXA_MAX;

	printStack = new unsigned int[2*TAXA_MAX];
	FIXED_MEM_SZ+= sizeof(unsigned int)*2*TAXA_MAX;

	ushrtStack = new unsigned short[2*TAXA_MAX];
	FIXED_MEM_SZ+= sizeof(unsigned short)*2*TAXA_MAX;

	ushrtStack2 = new unsigned short[2*TAXA_MAX];
	FIXED_MEM_SZ+= sizeof(unsigned short)*2*TAXA_MAX;
	
	boolStack = new bool[2*TAXA_MAX];
	FIXED_MEM_SZ+= sizeof(bool)*2*TAXA_MAX;

	boolStack2 = new bool[2*TAXA_MAX];
	FIXED_MEM_SZ+= sizeof(bool)*2*TAXA_MAX;
	
	chBuf = new char[2*TAXA_MAX*TAXA_LEN_MAX];
	FIXED_MEM_SZ+= sizeof(char)*2*TAXA_MAX*TAXA_LEN_MAX;

}

void MemBuf::reset()
{
	for(unsigned int i=0; i<nodeBufInd; ++i){
		nodeBuf[i].reset();
	}
	nodeBufInd =0;
	treeBufInd =0;

	taxaBufInd =0;
	taxaIdSet.clear();
}

MemBuf::~MemBuf()
{
	delete []taxaBuf;
	delete []taxaLenBuf;
	delete []isTaxaValidArr;

	//tmp bufs
	delete []uintStack;
	delete []uintStack2;
	delete []uintStack3;
	delete []printStack;
	delete []ushrtStack;
	delete []ushrtStack2;
	delete []boolStack;
	delete []boolStack2;
	delete []chBuf;
}
