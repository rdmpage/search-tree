#ifndef __GTD_UTIL__
#define __GTD_UTIL__

#include"GTD_common.h"

template <class T>
class CustArrContainer{
	unsigned int blkI,vecI;//tmp indices, used in [] operator overloading
public:
	CustArrContainer():vecSz(0)
	{
		arr2d = new T*[1024];		
	}
	static const unsigned int BLOCK_SIZE = 1000000;
	//unsigned int ind;
	T **arr2d;
	unsigned short vecSz;
	//overload [] operator
	T& operator[] (const unsigned int nIndex){
		vecI = nIndex/BLOCK_SIZE;
		blkI = nIndex%BLOCK_SIZE;
		if(vecI < vecSz){
			return arr2d[vecI][blkI];
		}else{
			arr2d[vecSz] = new T[BLOCK_SIZE];			
			++vecSz;	
			CUST_MEM_SZ+= sizeof(T)*BLOCK_SIZE;
			return arr2d[vecI][blkI];
		}
	}

	void reset()
	{
		for(unsigned short i=0; i<vecSz; ++i){
			delete []arr2d[i];
			CUST_MEM_SZ-= sizeof(T)*BLOCK_SIZE;
		}		
		vecSz = 0;
	}

	~CustArrContainer(){
		for(unsigned short i=0; i<vecSz; ++i){
			delete []arr2d[i];
			CUST_MEM_SZ-= sizeof(T)*BLOCK_SIZE;
		}
		delete arr2d;
	}
};

class StatArr
{
public:
	StatArr(){
		taxaLostFractionArr = NULL;
		ipMulTaxaFractionArr = NULL;
		ipNodeCArr = NULL;
		singNodeCArr = NULL;
		txCFreqArr_NdCAnalysis = NULL;
		leafCFreqArr = NULL;
	}

	double *taxaLostFractionArr;
	double *ipMulTaxaFractionArr;
	unsigned int *ipNodeCArr;
	unsigned int *singNodeCArr;
	unsigned short *txCFreqArr_NdCAnalysis;
	unsigned short *leafCFreqArr;
	
	~StatArr(){
		if(NULL != taxaLostFractionArr){
			delete [] taxaLostFractionArr;
		}
		if(NULL != ipMulTaxaFractionArr){
			delete [] ipMulTaxaFractionArr;
		}
		if(NULL != ipNodeCArr){
			delete [] ipNodeCArr;
		}
		if(NULL != singNodeCArr){
			delete [] singNodeCArr;
		}		
		if(NULL != txCFreqArr_NdCAnalysis){
			delete [] txCFreqArr_NdCAnalysis;
		}
		if(NULL != leafCFreqArr){
			delete [] leafCFreqArr;
		}
	}
};

class IpOpStat
{
public:
	unsigned int leavesIp;
	unsigned int internalNodesIp;
	unsigned int informativeEdgesIp;
	unsigned int taxaIp;
	unsigned int leavesOp;
	unsigned int internalNodesOp;
	unsigned int taxaOp;
	unsigned int mulTaxaIp;
	unsigned int mulTaxaOp;

	void reset()
	{
		leavesIp = leavesOp =0;
		internalNodesIp = internalNodesOp =0;
		informativeEdgesIp =0;
		taxaIp = taxaOp =0;
		mulTaxaIp = mulTaxaOp =0;
	}
	IpOpStat()
	{
		reset();
	}
};

class UL_HashTable
{//has table to store unsigned longs
public:
	UL_HashTable(unsigned int tblSz):bin(NULL),buf(NULL),bufInd(0),clashCount(0),maxBinSize(0),
		tableSize(tblSz),bufSize(tblSz/2),divFactor(10)
	{
		bin = new unsigned int[tableSize][2];
		for(unsigned int i=0; i<tableSize; ++i)
		{
			bin[i][0]= UINT_MAX;//it is a pointer (index) to the buf, needs to be initialized to null(UINT_MAX) as its defualt value 
			bin[i][1]= 0;
		}
		buf = new unsigned int[bufSize][3];
		for(unsigned int i=0; i<bufSize; ++i)
		{
			buf[i][0]= buf[i][1]= 0;
			buf[i][2]= UINT_MAX;//it is a pointer (index) to the (next node in)buf, needs to be initialized to null(UINT_MAX) as its defualt value 
		}
		genPrime();
	}

	unsigned int tableSize; //= 12000017 first prime greater than 12 million for sequence, 
	//1000003 = first prime greater than 1 milllion for taxa
	unsigned int bufSize;//size of the buffer which stores the elements = tableSize/2
	unsigned int divFactor;// = 100 for sequence, 10 for taxa;

	unsigned int clashCount;//number of elements not being the first element in any bin
	unsigned int maxBinSize;
	
	unsigned int (*bin)[2];//the hashTable bin,[0] stores index in buf and [1] stores numOfElements in bin
	unsigned int (*buf)[3];//buf to store number of elements
	//[0] has key, [1] has value and [2] points to next node, i.e. an index into buf
	static const unsigned char numPrimes =10; //considering UINT_MAX, we should never need more than 10 primes
	unsigned int prArr[numPrimes];//array to store the generated prime numbers
	unsigned int bufInd;//points to next available free slot, also a count of number of elements in the hash table

	void genPrime();
	void insert(unsigned int key,unsigned val);//inserts value into hash table
	unsigned int get(unsigned int key);//returns val for this key
	unsigned int * getPtr(unsigned int key);//returns ptr to the val for this key

	
	~UL_HashTable()
	{
		if(bin != NULL)
		{
			delete []bin;
			bin = NULL;
		}
		if(buf != NULL)
		{
			delete []buf;
			buf = NULL;
		}
	}
};

#endif //__GTD_UTIL__
