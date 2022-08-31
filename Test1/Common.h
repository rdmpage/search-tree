#ifndef __COMMON_H__
#define __COMMON_H__

#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<list>
#include <sstream>
#include<bitset>
#include<set>
#include<limits.h>
#include <cstdlib>

#ifdef WIN32
#pragma warning(default : 4242)
#endif


using namespace std;

const unsigned int thousand = 1000;
const unsigned int twoThousand = 2000;
const unsigned int fiveThousand = 5000;
const unsigned int tenThousand = 10000;
const unsigned int fiftyThousand = 50000;
const unsigned int hundredThousand = 100000;
const unsigned int twoHundredFiftyThousand = 250000;
const unsigned int fiveHundredThousand = 500000;
const unsigned int oneMillion =1000000;
const unsigned int fiveMillion = 5000000;
const unsigned int maxSeqCount = 7000000;//seven million
const unsigned int largePrimeForSeq = 14000029;//smallest prime greater than 14 million
const unsigned int maxTaxaCount = 2000000;
const unsigned int largePrimeForTaxa = 2000003;//smallest prime greater than 2 million
const unsigned int tenMillion = 10000000;
const unsigned int twentyMillion = 20000000;
const unsigned int fifteenMillion = 15000000;
const unsigned int fiftyMillion = 50000000;
const unsigned int hundredMillion = 100000000;
const unsigned int twoHundredMillion = 200000000;
const unsigned int twoFiftyMillion = 250000000;
const unsigned int maxFileSize = 1800000000;//1800 mb
//const unsigned int maxFileSize = 500000000;//500 mb
//const unsigned int numOfClusFiles =2;
const unsigned int maxClusOverlap = thousand;
const unsigned int twiceMaxClusOverlap = 2*maxClusOverlap;
const unsigned int maxTreesInClus = thousand;

extern ofstream T_LOG_FILE;
// all #defines should come here, for the ease of reference

#ifdef DO_LOG
#define T_LOG T_LOG_FILE
#else
#define T_LOG std::cout
#endif


//#define PRINT_SUP_FALSE //decides if sup % should be printed in the tree output


//file names
extern const char * genSeqOrgIdArrFileName;
extern const char * genTaxaOrgIdFileName;
extern const char * genSeqGenTaxaArrFileName;
extern const char * treeStrIdMapFileName;
extern const char * clusQueryTaxaCountFileName;
extern const char * strTaxaDBFileName;
extern const char * queryTaxaNameMapFileName;

//utility classes

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

	unsigned int tableSize; //= largePrimeForSeq for sequence
	//= largePrimeForTaxa for taxa
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

//stores memory related structures printing a tree, used by other operations too
struct PrintBuf
{
	PrintBuf():nodeIdName_ht(largePrimeForTaxa)
{}
	//strTaxa related structures
	unsigned char *strTaxaLenArr;//to store lenght of each string taxa, NOTE: it is indexed based on the genTaxaId to which this strTaxaId maps to
	char * strTaxaBuf;//to store all the string taxa
	unsigned int *genTaxaIndStrTaxaBufArr;//maps genTaxaId to the index of its corresponding strTaxa in strTaxaBuf

	//rest
	unsigned int *taxaUniqIdArr;
	char * charBuf;//used to print the tree in this buf
	unsigned int *queryTaxaUniqIdArr ;//maps uniqQueryTaxaIds to queryTaxaIds
	unsigned int *seqIdUniqQueryIdArr;//array which maps uniqSeqIds to their corresponding uniqQueryTaxaIds 	

	char *nodeIdNameBuf;//to store the node names
	UL_HashTable nodeIdName_ht;//maps nodeId to ptr in nodeIdNameBuf corresponding to node name



	~PrintBuf(){
		//strTaxa related structures
		delete strTaxaLenArr;
		delete strTaxaBuf;
		delete genTaxaIndStrTaxaBufArr;

		//rest
		delete taxaUniqIdArr;
		delete charBuf;
		delete queryTaxaUniqIdArr;
		delete seqIdUniqQueryIdArr;
		delete nodeIdNameBuf;
	}
};

//BELOW ARE TEST FUNCTIONS
void genAllMRT(unsigned short minLfSz, unsigned short maxLfSz);
void testConvertSeqIdToTaxonNames();

#endif //__COMMON_H__
