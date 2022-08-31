#ifndef __STRING_HASH__
#define __STRING_HASH__


#include<iostream>
#include<string>
#include <set>
#include <algorithm>
#include<fstream>
#include <cstring>
#include <climits>
#include <cstdlib>

using namespace std;

//a member of the string hash table which actually contains the data
class StrHashElem 
{
public:
	unsigned int key_hashCode;//hash code of the input string.
	char * key_strInd;//ptr in charBuf for the input string.
	unsigned char strLen;//length of the input string
	unsigned int val;//value pvoided by user while insertion.
	unsigned int next;//index of the next element in this bucket
};

//represents a bin in the hash table
class HashBin
{
public:
	HashBin():count(0)
	{}
	unsigned char count;//num of elements
	unsigned int ind;//buf ind
};

class String_Hash
{//has table to store unsigned int s
public:
	String_Hash(unsigned int tblSz,unsigned int charBubSz):bufInd(0),clashCount(0),
		tableSize(tblSz),bufSize(tblSz/2),divFactor(10),maxBinSz(0),binCount(0)
	{
		bin = new HashBin[tableSize];
		buf = new StrHashElem[bufSize];
		charBuf = new char[charBubSz];
		charBufInd = charBuf;
		genPrime();
	}

	unsigned int tableSize; //= largePrimeForSeq 
	//= largePrimeForTaxa
	unsigned int bufSize;//size of the buffer which stores the elements = tableSize/2
	unsigned int divFactor;// = 100 for sequence, 10 for taxa;

	unsigned int clashCount;//number of elements not being the first element in any bin
	unsigned char maxBinSz;//max size of a bin
	unsigned int binCount;//num of bins used

	
	HashBin *bin;//the hashTable bin
	StrHashElem (*buf);//buf to store elements of hash table
	unsigned int bufInd;//points to next available free slot, also a count of number of elements in the hash table
	char * charBuf;//to store all strings
	char * charBufInd;//index for charBuf
	
	static const unsigned char numPrimes = 10;
	unsigned int prArr[numPrimes];//array to store the generated prime numbers
	

	void genPrime();
	void insert(char *keyStr,unsigned int val);//inserts value into hash table
	unsigned int get(const char *keyStr);//returns val for this key
	unsigned int getHashCode(const char *keyStr);
	bool compareStr(const char *str1,const char *str2);
		
	~String_Hash()
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

		if(charBuf != NULL)
		{
			delete []charBuf;
			charBuf = NULL;
		}		
	}
};




#endif //__STRING_HASH__
