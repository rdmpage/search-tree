#ifndef __MRT_HEADER__
#define __MRT_HEADER__
#include"ClusterEntities.h"

//const defs
const unsigned int PRIME_S = 1000003;//prime number above million



//class defs



class BiPart
{
public:

	static unsigned char len;//num of 64 bit spaces needed to store a bi-partition. why 64 - this pc has a 64 bit processor and runs on 64 bit OS :=)
	unsigned short card;//cardinality of this bi-partition
	unsigned int hashCode;//hash code of this bi-paritition
	
	unsigned long long part[1+maxClusOverlap/64];//to store sub-tree part of this bi-partition
	unsigned int binNdInd;//points to index of the bin-node in the bin-nodeArr in which this bi-part is stored-
	//-this helps speed up, as the bi-part are not required to be fetched again from the hash tb in the Step 2 

	static unsigned long long lastOfset;//the last 64 bit space with bits representing extra taxa set to 0 - needed by reverse()
	

	void reset();//resets all 64-bit spaces for next use
	bool comp(BiPart &biPart);//compares two bi-parititions
	void addPart(BiPart &biPart); //add parition from biPart
	unsigned short diff(BiPart &biPart, ArrayBufs &arrBufs);//computes the symmetric difference between the biparitions
	void reverse();//reverses the biaprition i.e. now the bits represent the other side of the bipartition


};

//Represents Node in a Bin of the hash table
//Nodes are stored in the bin inform of linked list
class BinNd
{
public:

	unsigned short freq;//store freq for biPart
	unsigned int ind;//ptr to node of this bi-parition in the consensus tree, points to index in array of type MajTreeNodePtr
	BiPart *biPart;//ptr to the bi-parition which is stored in this bin
	unsigned int nxt;//ptr to the next node in the linked list for the Bin in which this node is contained
};

//Rrepresetns a bin or a bucket for the hash table
class Bin
{
public:
	Bin():ind(UINT_MAX)
	{}
	unsigned short len; //num of elements in this bin
	unsigned int ind;//index into some BinNd array where its first element start
	unsigned int add(BiPart *biPart,BinNd * binNdArr,unsigned int &binNdArrInd);//adds a Bi-partition to this bucket
	unsigned int getBiPar(BiPart *biPart, BinNd * binNdArr);

};

//when a maj bi-parition is added to the majTree as a node, then the corresponding node in the bin of the hash table
//contains a ptr(of type MajTreeNodePtr)  to this node in the majTree.  
class MajTreeNodePtr
{
public:
	unsigned short supNum;//to store support num of this maj bipartition
	unsigned short numChild;//to store numOf chilren it has - used as a temp variable
	unsigned int ind;//ptr to node in majTree
	unsigned int chInd;//ptr to its children in some childBuf
	unsigned int parInd;//ptr to parent of this bi-parition in the consensus tree
	unsigned int posInParChBuf;//position of this child in the childBuf of its parent, such arrangement is needed for parent update-
	//- operation while constructing the majority tree
	unsigned short bl;//to store the cumalative bl value 
	
	
};


//func defs
unsigned char genMRT(ClusterInfo &clus,BiPart * biPartArr,Bin *hashTb,BinNd * binNdArr
			,ArrayBufs &arrBufs,MajTreeNodePtr * ndPtrBuf,unsigned short &maxSupportNum, unsigned short &numOfUniqMajBiPart);
void genPrime(unsigned int *prArr,unsigned int numPr);

bool rootMRT(PhyloTree &refTree, PhyloTree &mrtTree, BiPart * biPartArr, ArrayBufs &arrBufs, unsigned int clusAvgSupInd);
void initBipartitions(PhyloTree &tree, BiPart * biPartArr, unsigned int biPartInd, ArrayBufs &arrBufs);
void initRefTreeRootBipartition(PhyloTree &refTree, BiPart * biPartArr, unsigned int biPartInd, ArrayBufs &arrBufs);


#endif //__MRT_HEADER__
