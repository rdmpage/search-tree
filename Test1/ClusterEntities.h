#ifndef __CLUSTER_ENTITIES__
#define __CLUSTER_ENTITIES__

#include"Common.h"


unsigned int NumberOfSetBits(unsigned int i);

//structure to store a node in the PhyloTree structure
class TreeNode
{
public:
	TreeNode():numOfChild(0),childInd(UINT_MAX)
	{}
	void reset() //re-init it for next use
	{
		numOfChild =0;
		childInd =UINT_MAX;		
	}
	unsigned short numOfChild; //num of children a node has; = node id if node is child, happens when childInd = -1
	//when it acts as a node id, it actually gives the index of the cluster's taxa array to which this tree belongs to
	unsigned short bl;//branch lenght of this node
	unsigned int childInd; //index into the GlobalBuf.childBuf from where the children start

	unsigned int sizeAsFile();//return total size to store this node in a file
	void writeToFile(ofstream &file);//write this node to the file
	void readFromFile(ifstream &file);//read this node from the file
};

//structure that stores a PhyloTree within a cluster
class PhyloTree
{
public:
	PhyloTree():numOfNodes(0),nodeInd(UINT_MAX)//intentionally init to -1 i.e max val
	{}
	void reset()//re-init for next use
	{
		numOfNodes=0;
		nodeInd=UINT_MAX;
	}
	unsigned short numOfNodes; //num of nodes a tree has = child nodes+intermediate nodes
	unsigned int nodeInd; //index into the GlobalBuf.nodeBuf from where the nodes start

	unsigned int sizeAsFile();//return total size to store this tree in a file
	void writeToFile(ofstream &file);//write this tree to the file
	void readFromFile(ifstream &file);//read this tree from the file
	void toString(unsigned short &len,unsigned int taxaInd,PrintBuf &printBuf,
		unsigned int (*stack)[4],unsigned char labelFormat,unsigned int avgSupInd);//converts it to a string
	//below is a test function
	void toSubtreeString(unsigned short &len,unsigned int taxaInd,char *buf
		,unsigned int *taxaUniqIdArr,unsigned int *seqIdUniqQueryIdArr,
		unsigned int *queryTaxaUniqIdArr,unsigned int (*stack)[4],unsigned char labelFormat,unsigned int avgSupInd);//test function: converts it to a string

};

//stats related to phylo-diversity value of the cluster
class PDStat
{
public:
	PDStat():mean(0),var(0),mrtPD(0)
	{}
	unsigned int mean;//mean pd values for the component trees
	unsigned int var;//variance of the pd values
	unsigned int mrtPD;//pd of the mrt
};

//structure that stores a cluster
class ClusterInfo
{
public:
	
	ClusterInfo():numOfTaxa(0),taxaInd(0),numOfTrees(0),treeInd(0),avgSupInd(0)
	{}

	void reset()//re-init for next use
	{
		numOfTaxa=0;
		taxaInd=0;
		numOfTrees=0;
		treeInd=0;
		avgSupInd=0;
		clusterAsTree.reset();
	}
	unsigned int clusterId; // evetually points to genUniqId for the cluster
	unsigned short numOfTaxa; //num of taxa this cluster consists of
	unsigned int taxaInd; //index into the GlobalBuf.taxaBuf from where the taxas start
	unsigned short numOfTrees; //num of trees a cluster has
	unsigned int treeInd;//index into the GlobalBuf.treeBuf from where the trees start
	unsigned int avgSupInd;//index into the GlobalBuf.avgSupBuf from where the avg support numbers for the nodes of the MRT start
	PhyloTree clusterAsTree; //to store this cluster as a tree
	PDStat pdStat;//to store stats related to pl values

	unsigned int sizeAsFile();//return total size to store this cluster in a file
	void writeToFile(ofstream &file);//write this cluster to the file
	void readFromFile(ifstream &file, unsigned short maxTreeC = USHRT_MAX);//read this cluster from the file
};


//structure that stores tmp buffers used during parsing of a newick tree
//they are needed to be passed to the parseNewick() to avoid their creation and deletion at every parse of a tree
//all such buffers are grouped together under this structure for ease of passing
class TmpBuff
{
public:
	TmpBuff():tmpStack(NULL),taxaArr(NULL),symArr(NULL),mapSeqIdGenId(largePrimeForSeq)		
	{}
	unsigned short * tmpStack;//stores the path from root to the current element at the top of the stack
	unsigned int * taxaArr;//to store the taxas
	unsigned short *blArr;//to store branch lenghts
	UL_HashTable mapSeqIdGenId;//stores mapping of taxaId to genTaxaId
	unsigned short * symArr;//store the symbols encountered in the in-order traversal
	~TmpBuff()
	{
		//delete the buffers
		if(NULL != tmpStack)
		{
			delete [] tmpStack;
		}
		if(NULL != taxaArr)
		{
			delete [] taxaArr;
		}
		if(NULL != symArr)
		{
			delete [] symArr;
		}
		if(NULL != blArr)
		{
			delete [] blArr;
		}
	}
};

class phyloTreeConstants{
public:
	const static string suffixString ;
	const static int maxTreeNAmeSize;
};

class GlobalBuf
{
public:
	GlobalBuf():childBuf(NULL),childBufInd(0),ll_childBuf(0),nodeBuf(NULL),nodeBufInd(0),ll_nodeBuf(0),
		treeBuf(NULL),treeBufInd(0),ll_treeBuf(0),taxaBuf(NULL),taxaBufInd(0),ll_taxaBuf(0),clusBuf(NULL),clusBufInd(0),ll_clusBuf(0),
		avgSupBuf(NULL),avgSupBufInd(0),ll_avgSupBuf(0)
	{
		childBuf = new unsigned short[CHILDBUFSIZE];
		nodeBuf = new TreeNode[NODEBUFSIZE]; 
		treeBuf = new PhyloTree[TREEBUFSIZE];
		taxaBuf = new unsigned int[TAXABUFSIZE];
		clusBuf = new ClusterInfo[CLUSBUFSIZE];
		avgSupBuf = new unsigned char[AVGSUPBUFSIZE];
	}

	void saveCurrBuf()//set the ll-lower limits of all the bufs so that cur data does not get reset
	{
		ll_childBuf = childBufInd;
		ll_nodeBuf = nodeBufInd;
		ll_treeBuf = treeBufInd;
		ll_taxaBuf = taxaBufInd;
		ll_clusBuf = clusBufInd;
		ll_avgSupBuf = avgSupBufInd;
	}

	void freeSavedBuf()//reset all the ll-lower limits of all the bufs to 0 
	{
		ll_childBuf = 0;
		ll_nodeBuf = 0;
		ll_treeBuf = 0;
		ll_taxaBuf = 0;
		ll_clusBuf = 0;
		ll_avgSupBuf =0;
		clearBuf();
	}

	void clearBuf()//clear all the bufs
	{
		childBufInd = ll_childBuf;
		
		for(unsigned int i=ll_nodeBuf; i<=nodeBufInd; ++i)
		{
			nodeBuf[i].reset();//reset nodes
		}
		nodeBufInd = ll_nodeBuf;

		for(unsigned int i=ll_treeBuf; i<=treeBufInd; ++i)
		{
			treeBuf[i].reset();//reset trees
		}
		treeBufInd = ll_treeBuf;

		taxaBufInd = ll_taxaBuf; 

		for(unsigned int i=ll_clusBuf; i<=clusBufInd; ++i)
		{
			clusBuf[i].reset();//reset clusters
		}
		clusBufInd = ll_clusBuf;
		avgSupBufInd = ll_avgSupBuf;
	}

	void printStatus()//prints status of current memory
	{
		T_LOG<<"ChildBuf: "<<childBufInd<<"("<<childBufInd/(CHILDBUFSIZE/100)<<"%)"<<"/"<<RESET_CHILDBUFSIZE<<"/"<<CHILDBUFSIZE<<"\n";
		T_LOG<<"NodeBuf: "<<nodeBufInd<<"("<<nodeBufInd/(NODEBUFSIZE/100)<<"%)"<<"/"<<RESET_NODEBUFSIZE<<"/"<<NODEBUFSIZE<<"\n";
		T_LOG<<"TreeBuf: "<<treeBufInd<<"("<<treeBufInd/(TREEBUFSIZE/100)<<"%)"<<"/"<<RESET_TREEBUFSIZE<<"/"<<TREEBUFSIZE<<"\n";
		T_LOG<<"TaxaBuf: "<<taxaBufInd<<"("<<taxaBufInd/(TAXABUFSIZE/100)<<"%)"<<"/"<<RESET_TAXABUFSIZE<<"/"<<TAXABUFSIZE<<"\n";
		T_LOG<<"ClusterBuf: "<<clusBufInd<<"("<<clusBufInd/(CLUSBUFSIZE/100)<<"%)"<<"/"<<RESET_CLUSBUFSIZE<<"/"<<CLUSBUFSIZE<<"\n";
		T_LOG<<"AvgSupBuf: "<<avgSupBufInd<<"("<<avgSupBufInd/(AVGSUPBUFSIZE/100)<<"%)"<<"/"<<RESET_AVGSUPBUFSIZE<<"/"<<AVGSUPBUFSIZE<<"\n";
	}

	bool memoryBreached()
	{
		if(childBufInd>BREACH_CHILDBUFSIZE)
		{
			return true;
		}
		
		if(nodeBufInd>BREACH_NODEBUFSIZE)
		{
			return true;
		}
		//if(treeBufInd>RESET_TREEBUFSIZE)
		//{
		//	return true;
		//}
		//if(taxaBufInd>RESET_TAXABUFSIZE)
		//{
		//	return true;
		//}
		//if(clusBufInd>RESET_CLUSBUFSIZE)
		//{
		//	return true;
		//}
		return false;
	}

	void check()//check if respective buffers have crossed threshold usage limit and need to be re-init
	{
		if(childBufInd>RESET_CHILDBUFSIZE)
		{
			childBufInd = ll_childBuf;
		}
		if(nodeBufInd>RESET_NODEBUFSIZE)
		{
			for(unsigned int i=ll_nodeBuf; i<=nodeBufInd; ++i)
			{
				nodeBuf[i].reset();//reset nodes
			}
			nodeBufInd = ll_nodeBuf;
		}
		if(treeBufInd>RESET_TREEBUFSIZE)
		{
			for(unsigned int i=ll_treeBuf; i<=treeBufInd; ++i)
			{
				treeBuf[i].reset();//reset trees
			}
			treeBufInd = ll_treeBuf;
		}
		if(taxaBufInd>RESET_TAXABUFSIZE)
		{
			taxaBufInd = ll_taxaBuf; 
		}
		if(clusBufInd>RESET_CLUSBUFSIZE)
		{
			for(unsigned int i=ll_clusBuf; i<=clusBufInd; ++i)
			{
				clusBuf[i].reset();//reset clusters
			}
			clusBufInd = ll_clusBuf;
		}
		if(avgSupBufInd > RESET_AVGSUPBUFSIZE)
		{
			avgSupBufInd = ll_avgSupBuf;
		}
	}

	unsigned short *childBuf; //buffer to store tree nodes' children's indexes
	unsigned int childBufInd; //index pointing to the next available free slot in childBuf
	const static unsigned int CHILDBUFSIZE = twoHundredMillion; //size of childBuf 
	const static unsigned int RESET_CHILDBUFSIZE = twoHundredMillion - twentyMillion;//size after breaching which, childBuf needs to be re-initialized
	const static unsigned int BREACH_CHILDBUFSIZE = twoHundredMillion - fiveMillion;//size after breaching which, current query needs to be aborted 
	//and re-tried or aborted permanently
	unsigned int ll_childBuf;

	
	TreeNode * nodeBuf;//buffer to store tree nodes
	unsigned int nodeBufInd;//index pointing to the next available free slot in nodeBuf
	const static unsigned int NODEBUFSIZE = twoHundredMillion; //size of nodeBuf 
	const static unsigned int RESET_NODEBUFSIZE = twoHundredMillion - twentyMillion;//size after breaching which, nodeBuf needs to be re-initialized
	const static unsigned int BREACH_NODEBUFSIZE = twoHundredMillion - fiveMillion;//size after breaching which, current query needs to be aborted 
	//and re-tried or aborted permanently
	unsigned int ll_nodeBuf;

	PhyloTree * treeBuf;//buffer to store trees
	unsigned int treeBufInd; //index pointing to the next available free slot in treeBuf
	const static unsigned int TREEBUFSIZE = fiveMillion; //size of treeBuf 
	const static unsigned int RESET_TREEBUFSIZE = fiveMillion - oneMillion;//size after breaching which, treeBuf needs to be re-initialized
	unsigned int ll_treeBuf;

	unsigned int *taxaBuf; //buffer to store a clusters' taxas' indexes
	unsigned int taxaBufInd; //index pointing to the next available free slot in taxaBuf
	const static unsigned int TAXABUFSIZE = tenMillion; //size of taxaBuf 
	const static unsigned int RESET_TAXABUFSIZE = tenMillion - oneMillion;//size after breaching which, taxaBuf needs to be re-initialized
	unsigned int ll_taxaBuf;

	ClusterInfo *clusBuf;//buffer to store the clusters
	unsigned int clusBufInd;//index pointing to the next available free slot in clusBuf
	const static unsigned int CLUSBUFSIZE = hundredThousand; //size of clusBuf 
	const static unsigned int RESET_CLUSBUFSIZE = hundredThousand - tenThousand;//size after breaching which, clusBuf needs to be re-initialized
	unsigned int ll_clusBuf;

	unsigned char *avgSupBuf;//buffer to store a avg support number for each of the majoirty bi-paritions
	unsigned int avgSupBufInd;//index pointing to the next available free slot in avgSupBuf
	const static unsigned int AVGSUPBUFSIZE = tenMillion; //size of avgSupBuf 
	const static unsigned int RESET_AVGSUPBUFSIZE = tenMillion - oneMillion;//size after breaching which, avgSupBuf needs to be re-initialized
	unsigned int ll_avgSupBuf;


	~GlobalBuf()
	{
		if(NULL != childBuf)
		{
			delete [] childBuf;
		}
		if(NULL != nodeBuf)
		{
			delete [] nodeBuf;
		}
		if(NULL != treeBuf)
		{
			delete [] treeBuf;
		}
		if(NULL != clusBuf)
		{
			delete [] clusBuf;
		}
		if(NULL != avgSupBuf)
		{
			delete [] avgSupBuf;
		}
	}
};

class ClusterPptInfo
{
public:
	//address related info
	unsigned short taxaCount;
	unsigned short fileId;		
	unsigned int addr;
	//attribute related info
	unsigned char clusType;//single-loci, biclique, quasi-biclique
	unsigned int phylotaId;//gives a way to map back to phylota DB

	//constants
	const static unsigned char SINGLE_LOCI =0;
	const static unsigned char BICLIQUE =1;
	const static unsigned char QUASI_BICLIQUE =2;

};

//the various array buffers needed for MRT and Subtree are grouped here for convenience
class ArrayBufs
{
public:
	ArrayBufs():stackBuf(NULL),taxaArr(NULL),tmpChildBuf(NULL),taxaArr2(NULL)
		,prArr(NULL),childBuf(NULL),pdBuf(NULL),bitCntArr(NULL)
	{
		stackBuf = new unsigned int [fiftyThousand][4];//needed by genSubtree
		taxaArr= new unsigned short[maxSeqCount];//needed by genSubtree
		tmpChildBuf = new unsigned short[fiftyThousand];//needed by genSubtree
		taxaArr2 = new unsigned short[maxClusOverlap][2];//needed by mrt
		for(unsigned short i=0; i<maxClusOverlap; ++i)
		{//init taxaArr
			taxaArr2[i][0] = (i%64);//[0] entry stores the bit position at which this taxa would be represented at
			taxaArr2[i][1] = i/64;//[1] entry stores the n'th 64-bit space at which this taxa would be stored for a bi-parition representation
		}
		prArr = new unsigned int[maxClusOverlap];//needed by mrt
		childBuf = new unsigned int[maxClusOverlap*maxClusOverlap];//100mb -needed by mrt
		pdBuf = new unsigned int[maxTreesInClus];// needed by mrt-to store the pd values for each of the clusters
		bitCntArr = new unsigned char[USHRT_MAX];//required to compare bipartitions
		for(unsigned int i=0; i<USHRT_MAX; ++i){
			bitCntArr[i] = (unsigned char)NumberOfSetBits(i);
		}
	}
	unsigned int (* stackBuf)[4]; //needed by genSubtree
	unsigned short *taxaArr;//needed by genSubtree
	unsigned short *tmpChildBuf;//needed by genSubtree
	unsigned short (*taxaArr2)[2];//needed by mrt
	unsigned int *prArr;//needed by mrt
	unsigned int *childBuf;//4mb needed by mrt
	unsigned int *pdBuf;// needed by mrt-to store the pd values for each of the clusters
	unsigned char *bitCntArr;//required to compare bipartitions

	~ArrayBufs()
	{
		if(NULL != stackBuf)
		{
			delete []stackBuf;
		}
		if(NULL != taxaArr)
		{
			delete []taxaArr;
		}
		if(NULL != tmpChildBuf)
		{		
			delete []tmpChildBuf;
		}
		if(NULL != taxaArr2)
		{
			delete []taxaArr2;
		}
		if(NULL != prArr)
		{
			delete []prArr;
		}
		if(NULL != childBuf)
		{
			delete []childBuf;
		}
		if(NULL != pdBuf) 
		{
			delete []pdBuf;
		}
		if(NULL != bitCntArr){
			delete []bitCntArr;
		}
	}
};


extern GlobalBuf gGlobalBuf;
//following are global pts to their respectives bufs - for easy ref
extern unsigned short* gChildBuf;
extern TreeNode* gNodeBuf;
extern PhyloTree* gTreeBuf;
extern unsigned int* gTaxaBuf;
extern ClusterInfo* gClusBuf;
extern unsigned char* gAvgSupBuf;
//Parses a newick tree
void parseNewick(string &strTree,ClusterInfo & clusterInfo, string &strClusterId,TmpBuff & tmpBuff);
void parseClusters(ClusterInfo * arrClusterInfo,TmpBuff &tmpBuff);
void generateClusterTrees();
void generateSPRTrees(ClusterInfo &clusterInfo,unsigned int *tmpQueue,unsigned int (*sNodeArr)[5],unsigned int (*tNodeArr)[2]);
void genSubtree(ClusterInfo &sClus,ClusterInfo &tClus,ArrayBufs &arrBufs, bool reRootWrtRefTree);
void reRootAtSmTx(PhyloTree &tree ,ArrayBufs &arrBufs);

//ProcessGeneratedClusters.cpp
void processGeneratedClusters();
void readClusterTrees(ClusterInfo &clus,TmpBuff &tmpBuff,string &bootStrapDir,string &clusFileName,unsigned short *taxaArr);
void parseClusterTreeNewick(string &strTree,PhyloTree * clusTree,TmpBuff & tmpBuff,unsigned short *taxaArr);
#endif
