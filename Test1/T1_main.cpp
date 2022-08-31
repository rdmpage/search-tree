
#include"Common.h"
#include"ClusterEntities.h"
#include "QueryEntities.h"
#include<math.h>

//below includes needed for test code
#include"MRT.h"
#include"T1_test.h"

ofstream T_LOG_FILE("./log/Test1.log",ios::app);

//file names
const char * genSeqOrgIdArrFileName = "./files/genSeqOrgIdArr.db";
const char * genTaxaOrgIdFileName = "./files/genTaxaOrgIdArr.db";
const char * genSeqGenTaxaArrFileName = "./files/genSeqGenTaxaArr.db";
const char * treeStrIdMapFileName = "./files/treeStrIdMap.db";
const char * clusQueryTaxaCountFileName = "./files/clusTaxaCount.db";
const char * strTaxaDBFileName ="./files/strTaxaMap.db";
const char * queryTaxaNameMapFileName ="./files/pb.dmp.ti_name.184";



void bloatTree();

struct TreeAddrInfo
{
unsigned short nodeCount;
unsigned short fileId;
unsigned int addr;
};




int main(int argc, char * argv[])
{
	if(!T_LOG_FILE.is_open())
	{
		std::cout<<"FATAL ERROR! T_LOG_FILE could not be opened!\n";
		exit(-1);
	}

	//below is a test function
	//testSubtree();
	//genAllMRT(150,150);
	//return 0;

	//processGeneratedClusters();
	//return 0;
	
	
	processQuery();
	return 0;
}



void bloatTree()
{
	unsigned int numOfTrees = 1000000;
	unsigned int numOfFiles = 10;
	string fileName = "./files/GB159.PB.newick.11.30.07";
	ifstream dataFile(fileName.c_str());
	string line ;
	string tmpString;
	size_t pos1 ;

	for(unsigned int j=0; j<numOfFiles; j++)
	{
		
		char resFileName[100] ;
		sprintf(resFileName,"./files/trees/Tree%d",j);
		ofstream resFile_create(resFileName);
		if(! resFile_create.is_open())
		{
			T_LOG<<"\nOuch!";
			return;
		}
		resFile_create.close();
		ofstream resFile(resFileName, ios::out);

		if(dataFile.is_open())
		{
			for(unsigned int i=0; i<numOfTrees; ++i)
			{

				if(dataFile.eof())
				{
					dataFile.clear();
					dataFile.close();
					dataFile.open(fileName.c_str());
					if(!dataFile.is_open())
					{
						T_LOG<<"ouch again";
					}
				}
				char num[10] ;
				sprintf(num,"%d",i+(j*1000000));

				getline (dataFile,line);
				if(line.length() == 0)
				{
					--i;
					continue;
				}
				pos1 = line.find_first_of(" \t");
				resFile<<num<<line.substr(pos1)<<"\n";
			}
		}
		resFile.close();
	}
	dataFile.close();

}


unsigned int UL_HashTable::get(unsigned int key)
{
	unsigned int rem,ind=0,count=0,quo=key;
	while(0 != quo)
	{//calculate key
		rem = quo%divFactor;
		ind = (rem*prArr[count++] +ind)%tableSize;
		quo = quo/divFactor;
	}

	//find the value in the chained list
	for(ind=bin[ind][0]; ind!= UINT_MAX; ind=buf[ind][2])
	{
		if(buf[ind][0] == key)
		{
			return buf[ind][1];
		}
	}
	return ind;//this will return UINT_MAX i.e. when the key was not found

}

unsigned int * UL_HashTable::getPtr(unsigned int key)
{
	unsigned int rem,ind=0,count=0,quo=key;
	while(0 != quo)
	{//calculate key
		rem = quo%divFactor;
		ind = (rem*prArr[count++] +ind)%tableSize;
		quo = quo/divFactor;
	}
	//find the value in the chained list
	for(ind=bin[ind][0]; ind!= UINT_MAX; ind=buf[ind][2])
	{
		if(buf[ind][0] == key)
		{
			return &buf[ind][1];
		}
	}
	return NULL;//this will return UINT_MAX i.e. when the key was not found
}

void UL_HashTable::insert(unsigned int key,unsigned val)
{
	unsigned int rem,ind=0,count=0,quo=key;

	while(0 != quo)
	{
		rem = quo%divFactor;
		ind = (rem*prArr[count++] +ind)%tableSize;
		quo = quo/divFactor;
	}
	if(0 == bin[ind][1])
	{//no element in this bin till now
		bin[ind][0] = bufInd;
		++bin[ind][1];
		buf[bufInd][0]=key;//store key
		buf[bufInd][1]=val;//store val
		buf[bufInd][2]=UINT_MAX;//making ptr to next node null
		++bufInd;
	}
	else
	{//find next free slot
		++clashCount;
		++bin[ind][1];
		if(maxBinSize<bin[ind][1])
		{
			maxBinSize=bin[ind][1];
		}
		for(ind=bin[ind][0]; buf[ind][2] != UINT_MAX;)
		{
			ind=buf[ind][2];
		}

		buf[ind][2] = bufInd;
		buf[bufInd][0] = key;
		buf[bufInd][1] = val;
		buf[bufInd][2] = UINT_MAX;//making prt to next node null
		++bufInd;
	}
}

void UL_HashTable::genPrime()
{
	unsigned int prC=0;
	bool validPr;
	unsigned int curPr,rem,quo;
	set<unsigned int> uniqPr;
	pair<set<unsigned int>::iterator,bool> ret;
	for(unsigned i=0; i<numPrimes; ++i)
	{
		validPr = false;
		while(!validPr)
		{
			curPr = rand()+1;
			quo = tableSize/curPr;
			rem = tableSize%curPr;
			quo = rand()%(quo+1);
			rem = rand()%(rem+1);
			curPr = curPr*quo + rem;
			if(curPr < tableSize)
			{
				ret = uniqPr.insert(curPr);
				validPr = ret.second;
			}

		}
		prArr[i] = curPr;
	}
}


void genAllMRT(unsigned short minLfSz, unsigned short maxLfSz)
{
	char fileName[1024];
	char opDir[] = "E:/MiningFrequentTrees/codeFromZhang/sampleTrees/MRTsFromTOL_Data/";
	char clusListFileName[] = "clusList150.txt";
	sprintf(fileName,"%s/%s",opDir,clusListFileName);
	
	ofstream clusListFile(fileName);
	if(!clusListFile.is_open()){
		cout<<"Error! Could not open File:" <<fileName<<endl;
		exit(-1);
	}

	unsigned int numOfClusFiles=0;
	while(true)
	{//count number of clusterFiles
		char clusterFileName[257] ;
		sprintf(clusterFileName,"./files/clusters/binary/Cluster%d.db",numOfClusFiles);
		ifstream clusFile(clusterFileName,ios::in | ios::binary);
		if(!clusFile.is_open())
		{
			T_LOG<<"\n Number of cluster files:"<<numOfClusFiles<<endl;
			break;
		}
		++numOfClusFiles;
		clusFile.close();
	}
	//init file handlers for cluster files
	ifstream *clusFileArr = new ifstream[numOfClusFiles];
	for(unsigned int i=0; i<numOfClusFiles; ++i)
	{
		char clusterFileName[257] ;
		sprintf(clusterFileName,"./files/clusters/binary/Cluster%d.db",i);
		clusFileArr[i].open(clusterFileName,ios::in | ios::binary);
	}

	ArrayBufs arrBufs;//group diff buffers needed for sub tree and MRT
	for(unsigned int i=0;i<maxSeqCount; ++i)
	{
		arrBufs.taxaArr[i] = 0;//init all values to true
	}
	genPrime(arrBufs.prArr,thousand);

	ClusterPptInfo *clusAddrArr = new ClusterPptInfo[tenMillion];
	unsigned int numClus = loadClusterMap(clusAddrArr);
		//assuming ti_gi are present
	unsigned int (*clusId_ti_giArr)[2] = new unsigned int[tenMillion][2];
	char (*clusStrIdArr)[21] = new char [tenMillion][21];
	loadClusStrIdMap(clusStrIdArr,clusId_ti_giArr);

	BiPart *biPartArr = new BiPart[twoThousand*(maxTreesInClus+1)];//15mb max taxa =1000, hence max bi-parition=200, max num of trees=1000 (<1001)
	Bin *hashTb = new Bin[PRIME_S]();//8mb
	BinNd * binNdArr = new BinNd[thousand*(maxTreesInClus+1)];
	MajTreeNodePtr * ndPtrBuf = new MajTreeNodePtr[twoThousand];

	PrintBuf PB;//some of the used arrays are grouped here
	PB.queryTaxaUniqIdArr = new unsigned int [maxTaxaCount];//maps uniqQueryTaxaIds to queryTaxaIds
	PB.seqIdUniqQueryIdArr = new unsigned int[maxSeqCount];//array which maps uniqSeqIds to their corresponding uniqQueryTaxaIds 
	loadUniqQueryTaxaIdArr(PB.queryTaxaUniqIdArr);
	loadSeqIdUniqQueryTaxaId(PB.seqIdUniqQueryIdArr);

	//loads strTaxa and genTaxaId-strTaxa mapping
	loadStrTaxaInfo(PB);//also allocates required memory to relevant structures
	loadNodeIdNameInfo(PB);//loads nodeId-nodeName mapping in ht, clusterName has a component nodeId which is needed-
	//- in string form

	//for result output
	PB.taxaUniqIdArr = new unsigned int[maxSeqCount];
	PB.charBuf = new char[oneMillion];//used in toString function of phyloTree
	loadGenSeqIdToOrgSeqIdArr(PB.taxaUniqIdArr);




	unsigned short fileId;
//	unsigned short maxSupNum;
	unsigned int clusInd;
//	unsigned short len;
	for(unsigned int i=0; i<numClus; ++i){//load clusters
		if(clusAddrArr[i].taxaCount < minLfSz || clusAddrArr[i].taxaCount > maxLfSz){
			continue;
		}

		clusInd = gGlobalBuf.clusBufInd;
		ClusterInfo &clus =gClusBuf[clusInd];
		fileId = clusAddrArr[i].fileId;//file id where this cluster is stored
		clus.clusterId = i;//update clusterId in clusterInfo
		clusFileArr[fileId].seekg(clusAddrArr[i].addr);//seek to the address in the file
		clus.readFromFile(clusFileArr[fileId]);//read the cluster from the file
		//genMRT(clus,biPartArr,hashTb,binNdArr,arrBufs,ndPtrBuf,maxSupNum);
		cout<<i<<endl;
		clusListFile<<clusAddrArr[clus.clusterId].taxaCount<<"\t"<<
			clusId_ti_giArr[clus.clusterId][0]<<"_"<<clusId_ti_giArr[clus.clusterId][1]<<"\t";
		//clus.clusterAsTree.toString(len,clus.taxaInd,PB,arrBufs.stackBuf,'S',clus.avgSupInd);
		//clusListFile<<PB.charBuf<<std::endl;		
		clusListFile<<std::endl;

		++gGlobalBuf.clusBufInd;
		gGlobalBuf.clearBuf();
	}



	for(unsigned int i=0; i<numOfClusFiles; ++i)
	{
		clusFileArr[i].close();
	}
	delete []clusFileArr;

	delete []ndPtrBuf;
	delete []binNdArr;
	delete []hashTb;
	delete []biPartArr;
	delete []clusAddrArr;
	clusListFile.close();
}
