#ifndef __QUERY_ENTITIES__
#define __QUERY_ENTITIES__
#include"ClusterEntities.h"

#ifndef WIN32
#include "./../TreeIndexing/server_socket.h"
#define SOCK_PATH_SERVER "./../sockets/server"
#endif

extern ofstream RESPONSE;

void processQuery();
void handleMemoryBreached(bool &preCleanUp);

class QueryInfo
{
public:
	QueryInfo():numOfQTaxa(0),numOfClus(0),sClusI(UINT_MAX),tClusI(UINT_MAX),useSeqId(false),labelFormat(0),
		cumQueryTaxa(false)
	{
		qTaxaArr = new unsigned int[maxSeqCount];
		clusIdArr = new unsigned int[tenThousand];
		seqIdOverlapArr = new unsigned int[tenThousand];
		taxaIdOverlapArr = new unsigned int[tenThousand];
		scoreArr = new unsigned short[tenThousand];
		maxSupportNumArr = new unsigned short[tenThousand];
		uniqNonTrivBpCArr = new unsigned short[tenThousand];
		treeFileName="";
		tableFileName="";
	}
	unsigned int numOfQTaxa;//sequence id count
	unsigned int *qTaxaArr;//sequence ids

	unsigned int numOfClus;//# of clusters in the result
	unsigned int *clusIdArr;//array storing above clusters
	unsigned int *seqIdOverlapArr;//over lap of each cluster with seq id
	unsigned int *taxaIdOverlapArr;//over lap of each cluster with query taxa id
	unsigned short *scoreArr;//tree score based on normalized avg support % of the maj-bi paritions
	unsigned short *maxSupportNumArr;//max support number for any maj-bipartition
	unsigned short *uniqNonTrivBpCArr;//number of unique non-trivial bipartitions in the MRT, --
	//-- the (clus.clusterAsTree.numOfNodes - clus.numOfTaxa -1) cannot be trusted as in reRooting subtrees before the MRT, -
	//- in reRooting the MRT, and in degree two check of the MRT (to come soon), numOfNodes is getting messed up.

	unsigned int sClusI;//index in gClusBuf to store the source clusters
	unsigned int tClusI;//index in gClusBuf to store the target clusters

	//command line option start
	bool useSeqId;//whether query was based on sequence ids
	bool cumQueryTaxa;//whether cumQueryTaxaShould be reported
	string treeFileName;//name of treeFile where trees in nexus format would be written
	string tableFileName;//name of tableFile where results in tabular form would be written
	unsigned char labelFormat;//whether taxon labels should be in ti#### or gi#### or both ti####_gi####
	//labelFormat='T' for Taxa, ='S' for sequence, =0 for both(both is also the default value)
	unsigned short maxBST_Count;//max bootstrap trees to be read
	//command line option end

	~QueryInfo()
	{
		if(NULL != qTaxaArr)
		{
			delete [] qTaxaArr;
		}
		if(NULL != clusIdArr)
		{
			delete [] clusIdArr;
		}
		if(NULL != seqIdOverlapArr)
		{
			delete [] seqIdOverlapArr;
		}
		if(NULL != taxaIdOverlapArr)
		{
			delete [] taxaIdOverlapArr;
		}
		if(NULL != scoreArr)
		{
			delete [] scoreArr;
		}
		if(NULL != maxSupportNumArr)
		{
			delete [] maxSupportNumArr;
		}
		if(NULL != uniqNonTrivBpCArr)
		{
			delete [] uniqNonTrivBpCArr;
		}
		
	}
};

unsigned int loadClusterMap(ClusterPptInfo *clusAddrArr);
void loadClusStrIdMap(char (*clusAddrArr)[21],unsigned int (*clusId_ti_giArr)[2]);
void loadUniqQueryTaxaIdArr(unsigned int *queryTaxaUniqIdArr);
void loadSeqIdUniqQueryTaxaId(unsigned int *seqIdUniqQueryIdArr);
void loadResFile(int i,QueryInfo &query);
void initTestQueryInfo(QueryInfo &queryInfo);
void loadGenSeqIdToOrgSeqIdArr(unsigned int *taxaUniqIdArr);
void loadClusQueryTaxaCountArr(unsigned int *clusQueryTaxaCountArr);
void loadStrTaxaInfo(PrintBuf &pb);
void loadNodeIdNameInfo(PrintBuf &pb);
#ifdef WIN32
void getQueryFileNumber(int &queryFileNum);
void sendResFileNum(int resCount);
#else
void getQueryFileNumber(int &queryFileNum,ServerSocket &serverSock);
void sendResFileNum(int resCount,ServerSocket &serverSock);
#endif

#endif //__QUERY_ENTITIES__
