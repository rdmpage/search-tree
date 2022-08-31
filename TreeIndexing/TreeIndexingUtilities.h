#ifndef _TREE_INDEXING_UTIL_
#define _TREE_INDEXING_UTIL_

#include "common.h"

#ifndef WIN32
#include "client_socket.h"
#include "server_socket.h"
#include <stdlib.h>
#include <sys/time.h>
#define SOCK_PATH_CLIENT "./../sockets/client"
#define SOCK_PATH_SERVER "./../sockets/server"
#endif

#include"StringHash.h"

using namespace std;
typedef pair<unsigned int,unsigned int> ul_ul;
typedef pair<unsigned short,unsigned int> us_ul;


//extern const char * dataFileName ;
extern const char * mappingFileName;
extern const char * treeFileName;
extern const char * invListFileName;
extern const char * vocabFileName;
extern const char * queryTaxaMapFileName;
extern const char * queryTaxaDirFileName;
//extern const char * taxaUniqIdArrDumpFileName;
extern const char * genSeqOrgIdArrFileName;
//extern const char * taxaIDMapDumpFileName;
extern const char * genTaxaOrgIdFileName;
extern const char * genSeqGenTaxaArrFileName;
extern const char * querySeqFreqArrG_FileName;
extern const char * invListDirFileName;
extern const char * invListBaseDir;
extern const char * treeDirMapFileName;
extern const char * treeStrIdMapFileName;
extern const char * clusQueryTaxaCountFileName;
extern const char * queryTaxaNameMapFileName;
extern const char * queryGenusNameMapFileName;
extern const char * strTaxaDBFileName;
extern const char * strGenusDBFileName;


extern ofstream TI_LOG_FILE;
extern ofstream RESPONSE;

#ifdef DO_LOG
#define T_LOG TI_LOG_FILE
#else
#define T_LOG std::cout
#endif

extern const unsigned int thousand;
extern const unsigned int tenThousand;
extern const unsigned int hundredThousand;
extern const unsigned int oneMillion;
extern const unsigned int fiveMillion;
extern const unsigned int maxSeqCount;
extern const unsigned int largePrimeForSeq;
extern const unsigned int maxTaxaCount;
extern const unsigned int largePrimeForTaxa;
extern const unsigned int largePrimeForTaxaPlusGenus;
extern const unsigned int tenMillion;
extern const unsigned int fifteenMillion;
extern const unsigned int fiftyMillion;
extern const unsigned int twoFiftyMillion;
extern const unsigned int bufArrSize;

extern vector <ifstream *> vecTreeHandler;
extern vector<string> vecTree;
extern vector <ifstream *> vecInvListHandler;
extern vector <ifstream *> vecInvListDirHandler;



extern const unsigned int treeDirCount;
extern const unsigned int treeSubDirCount;
//extern const unsigned int invListCount;
extern const char * treeDirNW;
extern const char * treeDirBN;

extern const unsigned int taxaBlockSize;
extern const unsigned int maxTaxaBlockCount; // needs to be updated with taxaBlockSize

class phyloTreeConstants{
public:
	const static string suffixString ;
	const static int maxTreeNAmeSize;
};

struct TreeAddrInfo
{
unsigned short nodeCount;
unsigned short fileId;
unsigned int addr;
};


struct VocabDB
//vocab for all taxas (both being used in some tree and ones only present in taxaid-seq id mapping file)
{
	unsigned int *vocabBuf;//array to store vocab-
	//- maps taxaId(as index)to ptr in inv list buffer
	unsigned int taxaCount;//# of taxas mapped by this vocab
};

struct InvListDB
//stores invlist for all taxas being used in all the trees
{
	unsigned int *invListBuf;//array to store invList
	//-invList file is constructed here before being written-
	//-also used to hold invList during run time
	unsigned int size;//size of IvListBuf
};

class UI_HashTable
{//has table to store unsigned int s
public:
	UI_HashTable(unsigned int tblSz):bin(NULL),buf(NULL),bufInd(0),clashCount(0),
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

	unsigned int tableSize; //= largePrimeForSeq
	//= largePrimeForTaxa
	unsigned int bufSize;//size of the buffer which stores the elements = tableSize/2
	unsigned int divFactor;// = 100 for sequence, 10 for taxa;

	unsigned int clashCount;//number of elements not being the first element in any bin

	
	unsigned int (*bin)[2];//the hashTable bin,[0] stores index in buf and [1] stores numOfElements in bin
	unsigned int (*buf)[3];//buf to store number of elements
	//[0] has key, [1] has value and [2] points to next node, i.e. an index into buf
	static const unsigned char numPrimes = 10;//considering UINT_MAX, we shall not need more than 10 prime numbers
	unsigned int prArr[numPrimes];//array to store the generated prime numbers
	unsigned int bufInd;//points to next available free slot, also a count of number of elements in the hash table

	void genPrime();
	void insert(unsigned int key,unsigned int val);//inserts value into hash table
	unsigned int get(unsigned int key);//returns val for this key
	unsigned int * getPtr(unsigned int key);//returns ptr to the val for this key

	
	~UI_HashTable()
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

struct AllMaps;//forward declaration

class CmdLineOpts
{
public:
	CmdLineOpts():taxaIds(NULL)		
	{
		taxaIds = new unsigned int[maxSeqCount];
		taxonFileOpt = "-f";
		treeFileOpt = "-t";
		tableFileOpt = "-s";
		taxaIdsOpt = "-L";
		strTaxaIdOpt = "-N";
		strTaxaFileOpt ="-n";
		overLapOpt = "-k";
		useSeqIdOpt ="-q";
		labelFormatOpt ="-z";
		cumQueryTaxaOpt ="-u";
		seqIdOverlapOpt ="-o";
		maxCCOpt ="-c";
		maxBSTOpt ="-b";
		reset();
	}

	void reset()//resets all values for consequtive use
	{
		taxonFileName ="";
		treeFileName ="";
		tableFileName ="";
		numTaxaId =0;
		overLap = 3;
		labelFormat =0;

		useSeqId = false;
		cumQueryTaxa = false;
		seqIdOverlap = false;
		maxCC = UINT_MAX;
		maxBST_Count = USHRT_MAX;
	}
	string taxonFileName;//file from where the taxa ids should be read
	string taxonFileOpt;
	string treeFileName;//output file for nexus tree file format
	string treeFileOpt;
	string tableFileName;//output file with tab delimited table with summary statistics
	string tableFileOpt;
	
	bool useSeqId;//if true, then ids point to sequence ids, by default it is false i.e taxa id
	string useSeqIdOpt;
	unsigned int *taxaIds;//to store query taxa ids or sequence ids as the case may be
	string taxaIdsOpt;
	string strTaxaIdOpt;
	string strTaxaFileOpt;
	unsigned int numTaxaId;//num of taxa ids
	unsigned int overLap;//min overlap with query taxa/seq ids.
	string overLapOpt;
	unsigned char labelFormat;//whether taxon labels should be in ti#### or gi#### or both ti####_gi####
	//labelFormat='T' for Taxa, ='S' for sequence, =0 for both(both is also the default value), ='N' if name is to be printed
	string labelFormatOpt;
	bool cumQueryTaxa;//if true then the first line in the table file contains cumalative query taxas
	string cumQueryTaxaOpt;
	bool seqIdOverlap;//if true then over lap is done with seq id else with taxa id, default its taxa id
	string seqIdOverlapOpt;
	string maxCCOpt;//max number of valid clusters to be reported
	unsigned int maxCC;

	string maxBSTOpt;
	unsigned short maxBST_Count;//max number of bootstrap trees to read

	


	bool parse(int argc, string *argv,  AllMaps &allMaps);
	
	~CmdLineOpts()
	{
		if(NULL != taxaIds)
		{
			delete [] taxaIds;
		}
	}

};

struct AllMaps
//combines all information required in memory to process input
{
	AllMaps():seqHashTable(largePrimeForSeq),queryTaxaDir(largePrimeForTaxa),strTaxaHashTb(NULL)
		//largePrimeForSeq = for sequence, 
		//largePrimeForTaxa = for taxa
	{}
	unsigned int * queryTaxMapBuf;//buf storing mapping of queryTaxaId-sequenceId 
	VocabDB vocabDB;
	InvListDB invListDB;
	TreeAddrInfo * arrAllTree;
	unsigned int * seqIdUniqQueryIdArr;//stores mapping from sequientially generated unique sequence ids to -
	//-the corresponding generated unique query taxa ids
		
	UI_HashTable seqHashTable;//stores mapping from orignal seq ids to gen seq ids
	void initSeqHashTable();
	bool *dupQueryTaxaFlag;//to check for dup query taxa
	
	UI_HashTable queryTaxaDir;//stores mapping from orgQueryTaxaId to the corresponding position in queryTaxaDirBuff
	void initQueryTaxaDirHashTable();

	String_Hash *strTaxaHashTb;
	unsigned int *strTaxaGenTaxaMapBuf;//buf storing mapping of strTaxaId-genQueryTaxaIds
	unsigned int *strTaxaOrgQueryTaxaMapIndArr;//stores ind for each strTaxa in strTaxaGenTaxaMapBuf
	unsigned char *strTaxaIdFreqArr;//for each strTaxaId this stores # of genTaxaIds this strTaxaId is mapped to
	void initStrTaxaHashTb();

	
	CmdLineOpts cmdLineOpts;

	~AllMaps()
	{
		delete [] arrAllTree;
		delete [] queryTaxMapBuf;
		delete [] vocabDB.vocabBuf;
		delete [] invListDB.invListBuf;
		delete [] dupQueryTaxaFlag;
		delete strTaxaHashTb;		
		delete [] strTaxaGenTaxaMapBuf;
		delete [] strTaxaOrgQueryTaxaMapIndArr;
		delete [] strTaxaIdFreqArr;
	}
};

//class to process input queries
class ProcessInput
{
	unsigned int *treeIdResBuf;
	//buf to store all the res trees which containing any of the taxa, includes duplicated trees too
	
	unsigned int *treeResFreqSeqIdBuf;
	//buf to store uniq trees which containing any of the sequence id,and their frequency
	//it is a treeId-freq mappin with index representing the treeIds
	
	unsigned int (*treeResFreqTaxaIdBuf)[2];
	//buf to store uniq trees which containing any of the taxa id,and their frequency
	//it is a treeId-freq mappin with index representing the treeIds
	//[0] stores the freq, [1] value is used as a flag that

	unsigned int *treeUniqIdResBuf;
	//buf to store the uniq treeIds

	unsigned int (*queryIdBuf)[2];
	//stores all the query taxa ids, and their corresponding starting point in taxaIdBuf

	
	string *argvBuf;//to sotre the parsed command line arguments
	char strBuf[257];//used as the commong buffer to send data to sockets
public:
	int respCount;//stores genId for response and result files
	AllMaps allMaps;


#ifndef WIN32
	ClientSocket sockClient;
	ServerSocket sockServer;
#endif

public:
	ProcessInput()
#ifndef WIN32
		:sockClient(SOCK_PATH_SERVER,T_LOG),
		sockServer(SOCK_PATH_CLIENT,T_LOG)
#endif
	{
		treeIdResBuf = new unsigned int[fiftyMillion];
		treeResFreqSeqIdBuf = new unsigned int[tenMillion];
		treeResFreqTaxaIdBuf = new unsigned int[tenMillion][2];
		treeUniqIdResBuf = new unsigned int[oneMillion];
		queryIdBuf = new unsigned int[oneMillion][2];

		for(unsigned int tc=0;tc<tenMillion;++tc)
		{//initialize all freqs to 0
			treeResFreqSeqIdBuf[tc]=0;
			treeResFreqTaxaIdBuf[tc][0]=0;
			treeResFreqTaxaIdBuf[tc][1]=0;//init flag to zero
		}

		argvBuf = new string [10100];//for max 10k taxa and for 100 other command arguments
		respCount = 1;//let query start with 1, it saves a lot of problem vs. starting with 0

		
	}

	~ProcessInput()
	{
		if(NULL != treeIdResBuf)
		{
			delete [] treeIdResBuf;
		}
		if(NULL != treeResFreqSeqIdBuf)
		{
			delete [] treeResFreqSeqIdBuf;
		}
		if(NULL != treeResFreqTaxaIdBuf)
		{
			delete [] treeResFreqTaxaIdBuf;
		}
		if(NULL != treeUniqIdResBuf)
		{
			delete [] treeUniqIdResBuf;
		}
		if(NULL != queryIdBuf)
		{
			delete [] queryIdBuf;
		}
		if(NULL != argvBuf)
		{
			delete []argvBuf;
		}
	}

	void getQuery(int &argc);
	int process(int argc);
	void sendQueryNumToClient(int num);//sends the queryNum to client, so that client can write its query in the corresponding file
	int sendQueryNumToServer(int num);//sends the queryNum to sever, server reads the corresponding intermediate result as query
	void signalResultReadyToClient(int num);//signals client that the result is ready to be picked up
	
};

void initialize_part1();
void readTreeDir();
void parseTree(string &strTree, unsigned int *buffArr, unsigned int &tbufArrCount, string &treeId);
void parseTaxa(string &strTaxa, unsigned int &queryId, unsigned int &taxaId);
void loadTrees();
void initializeQueryTaxaMap();
void initQueryTaxaNameMap();//to init the hash table for string taxa ids 
void initGenusTaxaNameMap();//to init the hash table for string genus ids to taxa ids


void initialize_part2();
void createVocabInitInvListMap(VocabDB &vocabDB,InvListDB &invListDB);
//void crateInvListsUpdateDir(VocabDB &vocabDB,InvListDB &invListDB);
void createInvList(VocabDB &vocabDB,InvListDB &invListDB);

void loadMaps(AllMaps &allMaps);
void initializeTreeFileHandlers();
void initializeInvListFileHandlers();
void finalizeFileHandlers();
void loadQueryTaxaBuf(unsigned int * &queryTaxMapBuf);
void loadVocab(VocabDB &vocabDB);
void loadInvList(InvListDB &invListDB);
void loadMapTreeDir(TreeAddrInfo* &treeArr);
void loadSeqIdUniqQueryTaxaId(unsigned int *&seqIdUniqQueryIdArr);





//For Testing purposes

void hashSeqId();



#endif

