#include "TreeIndexingUtilities.h"
//const char * dataFileName = "./files/GB159.PB.newick.11.30.07";
const char * mappingFileName = "./files/seqToTaxaMapping.txt";
const char * treeFileName = "./files/tree.db";
const char * invListFileName = "./files/invList.db";
const char * vocabFileName = "./files/vocab.db";
const char * queryTaxaMapFileName = "./files/queryTaxaMap.db";
const char * queryTaxaDirFileName = "./files/queryTaxaDir.db";
//const char * taxaUniqIdArrDumpFileName = "./files/genSeqOrgIdArr.db";
const char * genSeqOrgIdArrFileName = "./files/genSeqOrgIdArr.db";
const char * genTaxaOrgIdFileName = "./files/genTaxaOrgIdArr.db";
const char * genSeqGenTaxaArrFileName = "./files/genSeqGenTaxaArr.db";
const char * treeDirNW = "./files/trees/newick";
const char * treeDirBN = "./files/trees/binary";
const char * querySeqFreqArrG_FileName = "./files/temp/fileCounts/querySeqFreqArrG.db";
const char * invListBaseDir = "./files/invLists";
const char * treeDirMapFileName = "./files/treeDirMap.db";
const char * treeStrIdMapFileName = "./files/treeStrIdMap.db";
const char * clusQueryTaxaCountFileName = "./files/clusTaxaCount.db";
const char * queryTaxaNameMapFileName ="./files/pb.dmp.ti_name.184";
const char * queryGenusNameMapFileName ="./files/pb.dmp.ti_genus.184";
const char * strTaxaDBFileName ="./files/strTaxaMap.db";
const char * strGenusDBFileName ="./files/strGenusMap.db";

ofstream TI_LOG_FILE("./log/TreeIndexing.log",ios::app);
ofstream RESPONSE;

const unsigned int thousand = 1000;
const unsigned int tenThousand = 10000;
const unsigned int hundredThousand = 100000;
const unsigned int oneMillion =1000000;
const unsigned int fiveMillion = 5000000;
const unsigned int maxSeqCount = 7000000;//six million
const unsigned int largePrimeForSeq = 14000029;//smallest prime greater than 14 million
const unsigned int maxTaxaCount = 1000000;//one milllion
const unsigned int largePrimeForTaxa = 2000003;//smallest prime greater than 2 million
const unsigned int largePrimeForTaxaPlusGenus = 3000017;//smallest prime greater than 3 million
const unsigned int tenMillion = 10000000;
const unsigned int fifteenMillion = 15000000;
const unsigned int fiftyMillion = 50000000;
const unsigned int twoFiftyMillion = 250000000;
const unsigned int bufArrSize = twoFiftyMillion;

const unsigned int treeDirCount = 1;
const unsigned int treeSubDirCount = 3;
//const unsigned int invListCount = 10;
vector <ifstream *> vecTreeHandler;
vector <ifstream *> vecInvListHandler;
vector <ifstream *> vecInvListDirHandler;
vector<string> vecTree;
const string phyloTreeConstants::suffixString = "                    ";


const int phyloTreeConstants::maxTreeNAmeSize =20;


