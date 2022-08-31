#ifndef __GTD_TEST__
#define __GTD_TEST__
#include"GTD_MemBuf.h"
#include"GTD_TreeProcessor.h"

class EvalST_MulTrees{
public:
	void evalAllClusters();
	void processTreeFile(char *clusFileName, char *clusPrefix, ofstream &resFile, UL_HashTable &mapOrgSeqOrgTaxa);
	void loadAllMaps(UL_HashTable &mapOrgSeqOrgTaxa);
	void parseMulNewickST(char *strTree,unsigned int treeInd, UL_HashTable &mapOrgSeqOrgTaxa);

	//specific to Feb 02, 2012 analysis on best trees
	void evalAllBestTrees();
	void processBestTreeFile(char *treeFileName, char *baseDir, char *dirPrefix, ofstream &resFile, UL_HashTable &mapOrgSeqOrgTaxa,
		StatArr &statArr);
	void parseMulNewickBestTree(char *strTree,unsigned int treeInd, UL_HashTable &mapOrgSeqOrgTaxa);
	void initMapsBestTree(UL_HashTable &mapOrgSeqOrgTaxa);
};





#endif //__GTD_TEST__
