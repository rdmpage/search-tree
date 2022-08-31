#ifndef __GTD_TREE_PROCESSOR__
#define __GTD_TREE_PROCESSOR__

#include"GTD_common.h"

class TreeProcessor
{
public:
	void parseMulNewick(char *strTree,unsigned int treeInd);	
	void parseMulNewickBQred(char *strTree,unsigned int treeInd);	
	void computeAllSuvValues(unsigned int treeInd);
	void computeNodeSuvValue(unsigned int treeInd, unsigned int curNdInd, unsigned int *parSt, unsigned short parStC);
	unsigned short collapseNonInformativeEdges(unsigned int treeInd);
	void removeNonParticipatingLeaves(unsigned int treeInd);
	void print(unsigned int treeInd);
	void printTableBQred(ofstream &resFile, unsigned int treeInd);
	void collapseDupEdges(unsigned int treeInd);
	void removeDupLeaves(unsigned int treeInd);
	void remAllRedundantLeaves(unsigned int treeInd);
	void removeRedundantLeaf(unsigned int treeInd, unsigned short taxaId);
	bool restrictToSingTaxa(unsigned int srcTreeInd, unsigned int destTreeInd);
	void subtree(unsigned int srcTreeI, unsigned int destTreeI, bool *taxaArr, bool unRooted);
	void degTwoCheck(unsigned int treeI);
	void computeStats(unsigned int treeInd);
	
};


#endif //__GTD_TREE_PROCESSOR__
