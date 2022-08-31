#ifndef __T1_TEST_H__
#define __T1_TEST_H__

#include"Common.h"
#include"ClusterEntities.h"
#include"MRT.h"

void testSubtree();
void parseTestNewick(char *strTree,PhyloTree &currTree, ArrayBufs &arrBufs);
void testPrint(PhyloTree &currTree, ArrayBufs &arrBufs, char *prBuf, unsigned int clusTxInd);

void testRootMRT();


#endif //__T1_TEST_H__
