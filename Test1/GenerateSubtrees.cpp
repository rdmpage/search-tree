#include"ClusterEntities.h"


void genSubtree(ClusterInfo &sClus, ClusterInfo &tClus, ArrayBufs &arrBufs, bool reRootWrtRefTree)
{
	unsigned short stackC;//counter for stackBuf
	unsigned short tmpChildBufC;//counter for tmpChildBuf
	unsigned int chI,ndI,chNdI,taxaI,j,k,l;//temp values-used repeatedly in for loops
	
	
	PhyloTree *sTree,*tTree;

	tClus.numOfTrees = reRootWrtRefTree? (sClus.numOfTrees -1):sClus.numOfTrees;//first tree is only for re-rooting

	tClus.treeInd = gGlobalBuf.treeBufInd;
	gGlobalBuf.treeBufInd+=tClus.numOfTrees;

	tClus.taxaInd = gGlobalBuf.taxaBufInd;
	unsigned short tc=0;//taxa count
	
	bool reRootAtSmTaxa;
	
	for(unsigned short i=0; i<sClus.numOfTaxa; ++i)
	{
		k = gTaxaBuf[sClus.taxaInd+i];
		if(USHRT_MAX != arrBufs.taxaArr[k])
		{			
			gTaxaBuf[tClus.taxaInd+tc] = k;//copy taxaId
			arrBufs.taxaArr[k] = tc; //make it point to the taxa index for target cluster
			++tc;//increment taxa count for next taxa
		}
	}
	tClus.numOfTaxa = (unsigned short)tc;
	gGlobalBuf.taxaBufInd+=tClus.numOfTaxa;//make it point to next free slot

	tClus.pdStat.mean =0; //init the mean val
	unsigned short ti;
	for(unsigned short si=0;si<sClus.numOfTrees;++si)
	{
		ti = reRootWrtRefTree?si-1:si;
		if(si != 0 || !reRootWrtRefTree){			
			arrBufs.pdBuf[ti] =0;//init the pd value to be 0
		}
		
		stackC=0;
		tmpChildBufC=0;

		sTree = &gTreeBuf[sClus.treeInd+si];
		if(si != 0 || !reRootWrtRefTree){
			tTree = &gTreeBuf[tClus.treeInd+ti];
		}else{
			tTree = &sClus.clusterAsTree;
		}

		++tTree->numOfNodes;//reserve first place for root as tTree->numOfNodes is being used as a counter		
		tTree->nodeInd=gGlobalBuf.nodeBufInd;//gGlobalBuf.nodeBufInd is not being init to next free slot, numOfNodes is not known yet


		//remember stackBuf always stores data pertaining to source tree except [stackC][3] which points to tmpChildBufC and 
		//tmpChildBufC stores gChildbuf values for the target tree
		//as root does not have any parent, it needs special treatement
		arrBufs.stackBuf[stackC][0] = sTree->nodeInd;//init to root
		arrBufs.stackBuf[stackC][1] = 0; //visited cildren  count
		arrBufs.stackBuf[stackC][2] = 0; //valid children count
		arrBufs.stackBuf[stackC][3] = tmpChildBufC;//to store valid children
		tmpChildBufC+= gNodeBuf[sTree->nodeInd].numOfChild;//max possible valid children


		//put the first child of the root on stack
		chI = gNodeBuf[sTree->nodeInd].childInd;//child index of root node
		arrBufs.stackBuf[++stackC][0] = sTree->nodeInd+gChildBuf[chI];//node index of first child

		while(USHRT_MAX != stackC)
		{
			ndI = arrBufs.stackBuf[stackC][0];//try not to alter ndI in the blocks below, as it is useful to have a ref to the node on the stack-
			//-esp while modifying the code, error tends to happen
			++arrBufs.stackBuf[stackC-1][1];//increment visited children count of this nodes parent
			if(UINT_MAX == gNodeBuf[ndI].childInd)
			{//a leaf node
				taxaI=gNodeBuf[ndI].numOfChild + sClus.taxaInd;//index of taxa in gTaxaBuf
				if(USHRT_MAX != arrBufs.taxaArr[gTaxaBuf[taxaI]])			
				{//i.e. a query taxa
					chI= arrBufs.stackBuf[stackC-1][3]+ arrBufs.stackBuf[stackC-1][2]++;//index in tmpChildBuf for this new child-
					//-also increment valid visited children node count
					arrBufs.tmpChildBuf[chI] = tTree->numOfNodes++;//update child node index-
					//-for this node in valid children arr for its parent
					chNdI=arrBufs.tmpChildBuf[chI]+tTree->nodeInd;//node index for this child in gNodeBuf
					gNodeBuf[chNdI].childInd = UINT_MAX;//a leaf node
					gNodeBuf[chNdI].numOfChild = arrBufs.taxaArr[gTaxaBuf[taxaI]];//update taxaId ind
					gNodeBuf[chNdI].bl = gNodeBuf[ndI].bl;//init the branch length for the new node
				}
				//this node has been visited, now find the next node in tree in pre-order i.e. node, left to right children
				--stackC;//go to parent of curr node
				while(true)
				{
					ndI = arrBufs.stackBuf[stackC][0];//index of the parent, ndI is being modified here since parent is the current node under consideration
					if(arrBufs.stackBuf[stackC][1] == gNodeBuf[ndI].numOfChild)
					{//all children of this node are visited
						if(1 == arrBufs.stackBuf[stackC][2])
						{//only one valid child
							if(0 == stackC)
							{//root reached, tree processing done, root would be processed separtely
								--stackC; //i.e.stackC= USHRT_MAX, termination condition for the outer loop
								break;
							}
							//this node will not exist in tree, instead its only valid child would be passed to its parent to handle
							chI= arrBufs.stackBuf[stackC-1][3]+ arrBufs.stackBuf[stackC-1][2]++;//index in tmpChildBuf for this child for its new parent-
							//-also increment valid visited children node count for its new parent
							arrBufs.tmpChildBuf[chI] =  arrBufs.tmpChildBuf[arrBufs.stackBuf[stackC][3]];//byspassing: this->parent->child = this->child
						}
						else if(arrBufs.stackBuf[stackC][2] > 1)
						{//more than one valid children, this node has to exist
							if(0 == stackC)
							{//root reached, tree processing done, root would be processed separtely
								--stackC; //i.e.stackC= USHRT_MAX, termination condition for the outer loop
								break;
							}
							chI= arrBufs.stackBuf[stackC-1][3]+ arrBufs.stackBuf[stackC-1][2]++;//index in tmpChildBuf for this new child-
							//-also increment valid visited children node count for its parent
							arrBufs.tmpChildBuf[chI] = tTree->numOfNodes++;//update child node index-
							//-for this node in valid children arr of its parent
							chNdI=arrBufs.tmpChildBuf[chI]+tTree->nodeInd;//node index for this child in gNodeBuf
							//this is only place where a node for the target tree is being created, of-course the root is created separately at the end
							gNodeBuf[chNdI].numOfChild = (unsigned short)arrBufs.stackBuf[stackC][2];//num of children = valid child count
							j= gNodeBuf[chNdI].childInd = gGlobalBuf.childBufInd;
							gNodeBuf[chNdI].bl= gNodeBuf[ndI].bl; //copy the branch length
							gGlobalBuf.childBufInd+= gNodeBuf[chNdI].numOfChild;//point to next free slot
							k= arrBufs.stackBuf[stackC][3];//starting point in tmpChildBuf
							if(si != 0 || !reRootWrtRefTree){
								for(l=0;l<gNodeBuf[chNdI].numOfChild;++l){
									gChildBuf[j+l]=arrBufs.tmpChildBuf[k+l];//copy each child
									arrBufs.pdBuf[ti]+= gNodeBuf[tTree->nodeInd + gChildBuf[j+l]].bl 
										- gNodeBuf[chNdI].bl; //for pd value, relative diff of the child w.r.t to its parent
								}
							}else{//reRootWrtRefTree: pd value calculation not needed for the first tree, which is used only for re-rooting
								for(l=0;l<gNodeBuf[chNdI].numOfChild;++l){
									gChildBuf[j+l]=arrBufs.tmpChildBuf[k+l];//copy each child									
								}
							}
						}
						--stackC;//this node visited,go to his parent continuing search for the next node in tree in pre-order
					}
					else
					{//some children are still left, put next child on top of stack
						chI = gNodeBuf[ndI].childInd + arrBufs.stackBuf[stackC][1];
						arrBufs.stackBuf[++stackC][0] = sTree->nodeInd+gChildBuf[chI];//node index of next child
						//no need to init this child, as it would be done when it picked from the stack, for the first time
						break;
					}
				}
			}
			else
			{//an intermediate node picked from the stack for the first time, init the ndoe and put its first child on the stack
				arrBufs.stackBuf[stackC][1] = 0; //visited cildren  count
				arrBufs.stackBuf[stackC][2] = 0;//valid children count
				arrBufs.stackBuf[stackC][3] = tmpChildBufC;//to store valid children
				tmpChildBufC+= gNodeBuf[ndI].numOfChild;//max possible valid children
				//now put its first child on stack
				++stackC;//increment stack
				arrBufs.stackBuf[stackC][0]= gChildBuf[gNodeBuf[ndI].childInd] + sTree->nodeInd;//put the first child on stack
			}
		}


		//process root node of target tree
		gGlobalBuf.nodeBufInd+= tTree->numOfNodes;//init to next free slot as now numOfNodes are calculated
		//if(sTree->numOfNodes < tTree->numOfNodes)
		//{
		//	T_LOG<<"\n Error!";
		//}		

		reRootAtSmTaxa = true;

		if(arrBufs.stackBuf[0][2]==2 && (si != 0 || !reRootWrtRefTree)){//root has only two children, not ok for unrooted tree
			//first check if one of the children is smallest taxa, then no need to re-root. Also the first tree in the cluster --
			//-- is already rooted in case rooting is required, thus no need to reroot it at the samllest taxa
			j= tTree->nodeInd + arrBufs.tmpChildBuf[arrBufs.stackBuf[0][3]];//first child index
			if(UINT_MAX == gNodeBuf[j].childInd && //first child is a leaf node --
				gNodeBuf[j].numOfChild == 0){//--and it is the smallest taxa
					reRootAtSmTaxa = false;
			}

			k= tTree->nodeInd + arrBufs.tmpChildBuf[arrBufs.stackBuf[0][3]+1];//second child index
			if(UINT_MAX == gNodeBuf[k].childInd && //second child is a leaf node --
				gNodeBuf[k].numOfChild == 0){//--and it is the smallest taxa
					reRootAtSmTaxa = false;
			}

			gNodeBuf[tTree->nodeInd].bl =0;//original root node exists for bl puposes, i.e. in the next statements -
			//- its children's bl will be calculated as cum bl of its children - bl of root, and bl of root = bl of the original root -
			//- the above holds true whether or not one of its children will collapse during re-rootin

			arrBufs.pdBuf[ti]+= gNodeBuf[j].bl 
				- gNodeBuf[tTree->nodeInd].bl; //for pd value, relative diff of the child w.r.t to its parent
			arrBufs.pdBuf[ti]+= gNodeBuf[k].bl 
				- gNodeBuf[tTree->nodeInd].bl; //for pd value, relative diff of the child w.r.t to its parent

			
			if(reRootAtSmTaxa){//i.e., this tree is going to be re-rooted, 
				//need to collapse the degree two root node				
				if(UINT_MAX != gNodeBuf[j].childInd){//if first child is not a leaf node
					//init root node with all children of the first child + second child as the last child					
					chI = gNodeBuf[tTree->nodeInd].childInd = gGlobalBuf.childBufInd;
					gGlobalBuf.childBufInd+= gNodeBuf[tTree->nodeInd].numOfChild = gNodeBuf[j].numOfChild +1;
					chNdI = gNodeBuf[j].childInd;
					for(l=0;l<gNodeBuf[j].numOfChild;++l){
						gChildBuf[chI + l] = gChildBuf[chNdI +l];//copy each child						
					}
					gChildBuf[chI + l] = k - tTree->nodeInd;//make second child as the last child
				}else{//if first child is a leaf node then second has to be non-leaf
					//init root node with all children of the second child + first child as the last child
					chI = gNodeBuf[tTree->nodeInd].childInd = gGlobalBuf.childBufInd;
					gGlobalBuf.childBufInd+= gNodeBuf[tTree->nodeInd].numOfChild = gNodeBuf[k].numOfChild +1;
					chNdI = gNodeBuf[k].childInd;
					for(l=0;l<gNodeBuf[k].numOfChild;++l){
						gChildBuf[chI + l] = gChildBuf[chNdI +l];//copy each child						
					}
					gChildBuf[chI + l] = j - tTree->nodeInd;//make first child as the last child
				}				
			}else{//init the root node
				gNodeBuf[tTree->nodeInd].bl =0;//root node exists as it is
				chI = gNodeBuf[tTree->nodeInd].childInd = gGlobalBuf.childBufInd;
				gGlobalBuf.childBufInd+= gNodeBuf[tTree->nodeInd].numOfChild = 2;
				gChildBuf[chI] = j - tTree->nodeInd;//add first child
				gChildBuf[chI + 1] = k - tTree->nodeInd;//add second child
			}
		}else if(arrBufs.stackBuf[0][2]==1){//i.e. root has only one valid child
			//in this case copy the the only child node in the position of root node
			ndI = tTree->nodeInd+arrBufs.tmpChildBuf[arrBufs.stackBuf[0][3]];//index to the only valid child
			//bl of the root node here is not 0 and should not be set to 0 either as it will change the bl values of all its children -
			//- as bl values are calculated based on diff of cum bl values
			gNodeBuf[tTree->nodeInd] = gNodeBuf[ndI];//copy node
			--tTree->numOfNodes;
			gNodeBuf[ndI].reset();
			--gGlobalBuf.nodeBufInd;//we can free the only valid child node becuase it was the -
			//- last node to be allocated space in nodeBuf becuse the way we have traversed the tree
			
			if(2 == gNodeBuf[tTree->nodeInd].numOfChild && (si != 0 || !reRootWrtRefTree)){
				//now we need to check if this new root has only two children and if it needs to be collapsed
				chI = gNodeBuf[tTree->nodeInd].childInd;
				j = tTree->nodeInd + gChildBuf[chI];
				if(UINT_MAX == gNodeBuf[j].childInd && //first child is a leaf
					0 == gNodeBuf[j].numOfChild){ //and it is the smallest taxa)
						reRootAtSmTaxa = false;
				}

				k = tTree->nodeInd + gChildBuf[chI+1];
				if(UINT_MAX == gNodeBuf[k].childInd && //second child is a leaf
					0 == gNodeBuf[k].numOfChild){ //and it is the smallest taxa)
						reRootAtSmTaxa = false;
				}

				if(reRootAtSmTaxa){//this is going to be reRooted as smallest taxa, need -
					//- to collapse degree two node before that
					if(UINT_MAX != gNodeBuf[j].childInd){//first child is not a leaf
						//init root node with all children of the first child + second child as the last child					
						chI = gNodeBuf[tTree->nodeInd].childInd = gGlobalBuf.childBufInd;
						gGlobalBuf.childBufInd+= gNodeBuf[tTree->nodeInd].numOfChild = gNodeBuf[j].numOfChild +1;
						chNdI = gNodeBuf[j].childInd;
						for(l=0;l<gNodeBuf[j].numOfChild;++l){
							gChildBuf[chI + l] = gChildBuf[chNdI +l];//copy each child						
						}
						gChildBuf[chI + l] = k - tTree->nodeInd;//make second child as the last child
					}else{//if first child is a leaf node then second has to be non-leaf
						//init root node with all children of the second child + first child as the last child
						chI = gNodeBuf[tTree->nodeInd].childInd = gGlobalBuf.childBufInd;
						gGlobalBuf.childBufInd+= gNodeBuf[tTree->nodeInd].numOfChild = gNodeBuf[k].numOfChild +1;
						chNdI = gNodeBuf[k].childInd;
						for(l=0;l<gNodeBuf[k].numOfChild;++l){
							gChildBuf[chI + l] = gChildBuf[chNdI +l];//copy each child						
						}
						gChildBuf[chI + l] = j - tTree->nodeInd;//make first child as the last child
					}				
				}
			}
		}else{
			ndI=tTree->nodeInd;
			gNodeBuf[ndI].bl =0;//original root exists as it is
			gNodeBuf[ndI].numOfChild = (unsigned short)arrBufs.stackBuf[0][2];//num of valid children
			j=gNodeBuf[ndI].childInd = gGlobalBuf.childBufInd;//init childInd for root node
			gGlobalBuf.childBufInd+= gNodeBuf[ndI].numOfChild;//point to next free slot
			k= arrBufs.stackBuf[0][3];//starting point in tmpChildBuf
			if(si != 0 || !reRootWrtRefTree){
				for(l=0;l<gNodeBuf[ndI].numOfChild;++l){
					gChildBuf[j+l]=arrBufs.tmpChildBuf[k+l];//copy each child
					arrBufs.pdBuf[ti]+= gNodeBuf[tTree->nodeInd + gChildBuf[j+l]].bl 
						- gNodeBuf[tTree->nodeInd].bl; //for pd value, relative diff of the child w.r.t to its parent
				}
			}else{//for reRootWrtRefTree - pd value calculation not needed for the first tree, which is used only for re-rooting
				for(l=0;l<gNodeBuf[ndI].numOfChild;++l){
					gChildBuf[j+l]=arrBufs.tmpChildBuf[k+l];//copy each child
				}
			}
		}


		if(reRootAtSmTaxa){
			if(si != 0 || !reRootWrtRefTree){
				reRootAtSmTx(*tTree,arrBufs);			
			}
		}

		if(si != 0 || !reRootWrtRefTree){//not req when reRootWrtRefTree is true and si==0
			tClus.pdStat.mean+=arrBufs.pdBuf[ti];//add the pd value of the current tree to sum, mean is currently acting as sum
		}
		//subtree found, thank God
	}

	
	//cal PD stats - start
	float tmpDl,varPD,meanPD;//needed to round off and used as temp vars

	//calculate mean
	meanPD = (float)tClus.pdStat.mean;
	meanPD = meanPD/tClus.numOfTrees;
	tClus.pdStat.mean = (unsigned int)meanPD;
	if(meanPD - tClus.pdStat.mean >= .5)
	{//round off
		++tClus.pdStat.mean;
	}

	//calculate var of pd values
	tClus.pdStat.var=0;	
	varPD=0;
	for(unsigned short i=0; i<tClus.numOfTrees; ++i)
	{
		tmpDl= arrBufs.pdBuf[i]-meanPD;//meanPD contains the mean value here
		varPD+= tmpDl*tmpDl;
	}
	varPD/=tClus.numOfTrees;
	tClus.pdStat.var = (unsigned int)varPD;
	if(varPD - tClus.pdStat.var >= .5)
	{//round off
		++tClus.pdStat.var;
	}

	//cal PD stats - end

}



void reRootAtSmTx(PhyloTree &tree ,ArrayBufs &arrBufs)
{
	unsigned short stackC;//counter for stackBuf
	unsigned int chI,ndI,j,k;//temp values-used repeatedly in for loops



	stackC=0;	
	//as root does not have any parent, it needs special treatement
	arrBufs.stackBuf[stackC][0] = tree.nodeInd;//init to root
	arrBufs.stackBuf[stackC][1] = 0; //visited child count	

	//put the first child of the root on stack
	chI = gNodeBuf[tree.nodeInd].childInd;//child index of root node
	arrBufs.stackBuf[++stackC][0] = tree.nodeInd+gChildBuf[chI];//node index of first child

	while(true)
	{
		ndI = arrBufs.stackBuf[stackC][0];//try not to alter ndI in the blocks below, as it is useful to have a ref to the node on the stack-
		//-esp while modifying the code, error tends to happen	
		++arrBufs.stackBuf[stackC-1][1];//increment visited child count of its parent
		if(UINT_MAX == gNodeBuf[ndI].childInd)
		{//a leaf node
			if(0 == gNodeBuf[ndI].numOfChild){//smallest taxa
				break;
			}			
			//this node has been visited, now find the next node in tree in pre-order i.e. node, left to right children
			--stackC;//go to parent of curr node
			while(true)
			{
				ndI = arrBufs.stackBuf[stackC][0];//index of the parent, ndI is being modified here since parent is the current node under consideration
				if(arrBufs.stackBuf[stackC][1] == gNodeBuf[ndI].numOfChild)
				{//all children of this node are visited					
					--stackC;//this node visited,go to his parent continuing search for the next node in tree in pre-order
				}
				else
				{//some children are still left, put next child on top of stack					
					chI = gNodeBuf[ndI].childInd + arrBufs.stackBuf[stackC][1];
					arrBufs.stackBuf[++stackC][0] = tree.nodeInd+gChildBuf[chI];//node index of next child
					//no need to init this child, as it would be done when it picked from the stack, for the first time
					break;
				}
			}
		}
		else
		{//an intermediate node picked from the stack for the first time, init the ndoe and put its first child on the stack
			arrBufs.stackBuf[stackC][1] = 0; //visited cildren  count
			//now put its first child on stack
			++stackC;//increment stack
			arrBufs.stackBuf[stackC][0]= gChildBuf[gNodeBuf[ndI].childInd] + tree.nodeInd;//put the first child on stack
		}
	}

	//first the existing root will have one less child, remove the child
	j = gNodeBuf[tree.nodeInd].childInd + arrBufs.stackBuf[0][1] -1;//index of the child to be removed
	k = gNodeBuf[tree.nodeInd].childInd + gNodeBuf[tree.nodeInd].numOfChild -1;//last child index
	gChildBuf[j] = gChildBuf[k];//overwrite the existing child with the last child
	--gNodeBuf[tree.nodeInd].numOfChild;//decrement the number of children


	//copy the existing root node to a new node so that tree.nodeInd is now free to become the new root
	gNodeBuf[gGlobalBuf.nodeBufInd] = gNodeBuf[tree.nodeInd];
	arrBufs.stackBuf[0][0]  = gGlobalBuf.nodeBufInd++;//replace tree.nodeInd with the newly allocated node and --
	//-- increment gGlobalBuf.nodeBufInd for the next slot
	//currently the smTx is on the top of stack, replace it with the new root
	gNodeBuf[tree.nodeInd].childInd = gGlobalBuf.childBufInd;
	gNodeBuf[tree.nodeInd].numOfChild = 2;
	gChildBuf[gGlobalBuf.childBufInd] = ndI - tree.nodeInd;//make smTx as the first child
	arrBufs.stackBuf[stackC][0] = tree.nodeInd;//put new root on the stack
	arrBufs.stackBuf[stackC][1] = 2;//mark that the new child will be added as second child
	gGlobalBuf.childBufInd+= 2;

	while(stackC){//now reverse the parents
		chI = gNodeBuf[arrBufs.stackBuf[stackC][0]].childInd;
		gChildBuf[chI+ arrBufs.stackBuf[stackC][1] -1] = arrBufs.stackBuf[stackC-1][0] - tree.nodeInd;
		--stackC;
	}	
}
