#include"GTD_TreeProcessor.h"
#include"GTD_common.h"
#include"GTD_MemBuf.h"

void TreeProcessor::parseMulNewick(char *strTree,unsigned int treeInd)
{
	char *pos1, *pos2;
	size_t len;
	char ch;//tmp char
	unsigned short count =0; //to store the bracket count, init to 0
	unsigned short taxaC=0;//to store taxa count, init to 0
	unsigned short stackC=0;//counter which points to the last stored element in tmpStack, init to 0
	unsigned short symC=0;// to index next avail pos in TmpBuff.symArr, init to 0
	unsigned int *stack = MB.uintStack;

	MulTree &currTree= MB.treeBuf[treeInd];	
	char * treeAsChars = strTree;
	pair<set<TaxaLabel>::iterator,bool> insRet;
	TaxaLabel tmpLabel;

	//Only one pass is required
	//NOTE: the way newick tree is written the traversal is going to be in pre-order i.e root, left child, right child
	pos1 = strchr(strTree,'(');//find where the root starts
	++pos1;//skip the root, as it does not have any parent hence needs to be dealt before
	unsigned int chI;//temp var	
	currTree.rootInd = MB.nodeBufInd++;//assign root a node
	stack[stackC] = currTree.rootInd;//init to root
	unsigned short lfC=0;//leaf count of current tree
	unsigned short curDepth=0;
	do
	{
		ch = *pos1;
		if('(' == ch)
		{//new intermediate node
			++curDepth;
			chI =  MB.nodeBufInd++;//assign a node to this new node
			MB.nodeBuf[stack[stackC]].addChild(stack[stackC],chI);//add this new node to its parent
			stack[++stackC] = chI;//put this new child on the stack
			++count;//bracket count
			++pos1;			
		}
		else if(')' == ch )
		{
			--count;//bracket count
			--stackC;//decrement stackC to go one level up
			pos1 = strpbrk(pos1+1,":,)");//skip bl or other nodal attributes
			--curDepth;
			++IO_Stat.internalNodesIp;
		}
		else if(',' == ch)
		{
			++pos1;
		}
		else if(':' == ch)//newick format can have this, needs to be ignored
		{
			pos1 = strpbrk(pos1,",)");//skip bl or other nodal attributes
		}
		else //taxaId found
		{//new leaf node
			chI =  MB.nodeBufInd++;//assign a node to this new leaf node
			MB.nodeBuf[stack[stackC]].addChild(stack[stackC],chI);//add this new node to its parent
			pos2 = strpbrk(pos1,":,)");//':' is to ignore other nodal attributes
			len = pos2-pos1;
			strncpy(MB.taxaBuf[MB.taxaBufInd],pos1,len);
			MB.taxaBuf[MB.taxaBufInd][len]='\0';//terminate the string
			MB.taxaLenBuf[MB.taxaBufInd]=len;//store the strlen

			tmpLabel.ind = MB.taxaBufInd;

			insRet = MB.taxaIdSet.insert(tmpLabel);
			if(insRet.second){
				++MB.taxaBufInd;
			}
			MB.nodeBuf[chI].numChild= insRet.first->ind;//taxaId
			pos1 = pos2;//skip taxa id for next pos	
			++IO_Stat.leavesIp;
		}
	}
	while(USHRT_MAX != count);//tree has been parsed, sine the '(' for root was skipped, hence-
	//- count will reach -1 (USHRT_MAX) only when the very last ')' is reached.
	--stackC;//not needed for the root node, but just to make it look neat and in synch with the above for loop
	currTree.numNodes= (unsigned short)(MB.nodeBufInd - currTree.rootInd);
	IO_Stat.taxaIp = MB.taxaBufInd;
}

void TreeProcessor::parseMulNewickBQred(char *strTree,unsigned int treeInd)
{
	char *pos1, *pos2;
	size_t len;
	char ch;//tmp char
	unsigned short count =0; //to store the bracket count, init to 0
	unsigned short taxaC=0;//to store taxa count, init to 0
	unsigned short stackC=0;//counter which points to the last stored element in tmpStack, init to 0
	unsigned short symC=0;// to index next avail pos in TmpBuff.symArr, init to 0
	unsigned int *stack = MB.uintStack;
	
	MulTree &currTree= MB.treeBuf[treeInd];	
	char * treeAsChars = strTree;
	pair<set<TaxaLabel>::iterator,bool> insRet;
	TaxaLabel tmpLabel;

	//Only one pass is required
	//NOTE: the way newick tree is written the traversal is going to be in pre-order i.e root, left child, right child
	pos1 = strchr(strTree,'(');//find where the root starts
	++pos1;//skip the root, as it does not have any parent hence needs to be dealt before
	unsigned int chI;//temp var	
	currTree.rootInd = MB.nodeBufInd++;//assign root a node
	stack[stackC] = currTree.rootInd;//init to root
	unsigned short lfC=0;//leaf count of current tree
	unsigned short curDepth=0;
	do
	{
		ch = *pos1;
		if('(' == ch)
		{//new intermediate node
			++curDepth;
			chI =  MB.nodeBufInd++;//assign a node to this new node
			MB.nodeBuf[stack[stackC]].addChild(stack[stackC],chI);//add this new node to its parent
			stack[++stackC] = chI;//put this new child on the stack
			++count;//bracket count
			++pos1;			
		}
		else if(')' == ch )
		{
			--count;//bracket count
			--stackC;//decrement stackC to go one level up
			pos1 = strpbrk(pos1+1,":,)");//skip bl or other nodal attributes
			--curDepth;
			++IO_Stat.internalNodesIp;
		}
		else if(',' == ch)
		{
			++pos1;
		}
		else if(':' == ch)//newick format can have this, needs to be ignored
		{
			pos1 = strpbrk(pos1,",)");//skip bl or other nodal attributes
		}
		else //taxaId found
		{//new leaf node
			chI =  MB.nodeBufInd++;//assign a node to this new leaf node
			MB.nodeBuf[stack[stackC]].addChild(stack[stackC],chI);//add this new node to its parent
			pos1+=2;//skip ti
			pos2 = strchr(pos1,'_');
			len = pos2-pos1;			
			strncpy(MB.taxaBuf[MB.taxaBufInd],pos1,len);
			MB.taxaBuf[MB.taxaBufInd][len]='\0';//terminate the string
			MB.taxaLenBuf[MB.taxaBufInd]=len;//store the strlen
			tmpLabel.ind = MB.taxaBufInd;

			insRet = MB.taxaIdSet.insert(tmpLabel);
			if(insRet.second){
				++MB.taxaBufInd;
			}
			MB.nodeBuf[chI].numChild= insRet.first->ind;//taxaId

			pos1 = pos2+3;//skip _gi
			pos2 = strpbrk(pos1,":,)");//':' is to ignore other nodal attributes
			ch = *pos2;//temporarily store
			*pos2 =0;//null char
			MB.nodeBuf[chI].lastChInd = atoi(pos1);	//assign gId to lastChInd			
			*pos2 = ch;//restore						
			pos1 = pos2;//skip taxa id for next pos	
			++IO_Stat.leavesIp;
		}
	}
	while(USHRT_MAX != count);//tree has been parsed, sine the '(' for root was skipped, hence-
	//- count will reach -1 (USHRT_MAX) only when the very last ')' is reached.
	--stackC;//not needed for the root node, but just to make it look neat and in synch with the above for loop
	currTree.numNodes= (unsigned short)(MB.nodeBufInd - currTree.rootInd);
	IO_Stat.taxaIp = MB.taxaBufInd;
}





//uses MB.uintStack,MB.uintStack2
void TreeProcessor::computeAllSuvValues(unsigned int treeInd)
{
	//OK1
	MulTree &currTree= MB.treeBuf[treeInd];	
	for(unsigned short i=0; i<MB.taxaBufInd; ++i){
		MB.isTaxaValidArr[i] = false;
	}


	unsigned int *nxtChSt = MB.uintStack;//next child to be visited  of the node on the stack	
	unsigned int *stack = MB.uintStack2;
	unsigned short stC=0;
	stack[stC] = MB.treeBuf[treeInd].rootInd;	
	nxtChSt[stC] = stack[stC];//to indicate that it is being visited the first time
	
	unsigned int ndI;//tmp var
	while(USHRT_MAX != stC){
		ndI = stack[stC];
		if(UINT_MAX == MB.nodeBuf[ndI].firstChInd){//a leaf node
			--stC;
		}else{//intermediate node			
			if(ndI == nxtChSt[stC]){//first time being visited
				nxtChSt[stC] = MB.nodeBuf[ndI].firstChInd;
			}
			if(UINT_MAX != nxtChSt[stC]){
				if(UINT_MAX != MB.nodeBuf[nxtChSt[stC]].firstChInd){//a non-leaf child
					computeNodeSuvValue(treeInd,nxtChSt[stC],stack,stC);
				}
				nxtChSt[stC+1] = stack[stC+1] = nxtChSt[stC];//put this child on stack
				nxtChSt[stC] = MB.nodeBuf[nxtChSt[stC]].nextSibInd;//make it piont to the next child to be visited				
				++stC;//keep increment separate as it is being referenced at quite a few places
			}else{//all childrent visited
				--stC;//remove it from the stack
			}
		}
	}
}

//uses MB.uintStack3, MB.boolStack, MB.boolStack2
void TreeProcessor::computeNodeSuvValue(unsigned int treeInd, unsigned int curNdInd, unsigned int *parSt, unsigned short parStC)
{
	//OK1
	// we need to travel in two directions for this node, one up and one down
	unsigned short stackC = 0;
	unsigned int *stack = MB.uintStack3;
	bool * SdnFlgArr = MB.boolStack;
	bool * SupFlgArr = MB.boolStack2;
	unsigned short numOfTaxa = MB.taxaBufInd;

	for(unsigned short i=0; i<numOfTaxa; ++i){
		SdnFlgArr[i] = SupFlgArr[i] = false;			
	}

	unsigned short taxaId;
	// first visit the down part, which is easier
	stack[stackC] = curNdInd;// put current node on stack
	unsigned int ndI,chI;
	while (USHRT_MAX != stackC) {
		ndI = stack[stackC];
		if (UINT_MAX == MB.nodeBuf[ndI].firstChInd) {// a leaf node
			taxaId = MB.nodeBuf[ndI].numChild;// taxa Id
			if (!SdnFlgArr[taxaId]) {// this taxa has not been encountered till now
				SdnFlgArr[taxaId] = true;
				++MB.nodeBuf[curNdInd].Sdn;
			}
			--stackC;
		} else {// an intermediate node
			--stackC;// remove it from the stack and put all its children on stack
			chI = MB.nodeBuf[ndI].firstChInd;
			while(UINT_MAX != chI){
				stack[++stackC] = chI;
				chI = MB.nodeBuf[chI].nextSibInd;
			}
		}
	}

	// now we travel upwards, which is the difficult part
	// travel from the current node till the root and put all the siblings encountered on the stack
	stackC = USHRT_MAX;
	ndI = curNdInd;
	unsigned int parChldInd;// child index of its parent
	unsigned int parInd;// nodeInd of parent

	while(USHRT_MAX != parStC){
		parInd = parSt[parStC--];
		parChldInd = MB.nodeBuf[parInd].firstChInd;
		while(UINT_MAX != parChldInd){
			if(ndI != parChldInd){
				// i.e. current child is a sibling and not the node itself
				stack[++stackC] = parChldInd;
			}
			parChldInd = MB.nodeBuf[parChldInd].nextSibInd;
		}
		ndI = parInd;// move up to its parent
	}


	// now traverse the induced subtrees at all the nodes in the stack
	while (USHRT_MAX != stackC) {
		ndI = stack[stackC];
		if (UINT_MAX == MB.nodeBuf[ndI].firstChInd) {// a leaf node
			taxaId = MB.nodeBuf[ndI].numChild;// taxa Id
			if (!SupFlgArr[taxaId]) {// this taxa has not been encountered till now
				SupFlgArr[taxaId] = true;
				++MB.nodeBuf[curNdInd].Sup;
			}
			--stackC;
		} else {// an intermediate node
			--stackC;// remove it from the stack and put all its children on stack
			chI = MB.nodeBuf[ndI].firstChInd;
			while(UINT_MAX != chI){
				stack[++stackC] = chI;
				chI = MB.nodeBuf[chI].nextSibInd;
			}
		}
	}

	// calculate the actual Sdn Sup values
	unsigned short tmpSwp = MB.nodeBuf[curNdInd].Sup;
	MB.nodeBuf[curNdInd].Sup = numOfTaxa -  MB.nodeBuf[curNdInd].Sdn;
	MB.nodeBuf[curNdInd].Sdn = numOfTaxa - tmpSwp;
	// as number of species exclusively on one side = total number of species - number of species on the other side
	if (MB.nodeBuf[curNdInd].isInformative()) {
		//an informative edge- OK to update the valid taxaFlgArr
		for(unsigned short i=0; i<numOfTaxa; ++i){
			if((SdnFlgArr[i] && !SupFlgArr[i]) || (!SdnFlgArr[i] && SupFlgArr[i])){
				//i.e. this taxa occurs only one side
				MB.isTaxaValidArr[i] = true;
			}
		}

	}
}



//uses MB.uintStack
unsigned short  TreeProcessor::collapseNonInformativeEdges(unsigned int treeInd)
{
	//OK1
	unsigned short stackC = 0, infoEdgeC=0;
	unsigned int *stack = MB.uintStack;
	unsigned int ndI,chI;

	stack[stackC] = MB.treeBuf[treeInd].rootInd;
	while (USHRT_MAX != stackC) {
		ndI = stack[stackC];
		// necessarily an intermediate node
		chI = MB.nodeBuf[ndI].firstChInd;
		while(UINT_MAX != chI){
			if(UINT_MAX != MB.nodeBuf[chI].firstChInd) {
				if(!MB.nodeBuf[chI].isInformative()){
					//not a leaf node and a non-informative edge					
					MB.nodeBuf[ndI].collapseChild(chI);
				}else{//not a leaf node and informative
					++infoEdgeC;
				}
			}
			chI = MB.nodeBuf[chI].nextSibInd;
		}
		--stackC;// remove it from the stack and put all its children on stack
		chI = MB.nodeBuf[ndI].firstChInd;
		while(UINT_MAX != chI ){
			if(UINT_MAX != MB.nodeBuf[chI].firstChInd){//not a leaf node
				stack[++stackC] = chI;
			}
			chI = MB.nodeBuf[chI].nextSibInd;
		}		
	}
	return infoEdgeC;
}


//uses MB.uintStack
void TreeProcessor::removeNonParticipatingLeaves(unsigned int treeInd)
{
	//OK1
	//it is assumed that non-informatie edges and dup-edges have already been collapsed
	unsigned short stackC = 0;
	unsigned int *stack = MB.uintStack;
	unsigned int ndI,chI;

	stack[stackC] = MB.treeBuf[treeInd].rootInd;
	while (USHRT_MAX != stackC) {
		ndI = stack[stackC];
		// necessarily an intermediate node		
		--stackC;// remove it from the stack and put all its children on stack
		chI = MB.nodeBuf[ndI].firstChInd;
		while(UINT_MAX != chI ){
			if(UINT_MAX != MB.nodeBuf[chI].firstChInd){//not a leaf node
				stack[++stackC] = chI;
			}else if(!MB.isTaxaValidArr[MB.nodeBuf[chI].numChild]){
				MB.nodeBuf[ndI].deleteChild(chI);
			}
			chI = MB.nodeBuf[chI].nextSibInd;
		}		
	}

	//the below is not required since the non-informatie and dup edges have already been collapsed
	//check if the root has only two children
	//ndI = MB.treeBuf[treeInd].rootInd;
	//if(2 == MB.nodeBuf[ndI].numChild){
	//	chI = MB.nodeBuf[ndI].firstChInd;
	//	if(UINT_MAX == MB.nodeBuf[chI].firstChInd){
	//		//first child is leaf, collapse the second child
	//		chI = MB.nodeBuf[ndI].lastChInd;
	//		cout<<"Error! first child is leaf, collapse the second child."<<endl;
	//	}
	//	MB.nodeBuf[ndI].collapseChild(chI);
	//}
}


//uses MB.uintStack, MB.uintStack2
void TreeProcessor::collapseDupEdges(unsigned int treeInd)
{
	//OK1
	//it is assumed that non-informative edges have been collapsed before calling this step
	unsigned short stackC = 0;
	unsigned int *stack = MB.uintStack;
	unsigned int rootInd = MB.treeBuf[treeInd].rootInd;
	stack[stackC] = rootInd;// put the root on stack
	unsigned int ndI,chI, infNb1, infNb2;
	int compRes;
	//unsigned int firstChInd;

	while (USHRT_MAX != stackC) {
		ndI = stack[stackC];
		// necessarily an intermediate node
		infNb1 = infNb2 = 0;
		//find 1st informative neighbor
		//check if this node itself is informative
		if(rootInd != ndI){//only if not the root node
			infNb1 = ndI;//necessarily informative as non-informative edges have been collapsed
		}
		
		//now check the children for the remaining informative nodes
		chI = MB.nodeBuf[ndI].firstChInd;
		while(UINT_MAX != chI){
			if(UINT_MAX != MB.nodeBuf[chI].firstChInd){//not a leaf node
				//necssarily informative as as non-informative edges have been collapsed
				if(!infNb1){
					infNb1 = chI;						
				}else{
					infNb2 = chI;
					break;
				}
			}
			chI = MB.nodeBuf[chI].nextSibInd;
		}			
		if(infNb1 && infNb2){//two informative edges found
			if(infNb1 == ndI){//parent-child relationship
				compRes = MB.nodeBuf[ndI].compInfWidChld(infNb2);
				if(MulTreeNode::INFO_EQUAL == compRes){
					//delete rest of the siblings
					chI = MB.nodeBuf[ndI].firstChInd;
					while(UINT_MAX != chI){
						if(chI != infNb2){							
							MB.nodeBuf[ndI].deleteChild(chI);
						}
						chI = MB.nodeBuf[chI].nextSibInd;
					}
					//collapse the child
					MB.nodeBuf[infNb2].markForCollapse();
				}else if(MulTreeNode::INFO_GREATER == compRes){
					MB.nodeBuf[infNb2].markForCollapse();
				}else if(MulTreeNode::INFO_LESSER == compRes){
					MB.nodeBuf[ndI].markForCollapse();
				}
			}else{//sibling relationship
				compRes = MB.nodeBuf[infNb1].compInfWidSib(infNb2);
				if(MulTreeNode::INFO_EQUAL == compRes){
					//only delete rest of the siblings. rootInd == ndI is necessarily the case as -
					//- all the non-informative edges have already been collapsed, thus there cannot exist a -
					//- non-informative subtree branching upwards
					chI = MB.nodeBuf[ndI].firstChInd;
					while(UINT_MAX != chI){
						if(chI != infNb1 && chI != infNb2){
							MB.nodeBuf[ndI].deleteChild(chI);
						}
						chI = MB.nodeBuf[chI].nextSibInd;
					}
					//collapse one of the children
					MB.nodeBuf[infNb1].markForCollapse();
				}else if(MulTreeNode::INFO_GREATER == compRes){
					MB.nodeBuf[infNb2].markForCollapse();
				}else if(MulTreeNode::INFO_LESSER == compRes){
					MB.nodeBuf[infNb1].markForCollapse();
				}
			}
		}

		--stackC;// remove it from the stack and put all its children on stack
		chI = MB.nodeBuf[ndI].firstChInd;
		while(UINT_MAX != chI){
			if(UINT_MAX != MB.nodeBuf[chI].firstChInd){//only put non-leaf children
				stack[++stackC] = chI;
			}
			chI = MB.nodeBuf[chI].nextSibInd;
		}
	}

	//now collapse the marked nodes
	stackC=0;
	stack[stackC]= rootInd;
	while (USHRT_MAX != stackC) {
		ndI = stack[stackC];
		//necessarily an intermediate node
		//collapse marked children
		chI = MB.nodeBuf[ndI].firstChInd;
		while(UINT_MAX != chI){
			MB.nodeBuf[ndI].chkN_Collapse(chI);
			chI = MB.nodeBuf[chI].nextSibInd;
		}
		--stackC;// remove it from the stack and put all its children on stack
		chI = MB.nodeBuf[ndI].firstChInd;
		while(UINT_MAX != chI){
			if(UINT_MAX != MB.nodeBuf[chI].firstChInd){//only put non-leaf children
				stack[++stackC] = chI;
			}
			chI = MB.nodeBuf[chI].nextSibInd;
		}
	}
}

//uses MB.uintStack, MB.boolStack, MB.ushrtStack
void TreeProcessor::removeDupLeaves(unsigned int treeInd)
{
	//OK1
	unsigned short stackC = 0;
	unsigned int *stack = MB.uintStack;
	unsigned int ndI,chI;
	bool * dupLfFlgArr = MB.boolStack;
	unsigned short * dupLfArr = MB.ushrtStack;//to keep count of leaves in dupLfFlgArr, so -
	//-that it can be reset
	unsigned short dupLfC;
	for(unsigned short i=0; i<MB.taxaBufInd; ++i){
		dupLfFlgArr[i] = false;
	}

	stack[stackC] = MB.treeBuf[treeInd].rootInd;
	while (USHRT_MAX != stackC) {
		ndI = stack[stackC];
		// necessarily an intermediate node		
		--stackC;// remove it from the stack and put all its children on stack
		chI = MB.nodeBuf[ndI].firstChInd;
		dupLfC =0;
		while(UINT_MAX != chI ){
			if(UINT_MAX != MB.nodeBuf[chI].firstChInd){//not a leaf node
				stack[++stackC] = chI;
			}else{//a leaf child
				if(!dupLfFlgArr[MB.nodeBuf[chI].numChild]){//first encounter
					dupLfFlgArr[MB.nodeBuf[chI].numChild] = true;
					dupLfArr[dupLfC++] = MB.nodeBuf[chI].numChild;//store this leaf for later reset of dupLfFlgArr
				}else{//a dup leaf, delete this child
					MB.nodeBuf[ndI].deleteChild(chI);
				}
			}
			chI = MB.nodeBuf[chI].nextSibInd;
		}		
		//reset dupLfFlgArr for the next use
		for(unsigned short i=0; i<dupLfC; ++i){
			dupLfFlgArr[dupLfArr[i]] = false;
		}
	}

	//the below is not required since the non-informatie edges have already been collapsed
	//check if the root has only two children
	//ndI = MB.treeBuf[treeInd].rootInd;
	//if(2 == MB.nodeBuf[ndI].numChild){
	//	chI = MB.nodeBuf[ndI].firstChInd;
	//	if(UINT_MAX == MB.nodeBuf[chI].firstChInd){
	//		//first child is leaf, collapse the second child
	//		chI = MB.nodeBuf[ndI].lastChInd;
	//	}
	//	MB.nodeBuf[ndI].collapseChild(chI);
	//}
}

//uses MB.uintStack, MB.ushrtStack
void TreeProcessor::remAllRedundantLeaves(unsigned int treeInd)
{
	//OK1
	MulTree &currTree= MB.treeBuf[treeInd];	
	unsigned int ndI;//temp var
	unsigned short stackC=0;
	unsigned int *stack = MB.uintStack;
	unsigned short *leafCntArr = MB.ushrtStack;
	for(unsigned short i=0; i<MB.taxaBufInd; ++i){
		leafCntArr[i] =0;
	}
	stack[stackC] = currTree.rootInd;
	while(USHRT_MAX != stackC)
	{
		ndI = stack[stackC];
		if(UINT_MAX == MB.nodeBuf[ndI].firstChInd)
		{//a leaf node
			++leafCntArr[MB.nodeBuf[ndI].numChild];
			--stackC;
		}else{//an intermediate node
			--stackC;//remove it from stack and put all its children on stack
			ndI = MB.nodeBuf[ndI].firstChInd;
			while(UINT_MAX != ndI){
				stack[++stackC] = ndI;//put the child on stack
				ndI =MB.nodeBuf[ndI].nextSibInd;
			}
		}
	}

	for(unsigned short i=0; i<MB.taxaBufInd; ++i){
		if(leafCntArr[i]>2){//at least 3 multiplicity is required for the possibility of this taxa having redundant leaves.-
			//- This is becuase duplicate and non-participating leaves have already been removed.
			removeRedundantLeaf(treeInd,i);
		}
	}
}

//uses MB.uintStack,MB.uintStack2, MB.boolStack
void TreeProcessor::removeRedundantLeaf(unsigned int treeInd,unsigned short taxaId)
{
	//OK1
	//WARNING: this function uses (modifies) Sdn and Sup values of nodes
	MulTree &currTree= MB.treeBuf[treeInd];	
	unsigned int ndI,chI;//temp var
	unsigned short stackC=0;
	unsigned int *stack = MB.uintStack;
	unsigned int *parStack = MB.uintStack2;
	bool *boolSt = MB.boolStack;
	stack[stackC] = currTree.rootInd;
	boolSt[stackC]= false;//mark root as not visited

	//first we will mark all the Sdn values of the nodes. Sdn =1 means taxaId is present in the down subtree of the edge
	while(USHRT_MAX != stackC)
	{
		ndI = stack[stackC];
		if(UINT_MAX == MB.nodeBuf[ndI].firstChInd)
		{//a leaf node
			if(taxaId == MB.nodeBuf[ndI].numChild){
				MB.nodeBuf[ndI].Sdn =1;//mark that the down-subtree corresponding to this leaf has taxa=taxaId
				MB.nodeBuf[parStack[stackC]].Sdn = 1;//mark that the down-subtree of its parent has taxa=taxaId
			}
			--stackC;
		}else
		{//an intermediate node
			if(boolSt[stackC])
			{//all its children have been visited,
				//if down subtree of this node contains taxaId, then mark the down-subtree of its parent as well
				if(1 == MB.nodeBuf[ndI].Sdn && 0 != stackC){//only if root node itself is not on the stack
					MB.nodeBuf[parStack[stackC]].Sdn = 1;
				}
				//remove it from the stack.	we need to keep it on the stack so that its children may refer to it
				--stackC;
			}
			else
			{//encountered first time
				boolSt[stackC] = true;//mark this node as visited
				//put all its children on stack and mark them as not-visited
				chI = MB.nodeBuf[ndI].firstChInd;
				while(UINT_MAX != chI){
					MB.nodeBuf[chI].Sdn =0;
					stack[++stackC] = chI;//put the child on stack
					parStack[stackC] = ndI;//put its parent on the parent stack
					boolSt[stackC] = false;//mark it as not-visited
					chI =MB.nodeBuf[chI].nextSibInd;
				}	
			}
		}
	}


	//now we will mark the Sup values
	stackC =0;
	stack[stackC] = currTree.rootInd;
	unsigned int valNd1, valNd2;
	unsigned short subGraphDegC;
	while(USHRT_MAX != stackC)
	{
		ndI = stack[stackC];
		//necessarliy an intermediate node
		valNd1 = valNd2 = 0;
		//now find two valid nodes in its neighbor having taxaId in their subtree
		//check if this node itself is valid
		if(ndI != currTree.rootInd && MB.nodeBuf[ndI].Sup){
			valNd1 = ndI;
		}
		//now check rest of its children for remaining valid nodes
		chI = MB.nodeBuf[ndI].firstChInd;
		while(UINT_MAX != chI){
			if(MB.nodeBuf[chI].Sdn){
				if(!valNd1){
					valNd1 =chI ;
				}else{
					valNd2 =chI ;
					break;
				}
			}
			chI = MB.nodeBuf[chI].nextSibInd;
		}

		if(valNd1 && valNd2){//both valid
			//mark all the children's Sup=1
			chI = MB.nodeBuf[ndI].firstChInd;
			while(UINT_MAX != chI){				
				MB.nodeBuf[chI].Sup = 1;
				chI = MB.nodeBuf[chI].nextSibInd;
			}
		}else{//only first valid
			if(valNd1 != ndI){//applicable only for a child node
				MB.nodeBuf[valNd1].Sup =0;
			}
			//mark rest of the children as valid
			chI = MB.nodeBuf[ndI].firstChInd;
			while(UINT_MAX != chI){
				if(valNd1 != chI){
					MB.nodeBuf[chI].Sup = 1;
				}
				chI = MB.nodeBuf[chI].nextSibInd;
			}
		}

		
		//now all the neighbors of this have their Sup and Sdn values set w.r.t. taxaId, now -
		//- we can delete its redundant leaf children
		subGraphDegC =0;
		//first calculate the sub graph deg count of this node
		if(ndI != currTree.rootInd){
			if(MB.nodeBuf[ndI].Sup){
				++subGraphDegC;
			}
		}

		chI = MB.nodeBuf[ndI].firstChInd;
		while(UINT_MAX != chI){
			if(MB.nodeBuf[chI].Sdn){
					++subGraphDegC;
					if(subGraphDegC>2){
						break;
					}
			}
			chI =MB.nodeBuf[chI].nextSibInd;
		}

		if(subGraphDegC>2){//check if this has any redundant leaf children, if so delete them
			chI = MB.nodeBuf[ndI].firstChInd;
			while(UINT_MAX != chI){
				if(UINT_MAX == MB.nodeBuf[chI].firstChInd &&
					taxaId == MB.nodeBuf[chI].numChild){
						MB.nodeBuf[ndI].deleteChild(chI);
				}
				chI =MB.nodeBuf[chI].nextSibInd;
			}
		}

		//remove it from the stack and put all its children on stack 
		--stackC;		
		chI = MB.nodeBuf[ndI].firstChInd;
		while(UINT_MAX != chI){
			if(UINT_MAX != MB.nodeBuf[chI].firstChInd){//not a leaf child
				stack[++stackC] = chI;//put the child on stack
			}
			chI =MB.nodeBuf[chI].nextSibInd;
		}	
	}
}


void TreeProcessor::print(unsigned int treeInd)
{
//#define PRINT_GEN_TAXA_ID
	MulTree &currTree= MB.treeBuf[treeInd];	
	unsigned int ndI;//temp var
	unsigned short stackC=0;
	char *ep=MB.chBuf;
	unsigned int *prSt = MB.printStack;
	bool *boolSt = MB.boolStack;
	prSt[stackC] = currTree.rootInd;
	boolSt[stackC]= false;//mark root as not visited

	while(USHRT_MAX != stackC)
	{
		ndI = prSt[stackC];
		if(UINT_MAX == MB.nodeBuf[ndI].firstChInd)
		{//a leaf node
#ifdef PRINT_GEN_TAXA_ID
			itoa(MB.nodeBuf[ndI].numChild,ep,10);
			ep = ep + strlen(ep);
#else
			strcpy(ep,MB.taxaBuf[MB.nodeBuf[ndI].numChild]);
			ep = ep + MB.taxaLenBuf[MB.nodeBuf[ndI].numChild];//contains strlen of the taxa
#endif
			*(ep++) = ',';
			--stackC;
			
		}
		else
		{//an intermediate node
			if(boolSt[stackC])
			{//all its children have been visited
				*(--ep) = ')';//overwirte the extra ','		
				*(++ep) =',';
				++ep;
				--stackC;
				
			}
			else
			{//encountered first time
				*(ep++)='(';
				boolSt[stackC] = true;//mark this node as visited
				//put all its children on stack and mark them as not-visited
				stackC+= MB.nodeBuf[ndI].numChild;				
				ndI = MB.nodeBuf[ndI].firstChInd;
				while(UINT_MAX != ndI){
					prSt[stackC] = ndI;//put the child on stack
					boolSt[stackC--] = false;//makr it as not-visited
					ndI =MB.nodeBuf[ndI].nextSibInd;
				}	
				stackC+= MB.nodeBuf[prSt[stackC]].numChild;
			}
		}
	}
	*(--ep) = '\0'; //overwrite the last ','	


}




void TreeProcessor::printTableBQred(ofstream &resFile, unsigned int treeInd)
{

//Print tree

	MulTree &currTree= MB.treeBuf[treeInd];	
	unsigned int ndI;//temp var
	unsigned short stackC=0;
	unsigned int *prSt = MB.printStack;
	bool *boolSt = MB.boolStack;
	prSt[stackC] = currTree.rootInd;
	boolSt[stackC]= false;//mark root as not visited	

	while(USHRT_MAX != stackC)
	{
		ndI = prSt[stackC];
		if(UINT_MAX == MB.nodeBuf[ndI].firstChInd)
		{//a leaf node
			resFile<<MB.taxaBuf[MB.nodeBuf[ndI].numChild]<<"\t"<<MB.nodeBuf[ndI].lastChInd<<endl;
			--stackC;
		}
		else
		{//an intermediate node
			if(boolSt[stackC])
			{//all its children have been visited
				--stackC;				
			}
			else
			{//encountered first time
				boolSt[stackC] = true;//mark this node as visited
				//put all its children on stack and mark them as not-visited
				stackC+= MB.nodeBuf[ndI].numChild;				
				ndI = MB.nodeBuf[ndI].firstChInd;
				while(UINT_MAX != ndI){
					prSt[stackC] = ndI;//put the child on stack
					boolSt[stackC--] = false;//makr it as not-visited
					ndI =MB.nodeBuf[ndI].nextSibInd;
				}	
				stackC+= MB.nodeBuf[prSt[stackC]].numChild;
			}
		}
	}

}




void TreeProcessor::computeStats(unsigned int treeInd)
{
	MulTree &currTree= MB.treeBuf[treeInd];	
	unsigned int ndI;//temp var
	unsigned short stackC=0;
	
	unsigned int *stack = MB.uintStack;
	bool *boolSt = MB.boolStack;
	unsigned int *isTaxaMul = MB.uintStack2;
	stack[stackC] = currTree.rootInd;
	boolSt[stackC]= false;//mark root as not visited

	for(unsigned short i=0; i<MB.taxaBufInd; ++i){
		isTaxaMul[i] = UINT_MAX;
	}


	while(USHRT_MAX != stackC)
	{
		ndI = stack[stackC];
		if(UINT_MAX == MB.nodeBuf[ndI].firstChInd)
		{//a leaf node
			if(UINT_MAX == isTaxaMul[MB.nodeBuf[ndI].numChild]){
				isTaxaMul[MB.nodeBuf[ndI].numChild] = 0;//encountered first time				
			}else if(0 == isTaxaMul[MB.nodeBuf[ndI].numChild]){//eccountered second time
				isTaxaMul[MB.nodeBuf[ndI].numChild] = 1;
				++IO_Stat.mulTaxaOp;	
			}//else encountered third time or more
			--stackC;
			++IO_Stat.leavesOp;			
		}
		else
		{//an intermediate node
			if(boolSt[stackC])
			{//all its children have been visited
				--stackC;
				++IO_Stat.internalNodesOp;
			}
			else
			{//encountered first time
				boolSt[stackC] = true;//mark this node as visited
				//put all its children on stack and mark them as not-visited
				stackC+= MB.nodeBuf[ndI].numChild;				
				ndI = MB.nodeBuf[ndI].firstChInd;
				while(UINT_MAX != ndI){
					stack[stackC] = ndI;//put the child on stack
					boolSt[stackC--] = false;//makr it as not-visited
					ndI =MB.nodeBuf[ndI].nextSibInd;
				}	
				stackC+= MB.nodeBuf[stack[stackC]].numChild;
			}
		}
	}
	
	for(unsigned short i=0; i<MB.taxaBufInd; ++i){
		if(MB.isTaxaValidArr[i]){
			++IO_Stat.taxaOp;
		}
	}

}

bool TreeProcessor::restrictToSingTaxa(unsigned int srcTreeInd, unsigned int destTreeInd)
{
	MulTree &currTree= MB.treeBuf[srcTreeInd];	
	unsigned int ndI;//temp var
	unsigned short stackC=0;
	
	unsigned int *stack = MB.uintStack;
	bool *boolSt = MB.boolStack;

	unsigned short singTaxaC = 0;
	
	unsigned int *isTaxaMul = MB.uintStack2;
	stack[stackC] = currTree.rootInd;
	boolSt[stackC]= false;//mark root as not visited

	for(unsigned short i=0; i<MB.taxaBufInd; ++i){
		isTaxaMul[i] = UINT_MAX;
	}


	while(USHRT_MAX != stackC)
	{
		ndI = stack[stackC];
		if(UINT_MAX == MB.nodeBuf[ndI].firstChInd)
		{//a leaf node
			if(UINT_MAX == isTaxaMul[MB.nodeBuf[ndI].numChild]){
				isTaxaMul[MB.nodeBuf[ndI].numChild] = 0;//encountered first time	
				++singTaxaC;
			}else if(0 == isTaxaMul[MB.nodeBuf[ndI].numChild]){//eccountered second time
				isTaxaMul[MB.nodeBuf[ndI].numChild] = 1;				
			}//else encountered third time or more
			--stackC;			
		}
		else
		{//an intermediate node
			if(boolSt[stackC])
			{//all its children have been visited
				--stackC;				
			}
			else
			{//encountered first time
				boolSt[stackC] = true;//mark this node as visited
				//put all its children on stack and mark them as not-visited
				stackC+= MB.nodeBuf[ndI].numChild;				
				ndI = MB.nodeBuf[ndI].firstChInd;
				while(UINT_MAX != ndI){
					stack[stackC] = ndI;//put the child on stack
					boolSt[stackC--] = false;//makr it as not-visited
					ndI =MB.nodeBuf[ndI].nextSibInd;
				}	
				stackC+= MB.nodeBuf[stack[stackC]].numChild;
			}
		}
	}

	if(singTaxaC<4){//assuming unrooted tree
		return false;
	}
	for(unsigned short i=0; i<MB.taxaBufInd; ++i){
		if(0 == isTaxaMul[i]){
			MB.isTaxaValidArr[i] = true;
		}else{
			MB.isTaxaValidArr[i] = false;		
		}
	}
	subtree(srcTreeInd,destTreeInd,MB.isTaxaValidArr,true);
	if((MB.treeBuf[destTreeInd].numNodes - singTaxaC) == 1){//a star-like unrooted tree
		return false;
	}
	
	return true;
}

void TreeProcessor::subtree(unsigned int srcTreeI, unsigned int destTreeI, bool *taxaArr, bool unRooted)
{
	MulTree &srcTree = MB.treeBuf[srcTreeI];
	MulTree &destTree = MB.treeBuf[destTreeI];

	unsigned int *stack = MB.uintStack;
	unsigned int *nxtChSt = MB.uintStack2;
	unsigned int *lastValidChSt = MB.uintStack3;//stores the last valid descendent of the current node in the dest tree-
	//- used to find if the current node would exists in the destination tree or not
	unsigned short *validChCSt = MB.ushrtStack;//to stroe count of valid children
	unsigned short stC =0;
	bool * boolSt = MB.boolStack;
	
	destTree.rootInd = MB.nodeBufInd++;
	destTree.numNodes=0;

	stack[stC] = srcTree.rootInd;
	boolSt[stC] = true;//mark root as not-visited

	unsigned int ndI, lastValidChInd;//tmp var
	while(USHRT_MAX != stC){
		ndI = stack[stC];
		if(UINT_MAX == MB.nodeBuf[ndI].firstChInd){//a leaf node
			if(taxaArr[MB.nodeBuf[ndI].numChild]){//exists in the query set
				lastValidChInd = MB.nodeBufInd;				
				MB.nodeBuf[MB.nodeBufInd].numChild = MB.nodeBuf[ndI].numChild;
				MB.nodeBuf[MB.nodeBufInd].lastChInd = MB.nodeBuf[ndI].lastChInd;//lastChInd acts as gi id for biclique red requirement
				++MB.nodeBufInd;//keeping the increment separate
			}else{
				lastValidChInd = UINT_MAX;//to indicate this child is not valid
			}
			--stC;//remove it from the stack
		}else{//intermediate node
			if(boolSt[stC]){//being visited first time
				boolSt[stC] = false;
				nxtChSt[stC]= MB.nodeBuf[ndI].firstChInd;//put the first child on nxtChSt				
				lastValidChInd = UINT_MAX;
				validChCSt[stC]= 0;
			}
			if(lastValidChInd != UINT_MAX){
				++validChCSt[stC];
				if(1 == validChCSt[stC]){
					lastValidChSt[stC] = lastValidChInd;
				}else if(2 == validChCSt[stC]){//this node will exist, allocate mem
					MB.nodeBuf[MB.nodeBufInd].addChild(MB.nodeBufInd,lastValidChSt[stC]);
					MB.nodeBuf[MB.nodeBufInd].addChild(MB.nodeBufInd,lastValidChInd);
					lastValidChSt[stC] = MB.nodeBufInd;//store its self in the stack so that ind can be used while adding further children
					++MB.nodeBufInd;
				}else{
					MB.nodeBuf[lastValidChSt[stC]].addChild(lastValidChSt[stC],lastValidChInd);					
				}
			}
			if(UINT_MAX != nxtChSt[stC]){
				stack[stC+1] = nxtChSt[stC];//put this child on stack
				boolSt[stC+1] = true;//mark this child as not visited
				nxtChSt[stC] = MB.nodeBuf[nxtChSt[stC]].nextSibInd;//make it piont to the next child to be visited
				++stC;//keep increment separate as it is being referenced at quite a few places				
			}else{//all childrent visited
				if(0 == validChCSt[stC]){
					lastValidChInd = UINT_MAX;//this node does not exist, nor does any of its children
				}else {
					lastValidChInd = lastValidChSt[stC];//either this node exist or one of its children do pass whatever is the case
				}
				if(stack[stC] == srcTree.rootInd){
					if(UINT_MAX != lastValidChInd){//i.e. at least one node exist in the subtree
						MB.nodeBuf[destTree.rootInd] = MB.nodeBuf[lastValidChInd];						
					}
				}
				--stC;//remove it from the stack
			}
		}	
	}
	destTree.numNodes = MB.nodeBufInd - destTree.rootInd -1; //-1 for the extra root node which we use

	if(unRooted){		
		degTwoCheck(destTreeI);
	}
}

void TreeProcessor::degTwoCheck(unsigned int treeI)
{
	MulTreeNode &root = MB.nodeBuf[MB.treeBuf[treeI].rootInd];
	if(2 != root.numChild){
		return;
	}

	if(UINT_MAX == MB.nodeBuf[root.firstChInd].firstChInd){//first child is a leaf
		root.collapseChild(root.lastChInd);//collapse the second child
	}else if(UINT_MAX == MB.nodeBuf[root.lastChInd].firstChInd){//second child is a leaf
		root.collapseChild(root.firstChInd);//collapse the first child
	}else{// - both children have at least two children - no need to collapse
		return;
	}
	--MB.treeBuf[treeI].numNodes;
}