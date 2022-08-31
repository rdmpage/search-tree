#include"ClusterEntities.h"

void generateClusterTrees()
{
	ClusterInfo * arrClusterInfo = new ClusterInfo[hundredThousand];
	TmpBuff tmpBuff;
	parseClusters(arrClusterInfo,tmpBuff);
	gGlobalBuf.saveCurrBuf();
	unsigned int *tmpQueue = new unsigned int[tenThousand];//tmp queue
	unsigned int (*sNodeArr)[5] = new unsigned int[tenThousand][5];//tmp source nodes array
	unsigned int (*tNodeArr)[2] = new unsigned int[tenThousand][2];//tmp target nodes array
	ClusterPptInfo* clusAddrArr = new ClusterPptInfo[hundredThousand];
	unsigned short clusterFileCount=0;
	unsigned int currFileSize=0;
	ofstream *clusFile = NULL;
	{//to localize visibity
		char clusterFileName[257] ;
		sprintf(clusterFileName,"./files/clusters/binary/Cluster%d.db",clusterFileCount);
		clusFile  = new ofstream(clusterFileName,ios::out | ios::binary);
		//open binary cluster file for output
		if(!clusFile->is_open())
		{
			T_LOG<<"\nError!clusFile could not be created";
			return;
		}
	}

	for(unsigned int i=0; i<hundredThousand; ++i)
	{
		if(currFileSize>maxFileSize)
		{//open new file
			clusFile->close();
			delete clusFile;
			char clusterFileName[257] ;
			sprintf(clusterFileName,"./files/clusters/binary/Cluster%d.db",++clusterFileCount);
			clusFile  = new ofstream(clusterFileName,ios::out | ios::binary);
			//open binary cluster file for output
			if(!clusFile->is_open())
			{
				T_LOG<<"\nError!clusFile could not be created";
				return;
			}
			currFileSize=0;
		}
		clusAddrArr[i].taxaCount = arrClusterInfo[i].numOfTaxa;
		clusAddrArr[i].fileId = clusterFileCount;
		clusAddrArr[i].addr = currFileSize;

		generateSPRTrees(arrClusterInfo[i],tmpQueue,sNodeArr,tNodeArr);
		arrClusterInfo[i].writeToFile(*clusFile);
		currFileSize+= arrClusterInfo[i].sizeAsFile();
	}
	clusFile->close();
	delete clusFile;

	char clusterAddrFileName[257] ;
	sprintf(clusterAddrFileName,"./files/ClusterMap.db");
	ofstream clusAddrFile(clusterAddrFileName,ios::out | ios::binary);
	//open binary cluster file for output
	if(!clusAddrFile.is_open())
	{
		T_LOG<<"\nError!clusFile could not be created";
		return;
	}
	
	unsigned int numOfClus = hundredThousand;
	clusAddrFile.write(reinterpret_cast<char *>(&numOfClus),sizeof(unsigned int));//write number of clusters
	for(unsigned int i=0; i<hundredThousand; ++i)
	{
		clusAddrFile.write(reinterpret_cast<char *>(&clusAddrArr[i]),sizeof(ClusterPptInfo)); 
	}
	clusAddrFile.close();

	delete [] tNodeArr;
	delete [] sNodeArr;
	delete [] tmpQueue;
	delete [] clusAddrArr;
	delete [] arrClusterInfo;
}

void generateSPRTrees(ClusterInfo &clusterInfo,unsigned int *tmpQueue,unsigned int (*sNodeArr)[5],unsigned int (*tNodeArr)[2])
{
	static const unsigned char NUM_OF_TREES =100;//num of trees this cluster has, init to default value of 100
	unsigned int sPrI_n,sNdI_n,sChI_n,tNdI_n; //node indices for randomly found nodes
	unsigned int sPrI_c,sNdI_c,tNdI_c;//child indices of above nodes	
	//s:source,t:target
	//Pr:Parent,Nd:Node,Ch:Child
	//I:index,_n:in gNodeBuf,_c:in gChildBuf
	//sPrI_n:source parent's index in gNodeBuf,
	//sPrI_c:soucr parent's child's index in gChildBuf

	unsigned int sNodeC=0,tNodeC=0;//keep count of # of elements in sNodeArr,tNodeArr
	unsigned int j,k,l;//used as index in for loop
	unsigned short queueF=0,queueB; //pointers to frong and back tmpQueue
	unsigned int tempI,tempI2,temp3;//temp indices/vars
	gGlobalBuf.check();//check gGlobalBuf for the threshold limits of its buffer
	clusterInfo.numOfTrees = NUM_OF_TREES;
	PhyloTree * clusTree = &clusterInfo.clusterAsTree;//for ease of access
	PhyloTree * currTree;//points to curr tree being generated

	clusterInfo.treeInd = gGlobalBuf.treeBufInd;//init treeInd
	gGlobalBuf.treeBufInd+=clusterInfo.numOfTrees;//update treeBufInd for the next available slot

	if(clusTree->numOfNodes <= clusterInfo.numOfTaxa+2)//trivial case of a star-like tree or having only one intermediate node
	{
		for(unsigned short i=0;i<clusterInfo.numOfTrees;++i)
		{
			gTreeBuf[clusterInfo.treeInd+i] = clusterInfo.clusterAsTree;
		}
		return;
	}
	
	tempI = clusTree->nodeInd+ clusTree->numOfNodes;//index to the one extra node at the end which will act as the new node
	gNodeBuf[tempI].numOfChild =2;//new node will have 2 children
	gNodeBuf[tempI].childInd = gGlobalBuf.childBufInd;//init childInd
	gGlobalBuf.childBufInd+= gNodeBuf[tempI].numOfChild;//to point to next free slot

	for(j=0; j< clusTree->numOfNodes; ++j)
	{//get all the candidate source and target nodes
		tNdI_n= sPrI_n = clusTree->nodeInd+j;
		if(UINT_MAX != gNodeBuf[sPrI_n].childInd)
		{//not a leaft node
			for(k=0;k<gNodeBuf[sPrI_n].numOfChild;++k)
			{
				//for source node
				sPrI_c =  gNodeBuf[sPrI_n].childInd+k;
				sNdI_n = clusTree->nodeInd+gChildBuf[sPrI_c];
				if(UINT_MAX != gNodeBuf[sNdI_n].childInd)
				{//not a leaf node
					for(l=0;l<gNodeBuf[sNdI_n].numOfChild;++l)
					{
						sNdI_c = gNodeBuf[sNdI_n].childInd + l;
						sChI_n = clusTree->nodeInd + gChildBuf[sNdI_c];
						sNodeArr[sNodeC][0] = sPrI_n;
						sNodeArr[sNodeC][1] = sPrI_c;
						sNodeArr[sNodeC][2] = sNdI_n;
						sNodeArr[sNodeC][3] = sNdI_c;
						sNodeArr[sNodeC][4] = sChI_n;
						++sNodeC;
					}
				}

				//for target node
				tNdI_c = gNodeBuf[tNdI_n].childInd + k;
				tNodeArr[tNodeC][0] = tNdI_n;
				tNodeArr[tNodeC][1] = tNdI_c;
				++tNodeC;				
			}
		}
	}
		

	for(unsigned short i=0;i<clusterInfo.numOfTrees;++i)
	{
		currTree = &gTreeBuf[clusterInfo.treeInd+i];//for ease of access
		//select a source node
		j= rand()%sNodeC;
		sPrI_n =sNodeArr[j][0];
		sPrI_c =sNodeArr[j][1];
		sNdI_n =sNodeArr[j][2];
		sNdI_c =sNodeArr[j][3];
		sChI_n =sNodeArr[j][4];
		
		//correponding to this some valid target node must exist
		//find a valid target node randomly
		while(true)
		{	
			j= rand()%tNodeC;
			tNdI_n = tNodeArr[j][0];
			tNdI_c = tNodeArr[j][1];
			//check if target node is not in the subtree of source node
			queueF=0;
			queueB=0;
			tmpQueue[queueB++] = sNdI_n;
			do
			{
				tempI = tmpQueue[queueF++];//pop element at front
				if(tempI == tNdI_n)
				{//same as target node
					break;
				}
				if(UINT_MAX !=gNodeBuf[tempI].childInd)
				{//not a leaf node
					j=gNodeBuf[tempI].childInd;
					k=j+gNodeBuf[tempI].numOfChild;
					for(;j<k;++j)
					{//put each child in queue
						tmpQueue[queueB++] = clusTree->nodeInd+gChildBuf[j];
					}
				}				
			}while(queueF !=queueB);//i.e. queue is empty
			if(tempI != tNdI_n)
			{//i.e. tNdI_n was never found in the subtree
				break;//successful target node found
			}
		}

		//valid source and target node found, now do the SPR operation
		//operation on source node side
		if(2 == gNodeBuf[sNdI_n].numOfChild)
		{//if only 2 children, then source node would not exist
			gChildBuf[sPrI_c] = gChildBuf[sNdI_c];//by passing the souce node as it is no longer needed
			tempI2 = gNodeBuf[sNdI_n].childInd + 1-(sNdI_c-gNodeBuf[sNdI_n].childInd);//the one which should not exist
			currTree->numOfNodes = clusTree->numOfNodes;//numOfNodes unchanged, one delted and one added
		}
		else//case of more than 2 children
		{//swap the child index with the last child index and decrement the child count by 1
			tempI2 = gNodeBuf[sNdI_n].childInd+gNodeBuf[sNdI_n].numOfChild - 1;//index to last child
			//swap now
			tempI = gChildBuf[sNdI_c];
			gChildBuf[sNdI_c] = gChildBuf[tempI2];
			gChildBuf[tempI2] = (unsigned short)tempI;
			--gNodeBuf[sNdI_n].numOfChild;//decrease num of children for the source node
			currTree->numOfNodes = clusTree->numOfNodes+1;//numOfNodes increased by one as none deleted and one added
		}
		
		//operation on target node side	
		tempI = clusTree->nodeInd+ clusTree->numOfNodes;//index to the one extra node at the end which will act as the new node
		tempI = gNodeBuf[tempI].childInd;
		gChildBuf[tempI] = gChildBuf[tNdI_c];//first child- targets node's child
		gChildBuf[tempI+1] = gChildBuf[tempI2];//second child - source node's child
		gChildBuf[tNdI_c] = clusTree->numOfNodes;//target node should point to the new node
		
		//SPR done, now copy the tree
		currTree->nodeInd = gGlobalBuf.nodeBufInd;//init nodeInd
		gGlobalBuf.nodeBufInd+= currTree->numOfNodes;//update nodeBufInd to point to next free slot
		queueF = queueB =0;
		currTree->numOfNodes=0;//would be used as counter
		tmpQueue[queueB++] = clusTree->nodeInd;//init to root
		tempI2 = currTree->nodeInd;//next slot for free node
		do
		{
			tempI = tmpQueue[queueF++];//front element in the queue
			//tempI points to the source-from-copy node and tempI2 points to the target-to-copy node
			gNodeBuf[tempI2].numOfChild = gNodeBuf[tempI].numOfChild;
			if(UINT_MAX == gNodeBuf[tempI].childInd)
			{//leaf node
				gNodeBuf[tempI2].childInd = UINT_MAX;//mark target node also as leaf node
			}
			else
			{
				temp3=gNodeBuf[tempI2].childInd = gGlobalBuf.childBufInd;//init childInd
				gGlobalBuf.childBufInd+= gNodeBuf[tempI2].numOfChild;//for next free slot
				j=gNodeBuf[tempI].childInd;
				k=j+gNodeBuf[tempI].numOfChild;
				for(;j<k;++j)
				{//copy each children of source node 
					gChildBuf[temp3++] = ++currTree->numOfNodes;//assign space for the child node
					tmpQueue[queueB++] = clusTree->nodeInd+gChildBuf[j];//store the children in the stack
				}
			}
			++tempI2;//to point to the node corresponding to the next element in front of the queue
		}while(queueF!=queueB);//i.e. empty queue
		++currTree->numOfNodes;//since it is used as a pre-increment counter which points to the next avail-free node, needs to be incremented by 1 
		//to reflect the actual count

		//tree created, now clusTree changes need to be reverted
		tempI = clusTree->nodeInd+ clusTree->numOfNodes;//index to the last node which is the new node
		tempI = gNodeBuf[tempI].childInd;
		gChildBuf[tNdI_c] = gChildBuf[tempI];//restore child of the target node
		if(currTree->numOfNodes == clusTree->numOfNodes)
		{//a node was deleted
			gChildBuf[sPrI_c] = (unsigned short)(sNdI_n - clusTree->nodeInd);//restore the bypass
		}
		else
		{
			tempI = gChildBuf[sNdI_c];
			tempI2 = gNodeBuf[sNdI_n].numOfChild +gNodeBuf[sNdI_n].childInd;
			//reverse the swap
			gChildBuf[sNdI_c] = gChildBuf[tempI2];
			gChildBuf[tempI2] = (unsigned short)tempI;
			++gNodeBuf[sNdI_n].numOfChild;//restore the count by increasing by 1
		}
	}
}

void parseClusters(ClusterInfo * arrClusterInfo,TmpBuff &tmpBuff)
{
	//TmpBuff tmpBuff;
	tmpBuff.symArr = new unsigned short[fiftyThousand];
	tmpBuff.taxaArr = new unsigned int[fiftyThousand];
	tmpBuff.tmpStack = new unsigned short[fiftyThousand];

	//Load mapSeqIdGenId from file
	ifstream genSeqOrgIdArrFile (genSeqOrgIdArrFileName,ios::in | ios::binary);
	//open seqId to uniqueSeqId map file to read
	if(genSeqOrgIdArrFile.is_open())
	{
		unsigned int seqId,numSeqId;
		genSeqOrgIdArrFile.read(reinterpret_cast<char *>(&numSeqId),sizeof(unsigned int));
		for(unsigned int i=0; i<numSeqId; ++i)
		{
			genSeqOrgIdArrFile.read(reinterpret_cast<char *>(&seqId),sizeof(unsigned int));
			tmpBuff.mapSeqIdGenId.insert(seqId,i);
		}
	}
	genSeqOrgIdArrFile.close();//seqId to uniqueSeqId map file read

	string fileName = "./files/GB159.PB.newick.11.30.07";
	ifstream dataFile(fileName.c_str());
	string line ;
	string clusterId;
	unsigned int uniqIdGen =0;
	unsigned int numOfClus = hundredThousand;
	if(dataFile.is_open())
	{
		for(unsigned int i=0; i<numOfClus; ++i)
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
			getline (dataFile,line);
			if(line.length() == 0)
				{
					--i;
					continue;
				}
			parseNewick(line,arrClusterInfo[i],clusterId,tmpBuff);
			arrClusterInfo[i].clusterId = uniqIdGen++;
			//gGlobalBuf.check(); -shoud not be needed, it it is actually needed, then it will overwirte the existing clusters.
		}
	}
}

void parseNewick(string &strTree,ClusterInfo & clusterInfo, string &strClusterId, TmpBuff & tmpBuff)
{
	static istringstream stm;
	static size_t pos1, pos2;
	static char ch;//tmp char
	static unsigned short count ; //to store the bracket count
	count =0;//init bracket count 
	static unsigned short stackCounter;//counter which points to the last stored element in tmpStack
	stackCounter =0; //init stackCounter
	static unsigned short symC;// to index next avail pos in TmpBuff.symArr
	symC=0;//init symC
	
	static unsigned int tmpIndex,tmpIndex2; //to tmp store the index value where the index expression gets longish

	PhyloTree *clusTree = &clusterInfo.clusterAsTree;//ptr for shorter reference
	
	pos1 = strTree.find_first_of(" \t");//skip clus type - single loci, biclique, quasi-biclique
	pos1 = strTree.find_first_of(" \t",pos1+1);//skip phylota-id
	pos2 = strTree.find_first_of(" \t",pos1+1);//skip base dir: 0-26	
	strClusterId = strTree.substr(pos1+1,pos2-pos1-1);//base dir
	pos1=pos2;
	pos2 = strTree.find_first_of(" \t",pos1+1);//skip clusFileName
	//strClusterId = phyloTreeConstants::suffixString.substr(0,20 - (pos1 + 1) ) + strTree.substr(0,pos1+1);//get string tree id
	strClusterId += "/"+ strTree.substr(pos1+1,pos2-pos1-1) ;
	//parse the strring based tree id
	pos1 = pos2;

	pos1 = strTree.find_first_not_of(" \t",pos1);

	static unsigned int maxNumOfNodes=0,maxSymC=0,maxTaxaC=0;
	//PASS ONE
	//first pass over the tree string is done to assert num of nodes it has and num of taxa and store the taxaIds
	//also the positions of relevant symbols is stored so that string need not be parsed again in the second pass
	//NOTE: the way newick tree is written the traversal is going to be in pre-order i.e root, left child, right child
	do
	{
		
		ch = strTree.at(pos1);
		if('(' == ch)
		{//new intermediate node
			++clusTree->numOfNodes;
			++count;//bracket count
			tmpBuff.symArr[symC++]= 0;//and the symbol it refers to, 0 means '('
			++pos1;
		}
		else if(')' == ch )
		{
			--count;//bracket count
			tmpBuff.symArr[symC++]= 1;//and the symbol it refers to, 1 means ')'
			++pos1;
		}
		else if(',' == ch)
		{
			//do nothing
			//tmpBuff.symArr[symC++]= 2;//and the symbol it refers to, 2 means ','
			++pos1;
		}
		else //taxaId found
		{//new leaf node
			tmpBuff.symArr[symC++]= 3;//and the symbol it refers to, 3 means a taxa Id
			
			++clusTree->numOfNodes;
			pos2 = strTree.find_first_not_of("0123456789",pos1);
			stm.str(strTree.substr(pos1,pos2 - pos1));
			stm >> tmpBuff.taxaArr[clusterInfo.numOfTaxa++];
			stm.clear();
			pos1 = pos2;//skip taxa id for next pos
		}
	}
	while(0 != count);//tree has been parsed
	if(fiftyThousand <= clusTree->numOfNodes || 
		fiftyThousand <= clusterInfo.numOfTaxa ||
		fiftyThousand <= symC)
		{
			T_LOG<<"\nERROR! numOfNodes needs to be increased!";
			T_LOG<<"\n"<<symC<<" "<<clusTree->numOfNodes<<" "<<clusterInfo.numOfTaxa;
		}

	bool maxValuesChanged = false;
	if(maxNumOfNodes<clusTree->numOfNodes)
	{
		maxNumOfNodes=clusTree->numOfNodes;
		maxValuesChanged = true;
	}

	if(maxSymC<symC)
	{
		maxSymC = symC;
		maxValuesChanged = true;
	}
	if(maxTaxaC < clusterInfo.numOfTaxa)
	{
		maxTaxaC = clusterInfo.numOfTaxa;
		maxValuesChanged = true;
	}
	if(maxValuesChanged)
	{
		T_LOG<<"\n"<<maxSymC<<" "<<maxNumOfNodes<<" "<<maxTaxaC;
	}

	clusTree->nodeInd = gGlobalBuf.nodeBufInd; //init nodeInd
	gGlobalBuf.nodeBufInd+= clusTree->numOfNodes+1;//update nodeBufInd for the next free slot
	//the +1 is kept extra as this extra node is required in SPR tree generation
	clusterInfo.taxaInd = gGlobalBuf.taxaBufInd;//init taxaInd
	gGlobalBuf.taxaBufInd+= clusterInfo.numOfTaxa;//update taxaBufInd for the next free slot
	
	unsigned int seqId;
	for(count =0; count<clusterInfo.numOfTaxa; ++count)
	{
		seqId = tmpBuff.mapSeqIdGenId.get(tmpBuff.taxaArr[count]);
		if(UINT_MAX != seqId)
		{
			gTaxaBuf[clusterInfo.taxaInd+count] = seqId; //update with the genTaxaId
		}
		else
		{
			T_LOG<<"\n ERROR! tmpBuff.mapSeqIdGenId.find failed for:\t"<<tmpBuff.taxaArr[count];
			gTaxaBuf[clusterInfo.taxaInd+count] = 0;
		}
	}

	//PASS TWO
	//second pass stores the children count for each node 
	clusTree->numOfNodes=0;//init to 0, would be used as pointer to the position of curr new node in clusTree->nodeArr
	stackCounter =0;//init stackCounter
	tmpBuff.tmpStack[stackCounter] = clusTree->numOfNodes++; //init to root
	for(count =1; count<symC; ++count)
		//count =1 to skip the opening bracket for root, as root can not have any parent hence need to be considered separately
	{
		switch(tmpBuff.symArr[count])
		{
		case 0: // '('
			//increment child count for the parent at the top of the stack
			++gNodeBuf[clusTree->nodeInd+tmpBuff.tmpStack[stackCounter]].numOfChild;
			tmpBuff.tmpStack[++stackCounter] = clusTree->numOfNodes++;//update tmpStack for the new node
			break;
		case 1: //')'
			--stackCounter;//move up one level as subtree for curr node has been parsed
			break;
		//case 2: //','
		//	//do nothing
		//	break;
		default: // a taxa id
			++clusTree->numOfNodes;//increment node count
			//increment child count for the parent at the top of the stack
			++gNodeBuf[clusTree->nodeInd+tmpBuff.tmpStack[stackCounter]].numOfChild;
		}
	}

	//PASS THREE
	//third assigns memory to childArr for each of the intermediate nodes as well updates pointers to those children in the childArr
	clusTree->numOfNodes=0;//init to 0, would be used as pointer to the position of curr new node in clusTree->nodeArr
	clusterInfo.numOfTaxa=0;//||rly as above, would point to pos for curr new taxa

	stackCounter =0;//init stackCounter
	tmpBuff.tmpStack[stackCounter] = clusTree->numOfNodes++; //init to root
	//as root can not have any parent hence need to be considered separately
	gNodeBuf[clusTree->nodeInd+0].childInd = gGlobalBuf.childBufInd; //init childInd
	gGlobalBuf.childBufInd+= gNodeBuf[clusTree->nodeInd+0].numOfChild;//update childBufInd to point to next free slot

	for(count =1; count<symC; ++count)
		//count =1 to skip the opening bracket for root
	{
		switch(tmpBuff.symArr[count])
		{
		case 0: // '('
			//a new child found, update it in the childArr for the its parent which is the node at the top of the stack
			tmpIndex = clusTree->nodeInd+tmpBuff.tmpStack[stackCounter];//index for the parent
			tmpIndex2 = gNodeBuf[tmpIndex].childInd;//childInd for the parent
			gChildBuf[tmpIndex2+ --gNodeBuf[tmpIndex].numOfChild] = clusTree->numOfNodes;
			gNodeBuf[clusTree->nodeInd + clusTree->numOfNodes].childInd = gGlobalBuf.childBufInd; //init childInd for the new node
			gGlobalBuf.childBufInd+=gNodeBuf[clusTree->nodeInd+clusTree->numOfNodes].numOfChild;//update childBufInd for the next avail slot
			tmpBuff.tmpStack[++stackCounter] = clusTree->numOfNodes++;//update tmpStack for the new node and increment numOfNodes for the next new node
			break;
		case 1: //')'
			--stackCounter;//move up one level as subtree for curr node has been parsed
			break;
		//case 2: //','
		//	//do nothing
		//	break;
		default: // a taxa id
			gNodeBuf[clusTree->nodeInd+clusTree->numOfNodes].numOfChild = clusterInfo.numOfTaxa++; //here numOfChild is acting as taxaId,
				//increment taxaCount for next taxa
			//a new child found, update it in the childArr for the its parent which is the node at the top of the stack
			tmpIndex = clusTree->nodeInd+tmpBuff.tmpStack[stackCounter];//index for the parent
			tmpIndex2 = gNodeBuf[tmpIndex].childInd;//childInd for the parent
			gChildBuf[tmpIndex2+ --gNodeBuf[tmpIndex].numOfChild] = clusTree->numOfNodes++;//also update numOfChild of parent for next child
		}
	}

	//PASS FOUR : to recount num of children for each node. its a bit sad that i had to do this fourth pass, but i have to do. with some tweak 
	//this fourth pass can be avoided but it will make a very little difference over all and i am not in a mood to use my brain for that tweak as i
	// feel this is already quite a bit complicated. At the same time I do want to introduce any container in tree structure as the intention is to keep 
	//its memory footprint to minimum
	//fourth pass stores the children count for each node
	clusTree->numOfNodes=0;//init to 0, would be used as pointer to the position of curr new node in clusTree->nodeArr
	stackCounter =0;//init stackCounter
	tmpBuff.tmpStack[stackCounter] = clusTree->numOfNodes++; //init to root
	for(count =1; count<symC; ++count)
		//count =1 to skip the opening bracket for root, as root can not have any parent hence need to be considered separately
	{
		switch(tmpBuff.symArr[count])
		{
		case 0: // '('
			//increment child count for the parent at the top of the stack
			++gNodeBuf[clusTree->nodeInd+tmpBuff.tmpStack[stackCounter]].numOfChild;
			tmpBuff.tmpStack[++stackCounter] = clusTree->numOfNodes++;//update tmpStack for the new node
			break;
		case 1: //')'
			--stackCounter;//move up one level as subtree for curr node has been parsed
			break;
		//case 2: //','
		//	//do nothing
		//	break;
		default: // a taxa id
			++clusTree->numOfNodes;//increment node count
			//increment child count for the parent at the top of the stack
			++gNodeBuf[clusTree->nodeInd + tmpBuff.tmpStack[stackCounter]].numOfChild;
		}
	}
	
}
