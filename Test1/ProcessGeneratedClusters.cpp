#include"ClusterEntities.h"


void processGeneratedClusters()
{
	ClusterInfo clus;
	TmpBuff tmpBuff;//needed to parse a newick tree and to store taxaId-sequence id mapping
	tmpBuff.symArr = new unsigned short[fiftyThousand];
	tmpBuff.taxaArr = new unsigned int[fiftyThousand];
	tmpBuff.blArr = new unsigned short[fiftyThousand];
	tmpBuff.tmpStack = new unsigned short[fiftyThousand];
	unsigned short * taxaArr = new unsigned short[maxSeqCount];//for every clusters, stores the the index position of every taxa in gTaxaBuf-
	// -relative ClusterInfo.taxaInd

	ClusterPptInfo* clusPptInfroArr = new ClusterPptInfo[tenMillion];//to store cluster Addresses


	unsigned int clusC=0;

#ifdef WIN32
	string bootStrapDir ="G:/SearchTree_Dec11_2010_Trees/live/parsed";//source dir form where the clusters are being read from 
#else
	string bootStrapDir ="./../../parsedblTrees";//source dir form where the clusters are being read from 
#endif


	string line ;//to read a line from file
	//string clusterId;//str based string id read from the cluster file
	unsigned int uniqIdGen =0;

	//Load mapSeqIdGenIds from file
	ifstream genSeqOrgIdArrFile (genSeqOrgIdArrFileName,ios::in | ios::binary);
	//open taxa id to unique id map file to read
	if(genSeqOrgIdArrFile.is_open())
	{
		unsigned int seqId,numSeq;
		genSeqOrgIdArrFile.read(reinterpret_cast<char *>(&numSeq),sizeof(unsigned int));
		for(unsigned int i=0; i<numSeq; ++i)
		{
			genSeqOrgIdArrFile.read(reinterpret_cast<char *>(&seqId),sizeof(unsigned int));
			tmpBuff.mapSeqIdGenId.insert(seqId,i);
		}
	}
	genSeqOrgIdArrFile.close();//mapSeqIdGenId file read

	unsigned short clusterFileCount=0;//count for the source(input) cluster files
	unsigned int currFileSize=0;
	ofstream *clusFile_w = NULL;//o/p binary cluster file where is the cluster are written
	{//to localize visibity
		char clusterFileName_w[257] ;
		sprintf(clusterFileName_w,"./files/clusters/binary/Cluster%d.db",clusterFileCount);
		clusFile_w  = new ofstream(clusterFileName_w,ios::out | ios::binary);
		//open binary cluster file for output
		if(!clusFile_w->is_open())
		{
			T_LOG<<"\nError!clusFile could not be created";
			return;
		}
	}

	unsigned short clusDbFileCount=0;
	size_t pos1, pos2;

	while(true)
	{
		ifstream clusFile;//input cluster file
		string clusId;//cluster id is fetched into this

		char clusterFileName[257] ;
		sprintf(clusterFileName,"./files/clusters/newick/Tree%d",clusDbFileCount);
		clusFile.open(clusterFileName);
		//open source cluster file for input
		if(!clusFile.is_open())
		{
			T_LOG<<"Cluster File Stopped at:"<<clusDbFileCount<<std::endl;
			break;
		}
		++clusDbFileCount;//increment count for the source(input) cluster files 
		while(!clusFile.eof())
		{
			getline (clusFile,line);
			if(line.length() == 0)
			{
				continue;
			}
			parseNewick(line,clus,clusId,tmpBuff);
			clus.clusterId = uniqIdGen++;
			//extract clus-type and phylotaId from clusId
			pos1= line.find_first_of(" \t");//skip clusType - sinle loci, biclique or quasi-biclique
			clusPptInfroArr[clusC].clusType = atoi(line.substr(0,pos1).c_str());
			pos2= line.find_first_of(" \t", pos1+1);
			clusPptInfroArr[clusC].phylotaId = atoi(line.substr(pos1+1,pos2-pos1-1).c_str());
			
			

			readClusterTrees(clus,tmpBuff,bootStrapDir,clusId,taxaArr);//reads the trees corresponding to the cluster and loads them into -
			//- memory in ClusterInfo-clus

			if(currFileSize>maxFileSize)
			{//open new file
				clusFile_w->close();
				delete clusFile_w;
				char clusterFileName_w[257] ;
				sprintf(clusterFileName_w,"./files/clusters/binary/Cluster%d.db",++clusterFileCount);
				clusFile_w  = new ofstream(clusterFileName_w,ios::out | ios::binary);
				//open binary cluster file for output
				if(!clusFile_w->is_open())
				{
					T_LOG<<"\nError!clusFile could not be created";
					return;
				}
				currFileSize=0;
			}

			//fill in the data req for the address of this cluster
			clusPptInfroArr[clusC].taxaCount = clus.numOfTaxa;
			clusPptInfroArr[clusC].fileId = clusterFileCount;
			clusPptInfroArr[clusC].addr = currFileSize;
			++clusC;

			clus.writeToFile(*clusFile_w);
			currFileSize+=clus.sizeAsFile(); 
			clus.reset();//init for next use

			gGlobalBuf.check();
		}
		clusFile.close();//close the input cluster file
	}

	clusFile_w->close();//close the ouput binary cluster file
	delete clusFile_w;

	//the write the address info for all the clusters
	char clusterAddrFileName[257] ;
	sprintf(clusterAddrFileName,"./files/ClusterMap.db");
	ofstream clusAddrFile(clusterAddrFileName,ios::out | ios::binary);
	//open binary cluster file for output
	if(!clusAddrFile.is_open())
	{
		T_LOG<<"\nError!clusFile could not be created";
		return;
	}

	clusAddrFile.write(reinterpret_cast<char *>(&clusC),sizeof(unsigned int)); //write number of clusters
	for(unsigned int i=0; i<clusC; ++i)
	{
		clusAddrFile.write(reinterpret_cast<char *>(&clusPptInfroArr[i]),sizeof(ClusterPptInfo)); 
	}
	clusAddrFile.close();


	delete []clusPptInfroArr;
	delete []taxaArr;
}

//for the ClusterInfo-clus, reads all its trees into memory
void readClusterTrees(ClusterInfo &clus,TmpBuff &tmpBuff,string &bootStrapDir,string &clusFileName,unsigned short *taxaArr)
{
	static unsigned int clusCount =0;
	//stuff start
	//below stuff is needed as we need to trim the string which has some unix based characters too
	unsigned int _i,_j;
	char clusterFileName[513] ;
	for(_i=0;_i<bootStrapDir.length();++_i)
	{
		clusterFileName[_i] = bootStrapDir[_i];
	}
	clusterFileName[_i++]='/';
	for(_j=0;_j<clusFileName.length();++_j)
	{
		if(isspace(clusFileName[_j])==0)
		{
			clusterFileName[_i++]=clusFileName[_j];
		}
	}
	clusterFileName[_i]='\0';
	//stuff ends

	ifstream clusFile(clusterFileName);

	++clusCount;
	if(clusCount%1000 == 0)
	{
		T_LOG<<"\n"<<clusterFileName<<"| Count:"<<clusCount<<"\n";
	}

	istringstream stm;

	unsigned int numOfTrees;//number 

	string line;

	if(clusFile.is_open())
	{
		getline (clusFile,line);
		stm.str(line);
		stm>>numOfTrees;//number of boot strapped trees in the cluster
	}

	unsigned int taxaId;
	for(unsigned short i=0; i<clus.numOfTaxa; ++i)
	{//taxaArr,for every clusters, stores the the index position of every taxa in gTaxaBuf-
		// -relative ClusterInfo.taxaInd
		taxaId=gTaxaBuf[clus.taxaInd+i];
		taxaArr[taxaId] = i;
	}
	clus.numOfTrees = (unsigned short)numOfTrees;
	clus.treeInd = gGlobalBuf.treeBufInd;
	gGlobalBuf.treeBufInd+=clus.numOfTrees;

	for(unsigned short i=0; i<numOfTrees; ++i)
	{
		getline (clusFile,line);
		parseClusterTreeNewick(line,&gTreeBuf[clus.treeInd+i],tmpBuff,taxaArr);
	}
}

void parseClusterTreeNewick(string &strTree,PhyloTree * clusTree,TmpBuff & tmpBuff,unsigned short *taxaArr)
{
	static istringstream stm;
	static size_t pos1, pos2;
	static char ch;//tmp char
	static unsigned short count ; //to store the bracket count
	count =0;//init bracket count 
	static unsigned short taxaC=0;//to store taxa count
	taxaC=0;//init taxa count
	static unsigned short stackCounter;//counter which points to the last stored element in tmpStack
	stackCounter =0; //init stackCounter
	static unsigned short symC;// to index next avail pos in TmpBuff.symArr
	symC=0;//inti symC
	
	unsigned short blInd=0;//indexes the branch length values in the order in which the respective nodes were encountered


	static unsigned int tmpIndex,tmpIndex2,tmpIndex3; //to store tmp store the index value where the index expression gets longish

	pos1 = 0;//init pos1 to start of the string, only tree is there no treeId

	static unsigned int maxNumOfNodes=0,maxSymC=0;
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
			//no need to store ')' as ':' the new delimiter
			//tmpBuff.symArr[symC++]= 1;//and the symbol it refers to, 1 means ')'
			++pos1;
		}
		else if(',' == ch)
		{
			//do nothing
			//tmpBuff.symArr[symC++]= 2;//and the symbol it refers to, 2 means ','
			++pos1;
		}
		else if(':' == ch)
		{//branch lenght found
			tmpBuff.symArr[symC++]= 4;//and the symbol it refers to, 4 means ':'
			++pos1;
			pos2 = strTree.find_first_not_of("0123456789",pos1);
			stm.str(strTree.substr(pos1,pos2 - pos1));
			stm >> tmpBuff.blArr[blInd++];
			stm.clear();
			pos1 = pos2;//skip branch length for next pos
		}
		else //taxaId found
		{//new leaf node
			tmpBuff.symArr[symC++]= 3;//and the symbol it refers to, 3 means a taxa Id

			++clusTree->numOfNodes;
			pos2 = strTree.find_first_not_of("0123456789",pos1);
			stm.str(strTree.substr(pos1,pos2 - pos1));
			stm >> tmpBuff.taxaArr[taxaC++];
			stm.clear();
			pos1 = pos2;//skip taxa id for next pos
		}
	}
	while(0 != count);//tree has been parsed

	if(fiftyThousand <= clusTree->numOfNodes ||
		fiftyThousand <= symC)
	{
		T_LOG<<"\nERROR! numOfNodes needs to be increased!";
		T_LOG<<"\n"<<symC<<" "<<clusTree->numOfNodes;
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

	if(maxValuesChanged){
		T_LOG<<"\n"<<maxSymC<<" "<<maxNumOfNodes;
	}
	



	clusTree->nodeInd = gGlobalBuf.nodeBufInd; //init nodeInd
	gGlobalBuf.nodeBufInd+= clusTree->numOfNodes;//update nodeBufInd for the next free slot

	unsigned int seqIdGenId;

	for(count =0; count<taxaC; ++count)
	{
		seqIdGenId = tmpBuff.mapSeqIdGenId.get(tmpBuff.taxaArr[count]);
		if(UINT_MAX != seqIdGenId)
		{
			tmpBuff.taxaArr[count] = taxaArr[seqIdGenId]; //update with the genTaxaId			
		}
		else
		{
			T_LOG<<"\n ERROR! tmpBuff.mapSeqIdGenId.find failed for:\t"<<tmpBuff.taxaArr[count];			
			tmpBuff.taxaArr[count] = 0;			
		}
	}

	//PASS TWO
	//second pass stores the children count for each node and assigns bl 
	blInd=0;//reset the counter for branch length values
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
		//case 1: //')', ':' is the new delimiter for a branch
		//	--stackCounter;//move up one level as subtree for curr node has been parsed
		//	break;
			//case 2: //','
			//	//do nothing
			//	break;
		case 4: // ':' is the new delimiter for a branch as well as for taxa
			gNodeBuf[clusTree->nodeInd+tmpBuff.tmpStack[stackCounter]].bl= tmpBuff.blArr[blInd++];
			//bl is being added above as bl is being stored as the distance from the root
			--stackCounter;//move up one level as subtree for curr node has been parsed
			break;
		default: // a taxa id
			//increment child count for the parent at the top of the stack
			++gNodeBuf[clusTree->nodeInd+tmpBuff.tmpStack[stackCounter]].numOfChild;
			tmpBuff.tmpStack[++stackCounter] = clusTree->numOfNodes++;//update tmpStack for the new node			
		}
	}
	gNodeBuf[clusTree->nodeInd].bl =0;//set bl for the root to 0
	--stackCounter;//not needed for the root node, but just to make it look neat and in synch with the above for loop



	//PASS THREE
	//third assigns memory to childArr for each of the intermediate nodes as well updates pointers to those children in the childArr
	//it also calculates cumatlative bl
	clusTree->numOfNodes=0;//init to 0, would be used as pointer to the position of curr new node in clusTree->nodeArr
	taxaC=0;//||rly as above, would point to pos for curr new taxa

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
			tmpIndex3 = clusTree->nodeInd + clusTree->numOfNodes;//index for this child
			gNodeBuf[tmpIndex3].childInd = gGlobalBuf.childBufInd; //init childInd for the new node
			gGlobalBuf.childBufInd+=gNodeBuf[tmpIndex3].numOfChild;//update childBufInd for the next avail slot
			gNodeBuf[tmpIndex3].bl+= gNodeBuf[tmpIndex].bl;//add bl of its parent to calulate cumalative bl


			tmpBuff.tmpStack[++stackCounter] = clusTree->numOfNodes++;//update tmpStack for the new node and increment numOfNodes for the next new node
			break;
			//case 1: //')', ':' is the new delimiter
			//	--stackCounter;//move up one level as subtree for curr node has been parsed
			//	break;
			//case 2: //','
			//	//do nothing
			//	break;
		case 4: //':' is the new delimiter
			--stackCounter;//move up one level as subtree for curr node has been parsed
			break;
		default: // a taxa id
			gNodeBuf[clusTree->nodeInd+clusTree->numOfNodes].numOfChild = (unsigned short)tmpBuff.taxaArr[taxaC++]; //here numOfChild is acting as taxaId,
			//increment taxaCount for next taxa
			//a new child found, update it in the childArr for the its parent which is the node at the top of the stack
			tmpIndex = clusTree->nodeInd+tmpBuff.tmpStack[stackCounter];//index for the parent
			tmpIndex2 = gNodeBuf[tmpIndex].childInd;//childInd for the parent
			tmpIndex3 = clusTree->nodeInd + clusTree->numOfNodes;//index for this child
			gChildBuf[tmpIndex2+ --gNodeBuf[tmpIndex].numOfChild] = clusTree->numOfNodes++;//also update numOfChild of parent for next child
			gNodeBuf[tmpIndex3].bl+= gNodeBuf[tmpIndex].bl;//add bl of its parent to calulate cumalative bl
			++stackCounter;//update stack counter, as the next symbol is ':' which will decrement the stack
			}
	}
	--stackCounter;//again, not needed for the root node, but just to make it look neat and in synch with the above for loop
	

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
		//case 1: //')' ,':' is the new delimiter
		//	--stackCounter;//move up one level as subtree for curr node has been parsed
		//	break;
			//case 2: //','
			//	//do nothing
			//	break;
		case 4: //':' is the new delimiter
			--stackCounter;//move up one level as subtree for curr node has been parsed
			break;
		default: // a taxa id
			++clusTree->numOfNodes;//increment node count
			//increment child count for the parent at the top of the stack
			++gNodeBuf[clusTree->nodeInd + tmpBuff.tmpStack[stackCounter]].numOfChild;
			++stackCounter;//update stack counter, as the next symbol is ':' which will decrement the stack
			}
	}
	--stackCounter;//again, not needed for the root node, but just to make it look neat and in synch with the above for loop


}
