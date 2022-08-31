#include"MRT.h"
#include<set>
#include<cmath>

unsigned char BiPart::len;
unsigned long long BiPart::lastOfset;

//returns normalized tree score in terms of avg support number
unsigned char genMRT(ClusterInfo &clus,BiPart * biPartArr,Bin *hashTb,BinNd * binNdArr
					 ,ArrayBufs &arrBufs,MajTreeNodePtr * ndPtrBuf,unsigned short &maxSupportNum,
					 unsigned short &numOfUniqMajBiPart)
{
	unsigned short stackC;//stack counter
	unsigned int ndI,k;//tmp vars
	//int sL;//tmp vars
	unsigned short taxaId,j;//tmp vars
	unsigned int binNdInd,binNdIndCh;
	unsigned int biPartC=0;
	unsigned int binNdArrInd=0;//counter for binNdArr
	PhyloTree * curTree;//ptr to current tree;
	PhyloTree & majTree =clus.clusterAsTree;//for convenience
	unsigned int _C_ptr;//"ptr to last node corresponding to a majoirty bi-partition which is ansestor of current node in the traversal"-
	//- as per given in paper
	unsigned int childBufC=0;//init counter for childBuf
	unsigned short ndPtrBufC=0;//init counter for MajTreeNodePtr
	unsigned int mtNdInd;//ptr to node of a maj bi-part in the MajTreeNodePtr
	unsigned int parInd, chInd;//temp var
	unsigned short MAJ_LIM = clus.numOfTrees/2;
	unsigned int treeScore=0;//to caculate tree score based on support number
	numOfUniqMajBiPart=0;//to calculate number of non-trivial bi-paritions, used to normalize the treeScore as well
	unsigned short rootBL;//to store bl of the root while traversing a tree
	float tmpDl;//needed to round off and used as temp vars
	unsigned short maxSupNumSecondBest;


	unsigned long long ONE =1;//used to obtain the actual bi-parition from [0] entry of taxaArr as described above

	BiPart::len = clus.numOfTaxa/64 + 1;// init BiPart::len
	BiPart::lastOfset =0;
	for(unsigned short i= (clus.numOfTaxa/64)*64; i<clus.numOfTaxa; ++i){
		BiPart::lastOfset |= ONE<<arrBufs.taxaArr2[i][0];
	}

	maxSupportNum=0;
	maxSupNumSecondBest =0;



	//Step1- Go through nodes of each tree updating the freq count for each bi-partition
	for(unsigned short i=0; i<clus.numOfTrees; ++i)
	{
		stackC=0;
		curTree = &gTreeBuf[clus.treeInd+i];
		arrBufs.stackBuf[stackC][0] = curTree->nodeInd;//put root node on the stack
		arrBufs.stackBuf[stackC][1] = gNodeBuf[curTree->nodeInd].numOfChild;
		while(USHRT_MAX != stackC)
		{
			ndI = arrBufs.stackBuf[stackC][0];
			if(UINT_MAX == gNodeBuf[ndI].childInd)
			{//a leaf node
				taxaId = gNodeBuf[ndI].numOfChild;//gives a way to calculate index (in biPartArr) for this bi-partition in this tree 	
				if(0 == i)
				{//singleton paritions representing leaf nodes only need to be init once
					biPartArr[taxaId].reset();//reset parition before putting values
					biPartArr[taxaId].part[arrBufs.taxaArr2[taxaId][1]]= ONE<<arrBufs.taxaArr2[taxaId][0];//add this node to the parition-
					//-since the [0] entry stores the bit position at which this taxa would be represented at, hence it needs to be left shifted those many places-
					//-to get the bi-parition representation
					biPartArr[taxaId].card =1;//just this leaf node
					biPartArr[taxaId].hashCode = arrBufs.prArr[taxaId];//hash code for this leaf
					biPartArr[taxaId].binNdInd= hashTb[biPartArr[taxaId].hashCode%PRIME_S].add(&biPartArr[taxaId],binNdArr,binNdArrInd);
					//leaf node has to be a maj bi-part, hence freq count not required to check
				}
				--stackC;
			}
			else
			{//an intermediate node
				if(0 == arrBufs.stackBuf[stackC][1])
				{//all its children have been visited
					if( ndI != curTree->nodeInd)//not root parition
					{//not a root parition
						binNdInd = ((i+1)*twiceMaxClusOverlap+ ndI - curTree->nodeInd);//gives a way to calculate index for this bi-partition in this tree
						//intermediate nodes start from twoThousand onwards, the first twiceMaxClusOverlap is reserved for the leaf nodes
						//above reservation had to be done as leaf nodes are indexed based on their position in taxaArr of the cluster and it is same for all the trees
						//while the intermediate nodes are indexed based on their relative position to curTree->nodeInd
						k = gNodeBuf[ndI].childInd;

						//init bi-part
						biPartArr[binNdInd].hashCode=0;
						biPartArr[binNdInd].card=0;
						biPartArr[binNdInd].reset();

						for(j=0; j<gNodeBuf[ndI].numOfChild; ++j)
						{
							//find child's bin-node index in biPartArr
							if(UINT_MAX == gNodeBuf[curTree->nodeInd+gChildBuf[k +j]].childInd)
							{//this child is a leaf node
								binNdIndCh= gNodeBuf[curTree->nodeInd+gChildBuf[k +j]].numOfChild;//way of indexing leaf nodes
							}
							else
							{//this child is an int node
								binNdIndCh = (i+1)*twiceMaxClusOverlap+ gChildBuf[k+j];//way of indexing intermediate nodes
							}

							biPartArr[binNdInd].addPart(biPartArr[binNdIndCh]);//add this child's parition
							biPartArr[binNdInd].hashCode+= biPartArr[binNdIndCh].hashCode;//add child's hash code value
							biPartArr[binNdInd].card+=biPartArr[binNdIndCh].card;//add child's card
						}
						biPartArr[binNdInd].binNdInd = hashTb[biPartArr[binNdInd].hashCode%PRIME_S].add(&biPartArr[binNdInd],binNdArr,binNdArrInd);
						--stackC;
					}
					else
					{//a root parition, need to be init only once when i==0, i.e the first time
						if(0 !=i)
						{
							binNdInd = ((i+1)*twiceMaxClusOverlap+ ndI - curTree->nodeInd);//gives a way to calculate index for this bi-partition in this tree
							biPartArr[binNdInd].binNdInd = biPartArr[twiceMaxClusOverlap].binNdInd;//make it point to first root partition					
						}
						else
						{
							binNdInd = twiceMaxClusOverlap;//const position for root
							biPartArr[binNdInd].hashCode=0;
							biPartArr[binNdInd].card=clus.numOfTaxa;
							biPartArr[binNdInd].reset();
							k = gNodeBuf[ndI].childInd;
							for(j=0; j<gNodeBuf[ndI].numOfChild; ++j)
							{
								//find child's bin-node index in biPartArr
								if(UINT_MAX == gNodeBuf[curTree->nodeInd+gChildBuf[k +j]].childInd)
								{//this child is a leaf node
									binNdIndCh= gNodeBuf[curTree->nodeInd+gChildBuf[k +j]].numOfChild;//way of indexing leaf nodes
								}
								else
								{//this child is an int node
									binNdIndCh = (i+1)*twiceMaxClusOverlap+ gChildBuf[k+j];//way of indexing intermediate nodes
								}

								biPartArr[binNdInd].addPart(biPartArr[binNdIndCh]);//add this child's parition
								biPartArr[binNdInd].hashCode+= biPartArr[binNdIndCh].hashCode;//add child's hash code value			
							}
							biPartArr[binNdInd].binNdInd = hashTb[biPartArr[binNdInd].hashCode%PRIME_S].add(&biPartArr[binNdInd],binNdArr,binNdArrInd);
							binNdArr[biPartArr[binNdInd].binNdInd].freq = clus.numOfTrees;//occurs in all trees
						}
						--stackC;
					}
				}
				else
				{//encountered first time
					arrBufs.stackBuf[stackC][1] =0;//set all children visited and put all its children on the stack
					k = gNodeBuf[ndI].childInd;
					for(j=0; j<gNodeBuf[ndI].numOfChild; ++j)
					{
						arrBufs.stackBuf[++stackC][0]= curTree->nodeInd + gChildBuf[k+j];//put a child and increment stackC
						arrBufs.stackBuf[stackC][1] = gNodeBuf[arrBufs.stackBuf[stackC][0]].numOfChild;//mark that children needs to be visited
					}
				}
			}
		}
	}



	//Step2- Construct majority bi-parition tree in binNdArr


	majTree.nodeInd= gGlobalBuf.nodeBufInd++;//init nodeInd of majTree
	_C_ptr = biPartArr[twiceMaxClusOverlap].binNdInd;//index to root, root needs special treatement
	binNdArr[_C_ptr].ind =mtNdInd= ndPtrBufC++;
	ndPtrBuf[mtNdInd].supNum = clus.numOfTrees;//root is present in all the trees
	ndPtrBuf[mtNdInd].bl =0;//init the bl of the root to 0, the actual value would be added the first time the root is encoutered
	ndPtrBuf[mtNdInd].ind=majTree.nodeInd;
	ndPtrBuf[mtNdInd].chInd=childBufC;
	ndPtrBuf[mtNdInd].numChild=0;//init num of child to 0
	childBufC+=maxClusOverlap;//init for next free slot
	ndPtrBuf[mtNdInd].parInd=_C_ptr;//root is parent of itself


	for(unsigned short i=0; i<clus.numOfTrees; ++i)
	{
		_C_ptr = biPartArr[twiceMaxClusOverlap].binNdInd;//index to root
		stackC=0;
		curTree = &gTreeBuf[clus.treeInd+i];
		arrBufs.stackBuf[stackC][0] = curTree->nodeInd;//put root node on the stack
		arrBufs.stackBuf[stackC][1] = gNodeBuf[curTree->nodeInd].numOfChild;
		rootBL = gNodeBuf[curTree->nodeInd].bl;//we need to store it separately, the root node itself will change
		arrBufs.stackBuf[stackC][2] = _C_ptr;//value for _C_ptr for this bi-part


		while(USHRT_MAX != stackC)
		{
			ndI = arrBufs.stackBuf[stackC][0];//try not to alter ndI in the blocks below, as it is useful to have a ref to the node on the stack-
			//-esp while modifying the code, error tends to happen
			if(UINT_MAX == gNodeBuf[ndI].childInd)
			{//a leaf node- has to be a majority bi-part
				taxaId = gNodeBuf[ndI].numOfChild;//gives a way to calculate index (in biPartArr) for this bi-partition in this tree 
				binNdInd = biPartArr[taxaId].binNdInd;//stored when fetched from the hash table in Step1
				if(0 == i)
				{//singleton maj bi-paritions representing leaf nodes need to be added to majTree only once
					binNdArr[binNdInd].ind =mtNdInd = ndPtrBufC++;//allocate space to this majority tree node in MajTreeNodePtr buf
					ndPtrBuf[mtNdInd].bl = gNodeBuf[ndI].bl - rootBL; //init its bl value					
					ndPtrBuf[mtNdInd].supNum = clus.numOfTrees;//support num for the leaf node, a leaf would be there in all the trees
					ndPtrBuf[mtNdInd].ind = gGlobalBuf.nodeBufInd++;//assign a node to this leaf node in majTree
					ndPtrBuf[mtNdInd].parInd = _C_ptr;//  init its most recent maj bi-part ancestor as its parent in the majTree
					parInd = binNdArr[_C_ptr].ind;//get its parent index in MajTreeNodePtr buf
					ndPtrBuf[mtNdInd].posInParChBuf= ndPtrBuf[parInd].chInd + ndPtrBuf[parInd].numChild++;//the position at which this child-
					//- is store in the child buf of its parent
					arrBufs.childBuf[ndPtrBuf[mtNdInd].posInParChBuf] = mtNdInd;//store index in MajTreeNodePtr buf of this child in childBuf of its parent childbuf
					ndPtrBuf[mtNdInd].numChild = taxaId;
					ndPtrBuf[mtNdInd].chInd = UINT_MAX;//to indicate it as a taxa

					hashTb[biPartArr[taxaId].hashCode%PRIME_S].ind = UINT_MAX;//reset for next use
					//hashTb[biPartArr[taxaId].hashCode%PRIME_S].len =0;//reset for next use

				}
				else
				{//by now this node already exists
					mtNdInd = binNdArr[binNdInd].ind;//index of this node in MajTreeNodePtr buf
					ndPtrBuf[mtNdInd].bl+= gNodeBuf[ndI].bl - rootBL;//add the bl value
					parInd = ndPtrBuf[mtNdInd].parInd;//index of its parent in binNdArr
					if(binNdArr[parInd].biPart->card > binNdArr[_C_ptr].biPart->card)
					{//change of parent required
						parInd = binNdArr[parInd].ind;
						--ndPtrBuf[parInd].numChild;//decrement numChild count for its parent
						chInd = ndPtrBuf[parInd].chInd + ndPtrBuf[parInd].numChild;//index of child at the end which is being put at the position
						//- of the child being replaced
						ndPtrBuf[arrBufs.childBuf[chInd]].posInParChBuf = ndPtrBuf[mtNdInd].posInParChBuf;//update posInParChBuf for the child being put at the new position
						arrBufs.childBuf[ndPtrBuf[mtNdInd].posInParChBuf] = arrBufs.childBuf[ndPtrBuf[parInd].chInd + ndPtrBuf[parInd].numChild];//put the child at the end-
						//- at the new position

						//old parent changes done, now update new parent for this node
						ndPtrBuf[mtNdInd].parInd = _C_ptr;//new parent - index of its parent in binNdArr
						parInd = binNdArr[_C_ptr].ind;//index of the new parent in the MajTreeNodePtr buf
						ndPtrBuf[mtNdInd].posInParChBuf = ndPtrBuf[parInd].chInd + ndPtrBuf[parInd].numChild++;//update posInParChBuf-
						//-for the new parent, also incrment numChild for the new parent
						arrBufs.childBuf[ndPtrBuf[mtNdInd].posInParChBuf] = mtNdInd;//add this node in the child buf of its new parent						
					}
				}
				--stackC;
			}
			else
			{//an intermediate node
				if(0 == arrBufs.stackBuf[stackC][1])
				{//all its children have been visited
					binNdInd =arrBufs.stackBuf[stackC][3];//stored in stack when this node was encoutered first time
					if(binNdArr[binNdInd].freq > MAJ_LIM)
					{//a maj bi-parition
						//this maj bi-parition already exists - was created the first time it was encoutered and found to be maj bi-partition
						_C_ptr = arrBufs.stackBuf[stackC][2];//restore _C_ptr

						mtNdInd = binNdArr[binNdInd].ind;//index of this node in MajTreeNodePtr buf
						parInd = ndPtrBuf[mtNdInd].parInd;//index of its parent in binNdArr
						if(binNdArr[parInd].biPart->card > binNdArr[_C_ptr].biPart->card)
						{//change of parent required
							parInd = binNdArr[parInd].ind;
							--ndPtrBuf[parInd].numChild;//decrement numChild count for its parent
							chInd = ndPtrBuf[parInd].chInd + ndPtrBuf[parInd].numChild;//index of child at the end which is being put at the position
							//- of the child being replaced
							ndPtrBuf[arrBufs.childBuf[chInd]].posInParChBuf = ndPtrBuf[mtNdInd].posInParChBuf;//update posInParChBuf for the child being put at the new position
							arrBufs.childBuf[ndPtrBuf[mtNdInd].posInParChBuf] = arrBufs.childBuf[ndPtrBuf[parInd].chInd + ndPtrBuf[parInd].numChild];//put the child at the end-
							//- at the new position

							//old parent changes done, now update new parent for this node
							ndPtrBuf[mtNdInd].parInd = _C_ptr;//new parent - index of its parent in binNdArr
							parInd = binNdArr[_C_ptr].ind;//index of the new parent in the MajTreeNodePtr buf
							ndPtrBuf[mtNdInd].posInParChBuf = ndPtrBuf[parInd].chInd + ndPtrBuf[parInd].numChild++;//update posInParChBuf-
							//-for the new parent, also incrment numChild for the new parent
							arrBufs.childBuf[ndPtrBuf[mtNdInd].posInParChBuf] = mtNdInd;//add this node in the child buf of its new parent						
						}						
					}
					--stackC;
				}
				else
				{//encountered first time
					arrBufs.stackBuf[stackC][1] =0;//set all children visited and put all its children on the stack
					k =((i+1)*twiceMaxClusOverlap+ ndI - curTree->nodeInd);//gives a way to calculate index for this bi-partition in this tree
					hashTb[biPartArr[k].hashCode%PRIME_S].ind = UINT_MAX;//reset for next use
					binNdInd = biPartArr[k].binNdInd;//stored previously when fetched from the hash table in Step1
					arrBufs.stackBuf[stackC][3] = binNdInd;//so that next time, it is not fetched from the hash table
					if(binNdArr[binNdInd].freq > MAJ_LIM)
					{// a majoirity node encountered first time, needs to be created if it already does not exist in the tree
						if(UINT_MAX == binNdArr[binNdInd].ind)
						{//does not exists yet, needs to be created

							//NOTE that root has already been created, so this condition will never be true for root-
							//-root node was created outside the for loop

							//as every non trivial maj-bi parition will come only once here, hence an ideal place to 
							treeScore+=  binNdArr[binNdInd].freq;//add the support number score 
							++numOfUniqMajBiPart;//and also increase count of unique maj bi-partiions

							if(maxSupportNum < binNdArr[binNdInd].freq){//update maxSupportNum if required
								maxSupNumSecondBest = maxSupportNum;
								maxSupportNum = binNdArr[binNdInd].freq;								
							}else if(maxSupNumSecondBest< binNdArr[binNdInd].freq){
								maxSupNumSecondBest = binNdArr[binNdInd].freq;
							}


							binNdArr[binNdInd].ind =mtNdInd = ndPtrBufC++;//allocate space to this majority tree node in MajTreeNodePtr buf
							ndPtrBuf[mtNdInd].ind = gGlobalBuf.nodeBufInd++;//assign a node to this leaf node in majTree
							//only for non-root nodes do we need to add them to their parent child buf but-
							//-root has already been created, so this condition will never be true for root-
							//-root node was created outside the for loop
							ndPtrBuf[mtNdInd].parInd = _C_ptr;//  init its most recent maj bi-part ancestor as its parent in the majTree
							parInd = binNdArr[_C_ptr].ind;//get its parent index in MajTreeNodePtr buf
							ndPtrBuf[mtNdInd].posInParChBuf= ndPtrBuf[parInd].chInd + ndPtrBuf[parInd].numChild++;//the position at which this child-
							//- is store in the child buf of its parent
							arrBufs.childBuf[ndPtrBuf[mtNdInd].posInParChBuf] = mtNdInd;//store index in MajTreeNodePtr buf of this child in childBuf of its parent

							ndPtrBuf[mtNdInd].numChild =0;//init numChild =0
							ndPtrBuf[mtNdInd].chInd = childBufC;//init chInd
							childBufC+= maxClusOverlap;//make it point to next free slot
							ndPtrBuf[mtNdInd].bl =0;//init bl to 0, the bl value would be added immediately
							ndPtrBuf[mtNdInd].supNum= binNdArr[binNdInd].freq;//update the support num
						}
						ndPtrBuf[binNdArr[binNdInd].ind].bl+= gNodeBuf[ndI].bl -rootBL;//add the bl value
						arrBufs.stackBuf[stackC][2] =_C_ptr; //save _C_ptr
						_C_ptr = binNdInd;//update with new _C_ptr as this node
					}

					k = gNodeBuf[ndI].childInd;
					for(j=0; j<gNodeBuf[ndI].numOfChild; ++j)
					{//put all its children on the stack
						arrBufs.stackBuf[++stackC][0]= curTree->nodeInd + gChildBuf[k+j];//put a child and increment stackC
						arrBufs.stackBuf[stackC][1] = gNodeBuf[arrBufs.stackBuf[stackC][0]].numOfChild;//mark that children needs to be visited
					}
				}
			}
		}
	}


	//check if root is going to have only two children, in which case one of the children will collapse later and -
	//- its support value will have to be subtraced from the tree score
	unsigned short deductableSupVal =0;
	if(2 == ndPtrBuf[0].numChild){//yes, root is going to have two children only
		mtNdInd = arrBufs.childBuf[ndPtrBuf[0].chInd+0];//first child
		if(UINT_MAX == ndPtrBuf[mtNdInd].chInd){//first child is a leaf, second child will definitely not be a leaf
			mtNdInd = arrBufs.childBuf[ndPtrBuf[0].chInd+1];//second child
		}
		deductableSupVal= ndPtrBuf[mtNdInd].supNum;//deduct this from treeScore
		treeScore -= deductableSupVal;
		--numOfUniqMajBiPart;
		if(maxSupportNum == deductableSupVal){//i.e. the best support value was of this node only, need to replace it -
			//- with second best value
			maxSupportNum = maxSupNumSecondBest;
		}
	}

	//step 3- now transform the maj bi-part tree constructed in binNdArr to one in gNodeBuf
	stackC =0;
	//root needs special treatment
	arrBufs.stackBuf[stackC][0] = 0;//load the root
	majTree.numOfNodes = ndPtrBufC;// init numOfNodes
	clus.avgSupInd = gGlobalBuf.avgSupBufInd;//allocate memory for avg support numbers for this cluster
	gGlobalBuf.avgSupBufInd+= majTree.numOfNodes;//make it point to the next available slot.
	double tmpSup;

	while(USHRT_MAX != stackC)
	{
		mtNdInd = arrBufs.stackBuf[stackC--][0];//pop from the stack
		if(UINT_MAX == ndPtrBuf[mtNdInd].chInd)
		{//a leaf node
			ndI = ndPtrBuf[mtNdInd].ind;
			tmpSup = ndPtrBuf[mtNdInd].supNum*100.0/clus.numOfTrees;//normalize the % value
			gAvgSupBuf[clus.avgSupInd+ndI-majTree.nodeInd] = (unsigned char)floor(0.5 + tmpSup);//indexing of supNum is done in a way similar to indexing of-
			//-nodes of the tree, so that retreival of supNum can happen similarly
			gNodeBuf[ndI].childInd = UINT_MAX;//to indicate it as a taxaId
			gNodeBuf[ndI].numOfChild = ndPtrBuf[mtNdInd].numChild;//taxaId
			tmpDl = ndPtrBuf[mtNdInd].bl;
			tmpDl = tmpDl/ndPtrBuf[mtNdInd].supNum;//take the avg bl value
			gNodeBuf[ndI].bl = (unsigned short)tmpDl;
			if(tmpDl - gNodeBuf[ndI].bl >= .5)
			{
				++gNodeBuf[ndI].bl;
			}
		}
		else
		{//an intermediate node
			ndI = ndPtrBuf[mtNdInd].ind;
			tmpSup = ndPtrBuf[mtNdInd].supNum*100.0/clus.numOfTrees;//normalize the % value
			gAvgSupBuf[clus.avgSupInd+ndI-majTree.nodeInd] = (unsigned char)floor(0.5 + tmpSup);//indexing of supNum is done in a way similar to indexing of-
			//-nodes of the tree, so that retreival of supNum can happen similarly
			tmpDl = ndPtrBuf[mtNdInd].bl;
			tmpDl = tmpDl/ndPtrBuf[mtNdInd].supNum;//take the avg bl value
			gNodeBuf[ndI].bl = (unsigned short) tmpDl;
			if(tmpDl - gNodeBuf[ndI].bl >= .5)
			{
				++gNodeBuf[ndI].bl;
			}
			gNodeBuf[ndI].numOfChild = ndPtrBuf[mtNdInd].numChild;
			gNodeBuf[ndI].childInd = gGlobalBuf.childBufInd;//init childInd
			gGlobalBuf.childBufInd+=gNodeBuf[ndI].numOfChild;//point to the next free slot
			for(j=0; j<ndPtrBuf[mtNdInd].numChild; ++j)
			{//put all its children on the stack
				arrBufs.stackBuf[++stackC][0] = arrBufs.childBuf[ndPtrBuf[mtNdInd].chInd+j];
				gChildBuf[gNodeBuf[ndI].childInd+j] = (unsigned short)(ndPtrBuf[arrBufs.stackBuf[stackC][0]].ind - majTree.nodeInd);//relative position of this child
			}
		}
	}

	//because of outgrouping, root of the MRT will definitely have two children - not OK for an unrooted tree
	chInd = gNodeBuf[majTree.nodeInd].childInd;
	unsigned int ch1 = majTree.nodeInd + gChildBuf[chInd];//first child
	unsigned int ch2 = majTree.nodeInd + gChildBuf[chInd+1];//second child

	if(UINT_MAX != gNodeBuf[ch1].childInd){//first child is not a leaf; add all children of the first child-
		//- to the root and then add the second child as the last child 
		chInd = gNodeBuf[majTree.nodeInd].childInd = gGlobalBuf.childBufInd;
		gGlobalBuf.childBufInd+= gNodeBuf[majTree.nodeInd].numOfChild = gNodeBuf[ch1].numOfChild+1;		
		for(unsigned short i=0; i<gNodeBuf[ch1].numOfChild; ++i){
			gChildBuf[chInd + i] = gChildBuf[gNodeBuf[ch1].childInd +i];
		}

		gChildBuf[chInd + gNodeBuf[ch1].numOfChild] = ch2 - majTree.nodeInd;//add the second child as the last child
	}else{//second child is not a leaf; add all children of the second child-
		//- to the root and then add the first child as the last child 
		chInd = gNodeBuf[majTree.nodeInd].childInd = gGlobalBuf.childBufInd;
		gGlobalBuf.childBufInd+= gNodeBuf[majTree.nodeInd].numOfChild = gNodeBuf[ch2].numOfChild+1;		
		for(unsigned short i=0; i<gNodeBuf[ch2].numOfChild; ++i){
			gChildBuf[chInd + i] = gChildBuf[gNodeBuf[ch2].childInd +i];
		}
		gChildBuf[chInd + gNodeBuf[ch2].numOfChild] = ch1 - majTree.nodeInd;//add the first child as the last child
	}


	//code to calculate the pd value of the MRT - no need as we are not using this pd value
	//stackC =0;
	clus.pdStat.mrtPD =0;
	//arrBufs.stackBuf[stackC][0]= majTree.nodeInd;//put root on the stack
	//while(stackC != USHRT_MAX)
	//{
	//	ndI = arrBufs.stackBuf[stackC--][0];
	//	if(gNodeBuf[ndI].childInd == UINT_MAX)
	//	{
	//		continue;
	//	}
	//	for(int i=0; i<gNodeBuf[ndI].numOfChild; ++i)
	//	{
	//		//put the child on the stack
	//		arrBufs.stackBuf[++stackC][0] = majTree.nodeInd + gChildBuf[gNodeBuf[ndI].childInd+i];
	//		sL = gNodeBuf[arrBufs.stackBuf[stackC][0]].bl - gNodeBuf[ndI].bl; //pd val diff w.r.t its parent			
	//		if(sL > 0)
	//		{//its possible that bl of the parent might be greater than bl of its child by a slight margin only, but still
	//			clus.pdStat.mrtPD+= sL ;
	//		}
	//	}
	//}


	tmpSup= maxSupportNum*100.0/clus.numOfTrees;//normalize wrt to num of trees
	maxSupportNum = (unsigned short)floor(tmpSup+ .5);//round off

	tmpSup= treeScore*100.0/clus.numOfTrees;//normalize wrt to num of trees
	
	
	if(numOfUniqMajBiPart !=0)
	{
		//tmpSup/=numOfUniqMajBiPart;//normalize wrt num of maj bi-paritions
		tmpSup/=(clus.numOfTaxa-3);//normalize by maximum possible non-trivial biparitions for an unrooted tree
	}
	else
	{
		tmpSup =0;
	}
	treeScore = (unsigned int)floor(tmpSup+ .5);

	return (unsigned char)treeScore;

}


unsigned int Bin::getBiPar(BiPart *biPart, BinNd * binNdArr)
{
	static unsigned int i;
	i = ind;
	if(1 == len)
	{
		return ind;
	}
	while(true)
	{
		if(binNdArr[i].biPart->comp(*biPart))
		{//biPart found
			return i;
		}
		i = binNdArr[i].nxt;
	}
	return UINT_MAX;//this should never happen
}


unsigned int Bin::add(BiPart *biPart,BinNd * binNdArr,unsigned int &binNdArrInd)
{
	static unsigned int i;
	if(UINT_MAX == ind)
	{//first element being inserted in this bucket / bin
		ind=binNdArrInd++;
		binNdArr[ind].biPart = biPart;
		binNdArr[ind].nxt = UINT_MAX;
		binNdArr[ind].freq=1;
		binNdArr[ind].ind = UINT_MAX;//to indicate that this node does not exist in majTree yet
		len=1;
		return ind;
	}
	else
	{
		i = ind;
		while(true)
		{
			if(binNdArr[i].biPart->comp(*biPart))
			{//this biPart already exists
				++binNdArr[i].freq;//increment its freq
				return i;
			}
			if(UINT_MAX == binNdArr[i].nxt)
			{
				break;
			}
			i =  binNdArr[i].nxt;
		}
		binNdArr[i].nxt =binNdArrInd;
		binNdArr[binNdArrInd].biPart = biPart;
		binNdArr[binNdArrInd].freq =1;
		binNdArr[binNdArrInd].nxt = UINT_MAX;
		binNdArr[binNdArrInd].ind = UINT_MAX;//to indicate that this node does not exist in majTree yet
		++binNdArrInd;
		++len;
		return binNdArr[i].nxt;
	}
}

void BiPart::addPart(BiPart &biPart) //add parition from biPart
{
	static unsigned char i;
	part[0]|= biPart.part[0];
	if(len>1)
	{//ok.. we have a longish bi-partition, OR rest of 64 bit parts
		for(i=1; i<len; ++i)
		{
			part[i]|= biPart.part[i];
		}
	}

}

void BiPart::reset()//resets all 64-bit spaces for next use
{
	static unsigned char i;
	part[0]=0;
	if(len>1)
	{//ok.. we have a longish bi-partition, rest rest of 64 bit parts
		for(i=1; i<len; ++i)
		{
			part[i]=0;
		}
	}
}

bool BiPart::comp(BiPart &biPart)//compares two bi-parititions
{
	static unsigned char i;
	if(card != biPart.card)//check if cardinality is same
	{
		return false;
	}
	if(part[0] != biPart.part[0])//check if first 64 bit is same. 99% bi-partitions are expected have cardinality <64
	{
		return false;
	}
	if(len>1)
	{//ok.. we have a longish bi-partition, check for rest of 64 bit parts
		for(i=1; i<len; ++i)
		{
			if(part[i] != biPart.part[i])
			{
				return false;
			}
		}
	}
	return true;
}

void BiPart::reverse()
{
	static unsigned char i;
	part[0]= ~part[0];
	if(len>1)
	{//ok.. we have a longish bi-partition, reverse rest of 64 bit parts
		for(i=1; i<len; ++i)
		{
			part[i]= ~part[i];
		}
	}
	part[len-1]= part[len-1]&lastOfset;//zero the extra taxa in the last 64-bit space
}

unsigned short BiPart::diff(BiPart &biPart, ArrayBufs &arrBufs)
{
	unsigned long long ull;
	unsigned short * us;
	unsigned short diff1=0;
	unsigned short diff2=0;

	for(unsigned char i=0; i<len; ++i){
		ull = part[i] ^ biPart.part[i];//exclusive OR to get bit difference
		us = (unsigned short *)&ull;
		for(unsigned char j=0; j<4; ++j){
			diff1+= arrBufs.bitCntArr[us[j]];			
		}
	}

	reverse();//reverse this bipartition to get the other diff

	for(unsigned char i=0; i<len; ++i){
		ull = part[i] ^ biPart.part[i];//exclusive OR to get bit difference
		us = (unsigned short *)&ull;
		for(unsigned char j=0; j<4; ++j){
			diff2+= arrBufs.bitCntArr[us[j]];			
		}
	}

	//no need to restore, as we need to compare both partitions anyway
	//reverse();//restore this biparition to its original state

	if(diff1 < diff2){
		return diff1;
	}else{
		return diff2;
	}
}

void genPrime(unsigned int *prArr,unsigned int numPr)
{
	unsigned int prC=0;
	bool validPr;
	unsigned int curPr,rem,quo;
	set<unsigned int> uniqPr;
	pair<set<unsigned int>::iterator,bool> ret;
	for(unsigned i=0; i<numPr; ++i)
	{
		validPr = false;
		while(!validPr)
		{
			curPr = rand()+1;
			quo = PRIME_S/curPr;
			rem = PRIME_S%curPr;
			quo = rand()%(quo+1);
			rem = rand()%(rem+1);
			curPr = curPr*quo + rem;
			if(curPr < PRIME_S)
			{
				ret = uniqPr.insert(curPr);
				validPr = ret.second;
			}

		}
		prArr[i] = curPr;
	}

}

bool rootMRT(PhyloTree &refTree, PhyloTree &mrtTree, BiPart * biPartArr, ArrayBufs &arrBufs, unsigned int clusAvgSupInd)
{
	unsigned int newRootInd;
	unsigned short smDiff = USHRT_MAX, curDiff;
	if(gNodeBuf[refTree.nodeInd].numOfChild >2) {
		return false;//reference tree has more than two children, cannot be used for rooting
	}


	initRefTreeRootBipartition(refTree,biPartArr,0,arrBufs);
	BiPart & refRootBp = biPartArr[0];
	initBipartitions(mrtTree,biPartArr,1,arrBufs);

	unsigned short stackC;//counter for stackBuf
	unsigned int chI,ndI,j,k;//temp values-used repeatedly in for loops
	unsigned int (*stack)[4] = arrBufs.stackBuf;

	//first find the new root i.e. the node with the smallest bipartition difference

	stackC=0;

	stack[stackC][0] = mrtTree.nodeInd;//put the root on stack

	while(USHRT_MAX != stackC)
	{
		ndI = stack[stackC][0];//try not to alter ndI in the blocks below, as it is useful to have a ref to the node on the stack-
		--stackC;//remove it from the stack
		for(unsigned short i=0; i<gNodeBuf[ndI].numOfChild; ++i){
			chI = gChildBuf[gNodeBuf[ndI].childInd + i];
			curDiff = refRootBp.diff(biPartArr[chI +1],arrBufs);//+1 because biparitions of mrtTree starts at 1 in biPartArr
			chI += mrtTree.nodeInd;//actual child index
			if(curDiff < smDiff){
				smDiff = curDiff;
				newRootInd = chI;
			}
			if(UINT_MAX != gNodeBuf[chI].childInd){//this child is not a leaf, put it on the stack
				stack[++stackC][0] = chI;
			}
		}		
	}

	//now traverse the tree again till the newRootInd
	stackC=0;
	//as root does not have any parent, it needs special treatement	
	stack[stackC][0] = mrtTree.nodeInd;//init to root
	stack[stackC][1] = 0; //visited child count	

	//put the first child of the root on stack
	chI = gNodeBuf[mrtTree.nodeInd].childInd;//child index of root node
	stack[++stackC][0] = mrtTree.nodeInd+gChildBuf[chI];//node index of first child

	while(true)
	{
		ndI = stack[stackC][0];//try not to alter ndI in the blocks below, as it is useful to have a ref to the node on the stack-
		//-esp while modifying the code, error tends to happen	
		++stack[stackC-1][1];//increment visited child count of its parent
		if(ndI == newRootInd){
			break;
		}
		if(UINT_MAX == gNodeBuf[ndI].childInd)
		{//a leaf node			
			//this node has been visited, now find the next node in tree in pre-order i.e. node, left to right children
			--stackC;//go to parent of curr node
			while(true)
			{
				ndI = stack[stackC][0];//index of the parent, ndI is being modified here since parent is the current node under consideration
				if(stack[stackC][1] == gNodeBuf[ndI].numOfChild)
				{//all children of this node are visited					
					--stackC;//this node visited,go to his parent continuing search for the next node in tree in pre-order
				}
				else
				{//some children are still left, put next child on top of stack					
					chI = gNodeBuf[ndI].childInd + stack[stackC][1];
					stack[++stackC][0] = mrtTree.nodeInd+gChildBuf[chI];//node index of next child
					//no need to init this child, as it would be done when it picked from the stack, for the first time
					break;
				}
			}
		}
		else
		{//an intermediate node picked from the stack for the first time, init the node and put its first child on the stack
			stack[stackC][1] = 0; //visited cildren  count
			//now put its first child on stack
			++stackC;//increment stack
			stack[stackC][0]= gChildBuf[gNodeBuf[ndI].childInd] + mrtTree.nodeInd;//put the first child on stack
		}
	}

	if(stack[stackC-1][0] == mrtTree.nodeInd && 2 == gNodeBuf[mrtTree.nodeInd].numOfChild){
		//if the root edge is a child of the existing root and the existing root has only two children, then no need to re-root
		return true;
	}


	//first the existing root will have one less child, remove the child	
	j = gNodeBuf[mrtTree.nodeInd].childInd + stack[0][1] -1;//index of the child to be removed
	k = gNodeBuf[mrtTree.nodeInd].childInd + gNodeBuf[mrtTree.nodeInd].numOfChild -1;//last child index
	gChildBuf[j] = gChildBuf[k];//overwrite the existing child with the last child
	--gNodeBuf[mrtTree.nodeInd].numOfChild;//decrement the number of children
	//--mrtTree.numOfNodes;//decrement the number of nodes
	if(1 == gNodeBuf[mrtTree.nodeInd].numOfChild){//i.e. existing root has only one child now, overwrite the root with the only child
		j = gNodeBuf[mrtTree.nodeInd].childInd;
		gNodeBuf[mrtTree.nodeInd] = gNodeBuf[mrtTree.nodeInd + gChildBuf[j]];
		//--mrtTree.numOfNodes;//decrement the number of nodes OK
	}

	//copy the existing root node to a new node so that mrtTree.nodeInd is now free to become the new root
	gNodeBuf[gGlobalBuf.nodeBufInd] = gNodeBuf[mrtTree.nodeInd];
	stack[0][0]  = gGlobalBuf.nodeBufInd++;//replace mrtTree.nodeInd with the newly allocated node and --
	//-- increment gGlobalBuf.nodeBufInd for the next slot, it is OK to use like this, this wont disrupt --
	//-- any ndI - mrtTree.nodeInd mapping as this function is being called immediately after genMRT where --
	//-- the mrtTree was alloted nodes from gNodeBuf
	//currently the newRootInd is on the top of stack, replace it with mrtTree.nodeInd
	gNodeBuf[mrtTree.nodeInd].childInd = gGlobalBuf.childBufInd;
	gNodeBuf[mrtTree.nodeInd].numOfChild = 2;
	gChildBuf[gGlobalBuf.childBufInd] = ndI - mrtTree.nodeInd;//make newRootInd as the first child
	stack[stackC][0] = mrtTree.nodeInd;//put new root on the stack
	stack[stackC][1] = 2;//mark that the new child will be added as second child
	gGlobalBuf.childBufInd+= 2;
	gAvgSupBuf[clusAvgSupInd] = gAvgSupBuf[clusAvgSupInd + ndI - mrtTree.nodeInd];//assign the avgSup value of newRootInd to the new root --
	//-- i.e. mrtTree.nodeInd, so that while reversing parents, avgSup values are also properly reversed. Note that after parent reversal, --
	//-- the root node (mrtTree.nodeInd) will have two children with the same avgSup value

	while(stackC){//now reverse the parents		
		j = stack[stackC][0];//cur node index
		k = stack[stackC-1][0] - mrtTree.nodeInd;//old parent's relative node ind
		gChildBuf[gNodeBuf[j].childInd+ stack[stackC][1] -1] = (unsigned short)k;
		gAvgSupBuf[clusAvgSupInd + k] = gAvgSupBuf[clusAvgSupInd + j - mrtTree.nodeInd];//since the root direction has changed, cluster direction --
		//-- has also changed, thus avgSup values need to be reassigned to the new clusters from the old clusters. 
		--stackC;
	}	

	gAvgSupBuf[clusAvgSupInd] =100;//assign root 100% support
	return true;
}


void initBipartitions(PhyloTree &tree, BiPart * biPartArr, unsigned int biPartInd, ArrayBufs &arrBufs)
{
	unsigned long long ONE =1;//used to obtain the actual bi-parition from [0] entry of taxaArr
	unsigned int ndI,k,bpArrInd;
	unsigned short stackC=0,j;

	unsigned int (*stack)[4] = arrBufs.stackBuf;

	stack[stackC][0] = tree.nodeInd;//put root on stack
	stack[stackC][1] = gNodeBuf[tree.nodeInd].numOfChild;//num of remaning children to be visited


	while(USHRT_MAX != stackC)
	{
		ndI = stack[stackC][0];
		bpArrInd = biPartInd + ndI - tree.nodeInd;
		if(UINT_MAX == gNodeBuf[ndI].childInd)
		{//a leaf node		
			j = gNodeBuf[ndI].numOfChild;//taxaId
			biPartArr[bpArrInd].reset();
			biPartArr[bpArrInd].part[arrBufs.taxaArr2[j][1]]= ONE<<arrBufs.taxaArr2[j][0];			
			--stackC;
		}
		else
		{//an intermediate node
			if(0 == stack[stackC][1])
			{//all its children have been visited, add their partitions
				biPartArr[bpArrInd].reset();
				k = gNodeBuf[ndI].childInd;
				for(j=0; j<gNodeBuf[ndI].numOfChild; ++j){
					biPartArr[bpArrInd].addPart(biPartArr[biPartInd+gChildBuf[k+j]]);					
				}
				--stackC;
			}
			else
			{//encountered first time
				stack[stackC][1] =0;//set all children visited and put all its children on the stack
				k = gNodeBuf[ndI].childInd + gNodeBuf[ndI].numOfChild -1;
				for(j=0; j<gNodeBuf[ndI].numOfChild; ++j)
				{
					stack[++stackC][0]= tree.nodeInd + gChildBuf[k-j];//put a child and increment stackC, k-j so that --
					//-- the first child is popped first from the stack
					stack[stackC][1] = gNodeBuf[stack[stackC][0]].numOfChild;//mark that children needs to be visited					
				}
			}
		}
	}
}



void initRefTreeRootBipartition(PhyloTree &refTree, BiPart * biPartArr, unsigned int biPartInd, ArrayBufs &arrBufs)
{
	unsigned long long ONE =1;//used to obtain the actual bi-parition from [0] entry of taxaArr
	unsigned int ndI,k;
	unsigned short stackC=0,j;

	unsigned int (*stack)[4] = arrBufs.stackBuf;



	biPartArr[biPartInd].reset();
	ndI = gNodeBuf[refTree.nodeInd].childInd;
	ndI = refTree.nodeInd + gChildBuf[ndI];
	stackC =0;
	stack[stackC][0] = ndI;//put 1st child on stack, as we need the cluster corresponding to the subtree at one of its children
	stack[stackC][1] = gNodeBuf[ndI].numOfChild;//num of remaning children to be visited
	while(USHRT_MAX != stackC)
	{
		ndI = stack[stackC][0];		
		if(UINT_MAX == gNodeBuf[ndI].childInd)
		{//a leaf node		
			j = gNodeBuf[ndI].numOfChild;//taxaId
			biPartArr[biPartInd].part[arrBufs.taxaArr2[j][1]]|= ONE<<arrBufs.taxaArr2[j][0];			
			--stackC;
		}
		else
		{//an intermediate node
			if(0 == stack[stackC][1]){
				--stackC;//remove it from the stack, go to its parent
			}
			else
			{//encountered first time
				stack[stackC][1] =0;//set all children visited and put all its children on the stack
				k = gNodeBuf[ndI].childInd + gNodeBuf[ndI].numOfChild -1;
				for(j=0; j<gNodeBuf[ndI].numOfChild; ++j)
				{
					stack[++stackC][0]= refTree.nodeInd + gChildBuf[k-j];//put a child and increment stackC, k-j so that --
					//-- the first child is popped first from the stack
					stack[stackC][1] = gNodeBuf[stack[stackC][0]].numOfChild;//mark that children needs to be visited					
				}
			}
		}
	}

}

