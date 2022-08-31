#include"T1_test.h"
#include<cmath>

void testRootMRT()
{
	ArrayBufs arrBufs;
	BiPart *biPartArr = new BiPart[twiceMaxClusOverlap*101];//100mb max taxa =maxClusOverlap, max num of trees=100 (<101)
	char * prBuf = new char[oneMillion];

	unsigned short numTaxa = 7;//numTaxa should be exact else BiPart::lastOfset would not be init properly
	unsigned int clusTxInd = gGlobalBuf.taxaBufInd;
	gGlobalBuf.taxaBufInd+= numTaxa;
	for(unsigned short i=0; i<numTaxa; ++i){
		gTaxaBuf[clusTxInd + i] = i;
	}

	unsigned long long ONE =1;//used to obtain the actual bi-parition from [0] entry of taxaArr as described above
	BiPart::len = numTaxa/64 + 1;// init BiPart::len
	BiPart::lastOfset =0;
	for(unsigned short i= (numTaxa/64)*64; i<numTaxa; ++i){
		BiPart::lastOfset |= ONE<<arrBufs.taxaArr2[i][0];
	}


	char strRefTree[] = "(2#(3#0,1,2),(4#3,4,5,6));";
	char strMrtTree[] = "(2#0,(4#1,2,(3#3,4,5),6));";
	//char strMrtTree[] = "(2#0,(3#1,(2#2,3),4));";
	//char strMrtTree[] = "(2#1,(2#0,2));";
	PhyloTree refTree, mrtTree;
	
	parseTestNewick(strRefTree,refTree,arrBufs);
	testPrint(refTree,arrBufs,prBuf,clusTxInd);
	cout<<prBuf<<endl;

	parseTestNewick(strMrtTree,mrtTree,arrBufs);
	testPrint(mrtTree,arrBufs,prBuf,clusTxInd);
	cout<<prBuf<<endl;

	rootMRT(refTree, mrtTree,biPartArr,arrBufs,0);
	testPrint(mrtTree,arrBufs,prBuf,clusTxInd);
	cout<<prBuf<<endl;
	
	delete []prBuf;
	delete []biPartArr;
}

void testSubtree()
{
	ArrayBufs arrBufs;
	for(unsigned int i=0;i<USHRT_MAX; ++i)
	{
		arrBufs.taxaArr[i] = USHRT_MAX;//init all values to false
	}
	char *prBuf = new char[oneMillion];

	ClusterInfo &sClus = gClusBuf[gGlobalBuf.clusBufInd++];
	ClusterInfo &tClus = gClusBuf[gGlobalBuf.clusBufInd++];
	sClus.numOfTaxa = 5;
	sClus.taxaInd = gGlobalBuf.taxaBufInd;
	gGlobalBuf.taxaBufInd+= sClus.numOfTaxa;
	for(unsigned short i=0; i<sClus.numOfTaxa; ++i){
		gTaxaBuf[sClus.taxaInd + i] = i;
	}
	//arrBufs.taxaArr[0] = 0;
	arrBufs.taxaArr[1] = 0;
	arrBufs.taxaArr[2] = 0;
	arrBufs.taxaArr[3] = 0;
	arrBufs.taxaArr[4] = 0;

	sClus.numOfTrees = 2;
	string *treeArr = new string[sClus.numOfTrees];	
	treeArr[0] = "(2#0,(2#1,(2#2,(2#3,4))));";
	treeArr[1] = "(2#0,(2#1,(2#2,(2#3,4))));";
	sClus.treeInd = gGlobalBuf.treeBufInd;
	gGlobalBuf.treeBufInd+= sClus.numOfTrees;
	for(unsigned short i=0; i<sClus.numOfTrees; ++i){
		parseTestNewick(const_cast<char *>(treeArr[i].c_str()),gTreeBuf[sClus.treeInd+i],arrBufs);
		testPrint(gTreeBuf[sClus.treeInd +i],arrBufs,prBuf,sClus.taxaInd);
		cout<<prBuf<<endl;
	}
	genSubtree(sClus,tClus,arrBufs,true);
	testPrint(sClus.clusterAsTree,arrBufs,prBuf,tClus.taxaInd);
	cout<<prBuf<<endl;
	for(unsigned short i=0; i<sClus.numOfTrees-1; ++i){
		testPrint(gTreeBuf[tClus.treeInd +i],arrBufs,prBuf,tClus.taxaInd);
		cout<<prBuf<<endl;
	}
	
	delete []prBuf;
	delete []treeArr;
}

void parseTestNewick(char *strTree,PhyloTree &currTree, ArrayBufs &arrBufs)
{
	char *pos1, *pos2;
	char ch;//tmp char
	unsigned short count =0; //to store the bracket count, init to 0	
	unsigned short stackC=0;//counter which points to the last stored element in tmpStack, init to 0
	unsigned int  chI;//temp var	
	
	
	//Only one pass is required
	//NOTE: the way newick tree is written the traversal is going to be in pre-order i.e root, left child, right child
	pos1 = strchr(strTree,'(');//find where the root starts
	++pos1;//skip the root, as it does not have any parent hence needs to be dealt before
	chI = currTree.nodeInd = gGlobalBuf.nodeBufInd++;
	unsigned (*stack)[4] = arrBufs.stackBuf;
	stack[stackC][0] = chI;//init to root	 
	stack[stackC][1] = gNodeBuf[chI].childInd = gGlobalBuf.childBufInd;//index where next child will be added

	pos2 = strchr(pos1,'#');
	*pos2 = 0;		
	gGlobalBuf.childBufInd+= gNodeBuf[chI].numOfChild = atoi(pos1);//numChild of root;			
	pos1 = pos2 +1;//skip #
	
	do
	{		
		ch = *pos1;
		if('(' == ch)
		{//new intermediate node		
			chI =  gGlobalBuf.nodeBufInd++;//assign a node to this new node
			gChildBuf[stack[stackC][1]++] = chI - currTree.nodeInd;//add this new node to its parent			
			stack[++stackC][0] = chI;//put this new child on the stack
			stack[stackC][1] = gNodeBuf[chI].childInd = gGlobalBuf.childBufInd;//index where next child will be added for this new node
			++pos1;
			pos2 = strchr(pos1,'#');
			*pos2 = 0;
			gGlobalBuf.childBufInd+= gNodeBuf[chI].numOfChild = atoi(pos1);//numChild of this new node;				
			pos1 = pos2 +1;//skip #
			++count;//bracket count
		}
		else if(')' == ch )
		{
			--count;//bracket count
			--stackC;//decrement stackC to go one level up
			++pos1;
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
			chI =  gGlobalBuf.nodeBufInd++;//assign a node to this new node
			gChildBuf[stack[stackC][1]++] = chI - currTree.nodeInd;//add this new node to its parent						
			pos2 = strpbrk(pos1,":,)");//':' is to ignore other nodal attributes
			ch = *pos2;//save - to restore later
			*pos2 = 0;
			gNodeBuf[chI].numOfChild = atoi(pos1);//taxaId
			*pos2 = ch;//restore
			pos1 = pos2;//skip taxa id for next pos			
		}
	}
	while(USHRT_MAX != count);//tree has been parsed, sine the '(' for root was skipped, hence-
	//- count will reach -1 (USHRT_MAX) only when the very last ')' is reached.
	--stackC;//not needed for the root node, but just to make it look neat and in synch with the above for loop
	currTree.numOfNodes= (unsigned short)(gGlobalBuf.nodeBufInd - currTree.nodeInd);
}


void testPrint(PhyloTree &currTree, ArrayBufs &arrBufs, char *prBuf, unsigned int clusTxInd)
{
	unsigned int ndI,k;
	double taxaId;
	unsigned short stackC=0,j;
	size_t ep=0;
	unsigned char sz;
	unsigned int (*stack)[4] = arrBufs.stackBuf;
	
	stack[stackC][0] = currTree.nodeInd;//put root on stack
	stack[stackC][1] = gNodeBuf[currTree.nodeInd].numOfChild;
	while(USHRT_MAX != stackC)
	{
		ndI = stack[stackC][0];
		if(UINT_MAX == gNodeBuf[ndI].childInd)
		{//a leaf node			
			k = clusTxInd + gNodeBuf[ndI].numOfChild;
			taxaId  = k = gTaxaBuf[k];//taxaId
			if(0 != taxaId)
			{
				ep+= sz = (unsigned char)log10(taxaId)+1;//size of taxaId
				while(k!=0)
				{
					prBuf[--ep] = '0' + k%10;
					k=k/10;
				}
				ep+=sz;
			}
			else
			{
				prBuf[ep++] = '0';
			}							
			prBuf[ep++] = ',';
			--stackC;
		}
		else
		{//an intermediate node
			if(0 == stack[stackC][1])
			{//all its children have been visited
				prBuf[--ep] = ')';//overwirte the extra ','				
				++ep;
				prBuf[ep++] =',';
				--stackC;
			}
			else
			{//encountered first time
				prBuf[ep++]='(';
				stack[stackC][1] =0;//set all children visited and put all its children on the stack
				k = gNodeBuf[ndI].childInd + gNodeBuf[ndI].numOfChild -1;
				for(j=0; j<gNodeBuf[ndI].numOfChild; ++j)
				{
					stack[++stackC][0]= currTree.nodeInd + gChildBuf[k-j];//put a child and increment stackC, k-j so that --
					//-- the first child is popped first from the stack
					stack[stackC][1] = gNodeBuf[stack[stackC][0]].numOfChild;//mark that children needs to be visited
				}
			}
		}
	}

	prBuf[--ep]='\0';	//overwrite the last ','
}