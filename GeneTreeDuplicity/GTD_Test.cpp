#include"GTD_Test.h"


void EvalST_MulTrees::evalAllClusters()
{

	UL_HashTable mapOrgSeqOrgTaxa(12000017);//= 12000017 first prime greater than 12 million for sequence
	loadAllMaps(mapOrgSeqOrgTaxa);

	char parsedTreeBaseDir[] ="E:/TreeIndexing/BicliqueTrees_June19_9million/parsedblTrees";
	
	ifstream clusFile("E:/TreeIndexing/BicliqueTrees_June19_9million/clusters/Clus_0.db");
	if(!clusFile.is_open()){
		cout<<"ERROR! Could not open clusFile"<<endl;
		exit(-1);
	}

	ofstream resFile("E:/TreeIndexing/BicliqueTrees_June19_9million/MulTreeAnalysis2.txt");
	if(!resFile.is_open()){
		cout<<"ERROR! Could not open resFile."<<endl;
		exit(-1);
	}

	char *pos, *pos2;
	char type;
	char clusPrefix[1024],clusFileName[1024];
	resFile.setf(ios::fixed,ios::floatfield);
	resFile.precision(1);
	
	
	unsigned int clusC =0;
	while(true){
		clusFile.getline(MB.chBuf,2*MB.TAXA_MAX*MB.TAXA_LEN_MAX);
		if(0 == clusFile.gcount() || strlen(MB.chBuf) == 0){
			break;
		}
		type = MB.chBuf[0];
		if(type != '3'){
			continue;
		}
		if(clusC%1000 ==0){
			cout<<clusC<<"\t";
		}
		++clusC;
		pos = strchr(MB.chBuf,'\t');
		++pos;
		pos2 = strchr(pos,'\t');
		*(pos2) = 0;
		sprintf(clusPrefix,"%s",pos);
		sprintf(clusFileName,"%s/%c/%s",parsedTreeBaseDir,type,pos);
		processTreeFile(clusFileName,clusPrefix,resFile,mapOrgSeqOrgTaxa);

	}

	resFile.close();

	clusFile.close();
}



void EvalST_MulTrees::processTreeFile(char *clusFileName, char *clusPrefix, ofstream &resFile, UL_HashTable &mapOrgSeqOrgTaxa)
{

	ifstream ipFile(clusFileName);	
	if(!ipFile.is_open()){
		T_LOG<<"Error! Could not open file."<<endl;
		exit(-1);
	}

	unsigned short ipFracNoIntEdLoss =0, ipFracNoTxLoss =0, ipFracNonMulC =0, ipFracMulC =0 ;
	double avgIntEdLoss =0, avgTxLost =0, avgIntEdg=0, tmpRes, avgMul =0, avgIpRes =0, avgOpRes=0;


	ipFile.getline(MB.chBuf,2*MB.TAXA_MAX*MB.TAXA_LEN_MAX);
	unsigned short numTrees = atoi(MB.chBuf);
	unsigned short invalidTreeC =0, starTreeIpC =0;

	for(unsigned char i=0;i<numTrees; ++i){
		MB.reset();
		IO_Stat.reset();

		ipFile.getline(MB.chBuf,2*MB.TAXA_MAX*MB.TAXA_LEN_MAX);
		parseMulNewickST(MB.chBuf,MB.treeBufInd,mapOrgSeqOrgTaxa);
		
		if(i == 0){
			if(IO_Stat.leavesIp == IO_Stat.taxaIp){//not a mul-tree, skip the whole cluster
				ipFile.close();
				return;
			}
			resFile<<clusPrefix<<"\t";
		}
		TP.computeAllSuvValues(MB.treeBufInd);
		unsigned short infoEdgeC = TP.collapseNonInformativeEdges(MB.treeBufInd);
		
		if(infoEdgeC >=1){
			TP.collapseDupEdges(MB.treeBufInd);
			TP.removeNonParticipatingLeaves(MB.treeBufInd);
			TP.removeDupLeaves(MB.treeBufInd);
			TP.remAllRedundantLeaves(MB.treeBufInd);
			TP.computeStats(MB.treeBufInd);
			avgIntEdg+= IO_Stat.internalNodesIp-1;
		}else{	
			if(!IO_Stat.internalNodesIp){
				++starTreeIpC;				
				continue;
			}else{
				avgIntEdg+= IO_Stat.internalNodesIp-1;
				avgIntEdLoss+= 100;
				++invalidTreeC;
				continue;
			}
		}

		avgIpRes+= 100.0*(IO_Stat.internalNodesIp -1)/(IO_Stat.leavesIp -2);
		avgOpRes+= 100.0*(IO_Stat.internalNodesOp -1)/(IO_Stat.leavesOp -2);

		if(IO_Stat.internalNodesIp == IO_Stat.internalNodesOp){
			++ipFracNoIntEdLoss;			
		}else{
			avgIntEdLoss += (100.0*(IO_Stat.internalNodesIp - IO_Stat.internalNodesOp))/(IO_Stat.internalNodesIp-1);//-1 for unrooted tree
		}

		if(IO_Stat.taxaIp == IO_Stat.taxaOp){
			++ipFracNoTxLoss;
		}else{
			avgTxLost += (100.0*(IO_Stat.taxaIp - IO_Stat.taxaOp))/IO_Stat.taxaIp;
		}

		if(IO_Stat.leavesOp == IO_Stat.taxaOp){
			++ipFracNonMulC;
		}else{
			++ipFracMulC;
			avgMul+= IO_Stat.leavesOp- IO_Stat.taxaOp;
		}

		

	}

	numTrees-= starTreeIpC;

	resFile<<"<";
	resFile<<IO_Stat.taxaIp<<"\t";
	if(numTrees && invalidTreeC != numTrees){		
		avgIntEdg/= numTrees;
		tmpRes = avgTxLost/(numTrees-invalidTreeC);
		resFile<<tmpRes<<"\t";
		tmpRes = ipFracNoTxLoss*100.0/(numTrees-invalidTreeC);
		resFile<<tmpRes<<">\t";		
		resFile<<"<";		
		resFile<<avgIntEdg<<"\t";
		tmpRes = avgIntEdLoss/numTrees;
		resFile<<tmpRes<<"\t";	
		tmpRes = ipFracNoIntEdLoss*100.0/numTrees;
		resFile<<tmpRes<<">\t";				
		resFile<<"<";
		resFile<<avgIpRes/(numTrees-invalidTreeC)<<"\t";
		resFile<<avgOpRes/(numTrees-invalidTreeC)<<"\t";
		resFile<<100.0*(avgOpRes-avgIpRes)/avgIpRes<<">\t";
		resFile<<"<";
		tmpRes = ipFracNonMulC*100.0/(numTrees-invalidTreeC);
		resFile<<tmpRes<<"\t";		
		if(ipFracMulC){
			avgMul/= ipFracMulC;
		}
		resFile<<avgMul<<"\t";				
	}else{
		resFile<<"-\t->\t<-\t-\t->\t<-\t-\t->\t<-\t-\t";
	}
		
	tmpRes = invalidTreeC*100.0/(numTrees+starTreeIpC);
	resFile<<tmpRes<<"\t";
	tmpRes = starTreeIpC*100.0/(numTrees+starTreeIpC);
	resFile<<tmpRes<<">\t";	

	resFile<<endl;
	
	ipFile.close();	
}

void EvalST_MulTrees::parseMulNewickST(char *strTree,unsigned int treeInd, UL_HashTable &mapOrgSeqOrgTaxa)
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
	unsigned int mappedTaxaId;
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

			//replace with taxaId - start
			mappedTaxaId = atoi(MB.taxaBuf[MB.taxaBufInd]);
			mappedTaxaId = mapOrgSeqOrgTaxa.get(mappedTaxaId);
			if(UINT_MAX == mappedTaxaId){
				cout<<"ERROR! taxaId not found."<<endl;
				cout<<MB.taxaBuf[MB.taxaBufInd]<<"\t"<<mappedTaxaId<<endl;
				std::cin.get();
				exit(-1);
			}
			sprintf(MB.taxaBuf[MB.taxaBufInd],"%d",mappedTaxaId);
			len = strlen(MB.taxaBuf[MB.taxaBufInd]);
			//replace with taxaId - end

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

void EvalST_MulTrees::loadAllMaps(UL_HashTable &mapOrgSeqOrgTaxa)
{
	ifstream genTaxaOrgIdFile("../TreeIndexing/files/genTaxaOrgIdArr.db",ios::in | ios::binary);
	unsigned int numQueryTaxa;
	genTaxaOrgIdFile.read(reinterpret_cast<char *>(&numQueryTaxa),sizeof(unsigned int));//read number of clusters
	unsigned int * queryTaxaUniqIdArr = new unsigned int[numQueryTaxa];
	for(unsigned int i=0; i<numQueryTaxa; ++i)
	{
		genTaxaOrgIdFile.read(reinterpret_cast<char *>(&queryTaxaUniqIdArr[i]),sizeof(unsigned int));		
	}
	genTaxaOrgIdFile.close();
	

	ifstream genSeqGenTaxaArrFile("../TreeIndexing/files/genSeqGenTaxaArr.db",ios::in | ios::binary);
	unsigned int numSeqIds;
	genSeqGenTaxaArrFile.read(reinterpret_cast<char *>(&numSeqIds),sizeof(unsigned int));//read number of clusters
	unsigned int *seqIdUniqQueryIdArr = new unsigned int[numSeqIds];
	unsigned int queryId;
	for(unsigned int i=0; i<numSeqIds; ++i)
	{
		genSeqGenTaxaArrFile.read(reinterpret_cast<char *>(&queryId),sizeof(unsigned int));
		seqIdUniqQueryIdArr[i] = queryTaxaUniqIdArr[queryId];
	}
	genSeqGenTaxaArrFile.close();
	

	ifstream genSeqOrgIdArrFile("../TreeIndexing/files/genSeqOrgIdArr.db",ios::in | ios::binary);
	unsigned int numTaxa;
	genSeqOrgIdArrFile.read(reinterpret_cast<char *>(&numTaxa),sizeof(unsigned int));//read # of taxa ids	

	unsigned int orgSeqId;
	for(unsigned int i=0; i<numTaxa; ++i)
	{
		genSeqOrgIdArrFile.read(reinterpret_cast<char *>(&orgSeqId),sizeof(unsigned int));//read each of the taxa ids
		mapOrgSeqOrgTaxa.insert(orgSeqId,seqIdUniqQueryIdArr[i]);
	}
	genSeqOrgIdArrFile.close();

	delete []queryTaxaUniqIdArr;
	delete []seqIdUniqQueryIdArr;
}


void EvalST_MulTrees::evalAllBestTrees()
{

	UL_HashTable mapOrgSeqOrgTaxa(12000017);//= 12000017 first prime greater than 12 million for sequence
	initMapsBestTree(mapOrgSeqOrgTaxa);

	unsigned int maxTaxa = 2000;
	StatArr statArr;
	statArr.taxaLostFractionArr = new double[maxTaxa];
	statArr.ipMulTaxaFractionArr = new double[maxTaxa];
	statArr.ipNodeCArr = new unsigned int[maxTaxa];
	statArr.singNodeCArr = new unsigned int[maxTaxa];
	statArr.txCFreqArr_NdCAnalysis = new unsigned short[maxTaxa];
	statArr.leafCFreqArr = new unsigned short[maxTaxa];
	

	for(unsigned int i=0; i<maxTaxa; ++i){
		statArr.taxaLostFractionArr[i] =0;
		statArr.ipMulTaxaFractionArr[i] = 0;
		statArr.ipNodeCArr[i] = 0;
		statArr.singNodeCArr[i] = 0;
		statArr.txCFreqArr_NdCAnalysis[i] = 0;
		statArr.leafCFreqArr[i] =0;
	}
	
	char parsedTreeBaseDir[] ="E:/TreeIndexing/b_trees_for_Multree_analysis_feb02/original";
	
	ifstream clusFile("E:/TreeIndexing/b_trees_for_Multree_analysis_feb02/clusters/Clus_0.db");
	if(!clusFile.is_open()){
		cout<<"ERROR! Could not open clusFile"<<endl;
		exit(-1);
	}

	ofstream resFile("E:/TreeIndexing/b_trees_for_Multree_analysis_feb02/MulTreeAnalysis.txt");
	if(!resFile.is_open()){
		cout<<"ERROR! Could not open resFile."<<endl;
		exit(-1);
	}

	ofstream graphFile("E:/TreeIndexing/b_trees_for_Multree_analysis_feb02/TaxaLostGraph.dat");
	if(!graphFile.is_open()){
		cout<<"ERROR! Could not open graphFile."<<endl;
		exit(-1);
	}

	ofstream nodeCountFile("E:/TreeIndexing/b_trees_for_Multree_analysis_feb02/NodeCountGraph.dat");
	if(!graphFile.is_open()){
		cout<<"ERROR! Could not open nodeCountFile."<<endl;
		exit(-1);
	}



	char *pos;
	resFile.setf(ios::fixed,ios::floatfield);
	resFile.precision(1);

	graphFile.setf(ios::fixed,ios::floatfield);
	graphFile.precision(1);

	nodeCountFile.setf(ios::fixed,ios::floatfield);
	nodeCountFile.precision(1);

	
	
	char dirPrefix[1024], treeFileName[1024];
	unsigned int clusC =0;
	while(true){
		clusFile.getline(MB.chBuf,2*MB.TAXA_MAX*MB.TAXA_LEN_MAX);
		if(0 == clusFile.gcount() || strlen(MB.chBuf) == 0){
			break;
		}
		
		if(clusC%1000 ==0){
			cout<<clusC<<"\t";
		}
		++clusC;
		pos = strchr(MB.chBuf,'\t');
		*(pos) =0;		
		sprintf(dirPrefix,"%s",MB.chBuf);
		++pos;
		sprintf(treeFileName,"%s",pos);
		processBestTreeFile(treeFileName,parsedTreeBaseDir,dirPrefix,resFile,mapOrgSeqOrgTaxa,
			statArr);
	}

	for(unsigned int i=0; i<maxTaxa; ++i){
		if(statArr.leafCFreqArr[i]){
			graphFile<<i<<"\t"<<statArr.leafCFreqArr[i]<<"\t"<<statArr.taxaLostFractionArr[i]/statArr.leafCFreqArr[i]
			<<"\t"<<statArr.ipMulTaxaFractionArr[i]/statArr.leafCFreqArr[i]<<endl;
		}

	}

	for(unsigned int i=0; i<maxTaxa; ++i){
		if(statArr.txCFreqArr_NdCAnalysis[i]){
			nodeCountFile<<i<<"\t"<<statArr.txCFreqArr_NdCAnalysis[i]<<"\t"<<((double)statArr.ipNodeCArr[i])/statArr.txCFreqArr_NdCAnalysis[i]
			<<"\t"<<((double)statArr.singNodeCArr[i])/statArr.txCFreqArr_NdCAnalysis[i]<<endl;
		}		
	}

	resFile.close();
	graphFile.close();
	nodeCountFile.close(); 

	clusFile.close();
}



void EvalST_MulTrees::processBestTreeFile(char *treeFileName, char *baseDir, char *dirPrefix, ofstream &resFile, UL_HashTable &mapOrgSeqOrgTaxa,
										  StatArr &statArr)
{

	char fileName[1024];
	sprintf(fileName,"%s/%s/%s",baseDir,dirPrefix,treeFileName);
	ifstream ipFile(fileName);	
	if(!ipFile.is_open()){
		T_LOG<<"Error! Could not open file."<<endl;
		exit(-1);
	}

	unsigned short ipFracNoIntEdLoss =0, ipFracNoTxLoss =0, ipFracNonMulC =0, ipFracMulC =0 ;
	double avgIntEdLoss =0, avgTxLost =0, avgIntEdg=0, tmpRes, avgMulOp =0, avgMulIp =0, avgIpRes =0, avgOpRes=0,
		avgInfoEdgeIp=0, avgIntNdDeg =0;


	unsigned short invalidTreeC =0, starTreeIpC =0;
	unsigned short numTrees =1;

	for(unsigned char i=0;i<numTrees; ++i){
		MB.reset();
		IO_Stat.reset();

		ipFile.getline(MB.chBuf,2*MB.TAXA_MAX*MB.TAXA_LEN_MAX);
		parseMulNewickBestTree(MB.chBuf,MB.treeBufInd,mapOrgSeqOrgTaxa);

		if(i == 0){
			if(IO_Stat.leavesIp == IO_Stat.taxaIp){//not a mul-tree, skip the whole cluster
				ipFile.close();
				return;
			}
			resFile<<dirPrefix<<"/"<<treeFileName<<"\t";
		}
		TP.computeAllSuvValues(MB.treeBufInd);
		unsigned short infoEdgeC = TP.collapseNonInformativeEdges(MB.treeBufInd);		
		if(infoEdgeC >=1){
			avgInfoEdgeIp+= infoEdgeC;
			TP.collapseDupEdges(MB.treeBufInd);
			TP.removeNonParticipatingLeaves(MB.treeBufInd);
			TP.removeDupLeaves(MB.treeBufInd);
			TP.remAllRedundantLeaves(MB.treeBufInd);
			TP.computeStats(MB.treeBufInd);
			avgIntEdg+= IO_Stat.internalNodesIp-1;
			avgIntNdDeg+= (IO_Stat.leavesIp+IO_Stat.internalNodesIp-1.0)/IO_Stat.internalNodesIp;
		}else{	
			if(!IO_Stat.internalNodesIp){
				++starTreeIpC;				
				continue;
			}else{
				avgIntEdg+= IO_Stat.internalNodesIp-1;
				avgIntNdDeg+= (IO_Stat.leavesIp+IO_Stat.internalNodesIp-1.0)/IO_Stat.internalNodesIp;
				avgIntEdLoss+= 100;
				++invalidTreeC;
				continue;
			}
		}

		avgIpRes+= 100.0*(IO_Stat.internalNodesIp -1)/(IO_Stat.leavesIp -2);
		avgOpRes+= 100.0*(IO_Stat.internalNodesOp -1)/(IO_Stat.leavesOp -2);

		if(IO_Stat.internalNodesIp == IO_Stat.internalNodesOp){
			++ipFracNoIntEdLoss;			
		}else{
			avgIntEdLoss += (100.0*(IO_Stat.internalNodesIp - IO_Stat.internalNodesOp))/(IO_Stat.internalNodesIp-1);//-1 for unrooted tree
		}

		if(IO_Stat.taxaIp == IO_Stat.taxaOp){
			++ipFracNoTxLoss;
		}else{
			avgTxLost += (100.0*(IO_Stat.taxaIp - IO_Stat.taxaOp))/IO_Stat.taxaIp;
		}

		if(IO_Stat.leavesOp == IO_Stat.taxaOp){
			++ipFracNonMulC;
		}else{
			++ipFracMulC;
			avgMulOp+= IO_Stat.mulTaxaOp;
		}
		avgMulIp+= IO_Stat.mulTaxaIp;
	}


	numTrees-= starTreeIpC;

	resFile<<"<";
	resFile<<IO_Stat.taxaIp<<"\t";
	if(numTrees && invalidTreeC != numTrees){		
		avgIntEdg/= numTrees;
		tmpRes = avgTxLost/(numTrees-invalidTreeC);		
		resFile<<tmpRes<<">\t";
		tmpRes = ipFracNoTxLoss*100.0/(numTrees-invalidTreeC);
		//resFile<<tmpRes<<">\t";		
		resFile<<"<";		
		resFile<<avgIntEdg<<"\t";
		tmpRes = avgIntEdg*avgIntEdLoss/(100*numTrees);//number of internal edges lost
		resFile<<tmpRes<<"\t";
		tmpRes = avgIntEdLoss/numTrees;//fraction of internal edges lost
		resFile<<tmpRes<<"\t";
		tmpRes = avgInfoEdgeIp/(numTrees-invalidTreeC);		
		resFile<<tmpRes<<"\t";	
		tmpRes = avgIntEdg*(100-avgIntEdLoss)/(100*numTrees);//number of internal edges in the reduced tree
		resFile<<tmpRes<<">\t";	
		tmpRes = ipFracNoIntEdLoss*100.0/numTrees;
		//resFile<<tmpRes<<">\t";				
		//resFile<<"<";
		//resFile<<avgIpRes/(numTrees-invalidTreeC)<<"\t";
		//resFile<<avgOpRes/(numTrees-invalidTreeC)<<"\t";
		//resFile<<100.0*(avgOpRes-avgIpRes)/avgIpRes<<">\t";		
		resFile<<"<";
		tmpRes = ipFracNonMulC*100.0/(numTrees-invalidTreeC);
		//resFile<<tmpRes<<"\t";		
		if(1 == ipFracNonMulC){
			resFile<<"Y\t";//yes, reduced to singly-labelled
		}else{
			resFile<<"N\t";//No, not reduced to singly-labelled
		}
		if(ipFracMulC){
			avgMulOp/= ipFracMulC;
		}
		
		tmpRes = (avgMulOp*100)/(IO_Stat.taxaOp);
		resFile<<tmpRes<<"\t";	
		if(0 == ipFracNonMulC){//not reduced to a sing-tree
			statArr.taxaLostFractionArr[IO_Stat.taxaIp]+= tmpRes + avgTxLost/(numTrees-invalidTreeC);	
			++statArr.leafCFreqArr[IO_Stat.taxaIp];
		}
		tmpRes = (avgMulIp*100)/(numTrees*IO_Stat.taxaIp);
		if(0 == ipFracNonMulC){//not reduced to a sing-tree
			statArr.ipMulTaxaFractionArr[IO_Stat.taxaIp]+= tmpRes;
		}
		resFile<<tmpRes<<">";

		if(TP.restrictToSingTaxa(MB.treeBufInd,MB.treeBufInd+1)){
			++statArr.txCFreqArr_NdCAnalysis[IO_Stat.taxaIp]; 
			statArr.ipNodeCArr[IO_Stat.taxaIp]+= IO_Stat.leavesIp + IO_Stat.internalNodesIp;
			statArr.singNodeCArr[IO_Stat.taxaIp]+= MB.treeBuf[MB.treeBufInd+1].numNodes;
		}

		
	}else{
		resFile<<"->\t<-\t-\t-\t-\t->\t<-\t-\t->";
	}

	if(numTrees){
		resFile<<"\t<"<<avgIntNdDeg<<">";
	}else{
		resFile<<"\t<->";
	}
		
	tmpRes = invalidTreeC*100.0/(numTrees+starTreeIpC);
	//resFile<<tmpRes<<"\t";
	tmpRes = starTreeIpC*100.0/(numTrees+starTreeIpC);
	//resFile<<tmpRes<<">\t";	
	//resFile<<">";	

	resFile<<endl;
	
	ipFile.close();	
}

void EvalST_MulTrees::parseMulNewickBestTree(char *strTree,unsigned int treeInd, UL_HashTable &mapOrgSeqOrgTaxa)
{
	char *pos1, *pos2;
	size_t len;
	char ch;//tmp char
	unsigned short count =0; //to store the bracket count, init to 0
	unsigned short taxaC=0;//to store taxa count, init to 0
	unsigned short stackC=0;//counter which points to the last stored element in tmpStack, init to 0
	unsigned short symC=0;// to index next avail pos in TmpBuff.symArr, init to 0
	unsigned int *stack = MB.uintStack;
	bool *mulTaxaArr = MB.boolStack;	

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
	unsigned int mappedTaxaId;
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
			pos1+= 2;//to ignore 'gi'
			pos2 = strpbrk(pos1,":,)");//':' is to ignore other nodal attributes
			len = pos2-pos1;
			strncpy(MB.taxaBuf[MB.taxaBufInd],pos1,len);
			MB.taxaBuf[MB.taxaBufInd][len]='\0';//terminate the string

			//replace with taxaId - start
			mappedTaxaId = atoi(MB.taxaBuf[MB.taxaBufInd]);
			mappedTaxaId = mapOrgSeqOrgTaxa.get(mappedTaxaId);
			if(UINT_MAX == mappedTaxaId){
				cout<<"ERROR! taxaId not found."<<endl;
				cout<<MB.taxaBuf[MB.taxaBufInd]<<"\t"<<mappedTaxaId<<endl;
				std::cin.get();
				exit(-1);
			}
			sprintf(MB.taxaBuf[MB.taxaBufInd],"%d",mappedTaxaId);
			len = strlen(MB.taxaBuf[MB.taxaBufInd]);
			//replace with taxaId - end

			MB.taxaLenBuf[MB.taxaBufInd]=len;//store the strlen

			tmpLabel.ind = MB.taxaBufInd;

			insRet = MB.taxaIdSet.insert(tmpLabel);
			if(insRet.second){
				++MB.taxaBufInd;
				mulTaxaArr[insRet.first->ind] = false;
			}else if(!mulTaxaArr[insRet.first->ind]){//this taxa has been encountered second time
				mulTaxaArr[insRet.first->ind] = true;
				++IO_Stat.mulTaxaIp;
			}//else encountered third or more time

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

void EvalST_MulTrees::initMapsBestTree(UL_HashTable &mapOrgSeqOrgTaxa)
{
	ifstream orgSeqOrgTaxaMappingFile("E:/TreeIndexing/b_trees_for_Multree_analysis_feb02/pb.dmp.giti.184");
	char buf[1024];
	char *pos;
	unsigned int orgSeqId, orgTaxaId;
	while(true){
		orgSeqOrgTaxaMappingFile.getline(buf,1024);
		if(0 == orgSeqOrgTaxaMappingFile.gcount() || strlen(buf) == 0){
			break;
		}
		pos = strchr(buf,'\t');
		*(pos) = 0;
		orgSeqId = atoi(buf);
		++pos;
		orgTaxaId = atoi(pos);
		mapOrgSeqOrgTaxa.insert(orgSeqId,orgTaxaId);
	}
	
}

