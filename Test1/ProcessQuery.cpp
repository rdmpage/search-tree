#include "QueryEntities.h"
#include"MRT.h"
#include"Common.h"


ofstream RESPONSE;

void initTestQueryInfo(QueryInfo &queryInfo)
{	
	queryInfo.numOfClus=1;
	queryInfo.clusIdArr[0]=5;
	queryInfo.numOfQTaxa=6;
	queryInfo.qTaxaArr[0]=616753;//t0
	queryInfo.qTaxaArr[1]=616737;//t1
	queryInfo.qTaxaArr[2]=1257591;//t2
	queryInfo.qTaxaArr[3]=1253557;//t3
	queryInfo.qTaxaArr[4]=1446375;//t4
	queryInfo.qTaxaArr[5]=958531;//t5
}

void processQuery()
{
	PrintBuf PB;//some of the used arrays are grouped here

	ClusterPptInfo *clusPptInfoArr = new ClusterPptInfo[tenMillion];
	char (*clusStrIdArr)[21] = new char [tenMillion][21];
	//assuming ti_gi are present
	unsigned int (*clusId_ti_giArr)[2] = new unsigned int[tenMillion][2];

	//array to store query taxa count for each cluster
	unsigned int *clusQueryTaxaCountArr = new unsigned int[tenMillion];


	loadClusterMap(clusPptInfoArr);
	loadClusQueryTaxaCountArr(clusQueryTaxaCountArr);

	loadClusStrIdMap(clusStrIdArr,clusId_ti_giArr);

	gGlobalBuf.freeSavedBuf();

	PB.queryTaxaUniqIdArr = new unsigned int [maxTaxaCount];//maps uniqQueryTaxaIds to queryTaxaIds
	PB.seqIdUniqQueryIdArr = new unsigned int[maxSeqCount];//array which maps uniqSeqIds to their corresponding uniqQueryTaxaIds 
	bool *queryTaxaFlag_G= new bool[maxTaxaCount];//flag whether a queryTaxa is in use
	unsigned int *queryTaxaOverlapArr_G = new unsigned int [tenThousand];//buf to store ovelap with query taxa of current cluster
	unsigned int queryTaxaOverlapArrC_G =0;//count of queryTaxaOverlapArr

	bool *queryTaxaFlag_L= new bool[maxTaxaCount];//flag whether a queryTaxa is in use
	unsigned int *queryTaxaOverlapArr_L = new unsigned int [tenThousand];//buf to store ovelap with query taxa of current cluster
	unsigned int queryTaxaOverlapArrC_L;//count of queryTaxaOverlapArr

	for(unsigned int i=0; i<maxTaxaCount; ++i)
	{
		queryTaxaFlag_G[i] = queryTaxaFlag_L[i] = false;//init queryTaxaFlags to false
	}


	loadUniqQueryTaxaIdArr(PB.queryTaxaUniqIdArr);
	loadSeqIdUniqQueryTaxaId(PB.seqIdUniqQueryIdArr);



	//loads strTaxa and genTaxaId-strTaxa mapping
	loadStrTaxaInfo(PB);//also allocates required memory to relevant structures

	loadNodeIdNameInfo(PB);//loads nodeId-nodeName mapping in ht, clusterName has a component nodeId which is needed-
	//- in string form

	unsigned int clusId, fileId;//tmp vars
	ArrayBufs arrBufs;//group diff buffers needed for sub tree and MRT
	
	for(unsigned int i=0;i<maxSeqCount; ++i)
	{
		arrBufs.taxaArr[i] = USHRT_MAX;//init all values to false
	}
	
	//lot of declarations below limit number of taxa over lap<= maxClusOverlap
	//for MRT
	BiPart *biPartArr = new BiPart[twiceMaxClusOverlap*(maxTreesInClus+1)];//288mb for twiceMaxClusOverlap=2000, -
	//- maxTreesInClus=1000 (<101), and sizeof(BiPart) = 144 bytes
	//twiceMaxClusOverlap becuase (ndI - curTree->nodeInd) can be in [1,twiceMaxClusOverlap]
	Bin *hashTb = new Bin[PRIME_S]();//8mb
	BinNd * binNdArr = new BinNd[maxClusOverlap*maxTreesInClus];
	
	MajTreeNodePtr * ndPtrBuf = new MajTreeNodePtr[maxClusOverlap*2];//there can be maxClusOverlap*2 max majority nodes or clusters

	//for result output
	PB.taxaUniqIdArr = new unsigned int[maxSeqCount];
	PB.charBuf = new char[oneMillion];//used in toString function of phyloTree
	loadGenSeqIdToOrgSeqIdArr(PB.taxaUniqIdArr);

	
	genPrime(arrBufs.prArr,maxClusOverlap);

	unsigned int numOfClusFiles=0;
	while(true)
	{//count number of clusterFiles
		char clusterFileName[257] ;
		sprintf(clusterFileName,"./files/clusters/binary/Cluster%d.db",numOfClusFiles);
		ifstream clusFile(clusterFileName,ios::in | ios::binary);
		if(!clusFile.is_open())
		{
			T_LOG<<"\n Number of cluster files:"<<numOfClusFiles<<endl;
			break;
		}
		++numOfClusFiles;
		clusFile.close();
	}
	//init file handlers for cluster files
	ifstream *clusFileArr = new ifstream[numOfClusFiles];
	for(unsigned int i=0; i<numOfClusFiles; ++i)
	{
		char clusterFileName[257] ;
		sprintf(clusterFileName,"./files/clusters/binary/Cluster%d.db",i);
		clusFileArr[i].open(clusterFileName,ios::in | ios::binary);
	}

	QueryInfo queryInfo;

#ifndef WIN32
	ServerSocket serverSock(SOCK_PATH_SERVER,T_LOG);
#endif

bool reRootWrtRefTree;//if re-rooting w.r.t. ref tree will be done

bool preCleanUp = false; // if buffers need to be reset, as they could not be done so due to retry/abort of previous query operation
int queryFileNum;

	while(true)
	{//this loop denotes operations only relevant to query processing excluding one time memory setup operations

		if(preCleanUp)
		{//check if pre-resetting of certain buffers is required
			for(unsigned int i=0; i< queryInfo.numOfQTaxa; ++i)
			{
				arrBufs.taxaArr[queryInfo.qTaxaArr[i]] = USHRT_MAX;//reset qTaxaArr after genSubtree
			}

			for(unsigned short j=0; j<queryTaxaOverlapArrC_G; ++j)
			{//reset the flags for next use
				queryTaxaFlag_G[queryTaxaOverlapArr_G[j]] = false;
			}
			preCleanUp = false;
		}

#ifdef WIN32
		getQueryFileNumber(queryFileNum);
#else
		getQueryFileNumber(queryFileNum,serverSock);
#endif
		char respFileName[257];
		sprintf(respFileName,"./../common/queryResult/Response2_%d",queryFileNum);
		RESPONSE.open(respFileName);
		if(!RESPONSE.is_open())
		{
			T_LOG<<"Error! Could not open Response"<<queryFileNum<<endl;
			exit(-1);
		}


		loadResFile(queryFileNum,queryInfo);
		//initTestQueryInfo(queryInfo);

		for(unsigned int i=0; i< queryInfo.numOfQTaxa; ++i)
		{
			arrBufs.taxaArr[queryInfo.qTaxaArr[i]] = 0;//init qTaxaArr for genSubtree
		}


		queryInfo.sClusI=gGlobalBuf.clusBufInd;//init source cluster index of queryinfo
		gGlobalBuf.clusBufInd+= queryInfo.numOfClus;//point to next free slot
		queryInfo.tClusI=gGlobalBuf.clusBufInd;//init target cluster index of queryinfo
		gGlobalBuf.clusBufInd+= queryInfo.numOfClus;

		for(unsigned int i=0; i<queryInfo.numOfClus; ++i)
		{//load clusters
			clusId = queryInfo.clusIdArr[i];//cluster id
			fileId = clusPptInfoArr[clusId].fileId;//file id where this cluster is stored
			gClusBuf[queryInfo.sClusI+i].clusterId = clusId;//update clusterId in clusterInfo
			clusFileArr[fileId].seekg(clusPptInfoArr[clusId].addr);//seek to the address in the file
			gClusBuf[queryInfo.sClusI+i].readFromFile(clusFileArr[fileId],queryInfo.maxBST_Count);//read the cluster from the file
			if(gGlobalBuf.memoryBreached())
			{
				break;
			}
		}

		if(gGlobalBuf.memoryBreached())
		{	
			handleMemoryBreached(preCleanUp);
			RESPONSE.close();
#ifdef WIN32
			sendResFileNum(-1);//error
#else
			sendResFileNum(-1,serverSock);//error
#endif
			
			continue;
		}		

		
		for(unsigned int i=0; i<queryInfo.numOfClus; ++i)
		{//generate subtrees
			clusId = gClusBuf[queryInfo.sClusI+i].clusterId;
			reRootWrtRefTree = (clusPptInfoArr[clusId].clusType == ClusterPptInfo::SINGLE_LOCI)?
				true:false;
			gClusBuf[queryInfo.tClusI+i].clusterId = clusId;//init clusterId for target cluster
			genSubtree(gClusBuf[queryInfo.sClusI+i],gClusBuf[queryInfo.tClusI+i],arrBufs,reRootWrtRefTree);
			if(gGlobalBuf.memoryBreached())
			{
				break;
			}

		}
		
		if(gGlobalBuf.memoryBreached())
		{	
			handleMemoryBreached(preCleanUp);
			RESPONSE.close();
#ifdef WIN32
			sendResFileNum(-1);//error
#else
			sendResFileNum(-1,serverSock);//error
#endif			
			continue;
		}

		for(unsigned int i=0; i<queryInfo.numOfClus; ++i)
		{//generate MRT
			clusId = gClusBuf[queryInfo.sClusI+i].clusterId;
			queryInfo.scoreArr[i] = genMRT(gClusBuf[queryInfo.tClusI+i],biPartArr,hashTb,binNdArr
				,arrBufs,ndPtrBuf,queryInfo.maxSupportNumArr[i],queryInfo.uniqNonTrivBpCArr[i]);
			
			if(clusPptInfoArr[clusId].clusType == ClusterPptInfo::SINGLE_LOCI){
				rootMRT(gClusBuf[queryInfo.sClusI+i].clusterAsTree,gClusBuf[queryInfo.tClusI+i].clusterAsTree,biPartArr
					,arrBufs,gClusBuf[queryInfo.tClusI+i].avgSupInd);
			}

			if(i%1 == 0)
			{
				if(gGlobalBuf.memoryBreached())
				{
					break;
				}
			}
		}

		if(gGlobalBuf.memoryBreached())
		{
			handleMemoryBreached(preCleanUp);
			RESPONSE.close();
#ifdef WIN32
			sendResFileNum(-1);//error condition
#else
			sendResFileNum(-1,serverSock);//error condition
#endif
			
			continue;
		}

				

		char resFileName[100] ;
		char treeFileName[257];
		char tableFileName[257];

		queryInfo.treeFileName = "TRUE";
		queryInfo.tableFileName = "TRUE";
		sprintf(resFileName,"./../common/otherResults/FinalResult%d",queryFileNum);
		sprintf(treeFileName,"./../common/queryResult/TreeFile%d",queryFileNum);
		sprintf(tableFileName,"./../common/queryResult/TableFile%d",queryFileNum);

		unsigned short len;
		unsigned int seqId,uniqQueryTaxaId,seqIdCount;
		

		//write table file start
		if(0 != queryInfo.tableFileName.length())
		{
			ofstream tableFile(tableFileName);
			if(!tableFile.is_open())
			{
				RESPONSE<<"ERROR in opening "<<tableFileName<<"!\n";
			}
			else
			{				
				if(queryInfo.cumQueryTaxa)
				{
					queryTaxaOverlapArrC_G =0;
					seqIdCount=0;
					//now find the overlap with original query taxa ids for over all result
					for(unsigned short i=0; i<queryInfo.numOfClus; ++i)
					{
						for(unsigned short j=0; j<gClusBuf[queryInfo.tClusI+i].numOfTaxa; ++j)
						{
							seqId = gTaxaBuf[gClusBuf[queryInfo.tClusI+i].taxaInd+j];
							uniqQueryTaxaId =PB.seqIdUniqQueryIdArr[seqId];
							if(!queryTaxaFlag_G[uniqQueryTaxaId])
							{
								queryTaxaFlag_G[uniqQueryTaxaId] = true;
								queryTaxaOverlapArr_G[queryTaxaOverlapArrC_G++] = uniqQueryTaxaId;
								//tableFile<<queryTaxaUniqIdArr[uniqQueryTaxaId]<<" ";
							}
							if(arrBufs.taxaArr[seqId] != USHRT_MAX)
							{
								arrBufs.taxaArr[seqId] = USHRT_MAX;
								++seqIdCount;
							}
						}
					}

					tableFile<<"Total taxa overlap"<<"\t"<<"Total sequence overlap"<<std::endl;
					tableFile<<queryTaxaOverlapArrC_G<<"\t"<<seqIdCount;
					tableFile<<std::endl;
					for(unsigned short j=0; j<queryTaxaOverlapArrC_G; ++j)
					{//reset the flags for next use in writing the common result file
						queryTaxaFlag_G[queryTaxaOverlapArr_G[j]] = false;
					}
				}


				unsigned int nodeBufI;
				//write table headings
				tableFile<<"Node ID";//1
				tableFile<<"\t"<<"Cluster ID";//2
				tableFile<<"\t"<<"Cluster Type";//3
				tableFile<<"\t"<<"#Sequences";//4
				tableFile<<"\t"<<"#Taxa";//5
				tableFile<<"\t"<<"Taxa Overlap";//6
				tableFile<<"\t"<<"Sequence Overlap";//7
				tableFile<<"\t"<<"Quality score of MRT";//8
				tableFile<<"\t"<<"Max Support in MRT";//9
				tableFile<<"\t"<<"#Non-trivial Splits";//10
				tableFile<<"\t"<<"PD - Mean";//11
				tableFile<<"\t"<<"PD - Variance";//12
				tableFile<<"\t"<<"Phylota ID";//13
				tableFile<<std::endl;
				
				for(unsigned short i=0; i<queryInfo.numOfClus; ++i)
				{
					ClusterInfo &clus =gClusBuf[queryInfo.tClusI+i];
					nodeBufI = PB.nodeIdName_ht.get(clusId_ti_giArr[clus.clusterId][0]);
					tableFile<<&PB.nodeIdNameBuf[nodeBufI];//col # 1 write the node id in name form					
					tableFile<<"\t"<<clusId_ti_giArr[clus.clusterId][1];//col # 2 write the cluster id					
					tableFile<<"\t"<<(unsigned short)clusPptInfoArr[clus.clusterId].clusType;//col # 3 cluster type
					tableFile<<"\t"<<clusPptInfoArr[clus.clusterId].taxaCount;//col # 4 num of sequence ids
					tableFile<<"\t"<<clusQueryTaxaCountArr[clus.clusterId];//col # 5 num of taxa ids

					//FOR TIME BEING LET BELOW REMAIN COMMENTED
					//AS THERE IS NO USE IN CALCULATING OVERLAPPING QUERY TAXA IDS COUNT WHEN QUERY WAS BASED ON SEQ ID
					//now find the overlap with original query taxa ids for each cluster in the result set
					//queryTaxaOverlapArrC_L=0;
					//if(!queryInfo.useSeqId)
					//{//only in case if query was based on query taxa ids
					//	for(unsigned short j=0; j<gClusBuf[queryInfo.tClusI+i].numOfTaxa; ++j)
					//	{
					//		seqId = gTaxaBuf[gClusBuf[queryInfo.tClusI+i].taxaInd+j];
					//		uniqQueryTaxaId =seqIdUniqQueryIdArr[seqId];
					//		if(!queryTaxaFlag_L[uniqQueryTaxaId])
					//		{
					//			queryTaxaFlag_L[uniqQueryTaxaId] = true;
					//			queryTaxaOverlapArr_L[queryTaxaOverlapArrC_L++] = uniqQueryTaxaId;
					//			//resFile<<queryTaxaUniqIdArr[uniqQueryTaxaId]<<" ";
					//		}
					//	}	

					//	for(unsigned short j=0; j<queryTaxaOverlapArrC_L; ++j)
					//	{//reset the flags for next use
					//		queryTaxaFlag_L[queryTaxaOverlapArr_L[j]] = false;
					//	}
					//	tableFile<<"\t"<<queryTaxaOverlapArrC_L;//col # 5 number of overlapping taxa ids =0 if query was based on seq ids
					//}

					tableFile<<"\t"<<queryInfo.taxaIdOverlapArr[i];//col # 6 number of overlapping taxa ids =0 if query was based on seq ids
					tableFile<<"\t"<<queryInfo.seqIdOverlapArr[i];//col # 7 number of overlapping sequence ids
					tableFile<<"\t"<<queryInfo.scoreArr[i];//col # 8 avg support num for the MRT
					tableFile<<"\t"<<queryInfo.maxSupportNumArr[i];//col #9 max support num
					tableFile<<"\t"<<queryInfo.uniqNonTrivBpCArr[i];//col #10 number of non-trivial splits above 50%
					//Below cannot be trusted becuase cannot be trusted as in reRooting subtrees before the MRT, -
					//- in reRooting the MRT, and in degree two check of the MRT (to come soon), numOfNodes is getting messed up.
					//tableFile<<"\t"<<clus.clusterAsTree.numOfNodes - clus.numOfTaxa -1;//col #10 number of splits above 50%
					//=numOfNodes(allSplits) - clus.numOfTaxa (trivial leaf splits) -1(root node as a split)
					tableFile<<"\t"<<clus.pdStat.mean;//col #11 mean pd value of participating subtrees
					tableFile<<"\t"<<clus.pdStat.var;//col #12 var of pd values of participating subtrees
					//tableFile<<"\t"<<clus.pdStat.mrtPD;//col #13 pd value of the MRT
					tableFile<<"\t"<<clusPptInfoArr[clus.clusterId].phylotaId;//col #13 phylota id

					tableFile<<std::endl;
				}
			}
			tableFile.close();
		}
		//write table file end

				

		//write tree file start
		if(0 != queryInfo.treeFileName.length())
		{
			ofstream treeFile(treeFileName);
			
			if(!treeFile.is_open())
			{
				RESPONSE<<"ERROR in opening "<<treeFileName<<".\n";
			}
			else
			{
				treeFile<<"#Nexus\n";
				treeFile<<"Begin trees;\n";
				for(unsigned short i=0; i<queryInfo.numOfClus; ++i)
				{					
					ClusterInfo &clus =gClusBuf[queryInfo.tClusI+i];					
					treeFile<<"Tree Tree_"<<clusId_ti_giArr[clus.clusterId][0]<<"_"<<clusId_ti_giArr[clus.clusterId][1]<<"_T"
						<<(unsigned short)clusPptInfoArr[clus.clusterId].clusType<<" = ";
					gClusBuf[queryInfo.tClusI+i].clusterAsTree.toString(len,gClusBuf[queryInfo.tClusI+i].taxaInd
						,PB,arrBufs.stackBuf,queryInfo.labelFormat,clus.avgSupInd);
					treeFile<<PB.charBuf<<";"<<std::endl;
				}
				treeFile<<";"<<std::endl;
				treeFile<<"End;"<<std::endl;
			}
			treeFile.close();
		}
		//write tree file end

		ofstream resFile(resFileName);
		//ofstream testSubTreeFile("./files/results/testSubTreeFile");//need to comment once testing is done
		queryTaxaOverlapArrC_G =0;

		//now find the overlap with original query taxa ids for over all result
		for(unsigned short i=0; i<queryInfo.numOfClus; ++i)
		{
			for(unsigned short j=0; j<gClusBuf[queryInfo.tClusI+i].numOfTaxa; ++j)
			{
				seqId = gTaxaBuf[gClusBuf[queryInfo.tClusI+i].taxaInd+j];
				uniqQueryTaxaId =PB.seqIdUniqQueryIdArr[seqId];
				if(!queryTaxaFlag_G[uniqQueryTaxaId])
				{
					queryTaxaFlag_G[uniqQueryTaxaId] = true;
					queryTaxaOverlapArr_G[queryTaxaOverlapArrC_G++] = uniqQueryTaxaId;
					resFile<<PB.queryTaxaUniqIdArr[uniqQueryTaxaId]<<" ";
				}
			}
		}
		resFile<<std::endl;

		for(unsigned short i=0; i<queryInfo.numOfClus; ++i)
		{
			resFile<<clusStrIdArr[queryInfo.clusIdArr[i]]<<" "<<queryInfo.seqIdOverlapArr[i]<<" "<<queryInfo.scoreArr[i]<<std::endl;

			//now find the overlap with original query taxa ids for each cluster in the result set
			queryTaxaOverlapArrC_L=0;
			for(unsigned short j=0; j<gClusBuf[queryInfo.tClusI+i].numOfTaxa; ++j)
			{
				seqId = gTaxaBuf[gClusBuf[queryInfo.tClusI+i].taxaInd+j];
				uniqQueryTaxaId =PB.seqIdUniqQueryIdArr[seqId];
				if(!queryTaxaFlag_L[uniqQueryTaxaId])
				{
					queryTaxaFlag_L[uniqQueryTaxaId] = true;
					queryTaxaOverlapArr_L[queryTaxaOverlapArrC_L++] = uniqQueryTaxaId;
					resFile<<PB.queryTaxaUniqIdArr[uniqQueryTaxaId]<<" ";
				}
			}	
			resFile<<std::endl;
			for(unsigned short j=0; j<queryTaxaOverlapArrC_L; ++j)
			{//reset the flags for next use
				queryTaxaFlag_L[queryTaxaOverlapArr_L[j]] = false;
			}

			gClusBuf[queryInfo.tClusI+i].clusterAsTree.toString(len,gClusBuf[queryInfo.tClusI+i].taxaInd
				,PB,arrBufs.stackBuf,queryInfo.labelFormat,gClusBuf[queryInfo.tClusI+i].avgSupInd);
			resFile<<PB.charBuf<<std::endl;

			//test code start
			//for(unsigned short j=0; j< gClusBuf[queryInfo.tClusI+i].numOfTrees; ++j)
			//{
			//	unsigned int ind = gClusBuf[queryInfo.tClusI+i].treeInd+j;
			//	gTreeBuf[ind].toSubtreeString(len,gClusBuf[queryInfo.tClusI+i].taxaInd
			//	,charBuf,taxaUniqIdArr,seqIdUniqQueryIdArr,queryTaxaUniqIdArr,arrBufs.stackBuf,0
			//	,gClusBuf[queryInfo.tClusI+i].avgSupInd);
			//	testSubTreeFile<<charBuf<<";"<<std::endl;
			//}
			//test code end
		}

		//testSubTreeFile.close();
		resFile.close();
		
		RESPONSE.flush();
		RESPONSE.close();//so that client can resume writing


#ifdef WIN32
		sendResFileNum(queryFileNum);
#else
		sendResFileNum(queryFileNum,serverSock);
#endif




		for(unsigned int i=0; i< queryInfo.numOfQTaxa; ++i)
		{
			arrBufs.taxaArr[queryInfo.qTaxaArr[i]] = USHRT_MAX;//reset qTaxaArr after genSubtree
		}
		gGlobalBuf.printStatus();
		gGlobalBuf.check();
		//below is experimental
		gGlobalBuf.clearBuf();//for now on let the buf be cleared
		T_LOG.flush();

		for(unsigned short j=0; j<queryTaxaOverlapArrC_G; ++j)
		{//reset the flags for next use
			queryTaxaFlag_G[queryTaxaOverlapArr_G[j]] = false;
		}


#ifdef WIN32
		break;
		//T_LOG<<"";
#else
		//	break;
#endif



	}

	for(unsigned int i=0; i<numOfClusFiles; ++i)
	{
		clusFileArr[i].close();
	}
	delete []clusFileArr;

	delete [] ndPtrBuf;
	delete [] binNdArr; 
	delete [] hashTb;
	delete [] biPartArr;
	delete [] queryTaxaFlag_G;
	delete [] queryTaxaOverlapArr_G;
	delete [] queryTaxaFlag_L;
	delete [] queryTaxaOverlapArr_L;
	delete [] clusStrIdArr;
	delete [] clusQueryTaxaCountArr;
	delete [] clusId_ti_giArr;
	delete [] clusPptInfoArr;
}

//load strIds for the clusters
void loadClusStrIdMap(char (*clusAddrArr)[21],unsigned int (*clusId_ti_giArr)[2])
{
	unsigned int numClus;
	ifstream treeStrIdMapFile (treeStrIdMapFileName,ios::in | ios::binary);
	//open GenTreeId - StrTreeId mapping file to read
	string line;
	size_t pos;
	istringstream stm_i;//needed to chek if user input being numeric
	if(treeStrIdMapFile.is_open())
	{
		treeStrIdMapFile.read(reinterpret_cast<char *>(&numClus),sizeof(unsigned int));//read numOfTrees
		for(unsigned int tc=0; tc<numClus; ++tc)
		{
			treeStrIdMapFile.read(clusAddrArr[tc],phyloTreeConstants::maxTreeNAmeSize * sizeof(char)); //string based tree id
			//uniq tree ids are represented by numbers starting from 0 and strTreeIds are stored in that order only
			//hence no need to store the uniq tree ids
			clusAddrArr[tc][20]='\0';//intit last char to terminating string char
			line = clusAddrArr[tc];
			pos = line.find("_");//this is not a smile
			stm_i.str(line.substr(0,pos));
			stm_i>>clusId_ti_giArr[tc][0];						
			stm_i.clear();
			stm_i.str(line.substr(pos+1));
			stm_i>>clusId_ti_giArr[tc][1];
			stm_i.clear();
		}
	}
	treeStrIdMapFile.close();//GenTreeId - StrTreeId mapping file written
}



//load geSeqId to uniq query taxa mapping
void loadSeqIdUniqQueryTaxaId(unsigned int *seqIdUniqQueryIdArr)
{
	ifstream genSeqGenTaxaArrFile(genSeqGenTaxaArrFileName,ios::in | ios::binary);
	unsigned int numSeqIds;
	genSeqGenTaxaArrFile.read(reinterpret_cast<char *>(&numSeqIds),sizeof(unsigned int));//read number of clusters
	for(unsigned int i=0; i<numSeqIds; ++i)
	{
		genSeqGenTaxaArrFile.read(reinterpret_cast<char *>(&seqIdUniqQueryIdArr[i]),sizeof(unsigned int));
	}
	genSeqGenTaxaArrFile.close();
}

//load genQueryTaxa to queryTaxa mapping
void loadUniqQueryTaxaIdArr(unsigned int *queryTaxaUniqIdArr)
{

	ifstream genTaxaOrgIdFile(genTaxaOrgIdFileName,ios::in | ios::binary);
	unsigned int numQueryTaxa;
	genTaxaOrgIdFile.read(reinterpret_cast<char *>(&numQueryTaxa),sizeof(unsigned int));//read number of clusters
	for(unsigned int i=0; i<numQueryTaxa; ++i)
	{
		genTaxaOrgIdFile.read(reinterpret_cast<char *>(&queryTaxaUniqIdArr[i]),sizeof(unsigned int));		
	}
	genTaxaOrgIdFile.close();
}

//loads nodeId names
void loadNodeIdNameInfo(PrintBuf &pb)
{
	char buf[1024];
	char *ch;
	unsigned int bufInd=0;
	unsigned int nodeId;
	ifstream queryTaxaNameMapFile (queryTaxaNameMapFileName);
		if(!queryTaxaNameMapFile.is_open()){
			T_LOG<<"Error! Could not open queryTaxaNameMapFile."<<endl;
			return;
		}
		pb.nodeIdNameBuf = new char[twentyMillion];//actual file size is 13.7 mb
		while(true){
			queryTaxaNameMapFile.getline(buf,1023);
			if(0 == queryTaxaNameMapFile.gcount()){
				break;
			}
			ch = strchr(buf,'\t');
			*ch ='\0';
			nodeId = atol(buf);
			++ch;
			strcpy(&pb.nodeIdNameBuf[bufInd],ch);
			pb.nodeIdName_ht.insert(nodeId,bufInd);	
			bufInd+= (unsigned int)strlen(ch)+1;//+1 for the terminating character			
		}								
		queryTaxaNameMapFile.close();
}

//loads strTaxa and genTaxaId-strTaxa mapping
void loadStrTaxaInfo(PrintBuf &pb)
{

	//first we need to read numQueryTaxa
	ifstream genTaxaOrgIdFile(genTaxaOrgIdFileName,ios::in | ios::binary);
	if(!genTaxaOrgIdFile.is_open()){
		T_LOG<<"Error! Could not open genTaxaOrgIdFileName."<<endl;
		return;

	}
	unsigned int genQueryTaxaC;
	genTaxaOrgIdFile.read(reinterpret_cast<char *>(&genQueryTaxaC),sizeof(unsigned int));//read number of query taxa ids
	genTaxaOrgIdFile.close();
	pb.genTaxaIndStrTaxaBufArr = new unsigned int[genQueryTaxaC];
	pb.strTaxaLenArr = new unsigned char[genQueryTaxaC];

	ifstream strTaxaDBFile (strTaxaDBFileName,ios::in | ios::binary);
	if(!strTaxaDBFile.is_open()){
		T_LOG<<"Error! Could not open strTaxaDBFileName."<<endl;
		return;
	}
	
	size_t charBufSz;
	strTaxaDBFile.read(reinterpret_cast<char *>(&charBufSz),sizeof(size_t));
	pb.strTaxaBuf = new char[charBufSz];
	unsigned int strTaxaIdC, genTaxaId, tmp, bufInd=0;
	unsigned char orgTaxaIdC, strLen;
	strTaxaDBFile.read(reinterpret_cast<char *>(&tmp),sizeof(unsigned int));//mapBufSz:ignore, not needed here
	strTaxaDBFile.read(reinterpret_cast<char *>(&strTaxaIdC),sizeof(unsigned int));//# of strIds
	
	for(unsigned int i=0; i<strTaxaIdC; ++i){
		strTaxaDBFile.read(reinterpret_cast<char *>(&strLen),sizeof(unsigned char));
		strTaxaDBFile.read(&pb.strTaxaBuf[bufInd],sizeof(char)*strLen);
		strTaxaDBFile.read(reinterpret_cast<char *>(&orgTaxaIdC),sizeof(unsigned char));
		for(unsigned char j=0; j<orgTaxaIdC; ++j){
			strTaxaDBFile.read(reinterpret_cast<char *>(&tmp),sizeof(unsigned int));//orgTaxaId: ignore
			strTaxaDBFile.read(reinterpret_cast<char *>(&genTaxaId),sizeof(unsigned int));
			pb.genTaxaIndStrTaxaBufArr[genTaxaId]= bufInd;
			pb.strTaxaLenArr[genTaxaId] = strLen;//NOTE: strTaxaLenArr is indexed based on the genTaxaId to which this strTaxaId maps to
		}		
		bufInd+= strLen;
		pb.strTaxaBuf[bufInd++] = '\0';//terminating character and increment for the next free slot
	}
	strTaxaDBFile.close();
}

void loadClusQueryTaxaCountArr(unsigned int *clusQueryTaxaCountArr)
{
	ifstream clusQueryTaxaCountFile(clusQueryTaxaCountFileName,ios::in | ios::binary);
	unsigned int numClus;
	clusQueryTaxaCountFile.read(reinterpret_cast<char *>(&numClus),sizeof(unsigned int));//read number of clusters
	for(unsigned int i=0; i<numClus; ++i)
	{
		clusQueryTaxaCountFile.read(reinterpret_cast<char *>(&clusQueryTaxaCountArr[i]),sizeof(unsigned int));		
	}
	clusQueryTaxaCountFile.close();
}

//load addreass mapping for cluster
unsigned int loadClusterMap(ClusterPptInfo *clusAddrArr)
{
	char clusterAddrFileName[257] ;
	sprintf(clusterAddrFileName,"./files/ClusterMap.db");
	ifstream clusAddrFile(clusterAddrFileName,ios::in | ios::binary);
	unsigned int numClus;
	clusAddrFile.read(reinterpret_cast<char *>(&numClus),sizeof(unsigned int));//read number of clusters
	for(unsigned int i=0; i<numClus; ++i)
	{
		clusAddrFile.read(reinterpret_cast<char *>(&clusAddrArr[i]),sizeof(ClusterPptInfo));
	}
	clusAddrFile.close();

	return numClus;
}

//load query file corresponding to i- which is the result file from part 1
void loadResFile(int i,QueryInfo &query)
{
	char resFileName[257] ;

	sprintf(resFileName,"./../common/otherResults/Result_B%d",i);



	ifstream resFile(resFileName,ios::in | ios::binary);
	//read command line options start
	resFile.read(reinterpret_cast<char *>(&query.useSeqId),sizeof(bool));//if query was based on seq id
	resFile.read(reinterpret_cast<char *>(&query.cumQueryTaxa),sizeof(bool));//if cum query taxa is used

	size_t len;
	resFile.read(reinterpret_cast<char *>(&len),sizeof(size_t));//treeFileName size
	if(0 != len)
	{
		resFile.read(resFileName,(int)len);//read treeFileName
		resFileName[len] ='\0';//terminating character
		query.treeFileName = resFileName;
	}

	resFile.read(reinterpret_cast<char *>(&len),sizeof(size_t));//tableFileName size
	if(0 != len)
	{
		resFile.read(resFileName,(int)len);//read tableFileName
		resFileName[len] ='\0';//terminating character
		query.tableFileName = resFileName;
	}

	resFile.read(reinterpret_cast<char *>(&query.labelFormat),sizeof(unsigned char));//label format

	resFile.read(reinterpret_cast<char *>(&query.maxBST_Count),sizeof(query.maxBST_Count));//# of boot-strap trees to be read
	//read command line options end


	resFile.read(reinterpret_cast<char *>(&query.numOfQTaxa),sizeof(unsigned int));//read # of sequence ids
	for(unsigned int i=0; i<query.numOfQTaxa; ++i)
	{
		resFile.read(reinterpret_cast<char *>(&query.qTaxaArr[i]),sizeof(unsigned int));//read each of sequence ids
	}
	resFile.read(reinterpret_cast<char *>(&query.numOfClus),sizeof(unsigned int));//read # of clusters	
	if(query.numOfClus > tenThousand)
	{
		RESPONSE<<"Number of clusters in the result set > 10,000 limit! Considering only the first 10,000!\n";
		query.numOfClus = tenThousand;
	}
	for(unsigned int i=0; i<query.numOfClus; ++i)
	{
		resFile.read(reinterpret_cast<char *>(&query.clusIdArr[i]),sizeof(unsigned int));//read each of cluster ids
		resFile.read(reinterpret_cast<char *>(&query.taxaIdOverlapArr[i]),sizeof(unsigned int));//# of overlapping taxa ids
		resFile.read(reinterpret_cast<char *>(&query.seqIdOverlapArr[i]),sizeof(unsigned int));//# of overlapping seq ids
		if(query.seqIdOverlapArr[i]>maxClusOverlap)
		{
			T_LOG<<"\nDropping cluster:"<<query.clusIdArr[i]<<" ;Ovelap with seq ids: "<<query.seqIdOverlapArr[i]<<
				" greater than "<<maxClusOverlap<<" limit!"<<endl;
			--i;
			--query.numOfClus;
		}
	}
	resFile.close();
}

void loadGenSeqIdToOrgSeqIdArr(unsigned int *taxaUniqIdArr)
{
	ifstream genSeqOrgIdArrFile(genSeqOrgIdArrFileName,ios::in | ios::binary);
	unsigned int numTaxa;
	genSeqOrgIdArrFile.read(reinterpret_cast<char *>(&numTaxa),sizeof(unsigned int));//read # of taxa ids
	for(unsigned int i=0; i<numTaxa; ++i)
	{
		genSeqOrgIdArrFile.read(reinterpret_cast<char *>(&taxaUniqIdArr[i]),sizeof(unsigned int));//read each of the taxa ids
	}
	genSeqOrgIdArrFile.close();
}

void handleMemoryBreached(bool &preCleanUp)
{
	gGlobalBuf.printStatus();
	gGlobalBuf.clearBuf();
	//if(!retry)
	//{//retry the query
	//	retry = true;//indicate retry
	//	T_LOG<<"ERROR! Memory breached, attempting retry!\n";
	//}
	//else
	//{//retrying failed for the current query
	//	retry = false;
	//	T_LOG<<"ERROR! Retry failed, aborting query permanently!\n";
	//	RESPONSE<<"Insufficient memory to prcoess current query!\n";
	//}

	//no retry as of now
	//T_LOG<<"ERROR! Retry failed, aborting query permanently!\n";
	T_LOG<<"Insufficient memory to prcoess current query!\n";
	RESPONSE<<"Insufficient memory to prcoess current query!"<<endl;
	RESPONSE<<"Retrying the query with a larger taxon overlap may help!"<<endl;
	preCleanUp = true;//indicate pre-resetting of buffers

}

#ifdef WIN32
void getQueryFileNumber(int &queryFileNum)
{
	queryFileNum = 1;
}

void sendResFileNum(int resCount)
{
}

#else
void getQueryFileNumber(int &queryFileNum,ServerSocket &serverSock)
{
	istringstream stm;
	string strLine;
	while(true)
	{
		serverSock.connectTo();
		if(!serverSock.receiveFrom(strLine))
		{
			continue;
		}
		break;
	}
	
	stm.str(strLine);
	stm>>queryFileNum;
	T_LOG<<"\nQueryFileNum: "<<queryFileNum;

	//	queryFileNum = 1;
}

void sendResFileNum(int resCount,ServerSocket &serverSock)
{
	char strBuff[257];
	sprintf(strBuff,"%d",resCount);
	serverSock.sendTo(strBuff);
}
#endif
