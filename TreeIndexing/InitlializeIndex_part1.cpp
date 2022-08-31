#include "TreeIndexingUtilities.h"
#include"StringHash.h"

void initialize_part1()
{
	initializeQueryTaxaMap(); //read taxa mapping, create and store taxa map
	initQueryTaxaNameMap();//create strTaxaId- orgQueryTaxaId mapping
	initGenusTaxaNameMap();//create strGenusId- orgQueryTaxaId mapping

	readTreeDir(); // read tree files
	loadTrees(); // load trees and write them to binary forms
}

void readTreeDir()
{
	for(unsigned int i = 0; i<treeDirCount; ++i)
	{
		char fileName[257] ;
		sprintf(fileName,"%s/Tree%d",treeDirNW,i);
		vecTree.push_back(fileName);
	}
}

void initGenusTaxaNameMap()
{
	size_t charBufSz=0;
	char buf[1024];
	char *ch;
	unsigned int genTaxaId,orgTaxaId,genStrGenusId;//tmp vars

	
	UI_HashTable queryTaxaUniqId_ht(largePrimeForTaxa);//stores mapping from org query taxa to gen query taxa mapping
	String_Hash genusName_ht(largePrimeForTaxa,20*maxTaxaCount);//store strQueryTaxa to org queryTaxId
	unsigned int uniqQueryTaxaC = 0;
	unsigned int strQueryGenusC =0;

	ifstream genTaxaOrgIdFile(genTaxaOrgIdFileName,ios::in | ios::binary);
	if(!genTaxaOrgIdFile.is_open()){
		T_LOG<<"Error! Could not open mappingFileName."<<endl;
		return;
	}

	genTaxaOrgIdFile.read(reinterpret_cast<char *>(&uniqQueryTaxaC),sizeof(unsigned int));//both contain numTaxa
	unsigned int * genQueryTaxaOrgIdArr = new unsigned int[uniqQueryTaxaC];
	unsigned char * strGenusFreqArr = new unsigned char[uniqQueryTaxaC];//though this stores freq of strGenusId, uniqQueryTaxaC gives-
	//- an upper bound on num of strGenusIds 
	unsigned int *strTaxaIndArr = new unsigned int[uniqQueryTaxaC];//gives index into strQueryTaxaOrgTaxaMapBuf (declared later)
	for(unsigned int i=0; i<uniqQueryTaxaC; ++i){
		genTaxaOrgIdFile.read(reinterpret_cast<char *>(&orgTaxaId),sizeof(unsigned int));//both contain numTaxa
		queryTaxaUniqId_ht.insert(orgTaxaId,i);	
		genQueryTaxaOrgIdArr[i] = orgTaxaId;
	}
	genTaxaOrgIdFile.close();


	{
		ifstream queryGenusNameMapFile (queryGenusNameMapFileName);
		if(!queryGenusNameMapFile.is_open()){
			T_LOG<<"Error! Could not open queryGenusNameMapFile."<<endl;
			return;
		}
		while(true){
			queryGenusNameMapFile.getline(buf,1023);
			if(0 == queryGenusNameMapFile.gcount()){
				break;
			}
			ch = strchr(buf,'\t');
			*ch ='\0';
			orgTaxaId = atol(buf);
			++ch;
			genTaxaId = queryTaxaUniqId_ht.get(orgTaxaId);
			if(UINT_MAX != genTaxaId){//this taxaId exists in the one of the trees
				genQueryTaxaOrgIdArr[genTaxaId] = UINT_MAX;
				genStrGenusId = genusName_ht.get(ch);
				if(UINT_MAX != genStrGenusId){//i.e. this strGenusId occurs multiple times
					++strGenusFreqArr[genStrGenusId];//increment its frequency
				}else {//encountered first time
					genusName_ht.insert(ch,strQueryGenusC);
					strGenusFreqArr[strQueryGenusC++]=1;//init its freq to 1
					charBufSz+= strlen(ch)+1;//+1 for '\0'
				}
			}						
		}
		++charBufSz;//as it started with 0
		queryGenusNameMapFile.close();
	}		

	////to take care of the taxa which are not in the mapping file, we write their numerical value as string
	//for(unsigned int i=0; i<uniqQueryTaxaC; ++i){
	//	if(UINT_MAX != genQueryTaxaOrgIdArr[i]){
	//		sprintf(buf,"%d",genQueryTaxaOrgIdArr[i]);
	//		genusName_ht.insert(buf,strQueryGenusC);
	//		strGenusFreqArr[strQueryGenusC++]=1;		
	//		charBufSz+= strlen(buf)+1;//+1 for '\0'}
	//	}
	//}

	unsigned int *strQueryTaxaOrgTaxaMapBuf = new unsigned int[uniqQueryTaxaC];//each strQueryTaxa maps to multiple org query taxa, -
	//- this buf stores that mapping
	unsigned int mapBufInd =0;// index to strQueryTaxaOrgTaxaMapBuf

	for(unsigned int i=0; i<genusName_ht.bufInd; ++i){
		genStrGenusId = genusName_ht.buf[i].val;
		strTaxaIndArr[genStrGenusId] = mapBufInd;//pos in the mapping buf
		mapBufInd+= strGenusFreqArr[genStrGenusId];//point to the next free slot
		strGenusFreqArr[genStrGenusId]=0;//it would incremented in the next read of the file.
	}


	{
		ifstream queryGenusNameMapFile (queryGenusNameMapFileName);
		if(!queryGenusNameMapFile.is_open()){
			T_LOG<<"Error! Could not open queryGenusNameMapFile."<<endl;
			return;
		}
		while(true){
			queryGenusNameMapFile.getline(buf,1023);
			if(0 == queryGenusNameMapFile.gcount()){
				break;
			}
			ch = strchr(buf,'\t');
			*ch ='\0';
			orgTaxaId = atol(buf);
			++ch;
			genTaxaId = queryTaxaUniqId_ht.get(orgTaxaId);
			if(UINT_MAX != genTaxaId){//this taxaId exists in one of the trees -
				//- thus the corresponding strGenusId must have been inserted
				genStrGenusId = genusName_ht.get(ch);
				strQueryTaxaOrgTaxaMapBuf[strTaxaIndArr[genStrGenusId] + strGenusFreqArr[genStrGenusId]++]=
					orgTaxaId;//orgTaxaId as we can always get genTaxaId from genusName_ht but not vice versa
			}						
		}
		queryGenusNameMapFile.close();
	}

	////to take care of the taxa which are not in the mapping file, we write their numerical value as string
	//for(unsigned int i=0; i<uniqQueryTaxaC; ++i){
	//	if(UINT_MAX != genQueryTaxaOrgIdArr[i]){
	//		sprintf(buf,"%d",genQueryTaxaOrgIdArr[i]);
	//		genStrGenusId = genusName_ht.get(buf);
	//		strQueryTaxaOrgTaxaMapBuf[strTaxaIndArr[genStrGenusId] + strGenusFreqArr[genStrGenusId]++]=
	//			genQueryTaxaOrgIdArr[i];//orgTaxaId(genQueryTaxaOrgIdArr[i]) as we can always get genTaxaId from genusName_ht -
	//		//- but not vice versa

	//	}
	//}

	ofstream strGenusDBFile (strGenusDBFileName,ios::out | ios::binary);
	if(!strGenusDBFile.is_open()){
		T_LOG<<"Error! Could not open strTaxaDBFileName."<<endl;
		return;
	}	
	strGenusDBFile.write(reinterpret_cast<char *>(&charBufSz),sizeof(size_t)); // write size of char buf which would be needed to store all taxaIds-
	//- so that the str hash table could be initialized optimally
	strGenusDBFile.write(reinterpret_cast<char *>(&mapBufInd),sizeof(unsigned int));//write size of the buf which will store mapping of strGenusId -
	//- to multiple queryTaxaIds
	strGenusDBFile.write(reinterpret_cast<char *>(&strQueryGenusC),sizeof(unsigned int)); // write number of strGenusIds
	for(unsigned int i=0; i<genusName_ht.bufInd; ++i){
		strGenusDBFile.write(reinterpret_cast<char *>(&genusName_ht.buf[i].strLen),sizeof(unsigned char)); // write strlen
		strGenusDBFile.write(genusName_ht.buf[i].key_strInd,sizeof(char)*genusName_ht.buf[i].strLen); // write strGenusId
		strGenusDBFile.write(reinterpret_cast<char *>(&strGenusFreqArr[i]),sizeof(unsigned char)); //i = genStrGenusId,write count of -
		//- queryTaxaIds, this strGenusId maps to
		for(unsigned char j=0; j<strGenusFreqArr[i]; ++j){
			orgTaxaId = strQueryTaxaOrgTaxaMapBuf[strTaxaIndArr[i]+j];
			strGenusDBFile.write(reinterpret_cast<char *>(&orgTaxaId),sizeof(unsigned int)); // write orgTaxaId
			genTaxaId = queryTaxaUniqId_ht.get(orgTaxaId);
			strGenusDBFile.write(reinterpret_cast<char *>(&genTaxaId),sizeof(unsigned int)); // write genTaxaId
		}
	}
	strGenusDBFile.close();

	delete []strTaxaIndArr;
	delete []strQueryTaxaOrgTaxaMapBuf;
	delete []strGenusFreqArr;
	delete []genQueryTaxaOrgIdArr;
}


void initQueryTaxaNameMap()
{
	size_t charBufSz=0;
	char buf[1024];
	char *ch;
	unsigned int genTaxaId,orgTaxaId,genStrTaxaId;//tmp vars

	
	UI_HashTable queryTaxaUniqId_ht(largePrimeForTaxa);//stores mapping from org query taxa to gen query taxa mapping
	String_Hash taxaName_ht(largePrimeForTaxa,20*maxTaxaCount);//store strQueryTaxa to org queryTaxId
	unsigned int uniqQueryTaxaC = 0;
	unsigned int strQueryTaxaC =0;

	ifstream genTaxaOrgIdFile(genTaxaOrgIdFileName,ios::in | ios::binary);
	if(!genTaxaOrgIdFile.is_open()){
		T_LOG<<"Error! Could not open mappingFileName."<<endl;
		return;
	}

	genTaxaOrgIdFile.read(reinterpret_cast<char *>(&uniqQueryTaxaC),sizeof(unsigned int));//both contain numTaxa
	unsigned int * genQueryTaxaOrgIdArr = new unsigned int[uniqQueryTaxaC];
	unsigned char * strTaxaFreqArr = new unsigned char[uniqQueryTaxaC];//though this stores freq of strTaxaId, uniqQueryTaxaC gives-
	//- an upper bound on num of strQueryTaxaIds 
	unsigned int *strTaxaIndArr = new unsigned int[uniqQueryTaxaC];//gives index into strQueryTaxaOrgTaxaMapBuf (declared later)
	for(unsigned int i=0; i<uniqQueryTaxaC; ++i){
		genTaxaOrgIdFile.read(reinterpret_cast<char *>(&orgTaxaId),sizeof(unsigned int));//both contain numTaxa
		queryTaxaUniqId_ht.insert(orgTaxaId,i);	
		genQueryTaxaOrgIdArr[i] = orgTaxaId;
	}
	genTaxaOrgIdFile.close();


	{
		ifstream queryTaxaNameMapFile (queryTaxaNameMapFileName);
		if(!queryTaxaNameMapFile.is_open()){
			T_LOG<<"Error! Could not open queryTaxaNameMapFile."<<endl;
			return;
		}
		while(true){
			queryTaxaNameMapFile.getline(buf,1023);
			if(0 == queryTaxaNameMapFile.gcount()){
				break;
			}
			ch = strchr(buf,'\t');
			*ch ='\0';
			orgTaxaId = atol(buf);
			++ch;
			genTaxaId = queryTaxaUniqId_ht.get(orgTaxaId);
			if(UINT_MAX != genTaxaId){//this taxaId exists in the one of the trees
				genQueryTaxaOrgIdArr[genTaxaId] = UINT_MAX;
				genStrTaxaId = taxaName_ht.get(ch);
				if(UINT_MAX != genStrTaxaId){//i.e. this strTaxaId occurs multiple times
					++strTaxaFreqArr[genStrTaxaId];//increment its frequency
				}else {//encountered first time
					taxaName_ht.insert(ch,strQueryTaxaC);
					strTaxaFreqArr[strQueryTaxaC++]=1;//init its freq to 1
					charBufSz+= strlen(ch)+1;//+1 for '\0'
				}
			}						
		}
		++charBufSz;//as it started with 0
		queryTaxaNameMapFile.close();
	}		

	//to take care of the taxa which are not in the mapping file, we write their numerical value as string
	for(unsigned int i=0; i<uniqQueryTaxaC; ++i){
		if(UINT_MAX != genQueryTaxaOrgIdArr[i]){
			sprintf(buf,"%d",genQueryTaxaOrgIdArr[i]);
			taxaName_ht.insert(buf,strQueryTaxaC);
			strTaxaFreqArr[strQueryTaxaC++]=1;		
			charBufSz+= strlen(buf)+1;//+1 for '\0'}
		}
	}

	unsigned int *strQueryTaxaOrgTaxaMapBuf = new unsigned int[uniqQueryTaxaC];//each strQueryTaxa maps to multiple org query taxa, -
	//- this buf stores that mapping
	unsigned int mapBufInd =0;// index to strQueryTaxaOrgTaxaMapBuf

	for(unsigned int i=0; i<taxaName_ht.bufInd; ++i){
		genStrTaxaId = taxaName_ht.buf[i].val;
		strTaxaIndArr[genStrTaxaId] = mapBufInd;//pos in the mapping buf
		mapBufInd+= strTaxaFreqArr[genStrTaxaId];//point to the next free slot
		strTaxaFreqArr[genStrTaxaId]=0;//it would incremented in the next read of the file.
	}


	{
		ifstream queryTaxaNameMapFile (queryTaxaNameMapFileName);
		if(!queryTaxaNameMapFile.is_open()){
			T_LOG<<"Error! Could not open queryTaxaNameMapFile."<<endl;
			return;
		}
		while(true){
			queryTaxaNameMapFile.getline(buf,1023);
			if(0 == queryTaxaNameMapFile.gcount()){
				break;
			}
			ch = strchr(buf,'\t');
			*ch ='\0';
			orgTaxaId = atol(buf);
			++ch;
			genTaxaId = queryTaxaUniqId_ht.get(orgTaxaId);
			if(UINT_MAX != genTaxaId){//this taxaId exists in one of the trees -
				//- thus the corresponding strTaxaId must have been inserted
				genStrTaxaId = taxaName_ht.get(ch);
				strQueryTaxaOrgTaxaMapBuf[strTaxaIndArr[genStrTaxaId] + strTaxaFreqArr[genStrTaxaId]++]=
					orgTaxaId;//orgTaxaId as we can always get genTaxaId from taxaName_ht but not vice versa
			}						
		}
		queryTaxaNameMapFile.close();
	}

	//to take care of the taxa which are not in the mapping file, we write their numerical value as string
	for(unsigned int i=0; i<uniqQueryTaxaC; ++i){
		if(UINT_MAX != genQueryTaxaOrgIdArr[i]){
			sprintf(buf,"%d",genQueryTaxaOrgIdArr[i]);
			genStrTaxaId = taxaName_ht.get(buf);
			strQueryTaxaOrgTaxaMapBuf[strTaxaIndArr[genStrTaxaId] + strTaxaFreqArr[genStrTaxaId]++]=
				genQueryTaxaOrgIdArr[i];//orgTaxaId(genQueryTaxaOrgIdArr[i]) as we can always get genTaxaId from taxaName_ht -
			//- but not vice versa

		}
	}

	ofstream strTaxaDBFile (strTaxaDBFileName,ios::out | ios::binary);
	if(!strTaxaDBFile.is_open()){
		T_LOG<<"Error! Could not open strTaxaDBFileName."<<endl;
		return;
	}	
	strTaxaDBFile.write(reinterpret_cast<char *>(&charBufSz),sizeof(size_t)); // write size of char buf which would be needed to store all taxaIds-
	//- so that the str hash table could be initialized optimally
	strTaxaDBFile.write(reinterpret_cast<char *>(&mapBufInd),sizeof(unsigned int));//write size of the buf which will store mapping of strTaxaId -
	//- to multiple queryTaxaIds
	strTaxaDBFile.write(reinterpret_cast<char *>(&strQueryTaxaC),sizeof(unsigned int)); // write number of strTaxaIds
	for(unsigned int i=0; i<taxaName_ht.bufInd; ++i){
		strTaxaDBFile.write(reinterpret_cast<char *>(&taxaName_ht.buf[i].strLen),sizeof(unsigned char)); // write strlen
		strTaxaDBFile.write(taxaName_ht.buf[i].key_strInd,sizeof(char)*taxaName_ht.buf[i].strLen); // write strTaxaId
		strTaxaDBFile.write(reinterpret_cast<char *>(&strTaxaFreqArr[i]),sizeof(unsigned char)); //i = genStrTaxaId,write count of -
		//- queryTaxaIds, this strTaxaId maps to
		for(unsigned char j=0; j<strTaxaFreqArr[i]; ++j){
			orgTaxaId = strQueryTaxaOrgTaxaMapBuf[strTaxaIndArr[i]+j];
			strTaxaDBFile.write(reinterpret_cast<char *>(&orgTaxaId),sizeof(unsigned int)); // write orgTaxaId
			genTaxaId = queryTaxaUniqId_ht.get(orgTaxaId);
			strTaxaDBFile.write(reinterpret_cast<char *>(&genTaxaId),sizeof(unsigned int)); // write genTaxaId
		}
	}
	strTaxaDBFile.close();

	delete []strTaxaIndArr;
	delete []strQueryTaxaOrgTaxaMapBuf;
	delete []strTaxaFreqArr;
	delete []genQueryTaxaOrgIdArr;
}

void loadTrees()
{
	unsigned int * buffArr = new unsigned int[twoFiftyMillion];
	//tree file is constructed in this buffer before it is written

	string *allTreeStrIdArr = new string[tenMillion];
	//string based treeIds are stored here as they are read from the file

	unsigned int *allTreeTaxaCountArr = new unsigned int[tenMillion];
	//number of distinct query taxa occurring in each tree

	bool *queryTaxaFlag = new bool[oneMillion];
	//flag array for each query taxas based on their sequentially generated ids
	for(unsigned int i=0; i<oneMillion; ++i)
	{
		queryTaxaFlag[i] =false;//init all values to false
	}

	unsigned int *queryTaxaArr = new unsigned int[maxTaxaCount];
	//array to store all unique query taxas (not seq ids) corresponding to the sequence id 
	//contained in the tree
	unsigned int queryTaxaC;//counter for queryTaxaArr

	TreeAddrInfo *allTreeArr = new TreeAddrInfo[tenMillion];
	//array to store info of all the trees

	UI_HashTable mapSeqUniqIds(largePrimeForSeq);
	//map mapping seqIds to uniqSeqIds
	unsigned int seqUniqId;//to store value returned from above hash table



	unsigned int *querySeqFreqArrG = new unsigned int[maxSeqCount];
	for(unsigned int i=0; i<maxSeqCount; ++i)
	{
		querySeqFreqArrG[i] = 0;//init freq value to 0
	}


	unsigned int treeIdGenerator =0;
	//sequentially generates unique treeIds

	unsigned int taxaIdGenerator =0;
	//sequentially generates unique taxaIds

	//stores mapping from uniqGenSeqIds to uniqGenQueryTaxaIds
	unsigned int *genSeqGenTaxa = new unsigned int[maxSeqCount];

	//load above mapping
	ifstream genSeqGenTaxaArrFile(genSeqGenTaxaArrFileName,ios::in | ios::binary);
	unsigned int numSeqIds;
	genSeqGenTaxaArrFile.read(reinterpret_cast<char *>(&numSeqIds),sizeof(unsigned int));//read number of clusters
	for(unsigned int i=0; i<numSeqIds; ++i)
	{
		genSeqGenTaxaArrFile.read(reinterpret_cast<char *>(&genSeqGenTaxa[i]),sizeof(unsigned int));
	}
	genSeqGenTaxaArrFile.close();


	//Load seqUniqIdArr as orgSeq-uniqSeqId
	ifstream genSeqOrgIdArrFile (genSeqOrgIdArrFileName,ios::out | ios::binary);
	if(genSeqOrgIdArrFile.is_open())
	{
		unsigned int numOfSeq;
		genSeqOrgIdArrFile.read(reinterpret_cast<char *>(&numOfSeq),sizeof(unsigned int));//read numSeq
		unsigned int seqId;
		for(unsigned int i=0; i<numOfSeq; ++i)
		{
			genSeqOrgIdArrFile.read(reinterpret_cast<char *>(&seqId),sizeof(unsigned int));
			mapSeqUniqIds.insert(seqId,i);
		}
	}
	genSeqOrgIdArrFile.close();//seqId-uniqSeqId map file read


	for(unsigned short i=0; i<vecTree.size(); i++)
	{

		string dataFileName = vecTree[i];
		ifstream dataFile(dataFileName.c_str());
		//open newick tree file for input

		char strTreeFileName[257] ;
		sprintf(strTreeFileName,"%s/Tree%d.db",treeDirBN,i);
		ofstream treeFile (strTreeFileName,ios::out | ios::binary);
		//open binary tree file for output

		string line;//each line from tree file is read into this

		if (dataFile.is_open())
		{
			unsigned int bufArrCount=0;//current pointer to buffer
			unsigned int treeCount = 0;//maintains tree count
			unsigned int prevBuffArrCount =0;//must ; pointer to buffer before it was passed to parsTree
			unsigned int queryTaxaId;//tmp var to store genQueryTaxaId for a sequence id

			while (! dataFile.eof() )
			{
				getline (dataFile,line);
				if(line.length() == 0)
				{
					continue;
				}

				allTreeTaxaCountArr[treeCount] =0;//init query taxa count for this tree to 0
				queryTaxaC =0;// int counter for queryTaxaArr for this tree

				bufArrCount+=2;//once for the unique tree id, then for the node count
				parseTree(line,buffArr,bufArrCount,allTreeStrIdArr[treeIdGenerator]);//parse the line containing tree
				buffArr[prevBuffArrCount] = treeIdGenerator++; //store the tree id, post increment intentional
				buffArr[++prevBuffArrCount] = bufArrCount-prevBuffArrCount-1;//store node count

				//do not delete this, needed for FreqPhyloMiner
				//if(20 == buffArr[prevBuffArrCount]){
				//	cout<<"Found it:"<<treeCount<<endl ;
				//}

				for(++prevBuffArrCount; prevBuffArrCount<bufArrCount; ++prevBuffArrCount)
				{//for each taxa in tree

					seqUniqId = mapSeqUniqIds.get(buffArr[prevBuffArrCount]);
					//seqId need to be replaced with generated seqId
					if(UINT_MAX != seqUniqId)
					{
						buffArr[prevBuffArrCount] = seqUniqId;//replace seqId with generated seqId
						queryTaxaId =genSeqGenTaxa[seqUniqId]; 
						if(!queryTaxaFlag[queryTaxaId]){//a new query taxa id foound for this tree
							++allTreeTaxaCountArr[treeCount];//increment taxa count
							queryTaxaArr[queryTaxaC++] = queryTaxaId;
							queryTaxaFlag[queryTaxaId] = true;//mark this queryTaxaId
						}
					}
					else
					{//uniq id not created yet, this means tree file contains an id which mapping file does not contain
						T_LOG<<std::endl<<"ERROR! UNIQ SEQ ID NOT CREATED YET FOR:"<<buffArr[prevBuffArrCount];
						buffArr[prevBuffArrCount] =0;//assigned 0 seqId
					}

					//update taxa freq count for current file
					++querySeqFreqArrG[buffArr[prevBuffArrCount]];
				}
				++treeCount;//increment the tree count

				for(unsigned int i=0; i<queryTaxaC; ++i)
				{
					queryTaxaFlag[queryTaxaArr[i]] = false;//reset flags for use in next tree
				}
			}
			dataFile.close();//current newick file read completed

			unsigned short nodeCount;//of current tree
			unsigned int uniqTreeID;//of current tree
			bufArrCount = 0;//pointer to file buffer
			unsigned int currTreePos = 0;//position of current write pointer in the tree

			unsigned int tc;//for loop var for tree count
			unsigned int nc;//for loop var for node count


			for(tc=0; tc<treeCount; ++tc)
			{//for each tree
				uniqTreeID = buffArr[bufArrCount++]; //numerical tree id
				allTreeArr[uniqTreeID].fileId = i; //current tree fileID
				allTreeArr[uniqTreeID].addr = currTreePos;//current tree fileAddr
				allTreeArr[uniqTreeID].nodeCount = nodeCount = (unsigned short) buffArr[bufArrCount++]; // node count
				treeFile.write(reinterpret_cast<char *>(&uniqTreeID),sizeof(unsigned int)); // write the numerical tree id
				treeFile.write(reinterpret_cast<char *>(&nodeCount),sizeof(unsigned short));//write the node count
				for(nc=0; nc<nodeCount; ++nc)//write each uniq node id
				{
					treeFile.write(reinterpret_cast<char *>(&buffArr[bufArrCount++]),sizeof(unsigned int));
				}
				//phyloTreeConstants::maxTreeNAmeSize * sizeof(char)+
				currTreePos+=sizeof(unsigned int)+
					sizeof(unsigned short)+ nodeCount*sizeof(unsigned int);//update currTreePos, still faster than tellp()
			}
			treeFile.close();//current tree file written in binary form			

		}

	}

	ofstream treeDirMapFile (treeDirMapFileName,ios::out | ios::binary);
	//open tree dir file to write, this stores where each tree is stored
	if(treeDirMapFile.is_open())
	{
		for(unsigned int tc=0; tc<treeIdGenerator; ++tc)
		{
			treeDirMapFile.write(reinterpret_cast<char *>(&allTreeArr[tc]),sizeof(TreeAddrInfo)); //tree addr info
			//uniq tree ids are represented by numbers starting from 0 and TreeAddrInfos are stored in that order only
			//hence no need to store the uniq tree ids
		}
	}
	treeDirMapFile.close();//tree dir file written

	ofstream treeStrIdMapFile (treeStrIdMapFileName,ios::out | ios::binary);
	//open GenTreeId - StrTreeId mapping file to write
	if(treeStrIdMapFile.is_open())
	{
		treeStrIdMapFile.write(reinterpret_cast<char *>(&treeIdGenerator),sizeof(unsigned int));//write numOfTrees
		for(unsigned int tc=0; tc<treeIdGenerator; ++tc)
		{
			treeStrIdMapFile.write(allTreeStrIdArr[tc].c_str(),phyloTreeConstants::maxTreeNAmeSize * sizeof(char)); //string based tree id
			//uniq tree ids are represented by numbers starting from 0 and strTreeIds are stored in that order only
			//hence no need to store the uniq tree ids
		}
	}
	treeStrIdMapFile.close();//GenTreeId - StrTreeId mapping file written

	ofstream querySeqFreqArrG_File (querySeqFreqArrG_FileName,ios::out | ios::binary);
	//open global taxa freq file for o/p
	if(querySeqFreqArrG_File.is_open())
	{
		querySeqFreqArrG_File.write(reinterpret_cast<char *>(&numSeqIds),sizeof(unsigned int)); //write numOfSeqIds
		for(unsigned int fc=0; fc<numSeqIds; ++fc)
		{
			querySeqFreqArrG_File.write(reinterpret_cast<char *>(&querySeqFreqArrG[fc]),sizeof(unsigned int)); //write frequency
		}
	}
	querySeqFreqArrG_File.close();//global taxa freq file written


	ofstream clusQueryTaxaCountFile(clusQueryTaxaCountFileName,ios::out | ios::binary);
	//open queryTaxaCount for each cluster file for o/p
	if(clusQueryTaxaCountFile.is_open())
	{
		clusQueryTaxaCountFile.write(reinterpret_cast<char *>(&treeIdGenerator),sizeof(unsigned int));//write numOfClusters
		for(unsigned int i=0; i<treeIdGenerator; ++i)
		{//now write query taxa count for each tree
			clusQueryTaxaCountFile.write(reinterpret_cast<char *>(&allTreeTaxaCountArr[i]),sizeof(unsigned int));
		}
	}
	clusQueryTaxaCountFile.close();
	//queryTaxaCount for each cluster file written


	delete []querySeqFreqArrG;
	delete []genSeqGenTaxa;
	delete []allTreeArr;
	delete []queryTaxaArr;
	delete []queryTaxaFlag;
	delete []allTreeTaxaCountArr;
	delete []allTreeStrIdArr;
	delete []buffArr;
}

//parses a tree from newick file
void parseTree(string &strTree, unsigned int *buffArr, unsigned int &tbufArrCount, string &treeId)
{
	static istringstream stm;
	static size_t pos1, pos2;
	static unsigned int taxaId;

	pos1 = strTree.find_first_of(" \t");//skip clus type single-loci, biclique or quasi-biclique
	pos1 = strTree.find_first_of(" \t",pos1+1);//skip phylota-id
	pos1 = strTree.find_first_of(" \t",pos1+1);//skip base dir
	pos2 = strTree.find_first_of(" \t",pos1+1);
	//treeId = phyloTreeConstants::suffixString.substr(0,20 - (pos1 + 1) ) + strTree.substr(0,pos1+1) ;
	treeId = phyloTreeConstants::suffixString.substr(0,20 - (pos2 - pos1 -1) ) + strTree.substr(pos1+1,pos2-pos1-1) ;
	//parse the strring based tree id
	pos1 = pos2;

	while(true)
	{
		pos1 = strTree.find_first_of("0123456789",pos1);
		if(string::npos == pos1){
			break;
		}
		pos2 = strTree.find_first_not_of("0123456789",pos1);
		stm.str(strTree.substr(pos1,pos2 - pos1));
		stm >> taxaId;
		stm.clear();
		buffArr[tbufArrCount++] = taxaId;
		//store the node id in the buffer
		pos1 = pos2;
	}

}

void initializeQueryTaxaMap()
{

	UI_HashTable mapQueryTaxaUniqId(largePrimeForTaxa);
	unsigned int uniqTaxaId ;//to tmp store the uniqTaxaId returned by above map
	//maps queryTaxaId to its uniqGenQueryTaxaId

	unsigned int *queryTaxaFreqArr = new unsigned int[maxSeqCount];
	//map to store # of sequence ids a query taxa is mapped to

	unsigned int *queryTaxaUniqIdArr = new unsigned int [maxTaxaCount];
	//maps uniqQueryTaxaIds to queryTaxaIds

	unsigned int *genSeqGenTaxa = new unsigned int[maxSeqCount];
	//array which maps uniqSeqIds to their corresponding uniqQueryTaxaIds

	unsigned int * queryTaxaColArr = new unsigned int[maxSeqCount];
	//array to store query taxa column in the table read from taxa-sequence id mapping file

	unsigned int * querySeqColArr = new unsigned int[maxSeqCount];
	//array to store sequence id column in the table read from taxa-sequence id mapping file

	unsigned int * fileBuff = new unsigned int[fifteenMillion];
	//binary mapping file,containing all sequence ids for every query taxa, is written in this buffer first

	unsigned int *seqUniqIdArr = new unsigned int[maxSeqCount];
	//array maps uniq seqIds to org seqIds

	unsigned int seqCount =0;//must ; keeps count of no of rows read
	unsigned int buffCount =0;//must ; pointer to file buffer

	ifstream mappingFile(mappingFileName);
	//open sequence id-query taxa id mapping file to read

	ofstream queryTaxaMapFile(queryTaxaMapFileName, ios::out | ios::binary);
	//open  query taxa-sequence id mapping file to write

	ofstream queryTaxaDirFile(queryTaxaDirFileName, ios::out | ios::binary);
	//open query taxa dir file to write

	string line;//each line from mapping file is read into this

	unsigned int tempSwap = 0;

	UI_HashTable mapSeqUniqIds(largePrimeForSeq);
	//maps querySeq Ids to uniquely generated seqIds
	unsigned int uniqSeqId;//to store uniqSeqId returned by above map

	unsigned int queryTaxaIdGen =0;//generates unix queryTaxaIds
	unsigned int seqIdGen = 0;//generates uniq seqIds

	if(queryTaxaMapFile.is_open() && mappingFile.is_open() && queryTaxaDirFile.is_open())
	{
		while (! mappingFile.eof() )
		{
			getline(mappingFile,line);
			if(line.length() == 0)
			{
				continue;
			}
			parseTaxa(line,queryTaxaColArr[seqCount],querySeqColArr[seqCount]);//parses a taxa,sequence id pair
			uniqTaxaId = mapQueryTaxaUniqId.get(queryTaxaColArr[seqCount]);
			if(UINT_MAX == uniqTaxaId)//update freq map
			{
				mapQueryTaxaUniqId.insert(queryTaxaColArr[seqCount],queryTaxaIdGen);//insert in taxa hashTable
				queryTaxaFreqArr[queryTaxaIdGen] = 1;//set its frequency to 1
				queryTaxaUniqIdArr[queryTaxaIdGen] = queryTaxaColArr[seqCount];//store the uniqQueryTaxaId to queryTaxaId Mapping
				queryTaxaColArr[seqCount] = queryTaxaIdGen;//replace the queryTaxaId with uqniGenQueryTaxaId, -
				//- proves helpful while filling in genSeqGenTaxa
				++queryTaxaIdGen;//increment for next taxa
			}
			else
			{
				queryTaxaColArr[seqCount] = uniqTaxaId;//replace the queryTaxaId with uqniGenQueryTaxaId, -
				//- proves helpful while filling in genSeqGenTaxa
				++queryTaxaFreqArr[uniqTaxaId];//increment frequency i.e. number of seqIds this taxa is mapped to
			}
			++seqCount;//increment for next loop
		}

		queryTaxaDirFile.write(reinterpret_cast<char *>(&queryTaxaIdGen),sizeof(unsigned int));//write number of taxa
		for(unsigned int i=0; i<queryTaxaIdGen; ++i)
		{
			fileBuff[buffCount] = 0;//freq of a query taxa, initialized to 0
			tempSwap = queryTaxaFreqArr[i];//temp store the freq
			queryTaxaFreqArr[i] = buffCount; //store ptr in the fileBuff for this taxa
			queryTaxaDirFile.write(reinterpret_cast<char *>(&queryTaxaFreqArr[i]),sizeof(unsigned int));
			buffCount+= tempSwap+1;//increment it for next query id
		}

		queryTaxaDirFile.close(); //query taxa dir file written
		mappingFile.close();//read file query taxa-sequence id mapping closed

		//write the query taxa-sequence id mapping file in buffer
		for(unsigned int tc=0; tc<seqCount; ++tc)
		{

			//replace querySeqId with newly genSeqId or already existing seqId
			uniqSeqId = mapSeqUniqIds.get(querySeqColArr[tc]);
			if(UINT_MAX == uniqSeqId)
			{//uniq id has not been created for this taxa
				mapSeqUniqIds.insert(querySeqColArr[tc],seqIdGen);//instert this seq in the seq hash table
				seqUniqIdArr[seqIdGen] = querySeqColArr[tc];//insert query taxa id, in the uniq taxaIds to taxaIds mapping array
				querySeqColArr[tc]= seqIdGen;//replace taxa id with the uniq id
				genSeqGenTaxa[seqIdGen]=queryTaxaColArr[tc];// map uniqSeqIds to their corresponding uniqQueryTaxaIds -
				//- queryTaxaColArr[tc] has already been replaced with uniqQueryTaxaId
				++seqIdGen;//increment it for the next new seq
			}
			else
			{//should not happen as this means, that a seq id is mapped to multiple taxa id
				T_LOG<<std::endl<<"WARNING! SEQID:"<<querySeqColArr[tc]<<
					" MAPPED TO MULTIPLE TAXAIDS!";
				querySeqColArr[tc]=uniqSeqId;//replace taxa id with the uniq id
			}

			unsigned int tmpInd = queryTaxaFreqArr[queryTaxaColArr[tc]];
			//- queryTaxaColArr[tc] has already been replaced with uniqQueryTaxaId
			//here queryTaxaFreqArr contains ind for a queryTaxa in queryTaxa-seqId mapping fileBuf
			++fileBuff[tmpInd];//increment freq			
			fileBuff[tmpInd+ fileBuff[tmpInd]]= querySeqColArr[tc];//write the uqniSeqId in the buff

		}

		//write the buffer as file
		queryTaxaMapFile.write(reinterpret_cast<char *>(&buffCount),sizeof(unsigned int));
		//write the size of buffer
		for(unsigned int tc=0; tc<buffCount; ++tc)
		{//now write the buffer
			queryTaxaMapFile.write(reinterpret_cast<char *>(&fileBuff[tc]),sizeof(unsigned int));
		}
		queryTaxaMapFile.close();// query taxa-sequence id mapping file written

		//ofstream taxaIDMapDumpFile (taxaIDMapDumpFileName,ios::out | ios::binary);
		////open taxa id -unique id map file for write
		//if(taxaIDMapDumpFile.is_open())
		//{
		//	for(itrMapU=mapSeqUniqIds.begin(); itrMapU != mapSeqUniqIds.end(); ++itrMapU)
		//	{
		//		taxaIDMapDumpFile.write(reinterpret_cast<char *>(&(*itrMapU)),sizeof(pair<unsigned int, unsigned int>)); 
		//	}
		//}
		//taxaIDMapDumpFile.close();//taxa id-uniqTaxaId map file written

		ofstream genSeqOrgIdArrFile (genSeqOrgIdArrFileName,ios::out | ios::binary);
		//open uniq seqId to orgSeqId array file for o/p
		if(genSeqOrgIdArrFile.is_open())
		{
			genSeqOrgIdArrFile.write(reinterpret_cast<char *>(&(seqIdGen)),sizeof(unsigned int)); //write # of taxaIds
			for(unsigned int tc =0; tc<seqIdGen; ++tc)
			{
				genSeqOrgIdArrFile.write(reinterpret_cast<char *>(&(seqUniqIdArr[tc])),sizeof(unsigned int)); //write taxaId
			}
		}
		genSeqOrgIdArrFile.close();//taxa uniqTaxaId array file written

		ofstream genTaxaOrgIdFile (genTaxaOrgIdFileName,ios::out | ios::binary);
		//open queryTaxaId to uniqQueryTaxaId array file for o/p
		if(genTaxaOrgIdFile.is_open())
		{
			genTaxaOrgIdFile.write(reinterpret_cast<char *>(&(queryTaxaIdGen)),sizeof(unsigned int)); //write # of queryTaxaIds
			for(unsigned int qtc=0; qtc<queryTaxaIdGen; ++qtc)
			{
				genTaxaOrgIdFile.write(reinterpret_cast<char *>(&(queryTaxaUniqIdArr[qtc])),sizeof(unsigned int)); //write taxaId
			}
		}
		genTaxaOrgIdFile.close();//queryTaxaId to uniqQueryTaxaId array file written

		ofstream genSeqGenTaxaArrFile (genSeqGenTaxaArrFileName,ios::out | ios::binary);
		//open seqId to uniqQueryTaxaId array file for o/p
		if(genSeqGenTaxaArrFile.is_open())
		{
			genSeqGenTaxaArrFile.write(reinterpret_cast<char *>(&(seqIdGen)),sizeof(unsigned int)); //write # of seqIds
			for(unsigned int tc=0; tc<seqIdGen; ++tc)
			{
				genSeqGenTaxaArrFile.write(reinterpret_cast<char *>(&(genSeqGenTaxa[tc])),sizeof(unsigned int)); //write taxaId
			}
		}
		genSeqGenTaxaArrFile.close();//seqId to uniqQueryTaxaId array file written

	}


	delete [] queryTaxaFreqArr;
	delete [] seqUniqIdArr;
	delete [] queryTaxaColArr;
	delete [] querySeqColArr;
	delete [] fileBuff;
	delete [] genSeqGenTaxa;
	delete [] queryTaxaUniqIdArr;
}

//parses a taxa,sequence id pair
void parseTaxa(string &strTaxa, unsigned int &queryId, unsigned int &taxaId)
{
	static size_t currPos;
	static size_t endPos;
	static istringstream stm;

	currPos = strTaxa.find_first_not_of("0123456789",0);

	stm.str(strTaxa.substr(0,currPos));
	stm >> taxaId ;
	stm.clear();
	currPos = strTaxa.find_first_of("0123456789",currPos);
	endPos = strTaxa.find_first_not_of("0123456789",currPos);
	if( string::npos == endPos )
	{ 
		endPos = strTaxa.length();
	}
	stm.str(strTaxa.substr(currPos,endPos-currPos));
	stm >> queryId ;
	stm.clear();
}
