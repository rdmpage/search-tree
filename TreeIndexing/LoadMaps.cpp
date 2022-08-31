#include "TreeIndexingUtilities.h"


void loadMaps(AllMaps &allMaps)
{
	allMaps.initQueryTaxaDirHashTable();
	loadQueryTaxaBuf(allMaps.queryTaxMapBuf);
	loadVocab(allMaps.vocabDB);
	loadInvList(allMaps.invListDB);
	loadMapTreeDir(allMaps.arrAllTree);
	loadSeqIdUniqQueryTaxaId(allMaps.seqIdUniqQueryIdArr);
	allMaps.initSeqHashTable();
	allMaps.initStrTaxaHashTb();

	//initializeTreeFileHandlers();
	//initializeInvListFileHandlers();
}


//load seqId to uniq query taxa mapping
void loadSeqIdUniqQueryTaxaId(unsigned int *&seqIdUniqQueryIdArr)
{
	ifstream genSeqGenTaxaArrFile(genSeqGenTaxaArrFileName,ios::in | ios::binary);
	unsigned int numSeqIds;
	genSeqGenTaxaArrFile.read(reinterpret_cast<char *>(&numSeqIds),sizeof(unsigned int));//read number of clusters
	seqIdUniqQueryIdArr = new unsigned int[numSeqIds];
	for(unsigned int i=0; i<numSeqIds; ++i)
	{
		genSeqGenTaxaArrFile.read(reinterpret_cast<char *>(&seqIdUniqQueryIdArr[i]),sizeof(unsigned int));
	}
	genSeqGenTaxaArrFile.close();
}


void loadQueryTaxaBuf(unsigned int *&queryTaxMapBuf)
{
	ifstream queryTaxaMapFile(queryTaxaMapFileName, ios::in | ios::binary);
	if(queryTaxaMapFile.is_open())
	{
		unsigned int bufSize=0;//to store buf size
		queryTaxaMapFile.read(reinterpret_cast<char *>(&bufSize),sizeof(unsigned int));//read buf size	
		queryTaxMapBuf = new unsigned int[bufSize];//allocate memory to buf
		for(unsigned int i=0; i<bufSize; ++i)
		{//read the buufer as array in memory
			queryTaxaMapFile.read(reinterpret_cast<char *>(&queryTaxMapBuf[i]),sizeof(unsigned int));	
		}
	}
	queryTaxaMapFile.close();
}

void loadVocab(VocabDB &vocabDB)
{
	ifstream vocabFile(vocabFileName, ios::in | ios::binary);
	if(vocabFile.is_open())
	{
		vocabDB.taxaCount=0;
		vocabFile.read(reinterpret_cast<char *>(&vocabDB.taxaCount),sizeof(unsigned int));//read taxaCount
		vocabDB.vocabBuf = new unsigned int[vocabDB.taxaCount];//assign memory to vocabBuf
		for(unsigned int i=0; i<vocabDB.taxaCount; ++i)
		{//read the buufer as array in memory
			vocabFile.read(reinterpret_cast<char *>(&vocabDB.vocabBuf[i]),sizeof(unsigned int));	
		}
	}
	vocabFile.close();
	
}

void loadInvList(InvListDB &invListDB)
{
	invListDB.size=0;//to store size of invListBuf
	unsigned int ind=0;//used as common index while reading both the files

	{// read invList1.db
		char invListFile1Name[257] ;
		
		sprintf(invListFile1Name,"%s/invList%d.db",invListBaseDir,1);
		ifstream invListFile1(invListFile1Name,ios::in | ios::binary);
		if(invListFile1.is_open())
		{
			invListFile1.read(reinterpret_cast<char *>(&invListDB.size),sizeof(unsigned int));//read taxaCount
			invListDB.invListBuf = new unsigned int[invListDB.size];//assign memory to vocabBuf
			unsigned int midInd = invListDB.size/2;//partitioning the invLists into halves
			for(ind=0; ind<midInd; ++ind)
			{//read the buufer as array in memory
				invListFile1.read(reinterpret_cast<char *>(&invListDB.invListBuf[ind]),sizeof(unsigned int));	
			}
		}
	
		invListFile1.close();
	}
	{// read invList2.db
		char invListFile2Name[257] ;

		sprintf(invListFile2Name,"%s/invList%d.db",invListBaseDir,2);
		ifstream invListFile2(invListFile2Name,ios::in | ios::binary);
		if(invListFile2.is_open())
		{
			for(ind=invListDB.size/2; ind<invListDB.size; ++ind)
			{//read the buufer as array in memory
				invListFile2.read(reinterpret_cast<char *>(&invListDB.invListBuf[ind]),sizeof(unsigned int));	
			}
		}
		invListFile2.close();
	}

}

void loadMapTreeDir(TreeAddrInfo* &treeArr)
{
	unsigned int tc=0;//to count # of trees fetched- for check
	ifstream treeDirMapFile(treeDirMapFileName, ios::in | ios::binary);
	if(treeDirMapFile.is_open())
	{
		treeArr = new TreeAddrInfo[tenMillion];//allocate memory to allTreeArr
		TreeAddrInfo treeAddrInfo;
		while (true)
		{
			treeDirMapFile.read(reinterpret_cast<char *>(&treeAddrInfo),sizeof(treeAddrInfo)); //uniqTreeId
			if(treeDirMapFile.eof())
			{
				break;
			}
			treeArr[tc++] = treeAddrInfo;//post increment intentional
			//uniq tree ids are represented by numbers starting from 0 and TreeAddrInfos are stored in that order only
			//hence no need to store the uniq tree ids, hence treeIds are represented here by tc
			
		}
	}
	T_LOG<<"Tree Count :"<<tc<<"\n";
	treeDirMapFile.close();
}

void finalizeFileHandlers()
{
	for(unsigned int tc=0; tc<treeDirCount ; tc++)
	{
		vecInvListHandler[tc]->close();
		vecTreeHandler[tc]->close();
	}
	for(unsigned int ic=0; ic<vecInvListDirHandler.size() ; ic++)
	{
		vecInvListDirHandler[ic]->close();
	}
}

void initializeInvListFileHandlers()
{
	for(unsigned int tc=0;; tc++)
	{
		char invListFileName[257] ;
		sprintf(invListFileName,"%s/invList%d.db",invListBaseDir,tc);
		ifstream *invListFile = new ifstream(invListFileName,ios::in | ios::binary);
		if(!invListFile->is_open())
		{
			T_LOG<<"\n invListFile count stopped at : "<<tc;
			return;
		}
		vecInvListHandler.push_back(invListFile);
	}
}

void initializeTreeFileHandlers()
{
	for(unsigned int tc=0; tc<treeDirCount ; tc++)
	{
		char treeFileName[257] ;
		sprintf(treeFileName,"%s/Tree%d.db",treeDirBN,tc);
		ifstream *treeFile =new ifstream(treeFileName, ios::in | ios::binary);
		if(!treeFile->is_open())
		{
			T_LOG<<"\ntree file could not be opened"<<vecTree[tc];
			exit(-1);
		}
		vecTreeHandler.push_back(treeFile);
	}
}


unsigned int UI_HashTable::get(unsigned int key)
{
	unsigned int rem,ind=0,count=0,quo=key;
	while(0 != quo)
	{//calculate key
		rem = quo%divFactor;
		ind = (rem*prArr[count++] +ind)%tableSize;
		quo = quo/divFactor;
	}

	//find the value in the chained list
	for(ind=bin[ind][0]; ind!= UINT_MAX; ind=buf[ind][2])
	{
		if(buf[ind][0] == key)
		{
			return buf[ind][1];
		}
	}
	return ind;//this will return UINT_MAX i.e. when the key was not found

}

unsigned int * UI_HashTable::getPtr(unsigned int key)
{
	unsigned int rem,ind=0,count=0,quo=key;
	while(0 != quo)
	{//calculate key
		rem = quo%divFactor;
		ind = (rem*prArr[count++] +ind)%tableSize;
		quo = quo/divFactor;
	}

	//find the value in the chained list
	for(ind=bin[ind][0]; ind!= UINT_MAX; ind=buf[ind][2])
	{
		if(buf[ind][0] == key)
		{
			return &buf[ind][1];
		}
	}
	return NULL;//this will return UINT_MAX i.e. when the key was not found
}

void UI_HashTable::insert(unsigned int key,unsigned int val)
{
	unsigned int rem,ind=0,count=0,quo=key;

	while(0 != quo)
	{
		rem = quo%divFactor;
		ind = (rem*prArr[count++] +ind)%tableSize;
		quo = quo/divFactor;
	}
	if(0 == bin[ind][1])
	{//no element in this bin till now
		bin[ind][0] = bufInd;
		++bin[ind][1];
		buf[bufInd][0]=key;//store key
		buf[bufInd][1]=val;//store val
		buf[bufInd][2]=UINT_MAX;//making ptr to next node null
		++bufInd;
	}
	else
	{//find next free slot
		++clashCount;
		++bin[ind][1];

		for(ind=bin[ind][0]; buf[ind][2] != UINT_MAX;)
		{
			ind=buf[ind][2];
		}

		buf[ind][2] = bufInd;
		buf[bufInd][0] = key;
		buf[bufInd][1] = val;
		buf[bufInd][2] = UINT_MAX;//making prt to next node null
		++bufInd;
	}
}

void UI_HashTable::genPrime()
{
	unsigned int prC=0;//we will never need more than 20 prime numbers
	bool validPr;
	unsigned int curPr,rem,quo;
	set<unsigned int> uniqPr;
	pair<set<unsigned int>::iterator,bool> ret;
	for(unsigned i=0; i<numPrimes; ++i)
	{
		validPr = false;
		while(!validPr)
		{
			curPr = rand()+1;
			quo = tableSize/curPr;
			rem = tableSize%curPr;
			quo = rand()%(quo+1);
			rem = rand()%(rem+1);
			curPr = curPr*quo + rem;
			if(curPr < tableSize)
			{
				ret = uniqPr.insert(curPr);
				validPr = ret.second;
			}

		}
		prArr[i] = curPr;
	}
}

void AllMaps::initQueryTaxaDirHashTable()
{
	ifstream queryTaxaDirFile(queryTaxaDirFileName, ios::in | ios::binary);
	ifstream genTaxaOrgIdFile(genTaxaOrgIdFileName, ios::in | ios::binary);
	
	if(queryTaxaDirFile.is_open())
	{
		unsigned int posInBuf, numTaxa, orgTaxaId;
		queryTaxaDirFile.read(reinterpret_cast<char *>(&numTaxa),sizeof(unsigned int));	
		genTaxaOrgIdFile.read(reinterpret_cast<char *>(&numTaxa),sizeof(unsigned int));//both contain numTaxa
		T_LOG<<"Num of query taxa: "<<numTaxa<<endl;
		for(unsigned int i=0; i<numTaxa; ++i)
		{
			queryTaxaDirFile.read(reinterpret_cast<char *>(&posInBuf),sizeof(unsigned int));	
			genTaxaOrgIdFile.read(reinterpret_cast<char *>(&orgTaxaId),sizeof(unsigned int));//read org taxa Id
			queryTaxaDir.insert(orgTaxaId,posInBuf);
		}
	}
	genTaxaOrgIdFile.close();
	queryTaxaDirFile.close();	
}

void AllMaps::initSeqHashTable()
{
	ifstream genSeqOrgIdArrFile (genSeqOrgIdArrFileName,ios::in | ios::binary);
	unsigned int numSeq,maxSeqId=0,seqId;
	
	//open genSeqOrgIdArrFile for o/p
	if(genSeqOrgIdArrFile.is_open())
	{
		genSeqOrgIdArrFile.read(reinterpret_cast<char *>(&(numSeq)),sizeof(unsigned int)); //read # of taxaIds
		T_LOG<<"NumOfSeqId: "<<numSeq<<"\n";
		dupQueryTaxaFlag = new bool[numSeq];
		//below loop runs separately, just to make it look neat
		for(unsigned int i =0; i<numSeq; ++i){
			dupQueryTaxaFlag[i] = false;
		}
		for(unsigned int i =0; i<numSeq; ++i)
		{
			genSeqOrgIdArrFile.read(reinterpret_cast<char *>(&seqId),sizeof(unsigned int)); //read org seq 
			if(maxSeqId < seqId)
			{
				maxSeqId =seqId;
			}
			seqHashTable.insert(seqId,i);
		}
	}
	T_LOG<<"Max SeqId: "<<maxSeqId<<"\n";
	genSeqOrgIdArrFile.close();//genSeqOrgIdArrFile written
}

void AllMaps::initStrTaxaHashTb()
{
	ifstream strTaxaDBFile (strTaxaDBFileName,ios::in | ios::binary);
	if(!strTaxaDBFile.is_open()){
		T_LOG<<"Error! Could not open strTaxaDBFileName."<<endl;
			return;
	}

	ifstream strGenusDBFile (strGenusDBFileName,ios::in | ios::binary);
	if(!strGenusDBFile.is_open()){
		T_LOG<<"Error! Could not open strGenusDBFileName."<<endl;
			return;
	}
	size_t charBufSzTaxa,charBufSzGenus;

	strTaxaDBFile.read(reinterpret_cast<char *>(&charBufSzTaxa),sizeof(size_t));
	
	strGenusDBFile.read(reinterpret_cast<char *>(&charBufSzGenus),sizeof(size_t));


	strTaxaHashTb = new String_Hash(largePrimeForTaxaPlusGenus,(unsigned int)(charBufSzTaxa+ charBufSzGenus));	
	unsigned int mapBufSzTaxa,mapBufSzGenus,strTaxaIdC,strGenusIdC,mapBufInd=0;

	
	strTaxaDBFile.read(reinterpret_cast<char *>(&mapBufSzTaxa),sizeof(unsigned int));
	strGenusDBFile.read(reinterpret_cast<char *>(&mapBufSzGenus),sizeof(unsigned int));
	strTaxaGenTaxaMapBuf = new unsigned int[mapBufSzTaxa+mapBufSzGenus];

	strTaxaDBFile.read(reinterpret_cast<char *>(&strTaxaIdC),sizeof(unsigned int));
	strGenusDBFile.read(reinterpret_cast<char *>(&strGenusIdC),sizeof(unsigned int));
	strTaxaOrgQueryTaxaMapIndArr = new unsigned int[strTaxaIdC+strGenusIdC];
	strTaxaIdFreqArr = new unsigned char[strTaxaIdC+strGenusIdC];
	unsigned char strLen, orgTaxaIdC;
	unsigned int genTaxaId;
	char buf[1024];//to read strTaxaId
	
	//first read taxaIds
	for(unsigned int i=0; i<strTaxaIdC; ++i){
		strTaxaDBFile.read(reinterpret_cast<char *>(&strLen),sizeof(unsigned char));
		strTaxaDBFile.read(buf,sizeof(char)*strLen);
		buf[strLen] ='\0';
		//if(0 == strcmp(buf,"Aeschynomene villosa")){
		//	cout<<"";
		//}
		strTaxaHashTb->insert(buf,i);
		strTaxaDBFile.read(reinterpret_cast<char *>(&orgTaxaIdC),sizeof(unsigned char));
		strTaxaOrgQueryTaxaMapIndArr[i]= mapBufInd;
		strTaxaIdFreqArr[i] = orgTaxaIdC;
		for(unsigned char j=0; j<orgTaxaIdC; ++j){
			strTaxaDBFile.read(reinterpret_cast<char *>(&strTaxaGenTaxaMapBuf[mapBufInd+j]),sizeof(unsigned int));
			strTaxaDBFile.read(reinterpret_cast<char *>(&genTaxaId),sizeof(unsigned int));//ignore
			//if(561489 == strTaxaGenTaxaMapBuf[mapBufInd+j]){
			//	cout<<"";
			//}
		}		
		mapBufInd+= orgTaxaIdC;
	}

	//first read genusIds
	for(unsigned int i=strTaxaIdC; i<strTaxaIdC+strGenusIdC; ++i){
		strGenusDBFile.read(reinterpret_cast<char *>(&strLen),sizeof(unsigned char));
		strGenusDBFile.read(buf,sizeof(char)*strLen);
		buf[strLen] ='\0';
		//if(0 == strcmp(buf,"Aeschynomene villosa")){
		//	cout<<"";
		//}
		strTaxaHashTb->insert(buf,i);
		strGenusDBFile.read(reinterpret_cast<char *>(&orgTaxaIdC),sizeof(unsigned char));
		strTaxaOrgQueryTaxaMapIndArr[i]= mapBufInd;
		strTaxaIdFreqArr[i] = orgTaxaIdC;
		for(unsigned char j=0; j<orgTaxaIdC; ++j){
			strGenusDBFile.read(reinterpret_cast<char *>(&strTaxaGenTaxaMapBuf[mapBufInd+j]),sizeof(unsigned int));
			strGenusDBFile.read(reinterpret_cast<char *>(&genTaxaId),sizeof(unsigned int));//ignore
			//if(561489 == strTaxaGenTaxaMapBuf[mapBufInd+j]){
			//	cout<<"";
			//}
		}		
		mapBufInd+= orgTaxaIdC;
	}

	strGenusDBFile.close();
	strTaxaDBFile.close();
}

bool CmdLineOpts::parse(int argc, string * argv, AllMaps &allMaps)
{
	reset();//reset all command line options to their default value
	bool retVal = true;

	if(argc ==1)
	{
		RESPONSE<<"No arguments provided.\n";
		return false;
	}
	istringstream stm_i;//needed to chek if user input being numeric
	ostringstream stm_o;
	for(int i=1; i<argc; ++i)//first arg is name of the prog, hence skip it
	{
		if(taxonFileOpt.compare(argv[i]) == 0)
		{//"-f"
			taxonFileName = argv[++i];
			ifstream taxonFile (taxonFileName.c_str());
			if(!taxonFile.is_open())
			{
				RESPONSE<<"Could not open "<<taxonFileName<<".\n";
				return false;
			}
			string line;
			numTaxaId=0;//becuase you do not want both file option and list option to be used simultaneously
			while (! taxonFile.eof() )
			{				
				getline (taxonFile,line);
				if(line.length() == 0)
				{
					continue;
				}
				stm_i.str(line);
				stm_i>>taxaIds[numTaxaId];
				stm_i.clear();
				stm_o.str("");
				stm_o<<taxaIds[numTaxaId++];;
				if(0 != stm_o.str().compare(line))
				{
					RESPONSE<<"Invalid taxa-id:"<<line<<" in file "<<taxonFileName<<".\n";
					taxonFile.close();
					return false;
				}

			}
			taxonFile.close();
			if(0 == numTaxaId)
			{
				RESPONSE<<"No valid taxa-ids found.\n";
				return false;
			}
			
		}
		else if(taxaIdsOpt.compare(argv[i]) == 0)
		{//"-L"
			numTaxaId=0;//becuase you do not want both file option and list option to be used simultaneously
			while(i<argc-1)
			{
				stm_i.str(argv[++i]);
				stm_i>>taxaIds[numTaxaId];
				stm_i.clear();
				stm_o.str("");
				stm_o<<taxaIds[numTaxaId++];;
				if(0 != stm_o.str().compare(argv[i]))
				{
					if(argv[i][0] != '-')
					{//not an option but invalid id
						RESPONSE<<"Invalid numeric value for option:"<<taxaIdsOpt<<" "<<argv[i]<<"\n";
						return false;
					}
					--numTaxaId;
					--i;
					//next option, reduce i to point to it
					break;
				}
			}
			if(0 == numTaxaId)
			{
				RESPONSE<<"No valid taxa-ids found!\n";
				return false;
			}
		}
		else if(strTaxaFileOpt.compare(argv[i]) == 0)
		{//"-n"
			taxonFileName = argv[++i];
			ifstream taxonFile (taxonFileName.c_str());
			if(!taxonFile.is_open())
			{
				RESPONSE<<"Could not open "<<taxonFileName<<".\n";
				return false;
			}
			string line;
			numTaxaId=0;//becuase you do not want both file option and list option to be used simultaneously
			while (! taxonFile.eof() )
			{				
				getline (taxonFile,line);
				if(line.length() == 0)
				{
					continue;
				}
				size_t pos=0;
				//need to replace all '`' with ' '
				while((pos = line.find('`',pos))!= string::npos)
				{
					line.replace(pos,1,1,' ');
					++pos;
				}

				//if(0 == line.compare("Aeschynomene villosa")){
				//	cout<<"";
				//}

				unsigned int strTaxaId = allMaps.strTaxaHashTb->get(line.c_str());
				if(UINT_MAX != strTaxaId){
					unsigned int mapInd = allMaps.strTaxaOrgQueryTaxaMapIndArr[strTaxaId];
					unsigned char count = allMaps.strTaxaIdFreqArr[strTaxaId];
					for(unsigned char j=0; j<count; ++j){
						taxaIds[numTaxaId++] = allMaps.strTaxaGenTaxaMapBuf[mapInd+j];						
					}					
				}else{
					continue;//ignore this strTaxaId as it did not match any in strTaxaHashTb
				}
			}
			taxonFile.close();
			if(0 == numTaxaId)
			{
				RESPONSE<<"No valid query taxa names found in file: "<<taxonFileName<<".\n";
				return false;
			}
		}
		else if(strTaxaIdOpt.compare(argv[i]) == 0)
		{//"-N"
			numTaxaId=0;//becuase you do not want both file option and list option or name option to be used simultaneously
			while(i<argc-1)
			{
				if(argv[++i][0] == '-'){//next option reached
					--i;
					break;
				}

				size_t pos=0;
				//need to replace all '`' with ' '
				while((pos = argv[i].find('`',pos))!= string::npos)
				{
					argv[i].replace(pos,1,1,' ');
					++pos;
				}					

				unsigned int strTaxaId = allMaps.strTaxaHashTb->get(argv[i].c_str());
				if(UINT_MAX != strTaxaId){
					unsigned int mapInd = allMaps.strTaxaOrgQueryTaxaMapIndArr[strTaxaId];
					unsigned char count = allMaps.strTaxaIdFreqArr[strTaxaId];
					for(unsigned char j=0; j<count; ++j){
						taxaIds[numTaxaId++] = allMaps.strTaxaGenTaxaMapBuf[mapInd+j];						
					}					
				}else{
					continue;//ignore this strTaxaId as it did not match any in strTaxaHashTb
				}				
			}
			if(0 == numTaxaId)
			{
				RESPONSE<<"No valid query taxa names found!\n";
				return false;
			}
		}
		else if(overLapOpt.compare(argv[i]) == 0)
		{//"-k"
			stm_i.str(argv[++i]);
			stm_i>>overLap;
			stm_i.clear();
			stm_o.str("");
			stm_o<<overLap;
			--overLap;//we use > instead of >=, hence overlap value should be one less
			if(0 != stm_o.str().compare(argv[i])|| overLap<0)
			{
				RESPONSE<<"Invalid numeric value for option:"<<overLapOpt<<" "<<argv[i]<<"\n";
				return false;//invalid option
			}
			if(overLap<3){//we use > instead of >=, hence overlap value should be one less - actual min overlap value is 4
				overLap =3;
			}
		}
		else if(useSeqIdOpt.compare(argv[i]) == 0)
		{//"-q"
			++i;
			if(argv[i][0] == 'S')
			{
				useSeqId = true;
			}
			else if(argv[i][0] == 'T')
			{
				useSeqId = false;
			}
			else
			{
				RESPONSE<<"Invalid option for "<<useSeqIdOpt<<" "<<argv[i]<<";Only (S|T) allowed!\n";
				return false;//invalid option
			}

		}
		else if(treeFileOpt.compare(argv[i]) == 0)
		{//"-t"
			treeFileName= argv[++i];
			ofstream checkFile (treeFileName.c_str());
			if(!checkFile.is_open())
			{
				RESPONSE<<"Could not open "<<treeFileName<<".\n";
				return false;
			}
			checkFile.close();
		}
		else if(labelFormatOpt.compare(argv[i]) == 0)
		{//"-z"
			++i;
			if((argv[i][0] == 'S' && argv[i][1] == 'T') || (argv[i][0] == 'T' && argv[i][1] == 'S'))
			{
				labelFormat = 0;
			}
			else if(argv[i][0] == 'S')
			{
				labelFormat = 'S';
			}
			else if(argv[i][0] == 'T')
			{
				labelFormat = 'T';
			}
			else if(argv[i][0] == 'N')
			{
				labelFormat = 'N';
			}
			else if(argv[i][0] == 'A')
			{
				labelFormat = 'A';
			}
			else
			{
				RESPONSE<<"Invalid option for "<<labelFormatOpt<<" "<<argv[i]<<";only (S|T|ST|TS|N|A) allowed.\n";
				return false;//invalid option
			}

		}
		else if(cumQueryTaxaOpt.compare(argv[i]) == 0)
		{//"-u"
			cumQueryTaxa = true;
		}
		else if(tableFileOpt.compare(argv[i]) == 0)
		{//"-s"			
			tableFileName= argv[++i];
			ofstream checkFile (tableFileName.c_str());
			if(!checkFile.is_open())
			{
				RESPONSE<<"Could not open "<<checkFile<<".\n";
				return false;
			}
			checkFile.close();
		}else if(seqIdOverlapOpt.compare(argv[i]) == 0)
		{//"-o"
			 ++i;
			if(argv[i][0] == 'S')
			{
				seqIdOverlap = true;
			}
			else if(argv[i][0] == 'T')
			{
				seqIdOverlap = false;
			}
			else
			{
				RESPONSE<<"Invalid option for "<<seqIdOverlapOpt<<" "<<argv[i]<<";Only (S|T) allowed!\n";
				return false;//invalid option
			}
		}
		else if(maxCCOpt.compare(argv[i]) == 0)
		{//"-c"
			stm_i.str(argv[++i]);
			stm_i>>maxCC;
			stm_i.clear();
			stm_o.str("");
			stm_o<<maxCC;
			if(0 != stm_o.str().compare(argv[i])|| maxCC<0)
			{
				RESPONSE<<"Invalid numeric value for option:"<<maxCCOpt<<" "<<argv[i]<<"\n";
				return false;//invalid option
			}
		}else if(maxBSTOpt.compare(argv[i]) == 0)
		{//"-b"
			stm_i.str(argv[++i]);
			stm_i>>maxBST_Count;
			stm_i.clear();
			stm_o.str("");
			stm_o<<maxBST_Count;
			if(0 != stm_o.str().compare(argv[i])|| maxBST_Count<0)
			{
				RESPONSE<<"Invalid numeric value for option:"<<maxBST_Count<<" "<<argv[i]<<"\n";
				return false;//invalid option
			}
		}
		else
		{
			RESPONSE<<"Invalid option: "<<argv[i]<<"\n";
			return false;
		}
	}

	return retVal;
}
