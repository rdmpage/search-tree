#include "TreeIndexingUtilities.h"
#include<set>


int ProcessInput::process(int argc)
{

	bool quitFlag = false;
	int retVal = -1;
	CmdLineOpts &cmdLineOpts= allMaps.cmdLineOpts;//for convenience


	if(!cmdLineOpts.parse(argc,argvBuf,allMaps))
	{
		return -1;
	}

	//buf to store all the query id given as input
	unsigned int *taxaIdBuf = cmdLineOpts.taxaIds;//no need to create another taxaId buf

	unsigned int queryIdCount =0;//# of queryIds
	unsigned int taxaIdCount =0;//# of sequenceIds
	int id;//a query id
	istringstream stm_i;//needed to chek if user input being numeric
	ostringstream stm_o;
	string strId;

	bool readFromConsole = true;

	if(cmdLineOpts.numTaxaId !=0)
	{
		readFromConsole = false;
	}

	for(unsigned int resCount =1 ; !quitFlag ;resCount++)
	{

		queryIdCount=0;//must
		taxaIdCount=0;//must

		if(readFromConsole)
		{
			T_LOG<<"Enter Ids (Max count 1000), with each followed by enter and finally -1 to see the result, -2 to exit"<<std::endl;
		}
		while(true)
		{
			if(readFromConsole)
			{
				//cin>> id;
				cin>>strId;
				stm_i.str(strId);
				stm_i>>id;
				stm_i.clear();
				stm_o.str("");
				stm_o<<id;
				if(0 != strId.compare(stm_o.str()) || id<-2)
				{
					T_LOG<<strId<<":Invalid id, ignoring!"<<std::endl;
					continue;
				}


			}
			else
			{
				id =-1;
			}



			if(id == -2)
			{
				quitFlag = true;
				break;
			}else if (id != -1)
			{
				if(queryIdCount == thousand-1 || taxaIdCount==thousand-1)
				{
					T_LOG<<"Max count reached for Ids, cannot accept more, enter -1 to see the result!"<<std::endl;
				}
				else
				{
					if(!cmdLineOpts.useSeqId)
					{//i.e. input is taxa id
						queryIdBuf[queryIdCount++][0]=id;
					}
					else
					{//input is seq id
						taxaIdBuf[cmdLineOpts.numTaxaId++] = id;
					}
				}

			}else
			{

				if(!readFromConsole && !cmdLineOpts.useSeqId)
				{//not from console and taxa ids for query
					if(cmdLineOpts.numTaxaId>tenThousand)
					{
						RESPONSE<<"Taxon ids more than 10,000. Only first 10,000 being considered."<<std::endl;
						cmdLineOpts.numTaxaId = tenThousand;
					}
					queryIdCount=0;
					for(unsigned int tc=0; tc<cmdLineOpts.numTaxaId; ++tc)
					{
						queryIdBuf[queryIdCount++][0]=cmdLineOpts.taxaIds[tc];
					}


				}

				T_LOG<<"Processing, please wait for the result..."<<std::endl;
#ifdef WIN32
#else
				struct timeval stTime,endTime;
				gettimeofday(&stTime, NULL);

#endif


				char resFileName[257] ;
				sprintf(resFileName,"./../common/otherResults/Result%d",respCount);
				ofstream resFile(resFileName);
				//open the resultFile to write
				if(!resFile.is_open())
				{
					T_LOG<<std::endl<<"ERROR! RESULT FILE COULD NOT BE CREATED!";
					exit(-1);
				}

				unsigned int indTreeUniqIdResBuf =0;//index in treeUniqIdResBuf

				unsigned int simpleTC =0;//counts all the trees including the duplicated ones

				//read all the seq ids corresponding to all the queryTaxaIds
				unsigned int seqIdCount =0; //to store count of seqIds a queryTaxaId is mapped to

				if(!cmdLineOpts.useSeqId)
				{//i.e. taxa ids are being used for query; readFromConsole option does not matter
					unsigned int posInBuf;
					for(unsigned int qc=0; qc<queryIdCount; ++qc)
					{
						posInBuf = allMaps.queryTaxaDir.get(queryIdBuf[qc][0]);
						if(posInBuf != UINT_MAX)
						{//read seqIds for the corresponding queryId
							seqIdCount = allMaps.queryTaxMapBuf[posInBuf+1];//get the first seq id, seqIdCount-
							//- here is a seqId
							if(!allMaps.dupQueryTaxaFlag[seqIdCount]){//check for duplicacy
								//first time
								allMaps.dupQueryTaxaFlag[seqIdCount]= true;//mark it as encountered
								queryIdBuf[qc][1] = taxaIdCount;//store the index in taxaIdBuf
								seqIdCount = allMaps.queryTaxMapBuf[posInBuf];
								for(;seqIdCount>0;--seqIdCount)
								{
									taxaIdBuf[taxaIdCount++]=allMaps.queryTaxMapBuf[posInBuf+seqIdCount];
								}
							}else{
								//a dup taxa causing a dup seq
								queryIdBuf[qc][0]=queryIdBuf[queryIdCount-1][0];//swap it with the last value
								--queryIdCount;
								--qc;//decrement so that the swapped last value gets visisted
							}
						}
						else
						{
							queryIdBuf[qc][0] = UINT_MAX;//mark that this query id does not have any taxa id
							queryIdBuf[qc][1] = taxaIdCount;//store the index in taxaIdBuf - demarkation should stay correct
						}
					}
					queryIdBuf[queryIdCount][1] = taxaIdCount;//a demarkation needed in next for loop for counting tree frequencies
					//reset dup flag array for next use
					for(unsigned int tc=0; tc<taxaIdCount; ++tc){
						allMaps.dupQueryTaxaFlag[taxaIdBuf[tc]] = false;
					}

				}
				else
				{//sequence ids in use
					taxaIdCount = cmdLineOpts.numTaxaId;
					queryIdCount =0;
					unsigned int taxaId;
					//sequence ids are present, but in original form hence need to be converted to get the sequential form
					for(unsigned int tc=0; tc<taxaIdCount; ++tc)
					{
						taxaId = allMaps.seqHashTable.get(taxaIdBuf[tc]);
						if(UINT_MAX != taxaId)
						{
							if(!allMaps.dupQueryTaxaFlag[taxaId]){//check for duplicacy
								//first time
								allMaps.dupQueryTaxaFlag[taxaId]= true;//mark it as encountered
								taxaIdBuf[tc]=taxaId;
							}else{
								taxaIdBuf[tc]= taxaIdBuf[taxaIdCount-1];//swap it with the last value
								--taxaIdCount;
								--tc;//decrement so that the swapped last value gets visisted
							}
						}
						else
						{//invalid seq id
							taxaIdBuf[tc] = taxaIdBuf[taxaIdCount-1];//swap it with last value
							--tc;
							--taxaIdCount;
						}						
					}
					//reset dup flag array for next use
					for(unsigned int tc=0; tc<taxaIdCount; ++tc)
					{
						allMaps.dupQueryTaxaFlag[taxaIdBuf[tc]]= false;//mark it as encountered
					}
				}

				//now read all the trees -this includes duplicated trees
				unsigned int taxaId =0,queryId;//to store curr taxaId(seq id) and queryId (taxa id)
				unsigned int invListInd =0;//index in invListBuf for a taxaId
				unsigned int taxaFreq = 0;//# of treeIds a taxa is mapped to
				unsigned int treeId = 0;//to store curr treeId
				if(!cmdLineOpts.useSeqId)
				{//taxa ids in use
					for(unsigned int qc=0; qc<queryIdCount; ++qc)
					{
						if(UINT_MAX != queryIdBuf[qc][0])
						{//i.e. this taxa id has non zero seq ids
							taxaId = taxaIdBuf[queryIdBuf[qc][1]];//first seq id for this query taxa id
							queryId = allMaps.seqIdUniqQueryIdArr[taxaId];//corresponding uniquer query taxa id
							for(unsigned int tc=queryIdBuf[qc][1]; tc<queryIdBuf[qc+1][1]; ++tc)
							{
								taxaId=taxaIdBuf[tc];
								invListInd = allMaps.vocabDB.vocabBuf[taxaId];
								if(UINT_MAX != invListInd) //i.e. taxa id is present in some tree
								{
									taxaFreq = allMaps.invListDB.invListBuf[invListInd];//# of trees the taxa Id is mapped to
									for(;taxaFreq>0;--taxaFreq)
									{
										treeId = 
											treeIdResBuf[simpleTC++] = allMaps.invListDB.invListBuf[invListInd+taxaFreq];
										if(treeResFreqTaxaIdBuf[treeId][1] != queryId+1)//+1 is required as 0 means no query id hence will work for queryId =0
										{//first time this tree is figuring for this taxa
											++treeResFreqTaxaIdBuf[treeId][0];
											treeResFreqTaxaIdBuf[treeId][1] = queryId+1;
										}
										++treeResFreqSeqIdBuf[treeId];//increase freq for resTreeId
										if(1 == treeResFreqSeqIdBuf[treeId])
										{//first time, uniq treeIdFound, insert in uniqTreeIdArr
											treeUniqIdResBuf[indTreeUniqIdResBuf++]=treeId;
										}
									}
								}
							}
						}
					}
				}
				else
				{//sequence ids in use
					for(unsigned int tc=0;tc<taxaIdCount; ++tc)
					{
						taxaId=taxaIdBuf[tc];
						invListInd = allMaps.vocabDB.vocabBuf[taxaId];
						if(UINT_MAX != invListInd) //i.e. taxa id is present in some tree
						{
							taxaFreq = allMaps.invListDB.invListBuf[invListInd];//# of trees the taxa Id is mapped to
							for(;taxaFreq>0;--taxaFreq)
							{
								treeId = 
									treeIdResBuf[simpleTC++] = allMaps.invListDB.invListBuf[invListInd+taxaFreq];
								++treeResFreqSeqIdBuf[treeId];//increase freq for resTreeId
								if(1 == treeResFreqSeqIdBuf[treeId])
								{//first time, uniq treeIdFound, insert in uniqTreeIdArr
									treeUniqIdResBuf[indTreeUniqIdResBuf++]=treeId;
								}
							}
						}
					}
				}

				//required for part 2 start
				unsigned int validTreeCount=0;//counts trees with overlap of >= 3
				//required for part 2 end

				T_LOG<<"Processing done. "<<std::endl<<"Flushing the result in the result file...."<<std::endl;
				resFile<<"# of query taxa:"<<queryIdCount<<"\t# of sequence ids:"<<taxaIdCount<<"\n";
				resFile<<"# of unique trees/# of trees: in total ="<<indTreeUniqIdResBuf<<"/"<<simpleTC<<"\n";


				if(!cmdLineOpts.useSeqId && !cmdLineOpts.seqIdOverlap)
				{//taxa id in use and overlap with taxon id only
					resFile<<"TreeId \t"<<"#of Taxa\t"<<"#of Query Taxa \t"<<"#of Query SeqId\n";
					for(unsigned int tc=0; tc<indTreeUniqIdResBuf; ++tc)
					{
						treeId = treeUniqIdResBuf[tc];
						resFile<<treeId<<"\t"<<allMaps.arrAllTree[treeId].nodeCount
							<<"\t\t"<<treeResFreqTaxaIdBuf[treeId][0]<<"\t"
							<<treeResFreqSeqIdBuf[treeId]<<"\n";
						//required for part 2 start
						//if(treeResFreqSeqIdBuf[treeId]>2 && treeId< hundredThousand)
						if(treeResFreqTaxaIdBuf[treeId][0]>cmdLineOpts.overLap)
						{
							++validTreeCount;
						}
						//required for part 2 end
					}
				}
				else
				{//seq id in use or overlap with seq id
					resFile<<"TreeId \t"<<"#of Taxa\t"<<"#of Query Taxa \n";
					for(unsigned int tc=0; tc<indTreeUniqIdResBuf; ++tc)
					{
						treeId = treeUniqIdResBuf[tc];
						resFile<<treeId<<"\t"<<allMaps.arrAllTree[treeId].nodeCount
							<<"\t\t"<<treeResFreqSeqIdBuf[treeId]<<"\n";
						//required for part 2 start
						//if(treeResFreqSeqIdBuf[treeId]>2 && treeId< hundredThousand)
						if(treeResFreqSeqIdBuf[treeId]>cmdLineOpts.overLap)
						{
							++validTreeCount;
						}
						//required for part 2 end
					}
				}
				resFile.close();
				T_LOG<<"Intermediate Result file(for tracing purpose only):" <<resFileName<<std::endl;

				//write binary result file start
				char resFileNameB[257] ;

				sprintf(resFileNameB,"./../common/otherResults/Result_B%d",respCount);

				ofstream resFileB(resFileNameB,ios::out | ios::binary);
				if(!resFileB.is_open())
				{
					T_LOG<<"\nError!resFileB.is_open() failed";
					exit(-1);					
				}

				//write all the command line options-start
				size_t len;

				resFileB.write(reinterpret_cast<char *>(&cmdLineOpts.useSeqId),sizeof(bool));//if seqId is used for query
				resFileB.write(reinterpret_cast<char *>(&cmdLineOpts.cumQueryTaxa),sizeof(bool));//if seqId is used for query

				len=cmdLineOpts.treeFileName.length();
				resFileB.write(reinterpret_cast<char *>(&len),sizeof(size_t));//treeFileName
				if(len !=0)
				{//write the taxon file name
					resFileB.write(cmdLineOpts.treeFileName.c_str(),(int)len);//if len>0 i.e. treeFileName exists
				}

				len=cmdLineOpts.tableFileName.length();
				resFileB.write(reinterpret_cast<char *>(&len),sizeof(size_t));//tableFileName
				if(len !=0)
				{//write the taxon file name
					resFileB.write(cmdLineOpts.tableFileName.c_str(),(int)len);//if len>0 i.e. tableFileName exists
				}

				resFileB.write(reinterpret_cast<char *>(&cmdLineOpts.labelFormat),sizeof(unsigned char));//labelFormat = S|T|ST|TS|N

				resFileB.write(reinterpret_cast<char *>(&cmdLineOpts.maxBST_Count),sizeof(cmdLineOpts.maxBST_Count));//max boot strap trees to be read

				//write all the command line options-end

				resFileB.write(reinterpret_cast<char *>(&taxaIdCount),sizeof(unsigned int));//write # of sequence Ids
				for(unsigned int tc=0;tc<taxaIdCount; ++tc)
				{
					resFileB.write(reinterpret_cast<char *>(&taxaIdBuf[tc]),sizeof(unsigned int));//write each sequence id
				}

				//testing specific
				if(validTreeCount > cmdLineOpts.maxCC)
				{
					validTreeCount = cmdLineOpts.maxCC;
				}

				resFileB.write(reinterpret_cast<char *>(&validTreeCount),sizeof(unsigned int));//write validTreeCount

				if(!cmdLineOpts.useSeqId && !cmdLineOpts.seqIdOverlap)
				{//query based on taxa id and overlap with taxon id only
					for(unsigned int tc=0; tc<indTreeUniqIdResBuf; ++tc)
					{
						treeId = treeUniqIdResBuf[tc];
						if(treeResFreqTaxaIdBuf[treeId][0]>cmdLineOpts.overLap)
						{
							resFileB.write(reinterpret_cast<char *>(&treeUniqIdResBuf[tc]),sizeof(unsigned int));//write cluster id
							resFileB.write(reinterpret_cast<char *>(&treeResFreqTaxaIdBuf[treeId][0]),sizeof(unsigned int));//# of overlapping taxa ids
							resFileB.write(reinterpret_cast<char *>(&treeResFreqSeqIdBuf[treeId]),sizeof(unsigned int));//# of overlapping seq ids
						}
					}
				}
				else
				{//query based on seq id or overlap with seq id
					for(unsigned int tc=0; tc<indTreeUniqIdResBuf; ++tc)
					{
						treeId = treeUniqIdResBuf[tc];
						if(treeResFreqSeqIdBuf[treeId]>cmdLineOpts.overLap)
						{
							resFileB.write(reinterpret_cast<char *>(&treeUniqIdResBuf[tc]),sizeof(unsigned int));//write cluster id
							resFileB.write(reinterpret_cast<char *>(&treeResFreqTaxaIdBuf[treeId][0]),sizeof(unsigned int));//# of overlapping taxa ids NA in case of seq id query
							resFileB.write(reinterpret_cast<char *>(&treeResFreqSeqIdBuf[treeId]),sizeof(unsigned int));//# of overlapping seq ids
						}
					}
				}


				resFileB.close();
				//write binary result file end


				T_LOG<<"# of query taxa:"<<queryIdCount<<"\t# of sequence ids:"<<taxaIdCount<<"\n";
				T_LOG<<"# of unique trees/# of trees: in total ="<<indTreeUniqIdResBuf<<"/"<<simpleTC<<std::endl;
				T_LOG<<"# valid trees "<<validTreeCount<<std::endl;

				for(unsigned int tc=0;tc<indTreeUniqIdResBuf;++tc)
				{
					treeResFreqSeqIdBuf[treeUniqIdResBuf[tc]]=0;//init freq =0 for next use
					treeResFreqTaxaIdBuf[treeUniqIdResBuf[tc]][0]=0;//init freq =0 for next use
					treeResFreqTaxaIdBuf[treeUniqIdResBuf[tc]][1]=0;//int flag for next use
				}
				
#ifdef WIN32
				RESPONSE.close();
				retVal = respCount;
#else
				T_LOG<<"\nCalling step 2...";
				retVal = sendQueryNumToServer(respCount);
				if(retVal != respCount){
					T_LOG<<"\nServer returned Eorro in Step 2 ! Query num:"<<respCount<<", Server resp:"<<retVal<<endl;
					break;
				}
				T_LOG<<"\nStep 2 Done!Result file: \"FinalResult"<<respCount<<"\"\n";
				gettimeofday(&endTime, NULL);
				if(stTime.tv_usec > endTime.tv_usec)
				{
					endTime.tv_usec = endTime.tv_usec + 1000000;
					endTime.tv_sec = endTime.tv_sec-1;
				}
				unsigned int diffSec = endTime.tv_sec - stTime.tv_sec;
				unsigned int diffMiliSec =  (endTime.tv_usec - stTime.tv_usec)/1000;
				char respFileName[257];
				sprintf(respFileName,"./../common/queryResult/Response1_%d",respCount);				
				//RESPONSE.open(respFileName);//append
				//RESPONSE<<"SUCCESS"<<std::endl;
				float diffInTime = diffSec + (float)diffMiliSec/1000.0;
				//RESPONSE<<"Time taken: "<<diffSec<<"."<<diffMiliSec<<" s\n";
				//T_LOG<<"Time taken: "<<diffSec<<"."<<diffMiliSec<<" s\n";
				RESPONSE<<"Time taken: "<<diffInTime<<" s\n";
				T_LOG<<"Time taken: "<<diffInTime<<" s\n";				
				//RESPONSE.close();//done with processing of the current query; no need to close here, will be closed in the main routine
#endif

				break;
			}
		}
		if(!readFromConsole)
		{
			break;
		}
	}
	//delete [] taxaIdBuf;
	taxaIdBuf = NULL;//it only points to buf allMaps.cmdLineOpts.taxaIds, hence make the ref null

	return retVal;
}


void ProcessInput::getQuery(int &argc)
{
	char queryFileName[257];  
	argc=1;//first option is name of the prog
	sendQueryNumToClient(respCount);
	sprintf(queryFileName,"./../common/queryInput/Query%d",respCount);
	ifstream queryFile(queryFileName);
	sprintf(queryFileName,"./../common/queryResult/Response1_%d",respCount);
	RESPONSE.open(queryFileName);
	sprintf(queryFileName,"./../common/queryResult/Response2_%d",respCount);
	ofstream response2(queryFileName);//so that in case of error Response2_%d file exists
	if(queryFile.is_open() && RESPONSE.is_open() && response2.is_open())
	{
		string line;
		char * tok;
		while(!queryFile.eof())
		{
			getline(queryFile,line);
			if(line.length() == 0)
			{
				continue;
			}
			tok = strtok (const_cast<char *>(line.c_str()), "\t ");//first call needs line
			while (tok != NULL)
			{
				argvBuf[argc++] = tok;
				tok = strtok (NULL, "\t ");//subsequent calls need NULL to use line 
			}
		}
	}
	else
	{
		T_LOG<<"Cannot open Query"<<respCount<<" or Response _1/_2"<<respCount<<endl;
		signalResultReadyToClient(-1);//error
	}
}

void ProcessInput::sendQueryNumToClient(int num)
{
#ifndef WIN32
	sprintf(strBuf,"%d",num);
	istringstream stm;
	string strLine;
	while(true)
	{
		sockServer.connectTo();
		if(!sockServer.sendTo(strBuf))
		{
			continue;
		}
		if(!sockServer.receiveFrom(strLine))
		{
			continue;
		}
		break;
	}
	stm.str(strLine);
	stm>>num;
	T_LOG<<"Num received from client: "<<num<<"\n";
#endif	
}

int ProcessInput::sendQueryNumToServer(int num)
{
#ifndef WIN32
	sprintf(strBuf,"%d",num);
	sockClient.connectTo();
	if(!sockClient.sendTo(strBuf))
	{
		exit(-1);
	}
	istringstream stm;
	string strLine;
	stm.clear();
	if(!sockClient.receiveFrom(strLine))
	{
		exit(-1);
	}
	stm.str(strLine);
	stm>>num;
	T_LOG<<"Num received from server: "<<num<<"\n";
#endif
	return num;
}

void ProcessInput::signalResultReadyToClient(int num)
{
#ifndef WIN32
	sprintf(strBuf,"%d",num);
	sockServer.sendTo(strBuf);
	sockServer.disconnect();
#endif
}
