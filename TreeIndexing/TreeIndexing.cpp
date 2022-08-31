// TreeIndexing.cpp : Defines the entry point for the console application.

#include "TreeIndexingUtilities.h"
#include"Test_TI.h"
#ifndef WIN32
#include <signal.h>
#endif 


void initialize();


int main(int argc, char * argv[])
{

	//testConvertSeqIdToTaxonNames();
	//return 0;
#ifndef WIN32
signal( SIGPIPE, SIG_IGN );
#endif 
	if(!TI_LOG_FILE.is_open())
	{
		std::cout<<"FATAL ERROR!TI_LOG_FILE could not be openened.\n";
		exit(-1);
	}
	//initialize();
	//return 0;
	ProcessInput processI;
	loadMaps(processI.allMaps);
	
	while(true)
	{
		processI.getQuery(argc);
		if(argc == 0)
		{//query could not be read successfully
			continue;
		}
		int respVal = processI.process(argc);
		RESPONSE.close();
		processI.signalResultReadyToClient(respVal);
		if(respVal != processI.respCount){
			T_LOG<<"Processing of query "<<processI.respCount<<" failed."<<endl;
		}
		++processI.respCount;//if we reached this point, the query was processed, success or not depends upon respVal
		if(INT_MAX == processI.respCount){
			processI.respCount = 1;//round off
		}

#ifdef WIN32
	break;
#endif 
	}

	//testRead(mapVocab,mapQueryTaxaDir, mapTreeDir);
	//finalizeFileHandlers();

	TI_LOG_FILE.close();
	return 0;
}

void initialize()
{
	initialize_part1();
	initialize_part2();
	
}

int _main(int argc, char * argv[])
{
	
	hashSeqId();
	return 0;
}

void hashSeqId()
{
	ifstream genSeqOrgIdArrFile (genSeqOrgIdArrFileName,ios::in | ios::binary);
	unsigned int numTaxaId,maxTaxaId=0;;
	unsigned int *taxaUniqIdArr = NULL;
	UI_HashTable hashTable(largePrimeForSeq);
	set<unsigned int> queryTaxaSet;
	//open uniq id to taxa id array file for o/p
	if(genSeqOrgIdArrFile.is_open())
	{
		genSeqOrgIdArrFile.read(reinterpret_cast<char *>(&(numTaxaId)),sizeof(unsigned int)); //read # of taxaIds
		T_LOG<<"NumOfTaxaId: "<<numTaxaId<<"\n";
		taxaUniqIdArr = new unsigned int [numTaxaId];
		for(unsigned int tc =0; tc<numTaxaId; ++tc)
		{
			genSeqOrgIdArrFile.read(reinterpret_cast<char *>(&(taxaUniqIdArr[tc])),sizeof(unsigned int)); //write taxaId
			if(maxTaxaId < taxaUniqIdArr[tc])
			{
				maxTaxaId=taxaUniqIdArr[tc];
			}
			hashTable.insert(taxaUniqIdArr[tc],tc);
			queryTaxaSet.insert(taxaUniqIdArr[tc]);

		}
		for(unsigned int tc =0; tc<numTaxaId; ++tc)
		{
			if(tc != hashTable.get(taxaUniqIdArr[tc]))
			{
				T_LOG<<"Value Mismatch!\n";
			}
		}
	}
	T_LOG<<"Max TaxaId: "<<maxTaxaId<<"\n";
	T_LOG<<queryTaxaSet.size()<<"\n";
	genSeqOrgIdArrFile.close();//taxa uniqId array file written
	if(NULL !=taxaUniqIdArr)
	{
		delete taxaUniqIdArr;
	}
}

