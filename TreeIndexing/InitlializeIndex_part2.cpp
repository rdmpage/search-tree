#include "TreeIndexingUtilities.h"

void initialize_part2()
{
	VocabDB vocabDB;
	InvListDB invListDB;

	createVocabInitInvListMap(vocabDB,invListDB);//creates Vocab file and inits invList for freqs of each taxa

	createInvList(vocabDB,invListDB);

	delete[] vocabDB.vocabBuf;
	delete[] invListDB.invListBuf;
}

//creates Vocab file and inits invList for freqs of each taxa
void createVocabInitInvListMap(VocabDB &vocabDB,InvListDB &invListDB)
{
	unsigned int *querySeqFreqArrG = new unsigned int[maxSeqCount];
	invListDB.size=0;

	ifstream querySeqFreqArrG_File (querySeqFreqArrG_FileName,ios::in | ios::binary);
	//open global taxa freq file to read
	if(querySeqFreqArrG_File.is_open())
	{
		unsigned int numTaxa;
		querySeqFreqArrG_File.read(reinterpret_cast<char *>(&numTaxa),sizeof(unsigned int));
		
		for(unsigned int i=0; i<numTaxa; ++i)
		{
			querySeqFreqArrG_File.read(reinterpret_cast<char *>(&querySeqFreqArrG[i]),sizeof(unsigned int));//read frequency
			if(0 != querySeqFreqArrG[i])
			{//only for a taxa which exists in some tree
				invListDB.size+= 1+ querySeqFreqArrG[i];//space to store the tree count and -
				//-space to store each of the trees
			}
		}

		T_LOG<<std::endl<<"Size of invListBuf:"<<invListDB.size;
		invListDB.invListBuf = new unsigned int[invListDB.size];//assign memory to invListbuf for total taxacount
		vocabDB.taxaCount = numTaxa;
		vocabDB.vocabBuf = new unsigned int[vocabDB.taxaCount];//assign memory to vocabBuf
		unsigned int currPos =0; //curr ind in invList buf
		for(unsigned int i=0; i<numTaxa; ++i)
		{
			if(0 != querySeqFreqArrG[i])//if the taxa does exist in any tree
			{
				vocabDB.vocabBuf[i] =  currPos;
				invListDB.invListBuf[currPos] =0; //intialize freq to 0
				currPos+= querySeqFreqArrG[i]+1;
				//currPos incremented for tree freq count and id for each tree
			}
			else
			{//if the taxa does not exist in any tree
				vocabDB.vocabBuf[i] = UINT_MAX;//else initialize its pointer in vocab as symbolic null
			}
		}	
		//vocab file construction complete in vocabBuf, ready to be written 

	}
	querySeqFreqArrG_File.close();//global taxa freq file read

	ofstream vocabFile (vocabFileName,ios::out | ios::binary);
	//open vocab file to write vocabBuf into
	if(vocabFile.is_open())
	{
		vocabFile.write(reinterpret_cast<char *>(&vocabDB.taxaCount),sizeof(unsigned int)); //write taxaCount
		for(unsigned int tc=0; tc<vocabDB.taxaCount; ++tc)
		{		
			vocabFile.write(reinterpret_cast<char *>(&vocabDB.vocabBuf[tc]),sizeof(unsigned int));
		}
	}
	vocabFile.close();//vocab file written

	delete [] querySeqFreqArrG;
}

void createInvList(VocabDB &vocabDB,InvListDB &invListDB)
{
	for(unsigned short tc=0;tc<treeDirCount;++tc)
	{
		unsigned short numOfTaxas = 0; //to read in # of taxas 
		unsigned int taxaId = 0; //to read in taxaId
		unsigned int uniqTreeId = 0;//to read in uniqTreeId
		unsigned int taxaInd=0; //stores index to the taxaId in InvListBuf
		
		char treeFileName[257] ;
		sprintf(treeFileName,"%s/Tree%d.db",treeDirBN,tc);
		ifstream treeFile(treeFileName, ios::in | ios::binary);
		//open tree file to read
		if(treeFile.is_open())
		{
			while(true)
			{
				treeFile.read(reinterpret_cast<char *>(&uniqTreeId),sizeof(unsigned int));
				if(treeFile.eof())
				{
					break;//eof file reached
				}
				treeFile.read(reinterpret_cast<char *>(&numOfTaxas),sizeof(unsigned short));
				for(; numOfTaxas >0; --numOfTaxas)
				{
					treeFile.read(reinterpret_cast<char *>(&taxaId),sizeof(unsigned int));
					taxaInd=vocabDB.vocabBuf[taxaId];//index in invListBuf for the taxaId
					taxaInd+=++invListDB.invListBuf[taxaInd];//pre increment intentional;increase the frequency and calculate the index position 
					//for the new treeId to be written into
					invListDB.invListBuf[taxaInd] = uniqTreeId; //write curr tree uniqId into invListBuff
				}
			}
			treeFile.close();
		}
	}

	//write the invList buf into files
	//two files are written each taking half of invList buf

	{//just to localize the scope
		char invListFile1Name[257] ;
		sprintf(invListFile1Name,"%s/invList%d.db",invListBaseDir,1);
		ofstream invListFile1(invListFile1Name,ios::out | ios::binary);
		//open 1st invlist file to write
		if(invListFile1.is_open())
		{
			invListFile1.write(reinterpret_cast<char *>(&invListDB.size),sizeof(unsigned int));
			//write the invListBuf size in the first file
			unsigned int midInd=invListDB.size/2;//index dividing invListBuf into two parts
			for(unsigned int i=0; i<midInd; ++i)
			{
				invListFile1.write(reinterpret_cast<char *>(&invListDB.invListBuf[i]),sizeof(unsigned int));
			}
		}
		invListFile1.close();//1st invList file written
	}
	
	{//just to localize the scope
		char invListFile2Name[257] ;
		sprintf(invListFile2Name,"%s/invList%d.db",invListBaseDir,2);
		ofstream invListFile2(invListFile2Name,ios::out | ios::binary);
		//open 2nd invlist file to write
		if(invListFile2.is_open())
		{
			unsigned int midInd=invListDB.size/2;//index dividing invListBuf into two parts
			for(unsigned int i=midInd; i<invListDB.size; ++i)
			{
				invListFile2.write(reinterpret_cast<char *>(&invListDB.invListBuf[i]),sizeof(unsigned int));
			}
		}
		invListFile2.close();//2nd invList file written
	}
}

