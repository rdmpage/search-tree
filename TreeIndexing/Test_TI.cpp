#include"Test_TI.h"
#include "TreeIndexingUtilities.h"

void testConvertSeqIdToTaxonNames()
{
	UI_HashTable genSeqOrgId_ht(largePrimeForSeq);//stores mapping from org query taxa to gen query taxa mapping
	ifstream genSeqOrgIdFile(genSeqOrgIdArrFileName,ios::in | ios::binary);
	if(!genSeqOrgIdFile.is_open()){
		T_LOG<<"Error! Could not open mappingFileName."<<endl;
		return;
	}
	unsigned int numSeqIds,orgSeqId;
	genSeqOrgIdFile.read(reinterpret_cast<char *>(&numSeqIds),sizeof(unsigned int));//num of seqIds
	for(unsigned int i=0; i<numSeqIds; ++i){
		genSeqOrgIdFile.read(reinterpret_cast<char *>(&orgSeqId),sizeof(unsigned int));
		genSeqOrgId_ht.insert(orgSeqId,i);			
	}
	genSeqOrgIdFile.close();

	ifstream genSeqGenTaxaArrFile(genSeqGenTaxaArrFileName,ios::in | ios::binary);
	if(!genSeqGenTaxaArrFile.is_open()){
		T_LOG<<"Error! Could not open mappingFileName."<<endl;
		return;
	}
	genSeqGenTaxaArrFile.read(reinterpret_cast<char *>(&numSeqIds),sizeof(unsigned int));//read number of clusters
	unsigned *genSeqGenTaxaArr = new unsigned int[numSeqIds];
	for(unsigned int i=0; i<numSeqIds; ++i)
	{
		genSeqGenTaxaArrFile.read(reinterpret_cast<char *>(&genSeqGenTaxaArr[i]),sizeof(unsigned int));
	}
	genSeqGenTaxaArrFile.close();


	unsigned int *genTaxaIndStrTaxaBufArr;
	unsigned char *strTaxaLenArr;
	char *strTaxaBuf;
	loadStrTaxaInfo(genTaxaIndStrTaxaBufArr,strTaxaLenArr,strTaxaBuf);


	char srCtreeFileName[] ="E:/MiningFrequentTrees/codeFromZhang/sampleTrees/MRTsFromTOL_Data/comparison.txt";
	char destTreeFileName[] ="E:/MiningFrequentTrees/codeFromZhang/sampleTrees/MRTsFromTOL_Data/comparisonWithTaxonNames.txt";
	char buf[1024],destBuf[3072];
	ifstream srcTreeFile(srCtreeFileName);
	if(!srcTreeFile.is_open()){
		exception E;
		cout<<"Error in openeing: "<<srCtreeFileName<<"."<<endl;
		throw E;
	}

	ofstream destTreeFile(destTreeFileName);
	if(!destTreeFile.is_open()){
		exception E;
		cout<<"Error in openeing: "<<destTreeFileName<<"."<<endl;
		throw E;
	}

	char *pos,*pos2,*dest,tmpCh;
	unsigned int genSeqId,genTaxaId;
	while(true)	{
		srcTreeFile.getline(buf,1023);
		if(0 == srcTreeFile.gcount()){
			break;
		}
		pos = strchr(buf,'\t');
		++pos;
		pos = strchr(pos,'\t');
		++pos;
		pos = strchr(pos,'\t');
		*pos = '\0';
		++pos;
		destTreeFile<<buf<<"\t";
				
		dest = destBuf;
		while(*pos){
			if(*pos!=')' && *pos!='(' && *pos!=','){
				pos2 = strpbrk(pos,"(),");
				tmpCh = *pos2;//save the character before overwriting
				*pos2= '\0';
				orgSeqId = atoi(pos);
				*pos2 = tmpCh;
				pos = pos2;
				genSeqId = genSeqOrgId_ht.get(orgSeqId);
				genTaxaId = genSeqGenTaxaArr[genSeqId];
				strcpy(dest,&strTaxaBuf[genTaxaIndStrTaxaBufArr[genTaxaId]]);
				dest+= strTaxaLenArr[genTaxaId];//NOTE: strTaxaLenArr is indexed based on the genTaxaId to which this strTaxaId maps to				
			}else{
			*dest= *pos;
			++dest;
			++pos;
			}
		}
		*dest = '\0';
		destTreeFile<<destBuf<<endl;
	}

	srcTreeFile.close();
	destTreeFile.close();
	delete []genTaxaIndStrTaxaBufArr;
	delete []strTaxaLenArr;
	delete []strTaxaBuf;
	delete []genSeqGenTaxaArr;
}


//loads strTaxa and genTaxaId-strTaxa mapping
void loadStrTaxaInfo(unsigned int *&genTaxaIndStrTaxaBufArr,unsigned char *&strTaxaLenArr,char *&strTaxaBuf)
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
	genTaxaIndStrTaxaBufArr = new unsigned int[genQueryTaxaC];
	strTaxaLenArr = new unsigned char[genQueryTaxaC];

	ifstream strTaxaDBFile (strTaxaDBFileName,ios::in | ios::binary);
	if(!strTaxaDBFile.is_open()){
		T_LOG<<"Error! Could not open strTaxaDBFileName."<<endl;
		return;
	}
	
	size_t charBufSz;
	strTaxaDBFile.read(reinterpret_cast<char *>(&charBufSz),sizeof(size_t));
	strTaxaBuf = new char[charBufSz];
	unsigned int strTaxaIdC, genTaxaId, tmp, bufInd=0;
	unsigned char orgTaxaIdC, strLen;
	strTaxaDBFile.read(reinterpret_cast<char *>(&tmp),sizeof(unsigned int));//mapBufSz:ignore, not needed here
	strTaxaDBFile.read(reinterpret_cast<char *>(&strTaxaIdC),sizeof(unsigned int));//# of strIds
	
	for(unsigned int i=0; i<strTaxaIdC; ++i){
		strTaxaDBFile.read(reinterpret_cast<char *>(&strLen),sizeof(unsigned char));
		strTaxaDBFile.read(&strTaxaBuf[bufInd],sizeof(char)*strLen);
		strTaxaDBFile.read(reinterpret_cast<char *>(&orgTaxaIdC),sizeof(unsigned char));
		for(unsigned char j=0; j<orgTaxaIdC; ++j){
			strTaxaDBFile.read(reinterpret_cast<char *>(&tmp),sizeof(unsigned int));//orgTaxaId: ignore
			strTaxaDBFile.read(reinterpret_cast<char *>(&genTaxaId),sizeof(unsigned int));
			genTaxaIndStrTaxaBufArr[genTaxaId]= bufInd;
			strTaxaLenArr[genTaxaId] = strLen;//NOTE: strTaxaLenArr is indexed based on the genTaxaId to which this strTaxaId maps to
		}		
		bufInd+= strLen;
		strTaxaBuf[bufInd++] = '\0';//terminating character and increment for the next free slot
	}
	strTaxaDBFile.close();
}
