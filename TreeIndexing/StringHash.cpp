#include"StringHash.h"

bool String_Hash::compareStr(const char *str1, const char *str2)
{
	//assumed that stringt lenght is the same
	while(*str1){
		if(toupper(*str1++) != toupper(*str2++)){
			return false;
		}
	}
	return true;
}

unsigned int String_Hash::getHashCode(const char *str)
{
	//uses sdbm:http://www.cse.yorku.ca/~oz/hash.html
	unsigned int hash = 0;
	int c;

	while (c = *str++)
	{
		hash = c + (hash << 6) + (hash << 16) - hash;
	}
	return hash;

}

void String_Hash::insert(char *keyStr,unsigned int val)
{
	unsigned int key = getHashCode(keyStr);
	unsigned int rem,ind=0,quo=key;
	unsigned char count =0;

	while(0 != quo)
	{
		rem = quo%divFactor;
		ind = (rem*prArr[count++] +ind)%tableSize;
		quo = quo/divFactor;
	}
	if(0 == bin[ind].count)
	{//no element in this bin till now
		bin[ind].ind = bufInd;
		++bin[ind].count;
		StrHashElem &elem = buf[bufInd++];
		elem.key_hashCode =key;//store key
		elem.val =val;//store val		
		elem.next =UINT_MAX;//making ptr to next node null
		elem.strLen = (unsigned char)strlen(keyStr);
		elem.key_strInd = charBufInd;
		strcpy(charBufInd,keyStr);
		charBufInd+= elem.strLen+1;//+1 for '\0'
		++binCount;
	}
	else
	{//find next free slot
		++clashCount;
		if(maxBinSz<bin[ind].count+1){
			maxBinSz=bin[ind].count+1 ;
		}
		for(count=bin[ind].count++,ind=bin[ind].ind; count>1;--count,ind=buf[ind].next);
		buf[ind].next = bufInd;
		StrHashElem &elem = buf[bufInd++];
		elem.key_hashCode =key;//store key
		elem.val =val;//store val		
		elem.next =UINT_MAX;//making ptr to next node null
		elem.strLen =(unsigned char) strlen(keyStr);
		elem.key_strInd = charBufInd;
		strcpy(charBufInd,keyStr);
		charBufInd+= elem.strLen+1;//+1 for '\0'	
	}
}



unsigned int String_Hash::get(const char *keyStr)
{
	unsigned int key = getHashCode(keyStr);
	unsigned int rem,ind=0,quo=key;
	unsigned char count =0;
	while(0 != quo)
	{//calculate key
		rem = quo%divFactor;
		ind = (rem*prArr[count++] +ind)%tableSize;
		quo = quo/divFactor;
	}

	//find the value in the chained list
	for(count=bin[ind].count, ind=bin[ind].ind; count>0; --count,ind=buf[ind].next){
		if(buf[ind].key_hashCode == key){
			if(compareStr(buf[ind].key_strInd,keyStr)){
				return buf[ind].val;
			}
		}
	}
	
	//return ind;//this will return UINT_MAX i.e. when the key was not found
	return UINT_MAX;
}

void String_Hash::genPrime()
{
	unsigned int prC=0;
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
