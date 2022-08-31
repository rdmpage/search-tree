#include"GTD_util.h"

unsigned int UL_HashTable::get(unsigned int key)
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

unsigned int * UL_HashTable::getPtr(unsigned int key)
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

void UL_HashTable::insert(unsigned int key,unsigned val)
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
		if(maxBinSize<bin[ind][1])
		{
			maxBinSize=bin[ind][1];
		}
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

void UL_HashTable::genPrime()
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

