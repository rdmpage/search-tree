#include"ClusterEntities.h"
#include<cmath>
#include<cstring>

GlobalBuf gGlobalBuf;
unsigned short* gChildBuf = gGlobalBuf.childBuf;
TreeNode* gNodeBuf = gGlobalBuf.nodeBuf;
PhyloTree* gTreeBuf =gGlobalBuf.treeBuf;
unsigned int* gTaxaBuf =gGlobalBuf.taxaBuf;
ClusterInfo* gClusBuf =gGlobalBuf.clusBuf;
unsigned char* gAvgSupBuf =gGlobalBuf.avgSupBuf;

unsigned int NumberOfSetBits(unsigned int i)
{
    i = i - ((i >> 1) & 0x55555555);
    i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
    return ((i + (i >> 4) & 0xF0F0F0F) * 0x1010101) >> 24;
}

unsigned int TreeNode::sizeAsFile()
{
	unsigned int size =0;
	size+= sizeof(unsigned char); //to store the node type 0=intermediate, 1=leaf
	size+= sizeof(unsigned short);//to store numOfChild
	size+= sizeof(unsigned short);//to store bl
	if(UINT_MAX != childInd)
	{
		size+= sizeof(unsigned short) * numOfChild;//to store index of each of the children
	}
	return size;
}

void TreeNode::writeToFile(ofstream &file)
{
	unsigned char nodeType = UINT_MAX!=childInd?0:1;//asses nodeType, 0=intermediate, 1=leaf
	file.write(reinterpret_cast<char *>(&nodeType),sizeof(unsigned char)); // write the nodeType
	file.write(reinterpret_cast<char *>(&numOfChild),sizeof(unsigned short)); // write numOfChild (=nodeId for a leaf node)
	file.write(reinterpret_cast<char *>(&bl),sizeof(unsigned short)); // write bl i.e. branch length
	if(UINT_MAX != childInd)
	{
		for(unsigned short i=0; i<numOfChild; ++i)
		{
			file.write(reinterpret_cast<char *>(&gChildBuf[childInd+i]),sizeof(unsigned short)); // write the indexPtr for the child
		}
	}
}

void TreeNode::readFromFile(ifstream &file)
{
	unsigned char nodeType;
	file.read(reinterpret_cast<char *>(&nodeType),sizeof(unsigned char)); //read nodeType
	file.read(reinterpret_cast<char *>(&numOfChild),sizeof(unsigned short));//read numOfChild (=nodeId for a leaf node)
	file.read(reinterpret_cast<char *>(&bl),sizeof(unsigned short));//read bl i.e. branch length
	if(0 == nodeType) //an intermediate node 
	{
		childInd = gGlobalBuf.childBufInd;//init childInd with starting index in childBufInd
		gGlobalBuf.childBufInd+=numOfChild;//update childBufInd to point to next free slot
		//childArr = new unsigned short[numOfChild];//assign memory to childArr
		for(unsigned short i=0; i<numOfChild; ++i)
		{
			file.read(reinterpret_cast<char *>(&gChildBuf[childInd+i]),sizeof(unsigned short));//read each of the children
		}
	}
}


unsigned int PhyloTree::sizeAsFile()
{
	unsigned int size=0;
	size+=sizeof(unsigned short);//to store numOfNodes
	for(unsigned short i=0; i<numOfNodes; ++i)
	{
		size+=gNodeBuf[nodeInd+i].sizeAsFile();//size of each node
	}
	return size;
}

void PhyloTree::writeToFile(ofstream &file)
{
	file.write(reinterpret_cast<char *>(&numOfNodes),sizeof(unsigned short)); // write the numOfNodes
	for(unsigned short i=0; i<numOfNodes; ++i)
	{
		gNodeBuf[nodeInd+i].writeToFile(file);//write each node
	}
}

void PhyloTree::readFromFile(ifstream &file)
{
	file.read(reinterpret_cast<char *>(&numOfNodes),sizeof(unsigned short));//read numOfNodes
	nodeInd = gGlobalBuf.nodeBufInd;//init nodeInd
	gGlobalBuf.nodeBufInd+= numOfNodes;//update nodeBufInd for next free slot
	for(unsigned short i=0; i<numOfNodes; ++i)
	{
		gNodeBuf[nodeInd+i].readFromFile(file);//read each of the node
	}
}

//below is a test function
void PhyloTree::toSubtreeString(unsigned short &len,unsigned int taxaInd,char *buf
						 ,unsigned int *taxaUniqIdArr,unsigned int *seqIdUniqQueryIdArr,
						 unsigned int *queryTaxaUniqIdArr,unsigned int (*stack)[4],unsigned char labelFormat,unsigned int avgSupInd)
{
	unsigned int ndI,k,l;
	double taxaId;
	unsigned short stackC=0,j;
	size_t ep=0;
	unsigned char sz;
	stack[stackC][0] = nodeInd;//put root on stack
	stack[stackC][1] = gNodeBuf[nodeInd].numOfChild;
	stack[stackC][2] = nodeInd;//parent of the current node, root is the parent of itself
	//strTree = "";

	while(USHRT_MAX != stackC)
	{
		ndI = stack[stackC][0];
		if(UINT_MAX == gNodeBuf[ndI].childInd)
		{//a leaf node
			k = gNodeBuf[ndI].numOfChild;//taxaId in clusters taxa buf
			k = gTaxaBuf[taxaInd+k];//taxaId
			if(labelFormat == 0)
			{//both taxa id and seq id
				l= k;
				//first taxa id
				buf[ep++] = 't';
				buf[ep++] = 'i';
				taxaId= k = queryTaxaUniqIdArr[seqIdUniqQueryIdArr[k]];//org taxaId - here taxon id
				if(0 != taxaId)
				{
					ep+= sz = (unsigned char)log10(taxaId)+1;//size of taxaId
					while(k!=0)
					{
						buf[--ep] = '0' + k%10;
						k=k/10;
					}
					ep+=sz;
				}
				else
				{
					buf[ep++] = '0';
				}			
				
				buf[ep++] = '_';


				//now seq id
				buf[ep++] = 'g';
				buf[ep++] = 'i';
				taxaId = k = taxaUniqIdArr[l];//org taxaId - here seq id
				if(0 != taxaId)
				{
					ep+= sz = (unsigned char)log10(taxaId)+1;//size of taxaId
					while(k!=0)
					{
						buf[--ep] = '0' + k%10;
						k=k/10;
					}
					ep+=sz;
				}
				else
				{
					buf[ep++] = '0';
				}			
				
			}
			else if(labelFormat == 'S')
			{
				taxaId = k = taxaUniqIdArr[k];//org taxaId - here seq id
				if(0 != taxaId)
				{
					ep+= sz = (unsigned char)log10(taxaId)+1;//size of taxaId
					while(k!=0)
					{
						buf[--ep] = '0' + k%10;
						k=k/10;
					}
					ep+=sz;
				}
				else
				{
					buf[ep++] = '0';
				}			
				
			}
			else
			{
				taxaId= k = queryTaxaUniqIdArr[seqIdUniqQueryIdArr[k]];//org taxaId - here taxon id
				if(0 != taxaId)
				{
					ep+= sz = (unsigned char)log10(taxaId)+1;//size of taxaId
					while(k!=0)
					{
						buf[--ep] = '0' + k%10;
						k=k/10;
					}
					ep+=sz;
				}
				else
				{
					buf[ep++] = '0';
				}			
				
			}

			//write the branch length
			buf[ep++] = '-';
			taxaId = k = gNodeBuf[ndI].bl-gNodeBuf[stack[stackC][2]].bl;//branch length- take the relative diff with its parent
			if(0!= taxaId && gNodeBuf[ndI].bl >= taxaId)
			{
				ep+= sz = (unsigned char)log10(taxaId)+1;//size of taxaId
				while(k!=0)
				{
					buf[--ep] = '0' + k%10;
					k=k/10;
				}
				ep+=sz;
			}
			else
			{
				buf[ep++] = '0';
			}			
			
			buf[ep++] = ',';
			--stackC;
		}
		else
		{//an intermediate node
			if(0 == stack[stackC][1])
			{//all its children have been visited
				buf[--ep] = ')';//overwirte the extra ','

				//write the branch length
				++ep;
				buf[ep++] = '-';
				taxaId = k = gNodeBuf[ndI].bl-gNodeBuf[stack[stackC][2]].bl;//branch length- take the relative diff with its parent
				if(0!= taxaId && gNodeBuf[ndI].bl >= taxaId)
				{
					ep+= sz = (unsigned char)log10(taxaId)+1;//size of taxaId
					while(k!=0)
					{
						buf[--ep] = '0' + k%10;
						k=k/10;
					}
					ep+=sz;
				}
				else
				{
					buf[ep++] = '0';
				}			
				

				buf[ep++] =',';
				--stackC;
			}
			else
			{//encountered first time
				buf[ep++]='(';
				stack[stackC][1] =0;//set all children visited and put all its children on the stack
				k = gNodeBuf[ndI].childInd;
				for(j=0; j<gNodeBuf[ndI].numOfChild; ++j)
				{
					stack[++stackC][0]= nodeInd + gChildBuf[k+j];//put a child and increment stackC
					stack[stackC][1] = gNodeBuf[stack[stackC][0]].numOfChild;//mark that children needs to be visited
					stack[stackC][2] = ndI;//index of current parent, would be needed while calculating bl
				}
			}
		}
	}
	ep-= 3;//overwrite the last ':0,'
	buf[ep]='\0';
	len = (unsigned short)(ep+1);
}


void PhyloTree::toString(unsigned short &len,unsigned int taxaInd,PrintBuf &pb,
							  unsigned int (*stack)[4],unsigned char labelFormat,unsigned int avgSupInd)//converts it to a string
{
	char *buf = pb.charBuf;
	unsigned int ndI,k,l;
	double taxaId;
	unsigned short stackC=0,j;
	size_t ep=0;
	unsigned char sz;
	//check for root node deg two - START
	if(true){
	//if(2 != gNodeBuf[nodeInd].numOfChild){		
		stack[stackC][0] = nodeInd;//put root on stack
		stack[stackC][1] = gNodeBuf[nodeInd].numOfChild;
	}else{
		//assume a pseudo root node whose children are the combined set of children of both of its children
		stack[stackC][0] = nodeInd;//put root on stack
		stack[stackC][1] = 0;//but mark it as visited
		buf[ep++]='(';
		//now put the combined set of children
		ndI = gNodeBuf[nodeInd].childInd;
		for(j=0; j<2; ++j){
			unsigned int chI = nodeInd + gChildBuf[ndI+j];
			if(UINT_MAX != gNodeBuf[chI].childInd){
				k = gNodeBuf[chI].childInd;
				for(j=0; j<gNodeBuf[chI].numOfChild; ++j)
				{
					stack[++stackC][0]= nodeInd + gChildBuf[k+j];//put a child and increment stackC
					stack[stackC][1] = gNodeBuf[stack[stackC][0]].numOfChild;//mark that children needs to be visited
				}
			}else{//a leaf node 
				stack[++stackC][0]= chI;//put this node itself on the stack
				stack[stackC][1] = gNodeBuf[stack[stackC][0]].numOfChild;//mark that children needs to be visited
			}
		}		
	}
	//check for root node deg two - END	

	while(USHRT_MAX != stackC)
	{
		ndI = stack[stackC][0];
		if(UINT_MAX == gNodeBuf[ndI].childInd)
		{//a leaf node
			k = gNodeBuf[ndI].numOfChild;//taxaId in clusters taxa buf
			k = gTaxaBuf[taxaInd+k];//taxaId (seqId)
			buf[ep++] = '\'';
			if(labelFormat == 'N'){//taxon name
				l = pb.seqIdUniqQueryIdArr[k];//genQueryTaxaId
				strcpy(&buf[ep],&pb.strTaxaBuf[pb.genTaxaIndStrTaxaBufArr[l]]);
				ep+= pb.strTaxaLenArr[l];//NOTE: strTaxaLenArr is indexed based on the genTaxaId to which this strTaxaId maps to				
			} else if(labelFormat == 'A'){//all i.e. taxonName_taxaId_SeqId
				l = pb.seqIdUniqQueryIdArr[k];//genQueryTaxaId
				strcpy(&buf[ep],&pb.strTaxaBuf[pb.genTaxaIndStrTaxaBufArr[l]]);
				ep+= pb.strTaxaLenArr[l];//NOTE: strTaxaLenArr is indexed based on the genTaxaId to which this strTaxaId maps to
				l= k;

				buf[ep++] = '-';
				//buf[ep++] = '_';
				//first taxa id
				buf[ep++] = 't';
				buf[ep++] = 'i';
				taxaId= k = pb.queryTaxaUniqIdArr[pb.seqIdUniqQueryIdArr[k]];//org taxaId - here taxon id
				if(0 != taxaId)
				{
					ep+= sz = (unsigned char)log10(taxaId)+1;//size of taxaId
					while(k!=0)
					{
						buf[--ep] = '0' + k%10;
						k=k/10;
					}
					ep+=sz;
				}
				else
				{
					buf[ep++] = '0';
				}			

				buf[ep++] = '_';


				//now seq id
				buf[ep++] = 'g';
				buf[ep++] = 'i';
				taxaId = k = pb.taxaUniqIdArr[l];//org taxaId - here seq id
				if(0 != taxaId)
				{
					ep+= sz = (unsigned char)log10(taxaId)+1;//size of taxaId
					while(k!=0)
					{
						buf[--ep] = '0' + k%10;
						k=k/10;
					}
					ep+=sz;
				}
				else
				{
					buf[ep++] = '0';
				}				
			}
			else if(labelFormat == 0)
			{//both taxa id and seq id
				l= k;
				//first taxa id
				buf[ep++] = 't';
				buf[ep++] = 'i';
				taxaId= k = pb.queryTaxaUniqIdArr[pb.seqIdUniqQueryIdArr[k]];//org taxaId - here taxon id
				if(0 != taxaId)
				{
					ep+= sz = (unsigned char)log10(taxaId)+1;//size of taxaId
					while(k!=0)
					{
						buf[--ep] = '0' + k%10;
						k=k/10;
					}
					ep+=sz;
				}
				else
				{
					buf[ep++] = '0';
				}			
				
				buf[ep++] = '_';


				//now seq id
				buf[ep++] = 'g';
				buf[ep++] = 'i';
				taxaId = k = pb.taxaUniqIdArr[l];//org taxaId - here seq id
				if(0 != taxaId)
				{
					ep+= sz = (unsigned char)log10(taxaId)+1;//size of taxaId
					while(k!=0)
					{
						buf[--ep] = '0' + k%10;
						k=k/10;
					}
					ep+=sz;
				}
				else
				{
					buf[ep++] = '0';
				}			
				
			}
			else if(labelFormat == 'S')
			{
				taxaId = k = pb.taxaUniqIdArr[k];//org taxaId - here seq id
				if(0 != taxaId)
				{
					ep+= sz = (unsigned char)log10(taxaId)+1;//size of taxaId
					while(k!=0)
					{
						buf[--ep] = '0' + k%10;
						k=k/10;
					}
					ep+=sz;
				}
				else
				{
					buf[ep++] = '0';
				}			
				
			}
			else 
			{
				taxaId= k = pb.queryTaxaUniqIdArr[pb.seqIdUniqQueryIdArr[k]];//org taxaId - here taxon id
				if(0 != taxaId)
				{
					ep+= sz = (unsigned char)log10(taxaId)+1;//size of taxaId
					while(k!=0)
					{
						buf[--ep] = '0' + k%10;
						k=k/10;
					}
					ep+=sz;
				}
				else
				{
					buf[ep++] = '0';
				}			
				
			}

			
			//not required for leaf node as it will always have support = 100%
			////write the avgSupNum
			//buf[ep++] = '<';
			//
			//taxaId = k = gAvgSupBuf[avgSupInd+ ndI - nodeInd];//indexing of supNum is done in a way similar to indexing of-
			////-nodes of the tree, so that retreival of supNum can happen similarly
			//if(0 != taxaId)
			//{
			//	ep+= sz = (unsigned char)log10(taxaId)+1;//size of taxaId
			//	while(k!=0)
			//	{
			//		buf[--ep] = '0' + k%10;
			//		k=k/10;
			//	}
			//	ep+=sz;
			//}
			//else
			//{
			//	buf[ep++] = '0';
			//}			

			//buf[ep++] = '%';
			//buf[ep++] = '>';
			buf[ep++] = '\'';
			buf[ep++] = ',';
			--stackC;
		}
		else
		{//an intermediate node
			if(0 == stack[stackC][1])
			{//all its children have been visited
				buf[--ep] = ')';//overwirte the extra ','				
				++ep;

				//buf[ep++] = '<';//commented as of now as it does not look good for intermediate nodes
				//write the avgSupNum
#ifndef PRINT_SUP_FALSE
				taxaId = k = gAvgSupBuf[avgSupInd+ ndI - nodeInd];//indexing of supNum is done in a way similar to indexing of-
				//-nodes of the tree, so that retreival of supNum can happen similarly
				if(0 != taxaId)
				{
					ep+= sz = (unsigned char)log10(taxaId)+1;//size of taxaId
					while(k!=0)
					{
						buf[--ep] = '0' + k%10;
						k=k/10;
					}
					ep+=sz;
				}
				else
				{
					buf[ep++] = '0';
				}			
				
				buf[ep++] = '%';
#endif
				//buf[ep++] = '>';//commented as of now as it does not look good for intermediate nodes
				buf[ep++] =',';
				--stackC;
			}
			else
			{//encountered first time
				buf[ep++]='(';
				stack[stackC][1] =0;//set all children visited and put all its children on the stack
				k = gNodeBuf[ndI].childInd;
				for(j=0; j<gNodeBuf[ndI].numOfChild; ++j)
				{
					stack[++stackC][0]= nodeInd + gChildBuf[k+j];//put a child and increment stackC
					stack[stackC][1] = gNodeBuf[stack[stackC][0]].numOfChild;//mark that children needs to be visited
				}
			}
		}
	}
#ifndef PRINT_SUP_FALSE
	ep-= 5;//overwrite the last ':100,'
#else 
	ep-= 1;//overwrite the last ','	
#endif

	buf[ep]='\0';
	len = (unsigned short)(ep+1);
	//strTree.append(buf,ep+1);
}

unsigned int ClusterInfo::sizeAsFile()
{
	unsigned int size=0;
	//cluster id is not stored here as it is not required, it will be stored separately as a clusterId-filePtr mapping
	size+= sizeof(numOfTaxa); // to store numOfTaxa
	size+= sizeof(unsigned int)*numOfTaxa; //to store each of the taxaId
	size+= sizeof(numOfTrees); //to store numOfTrees
	for(unsigned short i=0; i<numOfTrees; ++i)
	{
		size+= gTreeBuf[treeInd+i].sizeAsFile();//add size of each tree
	}
	return size;
}

void ClusterInfo::writeToFile(ofstream &file)
{
	file.write(reinterpret_cast<char *>(&numOfTaxa),sizeof(numOfTaxa)); // write the numOfTaxa
	for(unsigned short i=0; i<numOfTaxa; ++i)
	{
		file.write(reinterpret_cast<char *>(&gTaxaBuf[taxaInd+i]),sizeof(unsigned int)); // write each of the taxaId
	}

	file.write(reinterpret_cast<char *>(&numOfTrees),sizeof(numOfTrees)); // write the numOfTrees
	for(unsigned short i=0; i<numOfTrees; ++i)
	{
		gTreeBuf[treeInd+i].writeToFile(file);//write each of the tree
	}
}

void ClusterInfo::readFromFile(ifstream &file, unsigned short maxTreeC)
{
	file.read(reinterpret_cast<char *>(&numOfTaxa),sizeof(unsigned short));//read numOfTaxa
	taxaInd= gGlobalBuf.taxaBufInd;//init taxaInd
	gGlobalBuf.taxaBufInd+=numOfTaxa;//update treeBufInd for the next free slot
	for(unsigned short i=0; i<numOfTaxa; ++i)
	{
		file.read(reinterpret_cast<char *>(&gTaxaBuf[taxaInd+i]),sizeof(unsigned int));//read each taxaId
	}

	file.read(reinterpret_cast<char *>(&numOfTrees),sizeof(numOfTrees));//read numOfTrees
	if(maxTreeC<numOfTrees){
		numOfTrees = maxTreeC;//i.e. read only the first maxTreeC trees
	}
	treeInd= gGlobalBuf.treeBufInd;//init treeInd
	gGlobalBuf.treeBufInd+= numOfTrees;//update treeBufInd for the next free slot
	for(unsigned short i=0; i<numOfTrees; ++i)
	{
		gTreeBuf[treeInd+i].readFromFile(file);//read each of the tree
	}
}



const string phyloTreeConstants::suffixString = "                    ";

const int phyloTreeConstants::maxTreeNAmeSize = 20;
