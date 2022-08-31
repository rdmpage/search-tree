#include"GTD_TreeEntities.h"
#include"GTD_MemBuf.h"

bool  operator <(const TaxaLabel & lhs, const TaxaLabel & rhs)
{
	return (strcmp(MB.taxaBuf[lhs.ind],MB.taxaBuf[rhs.ind])<0);
}

bool  operator >(const TaxaLabel & lhs, const TaxaLabel & rhs)
{
	return (strcmp(MB.taxaBuf[lhs.ind],MB.taxaBuf[rhs.ind])>0);
}

bool  operator ==(const TaxaLabel & lhs, const TaxaLabel & rhs)
{
	return (strcmp(MB.taxaBuf[lhs.ind],MB.taxaBuf[rhs.ind])==0);
}

void MulTreeNode::collapseChild(unsigned int chInd)
{//works only when chInd is not a leaf node
	//OK1
	if(firstChInd != chInd){
		MB.nodeBuf[MB.nodeBuf[chInd].prevSibInd].nextSibInd =
			MB.nodeBuf[chInd].firstChInd;
		MB.nodeBuf[MB.nodeBuf[chInd].firstChInd].prevSibInd = 
			MB.nodeBuf[chInd].prevSibInd;
	}else{
		firstChInd = MB.nodeBuf[chInd].firstChInd;
		MB.nodeBuf[MB.nodeBuf[chInd].firstChInd].prevSibInd = 
			UINT_MAX;
	}

	if(lastChInd != chInd){
		MB.nodeBuf[MB.nodeBuf[chInd].nextSibInd].prevSibInd =
			MB.nodeBuf[chInd].lastChInd;
		MB.nodeBuf[MB.nodeBuf[chInd].lastChInd].nextSibInd =
			MB.nodeBuf[chInd].nextSibInd;
	}else{
		lastChInd = MB.nodeBuf[chInd].lastChInd;
		MB.nodeBuf[MB.nodeBuf[chInd].lastChInd].nextSibInd = 
			UINT_MAX;
	}
	numChild+= MB.nodeBuf[chInd].numChild -1;
	//though chInd is removed but to maintain the flow of iterating children, make its nextSibInd point to the next child of its current-
	//-parent after the collapsing operation- which is its first child
	MB.nodeBuf[chInd].nextSibInd = MB.nodeBuf[chInd].firstChInd;
}

void MulTreeNode::deleteChild(unsigned int chInd)
{
	//OK1
	if(firstChInd != chInd){
		MB.nodeBuf[MB.nodeBuf[chInd].prevSibInd].nextSibInd =
			MB.nodeBuf[chInd].nextSibInd;
	}else{
		firstChInd = MB.nodeBuf[chInd].nextSibInd;
	}

	if(lastChInd != chInd){
		MB.nodeBuf[MB.nodeBuf[chInd].nextSibInd].prevSibInd =
			MB.nodeBuf[chInd].prevSibInd;
	}else{
		lastChInd = MB.nodeBuf[chInd].prevSibInd;
	}
	--numChild;
}

void MulTreeNode::addChild(unsigned int parentI, unsigned int childInd)
{
		if(UINT_MAX != firstChInd){
		MB.nodeBuf[lastChInd].nextSibInd = childInd;//adjust the nextSibInd of the current last child to the new last child
		MB.nodeBuf[childInd].prevSibInd = lastChInd;
		lastChInd = childInd;
		MB.nodeBuf[childInd].nextSibInd = UINT_MAX;//init the nextSib of the child to null
		++numChild;
	}else{//first child being inserted
		lastChInd = firstChInd = childInd;
		MB.nodeBuf[childInd].nextSibInd = UINT_MAX; //init the nextSib of the child to null
		MB.nodeBuf[childInd].prevSibInd = UINT_MAX;
		numChild =1;
	}
}
void MulTreeNode::markForCollapse()
{
	//OK1
	tmpVal|= DO_COLLAPSE;
}

void MulTreeNode::chkN_Collapse(unsigned int chInd)
{
	//OK1
	if(MB.nodeBuf[chInd].tmpVal & DO_COLLAPSE){
		collapseChild(chInd);
	}
}

bool MulTreeNode::isInformative()
{
	//OK1
	if(Sdn <2 || Sup <2){
		return false;
	}else{
		return true;
	}
}

int MulTreeNode::compInfWidChld(unsigned int ndI)
{
	//OK1
	//it is assumed that Sdn of this node and Sup of ndI refer to the same node i.e. this node is parent of ndI
	if(Sdn != MB.nodeBuf[ndI].Sdn){
		if(Sup != MB.nodeBuf[ndI].Sup){
			return INFO_NA;
		}else{
			return INFO_GREATER;//this > ndI
		}
	}else{
		if(Sup != MB.nodeBuf[ndI].Sup){
			return INFO_LESSER;//this < ndI
		}else{
			return INFO_EQUAL;//this == ndI
		}
	}
}

int MulTreeNode::compInfWidSib(unsigned int ndI)
{
	//OK1
	//it is assumed that Sup of this node and Sup of ndI refer to the same node i.e. their common parent node
	if(Sdn != MB.nodeBuf[ndI].Sup){
		if(Sup != MB.nodeBuf[ndI].Sdn){
			return INFO_NA;
		}else{
			return INFO_LESSER;
		}
	}else{
		if(Sup != MB.nodeBuf[ndI].Sdn){
			return INFO_GREATER;
		}else{
			return INFO_EQUAL;
		}
	}
}

