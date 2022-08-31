#include"GTD_MemBuf.h"
#include"GTD_TreeProcessor.h"

#include"GTD_Test.h"

unsigned int CUST_MEM_SZ =0;
unsigned int FIXED_MEM_SZ =0;
TreeProcessor TP;
MemBuf MB;
exception E; 
IpOpStat IO_Stat;
ofstream T_LOG_FILE("./GTD.log",ios::app);

#ifndef WIN32
	ServerSocket serverSock(SOCK_PATH_SERVER,T_LOG);
#endif
void processQuery(char *ipFileName, char *resFileName, char *statFileName);
void getQueryFileNumber(unsigned int &queryFileNum);
void sendResFileNum(unsigned int resCount);



int _main(int argc, char* argv[])
{
	EvalST_MulTrees em;
	//em.evalAllClusters();
	em.evalAllBestTrees();
	
	return 0;
}



int st_main(int argc, char* argv[])//for searchTree interface
{
	if(!T_LOG_FILE.is_open())
	{
		std::cout<<"FATAL ERROR! T_LOG_FILE could not be opened!\n";
		exit(-1);
	}

	char ipFileName[1024];
	char opFileName[1024];
	char statFileName[1024];
unsigned int queryFileNum;
	while(true){
		getQueryFileNumber(queryFileNum);
		T_LOG<<"Query received: "<<queryFileNum<<endl;
		sprintf(ipFileName,"../common/queryInput/gtd/Query%d",queryFileNum);
		sprintf(opFileName,"../common/queryResult/gtd/Result%d",queryFileNum);
		sprintf(statFileName,"../common/queryResult/gtd/Stat%d",queryFileNum);
		
		processQuery(ipFileName,opFileName,statFileName);
		sendResFileNum(queryFileNum);
		T_LOG<<"Query processed: "<<queryFileNum<<endl;
#ifdef WIN32
		break;
#endif
	}
	return 0;	
}

void processQuery(char *ipFileName, char *resFileName, char *statFileName)
{
	ifstream ipFile(ipFileName);
	ofstream resFile(resFileName);
	ofstream statFile(statFileName);
	
	if(!ipFile.is_open() || !resFile.is_open() || !statFile.is_open()){
		T_LOG<<"Error! Could not open file."<<endl;
		exit(-1);
	}
	//resFile<<"#Nexus\n";
	//resFile<<"Begin trees;\n";
	char *prefixPos;
	while(true){
		ipFile.getline(MB.chBuf,2*MB.TAXA_MAX*MB.TAXA_LEN_MAX);
		if(0 == ipFile.gcount() || strlen(MB.chBuf) == 0){
			break;
		}
		prefixPos = strchr(MB.chBuf,'(');
		*prefixPos = 0;
		//resFile<<MB.chBuf;
		*prefixPos = '(';
		
		TP.parseMulNewick(MB.chBuf,MB.treeBufInd);
		TP.computeAllSuvValues(MB.treeBufInd);
		unsigned short infoEdgeC = TP.collapseNonInformativeEdges(MB.treeBufInd);
		if(infoEdgeC <1){
			resFile<<"(No informative edge.)"<<";"<<std::endl;
			statFile<<"No informative edge found in the reduced tree."<<std::endl;
			MB.reset();
			IO_Stat.reset();
			continue;
		}
		TP.collapseDupEdges(MB.treeBufInd);
		TP.removeNonParticipatingLeaves(MB.treeBufInd);
		TP.removeDupLeaves(MB.treeBufInd);
		TP.remAllRedundantLeaves(MB.treeBufInd);		
		if(!TP.restrictToSingTaxa(MB.treeBufInd, MB.treeBufInd+1)){
			resFile<<"(A star-like tree.)"<<";"<<std::endl;
			statFile<<"Reduced tree is star-like."<<std::endl;
			MB.reset();
			IO_Stat.reset();
			continue;
		}
		++MB.treeBufInd;
		TP.print(MB.treeBufInd);
		resFile<<MB.chBuf<<";"<<std::endl;
		TP.computeStats(MB.treeBufInd);
		unsigned int ipEff = 2*100*(IO_Stat.internalNodesOp-1)/(IO_Stat.internalNodesIp -1 + IO_Stat.leavesIp);
		unsigned int opEff = 2*100*(IO_Stat.internalNodesOp-1)/(IO_Stat.internalNodesOp -1 + IO_Stat.leavesOp);
		//statFile<<"Efficiency gain: "<<ipEff<<"% => "<<opEff<<"%; ";
		statFile<<"Taxa lost: "<<IO_Stat.taxaIp-IO_Stat.taxaOp<<"; ";
		statFile<<"Nodes Removed: "<<IO_Stat.internalNodesIp + IO_Stat.leavesIp - IO_Stat.internalNodesOp - IO_Stat.leavesOp<<";"<<endl;
		MB.reset();
		IO_Stat.reset();
	}
	//resFile<<";"<<std::endl;
	//resFile<<"End;"<<std::endl;
	ipFile.close();
	resFile.close();
	statFile.close();
}


void getQueryFileNumber(unsigned int &queryFileNum)
{
#ifndef WIN32
	istringstream stm;
	string strLine;
	while(true)
	{
		serverSock.connectTo();
		if(!serverSock.receiveFrom(strLine))
		{
			continue;
		}
		break;
	}
	
	stm.str(strLine);
	stm>>queryFileNum;
	T_LOG<<"\nQueryFileNum: "<<queryFileNum;
#else
	queryFileNum = 1;
#endif
}

void sendResFileNum(unsigned int resCount)
{
#ifndef WIN32
	char strBuff[257];
	sprintf(strBuff,"%d",resCount);
	serverSock.sendTo(strBuff);
#endif
}
