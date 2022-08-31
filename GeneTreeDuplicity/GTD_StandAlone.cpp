#include"GTD_MemBuf.h"
#include"GTD_TreeProcessor.h"

#include"GTD_Test.h"

#ifdef STAND_ALONE
unsigned int CUST_MEM_SZ =0;
unsigned int FIXED_MEM_SZ =0;
TreeProcessor TP;
MemBuf MB;
exception E; 
IpOpStat IO_Stat;
ofstream T_LOG_FILE("./GTD.log",ios::app);
#endif //STAND_ALONE


void processCommandLineQuery(int argc, char* argv[]);
void processCommandLineQueryBQred(int argc, char* argv[]);//tailored for reducing biclique clusters - sp req of Mike

int main(int argc, char* argv[])
{
	processCommandLineQuery(argc, argv);
	//processCommandLineQueryBQred(argc, argv);
	return 0;
}


void processCommandLineQuery(int argc, char* argv[])
{
	if(argc != 3){
		T_LOG<<"Usage:"<<endl;
		T_LOG<<argv[0]<<" TreeFile ResultFile"<<endl;
		exit(-1);
		return;
	}
  
	char *ipFileName = argv[1];
	char *resFileName = argv[2];

	ifstream ipFile(ipFileName);
	ofstream resFile(resFileName);

	if(!ipFile.is_open()){
		T_LOG<<"Error! Could not open: "<<ipFileName<<endl;
		exit(-1);
	}
	
	if(!resFile.is_open()){
		T_LOG<<"Error! Could not open: "<<resFileName<<endl;
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
		resFile<<">>"<<std::endl;
		prefixPos = strchr(MB.chBuf,'(');
		//*prefixPos = 0;
		//resFile<<MB.chBuf;
		*prefixPos = '(';

		TP.parseMulNewick(MB.chBuf,MB.treeBufInd);
		TP.computeAllSuvValues(MB.treeBufInd);
		unsigned short infoEdgeC = TP.collapseNonInformativeEdges(MB.treeBufInd);
		if(infoEdgeC <1){
			resFile<<"(No informative edge.)"<<";"<<std::endl;			
			MB.reset();
			IO_Stat.reset();
			continue;
		}
		TP.collapseDupEdges(MB.treeBufInd);
		TP.removeNonParticipatingLeaves(MB.treeBufInd);
		TP.removeDupLeaves(MB.treeBufInd);
		TP.remAllRedundantLeaves(MB.treeBufInd);
		TP.print(MB.treeBufInd);
		resFile<<MB.chBuf<<";"<<std::endl;//print MRF
		TP.computeStats(MB.treeBufInd);
		if(IO_Stat.leavesOp != IO_Stat.taxaOp){
			if(!TP.restrictToSingTaxa(MB.treeBufInd, MB.treeBufInd+1)){
				resFile<<"(A star-like tree.)"<<";"<<std::endl;			
				MB.reset();
				IO_Stat.reset();
				continue;
			}
			++MB.treeBufInd;
			TP.print(MB.treeBufInd);
			resFile<<MB.chBuf<<";"<<std::endl;		
		}
		
		MB.reset();
		IO_Stat.reset();
	}
	//resFile<<";"<<std::endl;
	//resFile<<"End;"<<std::endl;
	ipFile.close();
	resFile.close();
}




void processCommandLineQueryBQred(int argc, char* argv[])
{
	if(argc != 3){
		T_LOG<<"Usage:"<<endl;
		T_LOG<<argv[0]<<" TreeFile ResultFile"<<endl;
		exit(-1);
		return;
	}
  
	char *ipFileName = argv[1];
	char *resFileName = argv[2];

	ifstream ipFile(ipFileName);
	ofstream resFile(resFileName);

	if(!ipFile.is_open()){
		T_LOG<<"Error! Could not open: "<<ipFileName<<endl;
		exit(-1);
	}
	
	if(!resFile.is_open()){
		T_LOG<<"Error! Could not open: "<<resFileName<<endl;
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
		//resFile<<">>"<<std::endl;
		prefixPos = strchr(MB.chBuf,'(');
		//*prefixPos = 0;
		//resFile<<MB.chBuf;
		*prefixPos = '(';

		TP.parseMulNewickBQred(MB.chBuf,MB.treeBufInd);
		TP.computeAllSuvValues(MB.treeBufInd);
		unsigned short infoEdgeC = TP.collapseNonInformativeEdges(MB.treeBufInd);
		if(infoEdgeC <1){
			//resFile<<"(No informative edge.)"<<";"<<std::endl;			
			MB.reset();
			IO_Stat.reset();
			continue;
		}
		TP.collapseDupEdges(MB.treeBufInd);
		TP.removeNonParticipatingLeaves(MB.treeBufInd);
		TP.removeDupLeaves(MB.treeBufInd);
		TP.remAllRedundantLeaves(MB.treeBufInd);
		TP.print(MB.treeBufInd);
		resFile<<MB.chBuf<<";"<<std::endl;//print MRF
		TP.computeStats(MB.treeBufInd);
		if(IO_Stat.leavesOp != IO_Stat.taxaOp){
			if(!TP.restrictToSingTaxa(MB.treeBufInd, MB.treeBufInd+1)){
				resFile<<"(A star-like tree.)"<<";"<<std::endl;			
				MB.reset();
				IO_Stat.reset();
				continue;
			}						
			++MB.treeBufInd;//increased only if the tree is restricted to singly-labelled taxa
			TP.print(MB.treeBufInd);
			resFile<<MB.chBuf<<";"<<std::endl;		
		}else{
			resFile<<"MRF is singly labelled!"<<endl;
		}
		TP.printTableBQred(resFile,MB.treeBufInd);
		
		
		MB.reset();
		IO_Stat.reset();
	}
	//resFile<<";"<<std::endl;
	//resFile<<"End;"<<std::endl;
	ipFile.close();
	resFile.close();
}


