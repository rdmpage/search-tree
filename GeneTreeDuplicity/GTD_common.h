#ifndef __GTD_COMMON__
#define __GTD_COMMON__
#include<iostream>
#include<set>
#include<ctime>
#include<fstream>
#include <vector>
#include <algorithm>
#include<limits.h>
#include <sstream>

using namespace std;
extern ofstream T_LOG_FILE;
// all #defines should come here, for the ease of reference

#ifdef DO_LOG
#define T_LOG T_LOG_FILE
#else
#define T_LOG std::cout
#endif

#ifndef WIN32
#include "./../TreeIndexing/server_socket.h"
#define SOCK_PATH_SERVER "./../sockets/gtd"
#endif

extern exception E; 

class TreeProcessor;
extern TreeProcessor TP;

class MemBuf;
extern MemBuf MB;

class IpOpStat;
extern IpOpStat IO_Stat;

extern unsigned int CUST_MEM_SZ;
extern unsigned int FIXED_MEM_SZ;

#endif //__GTD_COMMON__
