#ifndef __SERVER_SOCKET_H__
#define __SERVER_SOCKET_H__

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <unistd.h>
#include <string>
#include <iostream>
#include<fstream>

using namespace std;
#define DATA_SIZE 100
class ServerSocket
{
	int s, s2;
	struct sockaddr_un local, remote;
	bool connected;
	char dataBuf[DATA_SIZE+1];
	ostream &logFile;
public:
	ServerSocket(string sockPath,ostream &log_file);
	~ServerSocket();
	void connectTo();
	bool receiveFrom(string &strLine);
	bool sendTo(char *sendBuf);
	void disconnect();
};



#endif //__SERVER_SOCKET_H__
