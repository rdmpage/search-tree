#ifndef __CLIENT_SOCKET__
#define __CLIENT_SOCKET__

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


using namespace std;

class ClientSocket
{
	int s, t, len;
	struct sockaddr_un remote;
	bool connected;
	string sockPath;
	ostream &logFile;
	public:
	ClientSocket(string sock_path,ostream &log_file);
	~ClientSocket();
	void connectTo();
	bool sendTo(const char *str);
	bool receiveFrom(string &strLine);
	void disconnect();
};



#endif //__CLIENT_SOCKET__
