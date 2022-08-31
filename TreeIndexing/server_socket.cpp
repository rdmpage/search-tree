#ifndef WIN32
#include "server_socket.h"


ServerSocket::ServerSocket(string sockPath,ostream &log_file):logFile(log_file)
{
	connected = false;
	if ((s = socket(AF_UNIX, SOCK_STREAM, 0)) == -1) {
		logFile<<"SERR:socket\n";
		exit(1);
	}

	local.sun_family = AF_UNIX;
	strcpy(local.sun_path, sockPath.c_str());
	unlink(local.sun_path);
	int len = strlen(local.sun_path) + sizeof(local.sun_family);
	if (bind(s, (struct sockaddr *)&local, len) == -1) {
		logFile<<"SERR:bind\n";
		exit(1);
	}

	if (listen(s, 5) == -1) {
		logFile<<"SERR:listen\n";
		exit(1);
	}

}

ServerSocket::~ServerSocket()
{
	disconnect();
}

void ServerSocket::disconnect()
{
	close(s2);
	connected = false;
}
void ServerSocket::connectTo()
{
	if(!connected)
	{
		logFile<<"SMSG:Waiting for a connection...\n";
		socklen_t t = sizeof(remote);
		if ((s2 = accept(s, (struct sockaddr *)&remote, &t)) == -1) {
			logFile<<"SERR:accept.\n";
			return;
		}
		logFile<<"SMSG:Connected.\n";
		connected = true;
	}
}

bool ServerSocket::receiveFrom(string &strLine)
{
	if(!connected)
	{
		logFile<<"SERR:Not connected.\n";
		return false;
	}
	int n;
	strLine ="";
	n = recv(s2, dataBuf, DATA_SIZE, 0);
	if (n > 0) {
		dataBuf[n] = '\0';
		strLine+=dataBuf;
		return true;
	}
	else{
		logFile<<"SERR:Receive error!\n";
		disconnect();
		return false;
	}
}

bool ServerSocket::sendTo(char *sendBuf)
{
	if(connected)
	{
		if (send(s2, sendBuf,strlen(sendBuf), 0) < 0)
		{
			logFile<<"SERR:send.\n";
			disconnect();
			return false;
		}	
		return true;
	}
	else
	{
		logFile<<"Not connected.\n";
		return false;
	}

}

#endif



