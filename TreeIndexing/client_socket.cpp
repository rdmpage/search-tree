#ifndef WIN32

#include"client_socket.h"

ClientSocket::ClientSocket(string sock_path,ostream &log_file):logFile(log_file)
{
	sockPath = sock_path;
	connected = false;
	if ((s = socket(AF_UNIX, SOCK_STREAM, 0)) == -1)
	{
		logFile<<"SERR:socket.\n";
		exit(1);
	}
}

ClientSocket::~ClientSocket()
{
	disconnect();
}

void ClientSocket::disconnect()
{
	close(s);
	connected = false;
}
void ClientSocket::connectTo()
{
	if(!connected)
	{
		logFile<<"SMSG:Trying to connect...\n";

		remote.sun_family = AF_UNIX;
		strcpy(remote.sun_path, sockPath.c_str());
		len = strlen(remote.sun_path) + sizeof(remote.sun_family);
		if (connect(s, (struct sockaddr *)&remote, len) == -1) {
			logFile<<"SERR:connect.\n";
			return;
		}
		logFile<<"SMSG:Connected.\n";
		connected = true;
	}
}

bool ClientSocket::sendTo(const char *str)
{
	if(!connected)
	{
		logFile<<"SERR:Not connected.\n";
		return false;
	}
	if (send(s, str, strlen(str), 0) == -1)
	{
		logFile<<"SERR:send.\n";
		disconnect();
		return false;
	}
	return true;
}

bool ClientSocket::receiveFrom(string &strLine)
{
	char str[100] ={0};
	strLine="";
	if(!connected)
	{
		logFile<<"SERR:not connected\n.";
		return false;
	}
	if ((t=recv(s, str, 100, 0)) > 0) {
		str[t] = '\0';
		strLine+= str;
		logFile<<"SMSG:echo> "<<str<<" .\n";
		return true;
	} else {
		if (t < 0) logFile<<"SERR:recv.\n";
		else logFile<<"SMSG:Server closed connection.\n";
		disconnect();
		return false;
	}
	
}

#endif
