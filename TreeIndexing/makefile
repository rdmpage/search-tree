default:
	rm -f *.o *.out
	bash ./d2u.sh
#	CC *.cpp -D DO_LOG -D ST_CGI -g -lsocket -lnsl -o Step1Server.out
#-O3
#	g++ *.cpp -g  -lsocket -lnsl -o Step1Server.out
	CC *.cpp -O3  -lsocket -lnsl -o Step1Server.out
	CC STQueryStandAlone.cpp client_socket.cpp -g  -D ST_QUERY -lsocket -lnsl -o SearchTree.out
#	CC STQueryStandAlone.cpp client_socket.cpp -O3  -D ST_QUERY -lsocket -lnsl -o SearchTree.out
