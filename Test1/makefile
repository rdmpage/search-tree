default:
	rm -f *.o *.out
	bash ./d2u.sh
#	C *.cpp ./../TreeIndexing/server_socket.cpp -D DO_LOG -g -lsocket  -lnsl -o  Test1.out
	CC  *.cpp ./../TreeIndexing/server_socket.cpp -g -D DO_LOG -lsocket  -lnsl -o  Test1.out
#	CC  *.cpp ./../TreeIndexing/server_socket.cpp -O3 -lsocket  -lnsl -o  Test1.out
