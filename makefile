all: main 

main: main.o Help.o DDCal.o  Graph.o UncertainGraph.o truncated_normal.o
	g++ -o main main.o Help.o DDCal.o  Graph.o UncertainGraph.o truncated_normal.o /usr/local/lib/libigraph.a  -lm 
main.o : main.cpp
	g++ -Wall -c main.cpp -std=c++11 -I /usr/local/include/igraph/

DDCal.o: DDCal.cpp
	g++ -Wall -c DDCal.cpp -std=c++11 -I /usr/local/include/igraph/
Help.o: Help.cpp
	g++ -Wall -c Help.cpp  -std=c++11 -I /usr/local/include/igraph/
truncated_normal.o: truncated_normal.cpp
	g++ -Wall -c truncated_normal.cpp -std=c++11 -I /usr/local/include/igraph/
Graph.o: Graph.cpp
	g++ -Wall -c Graph.cpp  -std=c++11 -I /usr/local/include/igraph/
UncertainGraph.o: UncertainGraph.cpp
	g++ -Wall -c UncertainGraph.cpp -std=c++11  -I /usr/local/include/igraph/
clean:
	rm *.o
	rm main 
	
