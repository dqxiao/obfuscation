all: main

main: main.o Help.o Configuration.o 
	clang++ -std=c++0x main.o Help.o Configuration.o 

DDCal.o: DDCal.cpp 
	g++ -Wall -c DDCal.cpp

Help.o: Help.cpp Help.hpp
	g++ -Wall -c Help.cpp 

Graph.o: Graph.cpp Graph.hpp
	g++ -Wall -c Graph.cpp
UncertainGraph.o: UncertainGraph.cpp UncertainGraph.hpp
	g++ -Wall -c UncertainGraph.cpp

truncated_normal.o: truncated_normal.cpp truncated_normal.hpp
	g++ -Wall -c truncated_normal.cpp

RandHelp.o: RandHelp.cpp RandHelp.hpp
	g++ -Wall -c RandHelp.cpp
clean:
	rm -rf *.o 