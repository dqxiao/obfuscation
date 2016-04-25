rm *.o
g++ -Wall -c DDCal.cpp  -I /usr/local/include/igraph/
g++ -Wall -c Help.cpp  -I /usr/local/include/igraph/
g++ -Wall -c truncated_normal.cpp  -I /usr/local/include/igraph/
g++ -Wall -c Graph.cpp  -I /usr/local/include/igraph/
g++ -Wall -c UncertainGraph.cpp  -I /usr/local/include/igraph/
g++ -Wall -c main.cpp -I /usr/local/include/igraph/
g++  main.o Help.o DDCal.o  Graph.o UncertainGraph.o truncated_normal.o -o main
