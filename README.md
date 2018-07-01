# Totally-induced-edge-sampling
Implementation of TIES graph sampling algorithm, sourced from https://docs.lib.purdue.edu/cgi/viewcontent.cgi?article=2743&context=cstech, which has been optimised by making use a compressed sparse row representation of the graph.

ties_sequential uses the TIES algorithm to find a sample from a given input graph.
Usage: ties_sequential inputFilename outputFilename samplingRatio
Compile: g++ -std=c++11 ties_sequential.cpp -o ties_sequential.o

ties_parallel works the same, but makes use of multiple threads to sample the graph and induce edges.
Usage: ties_parallel inputFilename outputFilename samplingRatio
Compile: g++ -std=c++11 -fopenmp ties_parallel.cpp -o ties_parallel.o

The code expects the following input file conventions:
1)The file should have one line for each edge
2)Each line should be a pair of numbers, both greater than or equal to 1, specifying a directed edge
3)An undirected graph with edges a-b, should be represented with a line a b, as well as b a, for all edges.

The code runs on the data and creates an output file with a similar edge notation, 2 node numbers per line to specify an edge.
The code works with the c++11 standard, and OpenMP for the parallel code.

The sample datasets used are sourced from http://konect.uni-koblenz.de/networks/
The file datasetsInfo.txt has some basic information about the sample data sets. The files have been slightly modified to adhere to the aformentioned format.
