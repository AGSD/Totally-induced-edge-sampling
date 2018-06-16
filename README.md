# Totally-induced-edge-sampling
Implementation of TIES graph sampling algorithm

ties_sequential uses the TIES algorithm to find a sample from a given input graph.
The code expects the input file to have edges represented with 2 numbers, specifying the 2 nodes between which the edge exists.
There must be 1 line for representing each edge, for example
1 3
4 5
2 5
. .
. .

The code runs on the data and creates an output file with a similar edge notation, 2 node numbers per line to specify an edge

The code works with c++11 std.

The datasets used are sourced from http://konect.uni-koblenz.de/networks/
The file datasetsInfo.txt has some basic information about the sample data sets. The files have been slightly modified to adhere
to the aformentioned format.