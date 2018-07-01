/*  Author: Aditya Gaur
    Usage: ties_parallel inputFilename outputFilename samplingRatio
    About: Code to create a sample graph from an input graph with final number of nodes based on a sampling ratio.
           Algorithm used is TIES (Totally induced edge sampling), based on https://docs.lib.purdue.edu/cgi/viewcontent.cgi?article=2743&context=cstech
    Compile with: g++ -std=c++11 -fopenmp ties_parallel.cpp -o ties_parallel.o
*/
#include<iostream>
#include<fstream>
#include<vector>
#include<random>
#include<string>
#include<cstdio>
#include<ctime>
#include<utility>

//Toggle for verbose comments
#define vbs(x) x
//Number of threads to use
#define NUM_THREADS 8
//Number of additions in a sampling thread, after which we make an update to the total
#define UPDATE_COUNTER 20
//Toggle Debugging output
#define dbg(x)

using namespace std;

class edge{
public:
    long u,v;
    edge(long a, long b){
        u = a;
        v = b;
    }
    edge(){
    	u=0;
    	v=0;
    }
};

struct graph{
	long *vi;
	long *ei;
};

long currentVsize=0;		//currently number of sampled nodes
long requiredVsize;		//Required subsample size(number of vertices)
graph g;			//csr version of graph

long inputGraph(string, vector<edge>&);					//extract graph from input file
void make_csr(vector<edge>,long,long);					//convert graph to compressed sparse row format
pair<edge*,long> sampleGraph(vector<edge>&, double, long);		//sample original graph
vector<long> sampleEdges(vector<edge>&, bool*, long, long, long);	//sample edges within a thread
vector<edge> induceEdges(long*, graph&, bool*, long, long, long);	//induce edges within a thread
void writeGraph(string,edge*, long);					//write output graph to a file 
void parseCommandlineArgs(int,char* [], string&, string&, double&);	//parse input parameters from command line

int main(int argc, char* argv[]) {

    //Parameters to be specified, can be taken from commandline using parseCommandlineArgs
    string inFilename;
    string outFilename;
    double fi;
    
    //Working variables
    vector<edge> inputEdgeList;	//extracted list of edges from input file
    edge *sampledEdgeList;	//sampled output list of edges from sampling the graph
    long sampledSize;		//holds size of sampled subgraph
    long n,m;			//number of vertices and edges in original graph
    pair<edge*,long> sampled;	//holds return edge list and list size
    
    parseCommandlineArgs(argc, argv, inFilename, outFilename, fi);	//read input file name, ouput file name and value of fi(sampling ratio)

    vbs(cout<<"Getting input"<<endl;)

    n = inputGraph(inFilename,inputEdgeList);	//read edge list from input file, returns number of vertices
    m = inputEdgeList.size();			//number of edges in original graph
    
    vbs(cout<<"Input succeeded, graph has "<<n<<" vertices and "<<m<<" edges"<<endl;)

    make_csr(inputEdgeList, n, m);		//convert graph to csr format(helps in induction)
	
    sampled = sampleGraph(inputEdgeList,fi,n);	//sample the graph using fi
    sampledEdgeList = sampled.first;
    sampledSize = sampled.second;
    
    writeGraph(outFilename, sampledEdgeList, sampledSize);	//write sampled graph to file

    vbs(cout<<"Output graph written to "<<outFilename<<endl;)
    return 0;
}

pair<edge*,long> sampleGraph(vector<edge> &edgeList, double fi, long n){
	
	//initialising global variables
	currentVsize = 0;
	requiredVsize = fi*n;
	
	long m = edgeList.size();					//number of edges in original graph
	bool *vExist = new bool[n+1]();					//set to true for each vertex which is sampled

	dbg(printf("required size = %ld\n",requiredVsize);)

    	long totalNodes=0;					//Number of sampled vertices
	vector<long> tmpNodes[NUM_THREADS];		//individual vectors for each thread
	long *tmpNodeSizes = new long[NUM_THREADS]();	//Number of sampled vertices by each thread
	long *tmpNodeStart = new long[NUM_THREADS]();	//Starting index for copying vertices from each thread
	long psize = (m+NUM_THREADS-1)/NUM_THREADS;	//partition size for each thread sampling edges
	
	//Sample edges in parallel
    	#pragma omp parallel shared(edgeList, vExist, n, m, fi, psize, tmpNodes, tmpNodeSizes) num_threads(NUM_THREADS)
    	{
    		#pragma omp for
    		for(long i=0; i<NUM_THREADS; ++i) {
    			tmpNodes[i] = sampleEdges(edgeList, vExist, m, i, psize);
			long size = tmpNodes[i].size();
			tmpNodeSizes[i] = size;
			__sync_fetch_and_add(&totalNodes,size);
    		}
    	}
   
    	const long VsSize = totalNodes;			//Size of Vs
    	long *Vs = new long[VsSize]();			//Array of sampled vertices
   	
   	//Calculating starting point for copying different vectors to Vs
   	for(long i=1; i<NUM_THREADS; ++i)
   		tmpNodeSizes[i] += tmpNodeSizes[i-1];
   	for(long i=1; i<NUM_THREADS; ++i)
   		tmpNodeStart[i] = tmpNodeSizes[i-1];
   	tmpNodeStart[0] = 0;
   	
   	dbg(printf("tmpNodeStartValues and sizes\n");
   	for(long i=0; i<NUM_THREADS; ++i)
   		printf("%ld %ld\n",tmpNodeStart[i],tmpNodes[i].size());)
   	
   	//Copy values to Vs in parallel from all vectors
   	#pragma omp parallel shared(Vs, tmpNodes, tmpNodeSizes)
   	{
   		#pragma omp for
   		for(long i=0; i<NUM_THREADS; ++i){
   			long start = tmpNodeStart[i];
   			long end = start+tmpNodes[i].size();
   			for(long j=start; j<end; ++j)
   				Vs[j] = tmpNodes[i][j-start];
   		}
   	}
   	
    	vbs(cout<<"Sampled edges, number of vertices="<<VsSize<<endl;)
		
	long totalEdges=0;				//number of edges induced
	vector<edge> tmpEdges[NUM_THREADS];		//vector for each thread inducing edges
	long *tmpEdgeSizes = new long[NUM_THREADS]();	//Number of induced edges for each thread
	long *tmpEdgeStart = new long[NUM_THREADS]();	//Starting index for copying edges from each thread
	psize = (VsSize+NUM_THREADS-1)/NUM_THREADS;	//partition size for each thread inducing edges
	
	//Induce edges in parallel
    	#pragma omp parallel shared(tmpEdges,tmpEdgeSizes,Vs,g,vExist,n,m,psize)
    	{
    		#pragma omp for
    		for(long i=0; i<NUM_THREADS; ++i){
    			tmpEdges[i] = induceEdges(Vs, g, vExist, VsSize, i, psize);
    			long size = tmpEdges[i].size();
			tmpEdgeSizes[i] = size;
			__sync_fetch_and_add(&totalEdges,size);
    		}
    	}
    
	const long EsSize = totalEdges;			//Size of Es
	edge *Es = new edge[totalEdges];		//Array of induced edges
	
	//Calculating starting point for copying different vectors to Es
   	for(long i=1; i<NUM_THREADS; ++i)
   		tmpEdgeSizes[i] += tmpEdgeSizes[i-1];
   	for(long i=1; i<NUM_THREADS; ++i)
   		tmpEdgeStart[i] = tmpEdgeSizes[i-1];
   	tmpEdgeStart[0] = 0;
   	
   	dbg(printf("tmpEdgeStartValues and sizes\n");
   	for(long i=0; i<NUM_THREADS; ++i)
   		printf("%ld %ld\n",tmpEdgeStart[i],tmpEdges[i].size());)
	 
    	//Copy values to Es in parallel from all vectors
   	#pragma omp parallel shared(Es, tmpEdges, tmpEdgeSizes)
   	{
   		#pragma omp for
   		for(long i=0; i<NUM_THREADS; ++i){
   			long start = tmpEdgeStart[i];
   			long end = start+tmpEdges[i].size();
   			dbg(printf("%ld %ld\n",start,end);)
   			for(long j=start; j<end; ++j)
   				Es[j] = tmpEdges[i][j-start];
   		}
   	}
   	    
	vbs(cout<<"Induced edges, final count is "<<VsSize<<" vertices, and "<<EsSize<<" edges"<<endl;)
    
    	//releasing memory before exiting function
    	delete vExist, tmpNodeSizes, tmpNodeStart, Vs, tmpEdgeSize, tmpEdgeStart;
    
	return pair<edge*,long>(Es,EsSize);
}

void parseCommandlineArgs(int argc, char* argv[], string &in, string &out, double &fi){
    if(argc != 4){
        cout<<"Usage: ties_sequential inputFilename outputFilename samplingRatio [-u (flag for undirected graph)]";
        exit(1);
    }
    in = argv[1];
    out = argv[2];
    fi = atof(argv[3]);
}

long inputGraph(string filename, vector<edge> &el){ //creates the edge list in the passed vector, and returns the number of nodes(maximum node found in any edge)

    ifstream in(filename, ifstream::in);
    string line;
    long a,b; //input nodes for the edge
    long maxNode = 0;
    long currentMax;

    while(in>>a>>b){
        el.push_back(edge(a,b));
        maxNode = max(maxNode,max(a,b));
    }
    in.close();
    
    dbg(printf("Last 2 elements in inputVector\n");
    long s = el.size();
    printf("%ld %ld\n%ld %ld\n",el[s-1].u,el[s-1].v,el[s-2].u,el[s-2].v);)
    
    return maxNode;
}

void make_csr(vector<edge> el, long n, long m){
    long i;
   	g.vi = new long[n+2]();
	g.ei = new long[m+1]();
    
    for(i=0; i<m; ++i){
        long u = el[i].u;
        long v = el[i].v;
        g.vi[u]++;
    }
    for(i=1; i<n+1; ++i){
        g.vi[i] += g.vi[i-1];
    }
    for(i=0; i<m; ++i){
        long u = el[i].u;
        long v = el[i].v;
        g.ei[g.vi[u]--] = v;
    }
    for(i=1; i<n; ++i)
        g.vi[i] += 1;
    g.vi[0] = 0;
    g.vi[n+1] = m+1;
    
    dbg(printf("Printing CSR\n");
    for(i=0; i<n; ++i){
		long start = g.vi[i];
		long stop  = g.vi[i+1];
		printf("For node %ld: start=%ld stop=%ld\n",i,start,stop);
		for(long j=start; j<stop; ++j){
			printf("%ld ",g.ei[j]);
		}
		printf("\n\n");
    })
}

vector<long> sampleEdges(vector<edge> &el, bool* vExist, long m, long pnum, long psize){
    default_random_engine generator((unsigned int)time(0));
    long rnd;
    long lastValid = el.size()-1;
    long start = pnum*psize;
    long end = start+psize-1;
    end = end>m-1 ? m-1: end;
    long realSize = end-start+1;  
	vector<long> Vtmp;
	long counter=0;
	long currentOld;
		
    while(start < end){
        uniform_int_distribution<long> distribution(start,end);
        rnd = distribution(generator);
        long u = el[rnd].u;
        long v = el[rnd].v;

        if(__sync_bool_compare_and_swap(vExist+u,false,true)){    
            Vtmp.push_back(u);
			++counter;
        }
   
        if(__sync_bool_compare_and_swap(vExist+v,false,true)){    
            Vtmp.push_back(v);
			++counter;
        }

        swap(el[end],el[rnd]);
        --end;
        
        if(counter >= UPDATE_COUNTER){
        	currentOld = __sync_fetch_and_add(&currentVsize,counter);
        	if(currentOld+counter > requiredVsize){
        		dbg(printf("killing thread %ld, with Vtmp size %ld\n",currentOld,Vtmp.size());)
        		break;
        	}
        	counter = 0;
        }
    }
    
    return Vtmp;
}

vector<edge> induceEdges(long *Vs, graph &g, bool* vExist, long VsSize, long pnum, long psize){
    long i,j;
    
    long start = pnum*psize;
    long end = start+psize-1;
    end = end>VsSize-1 ? VsSize-1: end;

    vector<edge> Etmp;
    
    for(i=start; i<end; ++i){
        long curV = Vs[i];
        long csrStart = g.vi[curV];
        long csrStop  = g.vi[curV+1];
        for(j=csrStart; j < csrStop; ++j){
            long desV = g.ei[j];
            if(vExist[desV])
                Etmp.push_back(edge(curV,desV));
        }
    }
    
    return Etmp;
}


void writeGraph(string filename,edge *Es, long EsSize){
    ofstream out(filename,ofstream::out);
    long i;
    dbg(printf("EsSize at end = %ld\n",EsSize);)
    for(i=0; i<EsSize; ++i){
        out<<Es[i].u<<" "<<Es[i].v<<'\n';
    }
    out.close();
}
