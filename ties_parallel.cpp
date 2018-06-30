/*  Author: Aditya Gaur
    Usage: ties_parallel inputFilename outputFilename samplingRatio [-u (flag for undirected graph)]
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

//For verbose comments
#define vbs(x) x
//Number of threads to use
#define NUM_THREADS 8
//Number of additions in a sampling thread, after which we make an update to the total
#define UPDATE_COUNTER 20
//Toggle Debugging
#define dbg(x) x

using namespace std;

class edge{
public:
    long u,v;
    edge(long a, long b){ //for consistency, smaller of the two will be put in u, larger in v
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

bool undirected = false;
long currentVsize=0;
long requiredVsize;

long inputGraph(string, vector<edge>&);
void make_csr(vector<edge>,graph&,long,long);
vector<long> sampleEdges(vector<edge>&, bool*, long, long, double, long, long);
vector<edge> induceEdges(long*, graph&, bool*, long, long, long);
void writeGraph(string,edge*, long);
void parseCommandlineArgs(int,char* [], string&, string&, double&);

int main(int argc, char* argv[]) {

    //Parameters to be specified
    string inFilename;
    string outFilename;
    double fi;
    parseCommandlineArgs(argc, argv, inFilename, outFilename, fi);

    vector<edge> edgeList;

    long n, m;

    vbs(cout<<"Getting input"<<endl;)

    n = inputGraph(inFilename,edgeList);
    m = edgeList.size();
    vbs(cout<<"Input succeeded, graph has "<<n<<" nodes and "<<m<<" edges"<<endl;)
    vbs(if(undirected)cout<<"Considering input as undirected graph"<<endl;)

	graph g;
	g.vi = new long[n+2]();
	if(undirected)
		g.ei = new long[2*m+1]();
	else
		g.ei = new long[m+1]();
		
    bool *vExist = new bool[n+1]();
    long totalNodes=0;
    requiredVsize = fi*n;


	make_csr(edgeList, g, n, m);
	
	long psize = m/(NUM_THREADS-1);
	vector<long> tmpNodes[NUM_THREADS];	//individual vectors for each thread
	long *tmpNodeSizes = new long[NUM_THREADS]();
	long *tmpNodeStart = new long[NUM_THREADS]();
	
	//Sample edges in parallel
    #pragma omp parallel shared(edgeList, vExist, n, m, fi, psize, tmpNodes, tmpNodeSizes) num_threads(NUM_THREADS)
    {
    	#pragma omp for
    	for(long i=0; i<NUM_THREADS; ++i) {
    		tmpNodes[i] = sampleEdges(edgeList, vExist, n, m, fi, i, psize);
			long size = tmpNodes[i].size();
			tmpNodeSizes[i] = size;
			__sync_fetch_and_add(&totalNodes,size);
    	}
    }
   
    const long VsSize = totalNodes;
    long *Vs = new long[VsSize]();
   	
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
   	
    vbs(cout<<"Sampled edges, #nodes="<<VsSize<<endl;)
	
	dbg(
	long count = 0;
	for(long i=0; i<n; ++i){
		if(vExist[i])
			count++;
	}
	if(count!=VsSize)
		printf("# of unique nodes = %ld, VsSize = %ld\n",count,VsSize);)
		
	long totalEdges=0;
	vector<edge> tmpEdges[NUM_THREADS];
	long *tmpEdgeSizes = new long[NUM_THREADS]();
	long *tmpEdgeStart = new long[NUM_THREADS]();
	psize = VsSize/(NUM_THREADS-1);
	
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
    
	const long EsSize = totalEdges;
	edge *Es = new edge[totalEdges];
	
	dbg(cout<<"Induced"<<endl;)
	
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
   	    
    vbs(cout<<"Induced edges, final count is "<<VsSize<<" nodes, and "<<EsSize<<" edges"<<endl;)

    writeGraph(outFilename, Es, EsSize);
    vbs(cout<<"Output graph written to "<<outFilename<<endl;)

    return 0;
}

void parseCommandlineArgs(int argc, char* argv[], string &in, string &out, double &fi){
    if(argc < 4 || argc > 5){
        cout<<"Usage: ties_sequential inputFilename outputFilename samplingRatio [-u (flag for undirected graph)]";
        exit(1);
    }
    in = argv[1];
    out = argv[2];
    fi = atof(argv[3]);
    if(argc == 5){
        string u = argv[4];
        if(u.compare("-u")==0)
            undirected = true;

    }
}

long inputGraph(string filename, vector<edge> &el){ //creates the edge list in the passed vector, and returns the number of nodes(maximum node found in any edge)

    ifstream in(filename, ifstream::in);
    string line;
    long a,b; //input nodes for the edge
    long maxNode = 0;
    long currentMax;

    while(!in.eof()){
        in>>a>>b;
        el.push_back(edge(a,b));
        maxNode = max(maxNode,max(a,b));
    }
    in.close();
    return maxNode;
}

void make_csr(vector<edge> el, graph &g, long n, long m){
    long i;
    for(i=0; i<m; ++i){
        long u = el[i].u;
        long v = el[i].v;
        g.vi[u]++;
        if(undirected)
            g.vi[v]++;
    }
    for(i=1; i<n+1; ++i){
        g.vi[i] += g.vi[i-1];
    }
    for(i=0; i<m; ++i){
        long u = el[i].u;
        long v = el[i].v;
        g.ei[g.vi[u]--] = v;
        if(undirected)
            g.ei[g.vi[v]--] = u;
    }
    for(i=1; i<n; ++i)
        g.vi[i] += 1;
    g.vi[0] = 0;
    if(undirected)
        g.vi[n+1] = 2*m+1;
    else
        g.vi[n+1] = m+1;
}

vector<long> sampleEdges(vector<edge> &el, bool* vExist, long n, long m, double fi, long pnum, long psize){
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
        bool expected = false;
        //if(vExist[u].compare_exchange_strong(expected,true)){
        if(__sync_bool_compare_and_swap(vExist+u,false,true)){    
            //vExist[u] = true;
            Vtmp.push_back(u);
        }
        expected = false;
        //if(vExist[v].compare_exchange_strong(expected,true)){
        if(__sync_bool_compare_and_swap(vExist+v,false,true)){    
            //vExist[v] = true;
            Vtmp.push_back(v);
        }

        swap(el[end],el[rnd]);
        --end;
        ++counter;
        
        if(counter >= UPDATE_COUNTER){
        	currentOld = __sync_fetch_and_add(&currentVsize,counter);
        	if(currentOld+counter > requiredVsize)
        		break;
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
        out<<Es[i].u<<" "<<Es[i].v<<endl;
    }
    out.close();
}
