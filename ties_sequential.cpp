/*  Author: Aditya Gaur
    Usage: ties_sequential inputFilename outputFilename samplingRatio
    About: Code to create a sample graph from an input graph with final number of nodes based on a sampling ratio.
           Algorithm used is TIES (Totally induced edge sampling), based on https://docs.lib.purdue.edu/cgi/viewcontent.cgi?article=2743&context=cstech
    Compile with: g++ -std=c++11 ties_sequential.cpp -o ties_sequential.o
*/
#include<iostream>
#include<fstream>
#include<vector>
#include<random>
#include<string>
#include<ctime>

//For verbose comments
#define vbs(x) x

using namespace std;

class edge{
public:
    int u,v;
    edge(int a, int b){ //for consistency, smaller of the two will be put in u, larger in v
        u = a;
        v = b;
    }
};

struct graph{
	int *vi;
	int *ei;
};

bool undirected = false;

long inputGraph(string, vector<edge>&);
void make_csr(vector<edge>,graph&,long,long);
void sampleEdges(vector<edge>&, vector<long>&, bool*, long, long, double);
void induceEdges(vector<long>&, graph&, bool*, vector<edge>&, long);
void writeGraph(string,vector<edge>&);
void parseCommandlineArgs(int,char* [], string&, string&, double&);

int main(int argc, char* argv[]) {

    //Parameters to be specified
    string inFilename;
    string outFilename;
    double fi;
    parseCommandlineArgs(argc, argv, inFilename, outFilename, fi);

    vector<edge> edgeList;
    vector<long> Vs;	//del

    long n, m;

    vbs(cout<<"Getting input"<<endl;)

    n = inputGraph(inFilename,edgeList);
    m = edgeList.size();
    vbs(cout<<"Input succeeded, graph has "<<n<<" nodes and "<<m<<" edges"<<endl;)
    vbs(if(undirected)cout<<"Considering input as undirected graph"<<endl;)

	graph g;
	g.vi = new int[n+2]();
	if(undirected)
		g.ei = new int[2*m+1]();
	else
		g.ei = new int[m+1]();
	
    bool *vExist = new bool[n+1]();

	make_csr(edgeList, g, n, m);

    sampleEdges(edgeList, Vs, vExist, n, m, fi);
    vbs(cout<<"Sampled edges"<<endl;)

    vector<edge> Es;

    induceEdges(Vs, g, vExist, Es, n);
    vbs(cout<<"Induced edges, final count is "<<Vs.size()<<" nodes, and "<<Es.size()<<" edges"<<endl;)

    writeGraph(outFilename, Es);
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

    while(in>>a>>b){
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

void sampleEdges(vector<edge> &el, vector<long> &Vs, bool* vExist, long n, long m, double fi){
    default_random_engine generator(time(0));
    long rnd;
    long lastValid = el.size()-1;

    while(Vs.size() < fi*n){
        uniform_int_distribution<long> distribution(0,lastValid);
        rnd = distribution(generator);
        long u = el[rnd].u;
        long v = el[rnd].v;
        if(!vExist[u]){
            vExist[u] = true;
            Vs.push_back(u);
        }
        if(!vExist[v]){
            vExist[v] = true;
            Vs.push_back(v);
        }

        swap(el[lastValid],el[rnd]);
        --lastValid;
    }
}

void induceEdges(vector<long> &Vs, graph &g, bool* vExist, vector<edge> &Es, long n){
    bool *vChecked = new bool[n+1]();
    long i,j;
    for(i=0; i<Vs.size(); ++i){
        long curV = Vs[i];
        long start = g.vi[curV];
        long stop  = g.vi[curV+1];
        for(j=start; j<stop; ++j){
            long desV = g.ei[j];
            if(!undirected)
                Es.push_back(edge(curV,desV));
            else if(!vChecked[desV] && vExist[desV])
                Es.push_back(edge(curV,desV));
        }
        vChecked[curV] = true;
    }
}


void writeGraph(string filename,vector<edge> &Es){
    ofstream out(filename,ofstream::out);
    int i;
    for(i=0; i<Es.size(); ++i){
        out<<Es[i].u<<" "<<Es[i].v<<'\n';
    }
    out.close();
}
