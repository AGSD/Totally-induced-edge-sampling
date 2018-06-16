/*  Author: Aditya Gaur
    Usage: ties_sequential inputFilename outputFilename samplingRatio
    About: Code to create a sample graph from an input graph with final number of nodes based on a sampling ratio.
           Algorithm used is TIES (Totally induced edge sampling), based on https://docs.lib.purdue.edu/cgi/viewcontent.cgi?article=2743&context=cstech
*/
#include<iostream>
#include<fstream>
#include<vector>
#include<algorithm>
#include<set>
#include<random>
#include<string>
#include<cstdlib>

//For verbose comments
#define vbs(x) x

using namespace std;

class edge{
public:
    int u,v;
    edge(int a, int b){ //for consistency, smaller of the two will be put in u, larger in v
        u = min(a,b);
        v = max(a,b);
    }
};

struct CustomCompare
{
    bool operator()(const edge& e1, const edge& e2)
    {
         if(e1.u == e2.u){
            if(e1.v < e2.v){
                return true;
            }
            else
                return false;
        }
        else if(e1.u < e2.u)
            return true;
        else
            return false;
    }
};

long inputGraph(string, vector<edge>&);
void sampleEdges(const vector<edge>&, set<long>&, long, long, double);
void induceEdges(const vector<edge>&, set<long>&, set<edge,CustomCompare>&, long);
void writeGraph(string,set<edge,CustomCompare>&);
void parseCommandlineArgs(int,char* [], string&, string&, double&);

int main(int argc, char* argv[]) {

    //Parameters to be specified
    string inFilename;
    string outFilename;
    double fi;
    parseCommandlineArgs(argc, argv, inFilename, outFilename, fi);


    vector<edge> edgeList;
    set<long> Vs;
    set<edge,CustomCompare> Es;

    long n, m;

    vbs(cout<<"Getting input"<<endl;)

    n = inputGraph(inFilename,edgeList);
    m = edgeList.size();
    vbs(cout<<"Input succeeded, graph has "<<n<<" nodes and "<<m<<" edges"<<endl;)

    sampleEdges(edgeList,Vs,n,m,fi);
    vbs(cout<<"Sampled edges"<<endl;)

    induceEdges(edgeList, Vs, Es, m);
    vbs(cout<<"Induced edges, final count is "<<Es.size()<<" edges"<<endl;)

    writeGraph(outFilename, Es);
    vbs(cout<<"Output graph written to "<<outFilename<<endl;)

    return 0;
}

void parseCommandlineArgs(int argc, char* argv[], string &in, string &out, double &fi){
    if(argc != 4){
        cout<<"Usage: ties_sequential inputFilename outputFilename samplingRatio";
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

    while(!in.eof()){
        in>>a>>b;
        el.push_back(edge(a,b));
        maxNode = max(maxNode,max(a,b));
    }
    in.close();
    return maxNode;
}

void sampleEdges(const vector<edge> &el, set<long> &Vs, long n, long m, double fi){
    default_random_engine generator;
    uniform_int_distribution<long> distribution(1,m);
    long rnd;

    while(Vs.size() < fi*n){
        rnd = distribution(generator);
        Vs.insert(el[rnd].u);
        Vs.insert(el[rnd].v);
    }
}

void induceEdges(const vector<edge>&el, set<long>&Vs, set<edge,CustomCompare>&Es, long m){
    long a,b;
    long i;

    for(i=0; i<m; ++i){
        a = el[i].u;
        b = el[i].v;
        if(Vs.find(a) != Vs.end() && Vs.find(b) != Vs.end()){   //edge is present
            Es.insert(edge(a,b));
        }
    }
}

void writeGraph(string filename,set<edge,CustomCompare> &Es){
    ofstream out(filename,ofstream::out);
    set<edge,CustomCompare>::iterator itr;
    for(itr = Es.begin(); itr!=Es.end(); ++itr){
        out<<itr->u<<" "<<itr->v<<endl;
    }
    out.close();
}

