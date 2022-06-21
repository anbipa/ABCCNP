
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <stack>
#include <set>
#include <map>
#include <queue>
#include <limits>
#include <chrono>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "algorismes_comuns.h"
using namespace std;

using namespace std::chrono;



#define maxn 30100


int vertex_num, num_edges;

bool deleted[maxn], visited[maxn];

vector<vector<int> > adjlist;

vector<int> nodes_degree;

vector<int> V;

int vertex_num_visited;

int vertex_list[maxn];

Conf C;

int weight[maxn];


bool articulation[maxn];
int artcost[maxn];

vector<int> S, bestS;

int SolCost, BestSolCost;


//parameters

int candSize = 5;
double prob = 0.5;
int K = 150;
string input_file = "BA5000";


void read_parameters(int argc, char **argv) {

    int iarg=1;
    while (iarg < argc) {
        if (strcmp(argv[iarg],"-i")==0) input_file = argv[++iarg];
        else if (strcmp(argv[iarg],"-K")==0) K = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-cs")==0) candSize = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-p")==0) prob = atof(argv[++iarg]);
        iarg++;
    }
}



int nonDetMaxDegree(const vector<int> &ndegree, const vector<bool> &cover){

    //random number between 0 and 1
    double r = (double)rand() / (RAND_MAX);
    cout<<"r: "<<r<<endl;

    int returnNode = -1;

    if (r < prob){

        int maxdeg = -1;
        int maxnode = 0;

        for(int i = 0; i < vertex_num; ++i){
            if(!cover[i] and ndegree[i] > maxdeg){
                maxdeg = ndegree[i];
                maxnode = i;
            }
        }

        returnNode = maxnode;
        
    } else {

        vector<int> maxdeg, maxnode;

        for(int i = 0; i < vertex_num; ++i){
            if(!cover[i]){
                if(maxdeg.size() < candSize){
                    maxdeg.push_back(ndegree[i]);
                    maxnode.push_back(i);
                    break;
                } else {
                    for(int k = 0; k < candSize; ++k){
                        if(maxdeg[k] < ndegree[i]){
                            maxdeg[k] = ndegree[i];
                            maxnode[k] = i;
                            break;
                        }
                    }
                }
            }
        }

        int d = rand() % maxnode.size();
        returnNode = maxnode[d];
    }

    return returnNode;
}


vector<bool> vertex_cover_sorted(){

    vector<bool> cover(vertex_num, false);

    int nedg = num_edges;
    cout<<"nedg: "<<nedg<<endl;

    vector<int> ndegree = nodes_degree;


    while(nedg > 0){

        int v = nonDetMaxDegree(ndegree, cover);
        cover[v] = true;
       
        for (int i = 0; i < adjlist[v].size(); i++)
        {
            if(!cover[adjlist[v][i]]) {
                nedg--;
            }
            ndegree[adjlist[v][i]]--;
        }
    }

    return cover;
}


void VCoverOrdered(){

    vector<bool> vcover = vertex_cover_sorted();

    int n = 0;

    fill(deleted, deleted + vertex_num, false);

    for(int i = 0; i < vertex_num; ++i){
        if (vcover[i]){
            ++n;
            deleted[i] = true;
            S.push_back(i);
        }  
    }

    int cost = computeConf();

    while (n > K)
    {
        int i = wnodei(S);
        int w = S[i];
        S.erase(S.begin() + i);
        restoreVertex(w);
        --n;
    }
}


int main( int argc, char **argv ){

    read_parameters(argc, argv);

    read_graph("../cnd/" + input_file + ".txt");

    srand(time(0));

    auto start = high_resolution_clock::now();
    VCoverOrdered();
    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(stop - start);
    cout.precision(5);
    cout << "Time taken by function: " <<float(duration.count())/1000 << " seconds" << endl;

    for(int i = 0; i < vertex_num; ++i){
        if(deleted[i]) cout<<i<<" ";
    }
    cout<<endl<<"---------------------------------"<<endl;
    cout<<"Cost :"<<compute_cost_func()<<endl;
    cout<<"size :"<<S.size()<<endl;

}