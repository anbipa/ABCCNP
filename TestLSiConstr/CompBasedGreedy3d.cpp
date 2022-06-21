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


string input_file;
double deltaK;
int K;
int lIters;
int iIters;
int epochs;


void read_parameters(int argc, char **argv) {

    int iarg=1;
    while (iarg < argc) {
        if (strcmp(argv[iarg],"-i")==0) input_file = argv[++iarg];
        else if (strcmp(argv[iarg],"-dk")==0) deltaK = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-K")==0) K = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-lIters")==0) lIters = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-iIters")==0) iIters = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-epochs")==0) epochs = atoi(argv[++iarg]);
        iarg++;
    }
}



bool comp(const pair<int,int> &p1, const pair<int,int> &p2){
    return p1.second > p2.second;
}

int nodeMaxDegree(const vector<int> &ndegree, const vector<bool> &cover){
    int maxdeg = -1;
    int maxnode = 0;

    for(int i = 0; i < vertex_num; ++i){
        if(!cover[i] and ndegree[i] > maxdeg){
            maxdeg = ndegree[i];
            maxnode = i;
        }
    }

    return maxnode;
}


vector<bool> vertex_cover_sorted(){

    vector<bool> cover(vertex_num, false);

    int nedg = num_edges;

    vector<int> ndegree = nodes_degree;


    while(nedg > 0){

        int v = nodeMaxDegree(ndegree,cover);

        cover[v] = true;

        for (int i = 0; i < adjlist[v].size(); i++)
        {
            if(!cover[adjlist[v][i]]){
                nedg--;
                ndegree[adjlist[v][i]]--;
            } 
        }
    }

    return cover;
}

void shuffle(vector<int> &vec){
    int l = vec.size();
    int temp = 0;
    int randomIdx = 0;
    for(int i = 0; i < l; ++i){
        randomIdx = rand()%l;
        temp = vec[i];
        vec[i] = vec[randomIdx];
        vec[randomIdx] = temp;
    }

}


vector<bool> vertex_cover_shuffled(){

    vector<bool> cover(vertex_num, false);

    vector<int> shuffled;

    for (int i = 0; i < vertex_num; i++)
    {
        shuffled.push_back(i);
    }
    
    shuffle(shuffled);

    int nedg = num_edges;

    int j = 0;
    while(nedg > 0){

        int v = shuffled[j];

        if(nodes_degree[v]==1) --nedg;
        
        else cover[v] = true;

        for (int i = 0; i < adjlist[v].size(); i++)
        {
            if(!cover[adjlist[v][i]] and nodes_degree[adjlist[v][i]] != 1) --nedg;

        }
        ++j;
    }

    return cover;
}

void greedy3d(){
    int n = 0;

    vector<bool> cover = vertex_cover_shuffled();

    S.clear();

    for(int i = 0; i < vertex_num; ++i){
        if(cover[i]){
            deleted[i] = true;
            S.push_back(i);
        }
        else deleted[i] = false;
    }

    SolCost = computeConf();
    BestSolCost = vertex_num*(vertex_num-1)/2;
    bestS = S;
    do
    {
        while(S.size() > K - K*deltaK){
            int iu = wnodei(S);
            int u = S[iu];
            restoreVertex(u);
            S.erase(S.begin()+iu);
            if(S.size() == K){
                if(SolCost <= BestSolCost){
                    BestSolCost = SolCost;
                    bestS = S;
                }
            }
        }
        ++n;
        while(S.size() < K + K*deltaK){
            
            
            int cc = randomLargeCC();
            //vector<int> & list = C.connectedComponents[cc];
            //int u = list[rand()%list.size()];
            int u, bestDecr;
            biggesDecrement(C.connectedComponents[cc][0], u, bestDecr);
            S.push_back(u);
            deleteVertex(u);
            if(S.size() == K){
                if(SolCost <= BestSolCost){
                    BestSolCost = SolCost;
                    bestS = S;
                }
            }
        }
        ++n;
    } while (n<=lIters);
    fill(deleted, deleted+vertex_num, false);
    for(int i = 0; i < bestS.size(); ++i){
        deleted[bestS[i]] = true;
    }
}


int main( int argc, char **argv ){

    read_parameters(argc, argv);

    string filename;

    filename = "../cnd/" + input_file + ".txt";
    
    read_graph(filename);

    srand(time(0));

    //repeat 30 times
    float total_time = 0;
    int total_cost = 0;
    int best_cost = numeric_limits<int>::max();
    int cost = 0;
    for(int i = 0; i < epochs; ++i){
        auto start = high_resolution_clock::now();
        greedy3d();
        auto stop = high_resolution_clock::now();

        auto duration = duration_cast<milliseconds>(stop - start);
        total_time += duration.count();
        cost = compute_cost_func();
        if (cost < best_cost) best_cost = cost;
        total_cost += cost;
        cout<<"EPOCH "<<i<<" cost: "<<cost<<endl;
    }

    //print avg time and cost
    cout<<"Average time: "<<(float(total_time)/1000)/epochs<<" seconds"<<endl;
    cout<<"Average cost: "<<float(total_cost)/epochs<<endl;

    //print best cost
    cout<<"Best cost: "<<best_cost<<endl;
    return best_cost;

}
