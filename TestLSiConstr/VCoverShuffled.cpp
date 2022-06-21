
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
#define MAXITERS 100
#define ELITESETSIZE 10
#define K 80

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

void VCoverShufled(){

    vector<bool> vcover = vertex_cover_shuffled();

    int n = 0;

    fill(deleted, deleted + vertex_num, false);
    S.clear();

    for(int i = 0; i < vertex_num; ++i){
        if (vcover[i]) {
            deleted[i] = true;
            S.push_back(i);
            ++n;
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


int main(){

    string filename, file;

    cout<<"specify network"<<endl;
    cin>>file;

    filename = "../cnd/" + file + ".txt";
    
    read_graph(filename);

    //construction();
    
    //show_adjlist();
    //show_critical_nodes();

    srand(time(0));

    //repeat 30 times
    float total_time = 0;
    int total_cost = 0;
    int best_cost = numeric_limits<int>::max();
    int cost = 0;
    for(int i = 0; i < 10; ++i){
        auto start = high_resolution_clock::now();
        VCoverShufled();
        auto stop = high_resolution_clock::now();

        auto duration = duration_cast<milliseconds>(stop - start);
        total_time += duration.count();
        cost = compute_cost_func();
        if (cost < best_cost) best_cost = cost;
        total_cost += cost;
    }

    //print avg time and cost
    cout<<"Average time: "<<(float(total_time)/1000)/10<<" seconds"<<endl;
    cout<<"Average cost: "<<float(total_cost)/10<<endl;

    //print best cost
    cout<<"Best cost: "<<best_cost<<endl;

   
    for(int i = 0; i < vertex_num; ++i){
        if(deleted[i]) cout<<i<<" ";
    }
    cout<<endl<<"---------------------------------"<<endl;

}