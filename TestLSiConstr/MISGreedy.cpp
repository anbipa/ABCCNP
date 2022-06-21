
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
#define K 200


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


int minDegree(const vector<bool> &labeled, const vector<int> &degree){
    int min = vertex_num+1;
    int mini = 0;
    for(int i = 0; i < vertex_num; ++i){
        if(!labeled[i] and degree[i] < min){
            min = degree[i];
            mini = i;
        }
    }

    return mini;
}

vector<int> MaxIndepSet(){
    vector<bool> labeled(vertex_num, false);
    int nlabeled = 0;
    vector<int> MIS;
    vector<int> degree = nodes_degree;

    while(nlabeled != vertex_num){
        int u = minDegree(labeled, degree);
        labeled[u] = true;

        MIS.push_back(u);

        ++nlabeled;

        for(auto& elem: adjlist[u]){
            if(!labeled[elem]) ++nlabeled;
            labeled[elem] = true;
            for(auto& n: adjlist[elem]){
                degree[n]--;
            }
        }
    }

    return MIS;

}

// que passa amb el node 0?? (millorar greedy MIS) 
// que passa quan K > |MIS|??
vector<int> MISGreedy(){

    vector<int> MIS = MaxIndepSet();
    fill(deleted, deleted + vertex_num, true);
    
    for(auto& elem: MIS){
        deleted[elem] = false;
    }

    S.clear();
    
    for(int i = 0; i < vertex_num; ++i){
        if(deleted[i]){
            S.push_back(i);
        }
    }

    int SolCost = computeConf();


    while(MIS.size() != vertex_num - K){
        int i = wnodei(S);
        int v = S[i];
        swap(S[i], S[S.size()-1]);

        S.pop_back();

        MIS.push_back(v);
        deleted[v] = false;

        restoreVertex(v);
    }


    return S;
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

    auto start = high_resolution_clock::now();
    MISGreedy();
    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(stop - start);
    cout.precision(5);
    cout << "Time taken by function: " <<float(duration.count())/1000 << " seconds" << endl;

    for(int i = 0; i < vertex_num; ++i){
        if(deleted[i]) cout<<i<<" ";
    }
    cout<<endl<<"---------------------------------"<<endl;
    cout<<"Cost :"<<compute_cost_func()<<endl;

}

