
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
#define MAXITERS 3000
#define ELITESETSIZE 10
#define K 150
#define EPOCH 30

int vertex_num, num_edges;

bool deleted[maxn], visited[maxn];

vector<vector<int> > adjlist;

vector<int> nodes_degree;

vector<int> V;




int selectNodeWithHigherWeigth(const vector<int> &comp, map<int,int> &weights){

    int max = -1;
    int idx = 0;
    for (int i = 0; i < comp.size(); i++)
    {
        if(weights[comp[i]] > max){
            idx = comp[i];
            max = weights[comp[i]];
        } else if (weights[comp[i]] == max and degree(comp[i]) > degree(idx)){
            idx = comp[i];
        }
    }
    
    return idx;

}



void compBasedNeighbSearch(){
    
    int bestSolCost, SolCost;
    int iter = 0;
    int count = 0;

    SolCost = compute_cost_func();
    bestSolCost = SolCost;

    map<int,int> weights;
    for (int i = 0; i < vertex_num; i++)
    {
        weights.insert(make_pair(i,0));
    }
    
    while(iter < MAXITERS){

        vector<vector<int> > listcc = obtainListConnectedComponents();

        int maxCC = 0, minCC = vertex_num, size;

        for(int i = 0; i < listcc.size(); ++i){
            size = listcc[i].size();
            if (size <= minCC) minCC = size;
            if (size >= maxCC) maxCC = size;
        }

        vector<vector<int> > nodesLargeCC;

        for(int i = 0; i < listcc.size(); ++i){
            if(listcc[i].size() >= (maxCC + minCC)/2){
                nodesLargeCC.push_back(listcc[i]);
            } 
        }

        vector<int> comp = nodesLargeCC[rand()%(nodesLargeCC.size())];

        int idx = selectNodeWithHigherWeigth(comp, weights);

        for (int i = 0; i < comp.size(); i++){
            if(comp[i] != idx) weights[comp[i]]++;
        }

       

        int w = fastWorstNode();

        deleted[idx] = true;

        deleted[w] = false;

        int incr = compute_incr(w);
        int decr = compute_decr(idx);

        SolCost += incr - decr;

        if(SolCost <= bestSolCost){
            iter = 0;
            weights[w] = 0;
            bestSolCost = SolCost;
        } else {
            deleted[idx] = false;
            deleted[w] = true; 
            iter++;
            SolCost -= incr - decr;
        }
    }
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



void feasibleSolution(){
    vector<int> shuffled;
    for(int i = 0; i < vertex_num; ++i){
        shuffled.push_back(i);
    }

    shuffle(shuffled);

    fill(deleted, deleted + vertex_num, false);

    for(int i = 0; i < K; ++i) deleted[shuffled[i]] = true;

}


void ComponentBased(){

    feasibleSolution();
    compBasedNeighbSearch();
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

    int costs = 0;
    double exectime = 0;

    int costiter;

    int bestCost = vertex_num*(vertex_num-1)/2;

    for(int j = 0; j < EPOCH ; ++j){
        
        auto start = high_resolution_clock::now();
        ComponentBased();
        auto stop = high_resolution_clock::now();

        auto duration = duration_cast<milliseconds>(stop - start);
        exectime += double(duration.count())/1000;
        costiter = compute_cost_func();
        if(costiter < bestCost) bestCost = costiter;
        costs += costiter;

        cout<<"epoch "<<j<<" done"<<endl;
    }

    cout << "avg Time taken by function: " <<exectime/EPOCH << " seconds" << endl;

    cout<<endl<<"---------------------------------"<<endl;
    cout<<"avg Cost :"<<costs/EPOCH<<endl;
    cout<<"best cost: "<<bestCost<<endl;


}