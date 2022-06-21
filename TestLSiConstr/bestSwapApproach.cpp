
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
#define MAXITERS 1000
#define ELITESETSIZE 10
#define K 150
#define EPOCH 10
#define I 2000

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
        swap(vec[i], vec[randomIdx]);
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


void BestSwapApproach(){

    //MISGreedy();
    feasibleSolution();
    computeConf();
    int iter = 0;
    while(bestSwap()){
       cout<<computeCost()<<endl;
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

    int costs = 0;
    double exectime = 0;

    int costiter;

    int bestCost = vertex_num*(vertex_num-1)/2;

    for(int j = 0; j < EPOCH ; ++j){
        
        auto start = high_resolution_clock::now();
        BestSwapApproach();
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